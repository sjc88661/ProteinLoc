"""FASTA、第三方原始结果和候选信号。所有残基位置采用 1-based 闭区间。"""
import ast
import csv
import io
import json
import math
import os
import re
import uuid
from collections import Counter
from pathlib import Path

STANDARD_AA = set('ACDEFGHIKLMNPQRSTVWY')
ALLOWED_AA = STANDARD_AA | set('XBZJUO')
RULE_VERSION = '2026-11-10.1'
MODULES = ('signalp', 'phobius', 'deeploc', 'targetp', 'predgpi', 'deeptmhmm', 'er', 'pts', 'pts1_svm')
PRO_LOCALIZATIONS = ('Cell wall & surface','Extracellular','Cytoplasmic',
                     'Cytoplasmic Membrane','Outer Membrane','Periplasmic')


def read_fasta(text):
    text = text.lstrip('\ufeff').strip()
    if not text:
        raise ValueError('请输入 FASTA 或蛋白序列')
    plain = not text.startswith('>')
    if plain:
        text = '>Sequence_1\n' + text
    entries, header, chunks = [], None, []
    for line in text.splitlines():
        if line.startswith('>'):
            if header is not None:
                entries.append((header, ''.join(chunks)))
            header, chunks = line[1:].strip(), []
        elif line.strip():
            chunks.append(''.join(line.split()))
    entries.append((header, ''.join(chunks)))
    counts = Counter((h.split()[0] if h else '') for h, _ in entries)
    records = []
    for index, (header, source_sequence) in enumerate(entries, 1):
        original_id = header.split()[0] if header else ''
        sequence = source_sequence.upper()
        warnings, errors = [], []
        if plain:
            warnings.append('未提供 FASTA 标题，分配 Sequence_1')
        if not header:
            errors.append('FASTA 标题为空')
        if not sequence:
            errors.append('序列为空')
        bad = set(sequence) - ALLOWED_AA
        if bad:
            errors.append('含非法字符：' + ''.join(sorted(bad)))
        if counts[original_id] > 1:
            warnings.append('原始 ID 重复，已分配不同内部 ID')
        if set(sequence) & set('XBZJUO'):
            warnings.append('含未知/特殊残基；各模型的实际提交序列和变更另行记录')
        records.append(dict(internal_id=f'P{index:08d}', original_id=original_id,
                            header=header, sequence=sequence, source_sequence=source_sequence,
                            sequence_length=len(sequence), warnings=warnings, errors=errors))
    return records


def _number(value):
    try:
        n = float(value)
        return n if math.isfinite(n) else None
    except (TypeError, ValueError):
        return None


def _table(path, delimiter='\t'):
    lines = Path(path).read_text(encoding='utf-8-sig').splitlines()
    header_at = None
    for i, line in enumerate(lines):
        candidate = line.lstrip('#').strip()
        fields = next(csv.reader([candidate], delimiter=delimiter))
        if any(field.strip() in ('ID','Protein_ID','Sequence_ID','ACC') for field in fields[:2]):
            header_at = i
            break
    if header_at is None:
        raise ValueError(f'无法识别结果表头：{Path(path).name}')
    lines[header_at] = lines[header_at].lstrip('#').strip()
    return csv.DictReader(io.StringIO('\n'.join(lines[header_at:])), delimiter=delimiter)


def _insert(output, key, data):
    key = str(key or '').strip().split()[0] if str(key or '').strip() else ''
    if not key or key in output:
        raise ValueError(f'模型输出 ID 为空或重复：{key!r}')
    output[key] = dict(status='completed', **data)


def _scores(row):
    return {k: n for k, v in row.items() if k is not None and (n := _number(v)) is not None}


def _cleavage(row):
    value = next((v for k, v in row.items() if k and ('CS' in k or 'Cleavage' in k)), '') or ''
    match = re.search(r'(?:pos:\s*)?(\d+)\s*-\s*(\d+)', value)
    return int(match[1]) if match else None


def parse_signalp(path):
    result = {}
    for row in _table(path):
        prediction = row.get('Prediction')
        if not prediction:
            raise ValueError('SignalP 缺少 Prediction')
        _insert(result, row.get('ID'), dict(prediction=prediction,
                signal_peptide=prediction.upper() != 'OTHER', cleavage_after=_cleavage(row),
                scores=_scores(row), raw=row))
    return result


def parse_targetp(path):
    result = {}
    for row in _table(path):
        if not row.get('Prediction'):
            raise ValueError('TargetP 缺少 Prediction')
        prediction = {'MT':'mTP','CH':'cTP','TH':'lTP'}.get(row['Prediction'],row['Prediction'])
        _insert(result, row.get('ID'), dict(prediction=prediction, raw_prediction=row['Prediction'],
                cleavage_after=_cleavage(row), scores=_scores(row), raw=row))
    return result


def _labels(value):
    value = (value or '').strip()
    if value.startswith('['):
        try:
            parsed = ast.literal_eval(value)
            if isinstance(parsed, list):
                return [str(x) for x in parsed]
        except (ValueError, SyntaxError):
            pass
    return [s.strip() for s in re.split('[;,|]', value) if s.strip()]


def parse_deeploc(path):
    result = {}
    for row in _table(path, ','):
        localizations = row.get('Localizations', row.get('Localization'))
        if not localizations:
            raise ValueError('DeepLoc 缺少定位标签')
        _insert(result, row.get('Protein_ID', row.get('ID')), dict(localizations=_labels(localizations),
                signals=_labels(row.get('Signals')), membrane_types=_labels(row.get('Membrane types',row.get('Membrane_types'))),
                scores=_scores(row), raw=row))
    return result


def parse_deeplocpro(path):
    result = {}
    for row in _table(path, ','):
        label = row.get('Localization')
        scores = {name: _number(row.get(name)) for name in PRO_LOCALIZATIONS}
        if label not in PRO_LOCALIZATIONS or any(v is None or not 0 <= v <= 1 for v in scores.values()):
            raise ValueError('DeepLocPro 缺少有效定位标签或六类概率')
        if abs(sum(scores.values()) - 1) > .002:
            raise ValueError('DeepLocPro 六类概率之和异常')
        _insert(result, row.get('ACC'), dict(tool='DeepLocPro',localizations=[label],
                prediction_type='multiclass',scores=scores,raw=row))
    return result


def parse_phobius(path):
    result = {}
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        if not line.strip() or line.startswith('#') or 'PREDICTION' in line:
            continue
        parts = line.split()
        if len(parts) != 4 or not parts[1].isdigit() or parts[2] not in ('Y', 'N', '0', '1'):
            raise ValueError('Phobius 行格式不合法：' + line[:120])
        key, count, sp, topology = parts
        regions = [[int(a), int(b)] for a, b in re.findall(r'[io](\d+)-(\d+)(?=[io]|$)', topology)]
        if len(regions) != int(count) or any(a < 1 or b < a for a, b in regions):
            raise ValueError(f'Phobius 拓扑不完整：{key}')
        tail = 'cytoplasmic' if topology.endswith('i') else 'non_cytoplasmic' if topology.endswith('o') else 'unknown'
        _insert(result, key, dict(tm_count=int(count), signal_peptide=sp in ('Y', '1'),
                topology=topology, tm_regions=regions, c_terminal_side=tail, raw={'line': line}))
    return result


def parse_deeptmhmm(path):
    lines=[line.strip() for line in Path(path).read_text(encoding='utf-8').splitlines() if line.strip()]
    if len(lines)%3:
        raise ValueError('DeepTMHMM 三行格式不完整')
    result={}
    for i in range(0,len(lines),3):
        header,sequence,topology=lines[i:i+3]
        if not header.startswith('>') or '|' not in header or len(sequence)!=len(topology) or set(topology)-set('IOMSBP'):
            raise ValueError('DeepTMHMM 标题、序列或逐残基拓扑不合法')
        key,prediction=header[1:].split('|',1)
        def regions(label):
            return [[m.start()+1,m.end()] for m in re.finditer(label+'+',topology)]
        tm,beta,signal=regions('M'),regions('B'),regions('S')
        _insert(result,key,dict(prediction=prediction.strip(),topology=topology,output_sequence=sequence,
                tm_count=len(tm),tm_regions=tm,beta_count=len(beta),beta_regions=beta,
                signal_peptide=bool(signal),signal_regions=signal,
                c_terminal_side={'I':'cytoplasmic','O':'non_cytoplasmic','P':'periplasmic'}.get(topology[-1:],'unknown'),
                raw=dict(header=header,sequence=sequence,topology=topology)))
    return result


def parse_predgpi(path):
    data = json.loads(Path(path).read_text(encoding='utf-8'))
    if not isinstance(data, list):
        raise ValueError('PredGPI JSON 应为列表')
    result = {}
    for row in data:
        if 'features' not in row or not isinstance(row['features'], list):
            raise ValueError('PredGPI 缺少 features')
        hits = [x for x in row['features'] if x.get('description') == 'GPI-anchor']
        hit = hits[0] if hits else {}
        omega = int(hit['begin']) if hit else None
        if omega is not None and omega < 1:
            raise ValueError('PredGPI omega 位点非法')
        _insert(result, row.get('accession'), dict(gpi_anchored=bool(hits), omega_site=omega,
                scores={'category_score': hit.get('score')} if hit else {}, raw=row))
    return result


def parse_pts1_svm(path):
    result = {}
    with Path(path).open(encoding='utf-8-sig',newline='') as stream:
        for row in csv.DictReader(stream):
            score = _number(row.get('decision_function'))
            prediction = row.get('prediction')
            tail = row.get('c_terminal_14','')
            if prediction not in ('0','1') or score is None or len(tail)!=14:
                raise ValueError('PTS1 SVM 缺少有效类别、决策值或 C 端序列')
            positive = prediction=='1'
            if positive != (score>0) or row.get('is_pts1') != str(positive):
                raise ValueError('PTS1 SVM 类别与决策值不一致')
            length = int(row['length'])
            if length<14:
                raise ValueError('PTS1 SVM 报告长度小于 14 aa')
            _insert(result,row.get('id'),dict(tool='oKoch PTS1 SVM (reproduced)',
                    pts1_candidate=positive,scores={'decision_function':score},
                    output_tail=tail,output_length=length,raw=row,
                    interpretation='PTS1 候选模型；决策值不是概率，不等同于已确认的过氧化物酶体定位'))
    return result


PARSERS = dict(signalp=parse_signalp, phobius=parse_phobius, deeploc=parse_deeploc,
               targetp=parse_targetp, predgpi=parse_predgpi, deeptmhmm=parse_deeptmhmm,
               deeplocpro=parse_deeplocpro, pts1_svm=parse_pts1_svm)


def scan_signals(record, predictions, rules=None):
    rules = rules or {}
    sequence = record['sequence']
    er, pts = [], []
    if record.get('errors'):
        return {name: dict(status='unsupported', reason='输入序列无效', matches=[], rule_version=RULE_VERSION)
                for name in ('er', 'pts')}
    phobius = predictions.get('phobius', {})
    sp = predictions.get('signalp', {})
    for motif in ('HDEL', 'KDEL'):
        if sequence.endswith(motif):
            supported = (sp.get('status') == 'completed' and sp.get('signal_peptide') is True
                         and phobius.get('status') == 'completed' and phobius.get('tm_count') == 0)
            er.append(dict(kind='ER_luminal_retrieval', sequence=motif, start=len(sequence)-3,
                           end=len(sequence), background_supported=bool(supported)))
    match = re.search(r'(KK..|K.K..)$', sequence)
    if match:
        supported = (phobius.get('status') == 'completed' and phobius.get('tm_count', 0) > 0
                     and phobius.get('c_terminal_side') == 'cytoplasmic')
        er.append(dict(kind='ER_dilysine_candidate', sequence=match[0], start=match.start()+1,
                       end=len(sequence), background_supported=bool(supported)))
    if sequence.endswith('SKL'):
        pts.append(dict(kind='PTS1_classic', sequence='SKL', start=len(sequence)-2, end=len(sequence)))
    elif rules.get('pts1_extended', True) and re.search(r'[STAGCN][RKH][LIVMAFY]$', sequence):
        pts.append(dict(kind='PTS1_extended_candidate', sequence=sequence[-3:], start=len(sequence)-2,
                        end=len(sequence)))
    window = max(9, int(rules.get('pts2_window', 30)))
    for match in re.finditer(r'(?=([RK][LVI].{5}[HQ][LA]))', sequence[:window]):
        pts.append(dict(kind='PTS2_candidate', sequence=match[1], start=match.start()+1, end=match.start()+9))
    result = {}
    for name, matches in (('er', er), ('pts', pts)):
        result[name] = dict(status='completed', matches=matches, rule_version=RULE_VERSION,
                            interpretation='候选信号注释；未命中不能排除该定位',
                            sources=(['https://pmc.ncbi.nlm.nih.gov/articles/PMC2120172/',
                                      'https://pubmed.ncbi.nlm.nih.gov/23481256/'] if name == 'er' else
                                     ['https://prosite.expasy.org/PDOC00299',
                                      'https://pubmed.ncbi.nlm.nih.gov/8125978/']))
    result['pts']['scan_window'] = window
    return result


def annotate(predictions):
    notes = []
    sp, tm, dl, tp, gpi = [predictions.get(k, {}) for k in ('signalp', 'phobius', 'deeploc', 'targetp', 'predgpi')]
    if sp.get('status') == tm.get('status') == 'completed' and sp.get('signal_peptide') != tm.get('signal_peptide'):
        notes.append('SignalP 与 Phobius 的信号肽预测不一致')
    dtm=predictions.get('deeptmhmm',{})
    if dtm.get('status') == tm.get('status') == 'completed':
        if dtm.get('tm_count')!=tm.get('tm_count') or dtm.get('signal_peptide')!=tm.get('signal_peptide'):
            notes.append('Phobius 与 DeepTMHMM 的跨膜/信号肽预测不一致；分别保留原始证据')
    if gpi.get('status') == 'completed' and gpi.get('gpi_anchored'):
        notes.append('GPI 锚定候选；结合分泌信号与拓扑评估表面锚定')
    if tp.get('status') == dl.get('status') == 'completed' and tp.get('prediction') == 'mTP':
        if any('mitochond' in x.lower() for x in dl.get('localizations', [])):
            notes.append('TargetP 与 DeepLoc 提供一致的线粒体定位证据')
    return notes


def summary_row(record):
    predictions = record.get('predictions', {})
    localization = 'deeplocpro' if 'deeplocpro' in predictions else 'deeploc'
    def label(name, key, boolean=False):
        p = predictions.get(name, {})
        if p.get('applicable') is False:
            return '不适用'
        if p.get('status') != 'completed':
            return p.get('status', 'not_run')
        value = p.get(key)
        if boolean:
            return 'Yes' if value else 'No'
        return '; '.join(map(str, value)) if isinstance(value, list) else value
    def motifs(name):
        p = predictions.get(name, {})
        if p.get('applicable') is False:
            return '不适用'
        if p.get('status') != 'completed':
            return p.get('status', 'not_run')
        return '; '.join(f"{x['kind']}:{x['sequence']}@{x['start']}-{x['end']}" for x in p.get('matches', [])) or '未检出所扫描基序'
    row = {'Protein ID': record['original_id'], 'Internal ID': record['internal_id'],
           'Length': record['sequence_length'], 'Status': record.get('status', 'pending'),
           'Localization': label(localization, 'localizations'),
           'Localization model': 'DeepLocPro' if localization=='deeplocpro' else 'DeepLoc 2.1',
           'SP': label('signalp', 'signal_peptide', True),
           'TM': label('phobius', 'tm_count'), 'DeepTMHMM': label('deeptmhmm','prediction'), 'TargetP': label('targetp', 'prediction'),
           'GPI': label('predgpi', 'gpi_anchored', True), 'ER signal': motifs('er'), 'PTS candidate': motifs('pts'),
           'PTS1 SVM': label('pts1_svm','pts1_candidate',True),
           'PTS1 SVM score': label('pts1_svm','scores').get('decision_function') if predictions.get('pts1_svm',{}).get('status')=='completed' else '',
           'Notes': '; '.join(record.get('annotations', [])),
           'Warnings': '; '.join(record.get('warnings', []) + record.get('errors', []))}
    row.update({name + ' status': predictions.get(name, {}).get('status', 'not_run') for name in MODULES})
    row['deeplocpro status'] = predictions.get('deeplocpro',{}).get('status','not_run')
    row['Header'], row['Sequence'] = record['header'], record['sequence']
    return row


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name('.' + path.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def _spreadsheet_safe(value):
    return "'" + value if isinstance(value, str) and value.startswith(('=', '+', '-', '@', '\t', '\r')) else value


def export_results(records, folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    atomic_json(folder / 'results.json', records)
    rows = [summary_row(r) for r in records]
    columns = list(rows[0]) if rows else ['Protein ID', 'Internal ID', 'Status']
    temporary = folder / '.summary.tsv.tmp'
    with temporary.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows({k: _spreadsheet_safe(v) for k, v in row.items()} for row in rows)
    os.replace(temporary, folder / 'summary.tsv')
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font, PatternFill
    except ImportError:
        return {'xlsx': False, 'reason': '缺少 openpyxl；JSON/TSV 已保存'}
    workbook = Workbook(write_only=True)
    sheet = workbook.create_sheet('Proteins')
    from openpyxl.cell import WriteOnlyCell
    header = []
    for key in columns:
        cell = WriteOnlyCell(sheet, value=key)
        cell.fill, cell.font = PatternFill('solid', fgColor='DCEAF5'), Font(bold=True)
        header.append(cell)
    sheet.append(header)
    for row in rows:
        values = []
        for key in columns:
            value = _spreadsheet_safe(row[key])
            if isinstance(value, str) and len(value) > 32767:
                value = value[:32600] + '… [Excel 单元格长度受限；完整内容见 results.json]'
            values.append(value)
        sheet.append(values)
    temporary = folder / '.summary.xlsx.tmp'
    workbook.save(temporary)
    os.replace(temporary, folder / 'summary.xlsx')
    return {'xlsx': True}
