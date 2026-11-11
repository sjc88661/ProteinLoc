"""本地分析及独立后台队列。GUI 和 CLI 共用这些函数；不上传序列。"""
import argparse
import copy
import hashlib
import importlib.metadata
import json
import os
import shutil
import signal
import subprocess
import sys
import time
import traceback
import uuid
from datetime import datetime, timezone
from pathlib import Path

from parsers import (MODULES, PARSERS, RULE_VERSION, STANDARD_AA, annotate, atomic_json,
                     export_results, read_fasta, scan_signals)

ROOT = Path(__file__).resolve().parent
CORE = ('signalp', 'phobius', 'deeploc')
EUKARYOTES = ('fungi', 'plant', 'animal', 'other_euk')
PROKARYOTES = {'gram_negative':'negative','gram_positive':'positive','archaea':'archaea','prokaryote':'any'}
TERMINAL = ('completed', 'partial', 'failed', 'cancelled', 'interrupted')
_DIGESTS = {}


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def local_path(value):
    path = Path(value).expanduser()
    return path if path.is_absolute() else ROOT / path


def _upgrade_config(config):
    """Add new tools and replace only recognized built-in legacy organism literals."""
    for name in ('deeplocpro','pts1_svm'):
        if name not in config['tools']:
            config['tools'][name] = copy.deepcopy(read_json(ROOT/'config.json')['tools'][name])
    commands = [(config.get('metaloc',{}).get('command',[]),('nextflow','nextflow.exe'),
                 ('--organism',),('euk','other'),'{metaloc_organism}'),
                (config['tools'].get('signalp',{}).get('command',[]),('signalp6','signalp6.py'),
                 ('--organism','-org'),('euk','eukarya','other'),'{signalp_organism}')]
    for command,programs,flags,values,template in commands:
        if command and Path(str(command[0])).name in programs:
            for index,part in enumerate(command[:-1]):
                if part in flags and command[index+1] in values:
                    command[index+1] = template
    return config


def load_config(path=None):
    path = Path(path) if path else ROOT / 'config.json'
    config = _upgrade_config(read_json(path))
    if config.get('core_engine') not in ('metaloc', 'commands'):
        raise ValueError('core_engine 必须是 metaloc 或 commands')
    for name in PARSERS:
        tool = config['tools'][name]
        if not isinstance(tool.get('command'), list) or not tool['command']:
            raise ValueError(f'{name} command 必须是参数列表')
        if not tool.get('output_glob'):
            raise ValueError(f'{name} 缺少 output_glob')
    return config


def _context(config, output=None, fasta=None, organism='fungi'):
    native_python = Path(sys.prefix).resolve()/'bin/python'
    context = dict(root=str(ROOT), python=str(native_python) if native_python.is_file() else sys.executable, output=str(output or ROOT),
                   input=str(fasta or ROOT / 'input.faa'), targetp_organism='pl' if organism == 'plant' else 'non-pl',
                   metaloc_organism='euk' if organism in EUKARYOTES else 'other',
                   signalp_organism='eukarya' if organism in EUKARYOTES else 'other',
                   deeplocpro_group=PROKARYOTES.get(organism,'any'))
    meta = config.get('metaloc', {})
    context['metaloc'] = str(local_path(meta.get('path', 'upstream/metaLoc')))
    for key, path in meta.get('assets', {}).items():
        context[key] = str(local_path(path))
    if output:
        context.update(work=str(Path(output).parents[1] / 'nextflow-work'), override=str(Path(output) / 'runtime.config'))
    return context


def _expand(parts, context):
    command=[str(part).format_map(context) for part in parts]
    binary=Path(command[0])
    if binary.is_absolute() and binary.parent.is_dir():
        # Resolve the directory, preserving a venv's python symlink itself.
        command[0]=str(binary.parent.resolve()/binary.name)
    return command


def _command_problem(command):
    binary = command[0]
    found = Path(binary).is_file() if '/' in binary else shutil.which(binary, path=_env()['PATH'])
    if not found:
        return '未找到程序：' + binary
    if '/' in binary and not os.access(binary, os.X_OK):
        return '程序不可执行：' + binary
    for arg in command[1:]:
        if arg.endswith(('.py', '.pl', '.jar')) and '/' in arg and not Path(arg).is_file():
            return '未找到程序文件：' + arg
    return None


def _env(extra=None):
    env = os.environ.copy()
    env['PATH'] = str(ROOT / '.runtime/bin') + os.pathsep + env.get('PATH', '')
    env['NXF_HOME'] = str(ROOT / '.runtime/nextflow')
    env['MAMBA_ROOT_PREFIX'] = str(ROOT / '.runtime/mamba')
    env['TORCH_HOME'] = str(ROOT / '.runtime/model-cache/torch')
    env['HF_HOME'] = str(ROOT / '.runtime/model-cache/huggingface')
    env['OMP_NUM_THREADS'] = env.get('OMP_NUM_THREADS', '2')
    env.update(extra or {})
    return env


def _meta_problems(config, organism='fungi'):
    meta = config['metaloc']
    context = _context(config,organism=organism)
    errors = []
    if not (Path(context['metaloc']) / 'main.nf').is_file():
        errors.append('未找到 metaLoc main.nf')
    else:
        try:
            actual = subprocess.check_output(['git', '-C', context['metaloc'], 'rev-parse', 'HEAD'], text=True, timeout=5).strip()
            if actual != meta.get('commit'):
                errors.append('metaLoc commit 与配置不一致')
        except (OSError, subprocess.SubprocessError):
            errors.append('无法核验 metaLoc commit')
    for key in (('signalp_tar','phobius_tar') if organism in PROKARYOTES else ('signalp_tar','deeploc_tar','phobius_tar')):
        path = Path(context.get(key, ''))
        if not path.is_file() or path.stat().st_size == 0:
            errors.append('缺少授权安装包：' + str(path))
    command = _expand(meta['command'], dict(context, work='work', override='runtime.config'))
    problem = _command_problem(command)
    if problem:
        errors.append(problem)
    if meta.get('use_micromamba') and not shutil.which('micromamba', path=_env()['PATH']):
        errors.append('未找到 micromamba')
    if meta.get('runtime_environment') and not (local_path(meta['runtime_environment'])/'conda-meta/history').is_file():
        errors.append('模型环境尚未安装；运行 bash scripts/python.sh scripts/install_core.py')
    return errors


def check_environment(config=None):
    config = config or load_config()
    report = []
    meta_errors = _meta_problems(config) if config['core_engine'] == 'metaloc' else []
    context = _context(config)
    for name, tool in config['tools'].items():
        use_meta = name in (*CORE,'deeplocpro') and config['core_engine'] == 'metaloc'
        fallback = use_meta and meta_errors and config['metaloc'].get('fallback_local_commands', False)
        errors = list(meta_errors) if use_meta and not fallback else []
        if not use_meta or fallback or name=='deeplocpro':
            problem = _command_problem(_expand(tool['command'], context))
            if problem:
                errors.append(problem)
        for value in tool.get('required_files',[]):
            if not local_path(value).is_file():
                errors.append('缺少本地模型文件：'+value)
        report.append(dict(module=name, ready=not errors, version=tool.get('version'), issues=errors,
                           engine='metaloc' if use_meta and not fallback else 'local_command',
                           notes=list(meta_errors) if fallback else []))
    report.extend(dict(module=name, ready=True, version=RULE_VERSION, issues=[]) for name in ('er', 'pts'))
    return report


def _process_identity(pid):
    try:
        stat = Path(f'/proc/{pid}/stat').read_text()
        fields = stat[stat.rfind(')')+2:].split()
        return None if fields[0] == 'Z' else fields[19]
    except (OSError, IndexError):
        return None


def _alive(info):
    return bool(info and info.get('token') and _process_identity(info.get('pid')) == info['token'])


def _group_members(info):
    """Find owned processes, including descendants that created a new session."""
    if not info:
        return {}
    original, tagged = {}, {}
    verified = _alive(info)
    birth = int(info.get('token') or 0)
    for directory in Path('/proc').iterdir():
        if not directory.name.isdigit():
            continue
        try:
            stat = (directory/'stat').read_text()
            fields = stat[stat.rfind(')')+2:].split()
            if fields[0] == 'Z' or directory.stat().st_uid != os.getuid():
                continue
            pid = int(directory.name)
            if int(fields[2]) == info['pid'] and int(fields[3]) == info['pid']:
                original[pid] = fields[19]
            if info.get('run_tag') and int(fields[19]) >= birth:
                marker = ('PROTEINLOC_RUN_TOKEN='+info['run_tag']).encode()
                if marker in (directory/'environ').read_bytes().split(b'\0'):
                    tagged[pid] = fields[19]
        except (OSError, ValueError, IndexError):
            continue
    if verified or set(original).intersection(tagged):
        tagged.update(original)
    return tagged


def _signal_owned(info, sig):
    for pid, token in _group_members(info).items():
        descriptor = None
        try:
            if hasattr(os,'pidfd_open') and hasattr(signal,'pidfd_send_signal'):
                descriptor = os.pidfd_open(pid)
                if _process_identity(pid) == token:
                    signal.pidfd_send_signal(descriptor,sig)
            elif _process_identity(pid) == token:
                os.kill(pid,sig)
        except ProcessLookupError:
            pass
        finally:
            if descriptor is not None:
                os.close(descriptor)


def _stop_group(info, grace=2):
    if not _group_members(info):
        return
    _signal_owned(info,signal.SIGTERM)
    deadline = time.monotonic()+grace
    while _group_members(info) and time.monotonic()<deadline:
        time.sleep(.08)
    if _group_members(info):
        _signal_owned(info,signal.SIGKILL)
        deadline = time.monotonic()+3
        while _group_members(info) and time.monotonic()<deadline:
            time.sleep(.08)
        if _group_members(info):
            raise RuntimeError('工具进程组尚未停止，保留任务记录以便恢复')


def _save(job, state):
    state['updated_at'] = now()
    atomic_json(Path(job) / 'job.json', state)
    if state['status']=='running' and state.get('worker',{}).get('pid')==os.getpid():
        atomic_json(ROOT/'.runtime/active.json',dict(job=str(Path(job).resolve()),
                    worker=state['worker'],process=state.get('process')))


def create_job(text, organism, config=None, runs_dir=None, species_name=''):
    config = _upgrade_config(copy.deepcopy(config or load_config()))
    records = read_fasta(text)
    runs_dir = Path(runs_dir) if runs_dir else local_path(config.get('runs_dir', 'runs'))
    runs_dir.mkdir(parents=True, exist_ok=True)
    job_id = datetime.now().strftime('%Y%m%d-%H%M%S') + '-' + uuid.uuid4().hex[:8]
    job = runs_dir / job_id
    job.mkdir()
    (job / 'original.fasta').write_text(text, encoding='utf-8')
    atomic_json(job / 'records.json', records)
    state = dict(id=job_id, status='queued', created_at=now(), updated_at=now(), organism=organism,
                 species_name=species_name.strip(), config=config, count=len(records), attempt=0,
                 input_sha256=file_hash(job / 'original.fasta'), records_sha256=file_hash(job / 'records.json'),
                 modules={}, message='等待分析', process=None, errors=[])
    _save(job, state)
    return job


def _lock(runs_dir):
    import fcntl
    runs_dir = Path(runs_dir)
    runs_dir.mkdir(parents=True, exist_ok=True)
    lock = (runs_dir / '.worker.lock').open('a+')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        return None
    return lock


def _execution_lock():
    import fcntl
    folder = ROOT / '.runtime'
    folder.mkdir(parents=True, exist_ok=True)
    lock = (folder/'execution.lock').open('a+')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        return None
    active = folder/'active.json'
    try:
        if active.exists():
            lease = read_json(active)
            previous_job = Path(lease['job'])
            previous = read_json(previous_job/'job.json') if (previous_job/'job.json').exists() else None
            if previous and previous['status']=='running':
                if _alive(previous.get('worker')):
                    raise RuntimeError('检测到仍活动的旧任务记录，暂不启动新的模型计算')
                _recover(previous_job,previous)
            else:
                _stop_group(lease.get('process'))
            active.unlink(missing_ok=True)
    except Exception:
        lock.close()
        raise
    return lock


def list_jobs(runs_dir=None):
    root = Path(runs_dir) if runs_dir else local_path(load_config().get('runs_dir', 'runs'))
    jobs = []
    for path in root.glob('*/job.json'):
        try:
            state = read_json(path)
            state['path'] = str(path.parent)
            jobs.append(state)
        except (OSError, ValueError):
            continue
    return sorted(jobs, key=lambda j: j['id'], reverse=True)


def ensure_worker(runs_dir):
    lock = _lock(runs_dir)
    if lock is None:
        return
    lock.close()
    if not any(j['status'] == 'queued' for j in list_jobs(runs_dir)):
        return
    log = (Path(runs_dir) / 'worker.log').open('ab')
    try:
        subprocess.Popen([_context(load_config())['python'], str(ROOT / 'pipeline.py'), 'worker', '--runs-dir', str(Path(runs_dir).resolve())],
                         cwd=ROOT, stdin=subprocess.DEVNULL, stdout=log, stderr=log,
                         start_new_session=True, close_fds=True, env=_env())
    finally:
        log.close()


def cancel_job(job):
    job = Path(job)
    (job / 'cancel.request').touch()
    state = read_json(job / 'job.json')
    if state['status'] == 'queued':
        state.update(status='cancelled', message='排队任务已取消')
        _save(job, state)


def resume_job(job, config=None):
    job = Path(job)
    state = read_json(job / 'job.json')
    if state['status'] == 'running' and _alive(state.get('worker')):
        raise ValueError('任务仍在运行')
    _stop_group(state.get('process'))
    _recover(job, state)
    if file_hash(job / 'original.fasta') != state['input_sha256'] or file_hash(job / 'records.json') != state['records_sha256']:
        raise ValueError('原始输入或 ID 映射已改变，请建立新任务')
    if config is not None:
        state['config'] = copy.deepcopy(config)
    state['config'] = _upgrade_config(state['config'])
    (job / 'cancel.request').unlink(missing_ok=True)
    state.update(status='queued', message='等待续跑；仅复用输入/版本/参数一致且完整的结果', errors=[], process=None)
    _save(job, state)


def _recover(job, state):
    if state['status'] == 'running' and not _alive(state.get('worker')):
        _stop_group(state.get('process'))
        state.update(status='interrupted', process=None, message='上次计算被中断，可手动续跑')
        _save(job, state)


def wait_job(job):
    job = Path(job)
    while True:
        state = read_json(job/'job.json')
        _recover(job, state)
        if state['status'] in TERMINAL:
            return state
        ensure_worker(job.parent)
        time.sleep(.3)


class Cancelled(Exception):
    def __init__(self, results=None):
        super().__init__('任务已取消')
        self.results = results or {}


def _run_process(job, state, command, cwd, env=None, stdout_file=None):
    problem = _command_problem(command)
    if problem:
        raise RuntimeError(problem)
    if (job / 'cancel.request').exists():
        raise Cancelled()
    with (job / 'run.log').open('ab') as log:
        log.write(('\n[' + now() + '] ' + json.dumps(command, ensure_ascii=False) + '\n').encode())
        log.flush()
        output = Path(stdout_file).open('wb') if stdout_file else log
        run_tag = uuid.uuid4().hex
        child_env = _env(env)
        child_env['PROTEINLOC_RUN_TOKEN'] = run_tag
        process = subprocess.Popen(command, cwd=cwd, env=child_env, stdin=subprocess.DEVNULL,
                                   stdout=output, stderr=log, start_new_session=True)
        state['process'] = dict(pid=process.pid, token=_process_identity(process.pid), run_tag=run_tag, command=command)
        _save(job, state)
        started = time.monotonic()
        try:
            while process.poll() is None:
                if (job / 'cancel.request').exists():
                    raise Cancelled()
                if time.monotonic() - started > float(state['config'].get('timeout_seconds', 86400)):
                    raise TimeoutError('工具运行超过配置时限')
                time.sleep(.15)
            return process.returncode
        finally:
            _stop_group(state['process'])
            process.wait(timeout=4)
            if stdout_file:
                output.close()
            state['process'] = None
            _save(job, state)


def _input_problem(record, tool):
    sequence = record['sequence']
    if record['errors']:
        return '; '.join(record['errors'])
    if len(sequence) < tool.get('min_length', 1):
        return f"序列长度小于本适配器验证下限 {tool['min_length']} aa"
    if tool.get('max_length') and len(sequence) > tool['max_length']:
        return f"序列长 {len(sequence)} aa，超过该模块设置的 {tool['max_length']} aa 上限"
    if tool.get('standard_only') and set(sequence) - STANDARD_AA:
        return '本适配器仅对标准 20 种氨基酸验证，特殊残基不作阴性处理'
    if tool.get('supported_alphabet'):
        analyzed = sequence[-int(tool.get('analysis_cterm', len(sequence))):]
        if set(analyzed) - set(tool['supported_alphabet']):
            return '模型实际分析区段含字母表以外的特殊残基，不作阴性处理'
    return ''


def _prepare(records, name, folder, config):
    folder.mkdir(parents=True, exist_ok=True)
    tool = config['tools'].get(name, {})
    eligible, unsupported, inputs = [], {}, {}
    for record in records:
        key, sequence = record['internal_id'], record['sequence']
        reason = _input_problem(record, tool)
        if reason:
            unsupported[key] = dict(status='unsupported', reason=reason)
            continue
        submitted = sequence if tool.get('preserve_residues') else ''.join(c if c in STANDARD_AA or c == 'X' else 'X' for c in sequence)
        changes = [[i, a, b] for i, (a, b) in enumerate(zip(sequence, submitted), 1) if a != b]
        inputs[key] = dict(submitted_sequence=submitted, submitted_length=len(submitted),
                           submitted_ranges=[[1, len(sequence)]], changed_residues=changes)
        eligible.append((key, submitted))
    fasta = folder / 'input.faa'
    content = ''.join(f'>{key}\n{seq}\n' for key, seq in eligible).encode('utf-8')
    if not fasta.is_file() or fasta.read_bytes() != content:
        fasta.write_bytes(content)
    atomic_json(folder / 'input_map.json', inputs)
    return fasta, inputs, unsupported


def _digest(path):
    stat = path.stat()
    identity = (stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    cached = _DIGESTS.get(str(path))
    if cached is None or cached[0] != identity:
        cached = (identity, file_hash(path))
        _DIGESTS[str(path)] = cached
    return cached[1]


def _provenance(config, name, organism='fungi'):
    settings = config['metaloc'] if name == 'core' else config['tools'][name]
    context = _context(config,organism=organism)
    paths = list(settings.get('provenance_paths', []))
    if name=='core' and organism in PROKARYOTES:
        paths.extend(path for path in config['tools']['deeplocpro'].get('provenance_paths',[]) if path not in paths)
    known_predgpi = name == 'predgpi' and any('/scripts/predgpi_compat.py' in str(x) for x in settings.get('command', []))
    if known_predgpi:
        paths = paths or settings.get('watch_paths', [])
    if name == 'core' and not paths:
        conda = ROOT/'.runtime/conda-envs'
        if conda.exists() and any(conda.glob('*/conda-meta')):
            paths = [str(conda)]
            if (ROOT/'.runtime/model-cache').exists():
                paths.append(str(ROOT/'.runtime/model-cache'))
    verified = bool(paths)
    data = dict(python=sys.version, files={})
    for value in paths:
        path = local_path(str(value).format_map(context)).resolve()
        if not path.exists():
            verified = False
            continue
        files = sorted(path.rglob('*')) if path.is_dir() else [path]
        for item in files:
            if item.is_file() and '__pycache__' not in item.parts and item.suffix not in ('.pyc','.lock'):
                data['files'][str(item)] = _digest(item)
    version_command = settings.get('version_command')
    if version_command:
        try:
            env = {k: str(v).format_map(context) for k,v in settings.get('env',{}).items()}
            probe = subprocess.run(_expand(version_command,context),capture_output=True,text=True,timeout=10,env=_env(env))
            data['actual_version'] = (probe.stdout+probe.stderr).strip()
            verified = verified and probe.returncode in settings.get('version_success_codes',[0]) and bool(data['actual_version'])
        except (OSError,subprocess.SubprocessError,ValueError,KeyError) as error:
            data['version_error'] = str(error)
            verified = False
    elif known_predgpi:
        try:
            data['numpy'] = importlib.metadata.version('numpy')
        except importlib.metadata.PackageNotFoundError:
            verified = False
    elif name != 'core':
        verified = False
    data['verified'] = verified
    return data


def _signature(config, name, state):
    data = dict(input=state['records_sha256'], organism=state['organism'],
                parser=file_hash(ROOT / 'parsers.py'), runner=file_hash(ROOT / 'pipeline.py'))
    data['settings'] = config['metaloc'] if name == 'core' else config['tools'][name]
    if name == 'core':
        names = ('signalp','phobius','deeplocpro') if state['organism'] in PROKARYOTES else CORE
        data['tools'] = {k: config['tools'][k] for k in names}
    context = _context(config,organism=state['organism'])
    candidates = []
    parts = data['settings'].get('command', [])
    for index, part in enumerate(parts):
        try:
            value = str(part).format_map(dict(context, work='work', override='runtime.config'))
            if index == 0:
                path = Path(value) if '/' in value else Path(shutil.which(value, path=_env()['PATH']) or value)
            elif '/' in value:
                path = Path(value)
            else:
                continue
            if path.is_file():
                candidates.append(path)
        except (KeyError, ValueError):
            continue
    for watched in data['settings'].get('watch_paths', []):
        path = local_path(watched).resolve()
        if path.is_dir():
            candidates.extend(p for p in sorted(path.rglob('*')) if p.is_file() and '__pycache__' not in p.parts)
        else:
            candidates.append(path)
    if name == 'core':
        candidates.extend(local_path(x) for x in config['metaloc'].get('assets', {}).values())
        candidates.extend(Path(context['metaloc']) / x for x in ('main.nf', 'nextflow.config'))
    for path in candidates:
        if path.is_file():
            data.setdefault('files', {})[str(path)] = _digest(path)
    provenance = _provenance(config,name,state['organism'])
    data['provenance'] = provenance
    return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest(), provenance


def _cached(job, name, fingerprint):
    try:
        data = read_json(job / 'cache' / (name + '.json'))
        if data['fingerprint'] == fingerprint and data['complete']:
            if all(Path(p).is_file() and file_hash(p) == digest for p, digest in data['files'].items()):
                return data['results']
    except (OSError, ValueError, KeyError):
        pass
    return None


def _cache(job, name, fingerprint, results, files, complete):
    atomic_json(job / 'cache' / (name + '.json'), dict(fingerprint=fingerprint, results=results,
                files={str(path): file_hash(path) for path in files}, complete=complete))


def _find_output(folder, pattern):
    files = sorted(p for p in folder.glob(pattern) if p.is_file())
    if len(files) != 1:
        raise ValueError(f'结果文件应恰有一个，找到 {len(files)} 个：{pattern}')
    return files[0]


def _align(name, path, records, inputs, unsupported, config, job, missing_status='failed'):
    parsed = PARSERS[name](path)
    unknown = set(parsed) - set(inputs)
    if unknown:
        raise ValueError('输出含未知 ID：' + ', '.join(sorted(unknown)[:5]))
    result = copy.deepcopy(unsupported)
    for record in records:
        key = record['internal_id']
        if key not in inputs:
            continue
        if key not in parsed:
            result[key] = dict(status=missing_status, reason='工具输出缺少该蛋白 ID')
            continue
        item = parsed[key]
        if item.get('output_sequence') is not None and item['output_sequence']!=inputs[key]['submitted_sequence']:
            raise ValueError('工具输出序列与该内部 ID 的提交序列不一致：'+key)
        if item.get('output_tail') is not None and item['output_tail']!=inputs[key]['submitted_sequence'][-len(item['output_tail']):]:
            raise ValueError('工具输出 C 端与该内部 ID 的提交序列不一致：'+key)
        if item.get('output_length') is not None and item['output_length']!=inputs[key]['submitted_length']:
            raise ValueError('工具报告长度与该内部 ID 的提交序列不一致：'+key)
        length = record['sequence_length']
        if any(a < 1 or b > length or b < a for a, b in item.get('tm_regions', [])):
            raise ValueError('跨膜坐标超出原始序列：' + key)
        for coordinate in ('cleavage_after', 'omega_site'):
            if item.get(coordinate) is not None and not (1 <= item[coordinate] <= length):
                raise ValueError('切割/加工坐标超出序列：' + key)
        max_length = config['tools'][name].get('max_length')
        if max_length and length > max_length:
            result[key] = dict(status='unsupported', reason='模型可能内部截断；本地版本的区段映射未验证，原始输出已保留')
            continue
        item.update(version=config['tools'][name].get('version', 'unknown'),
                    raw_file=str(path.relative_to(job)), input_map=str(Path(inputs[key]['map_file']).relative_to(job)),
                    submitted_length=inputs[key]['submitted_length'], changed_residues=inputs[key]['changed_residues'])
        if config['tools'][name].get('analysis_nterm'):
            window=int(config['tools'][name]['analysis_nterm'])
            item.update(analyzed_ranges=[[1,min(length,window)]],
                        analysis_note=f'本地工具内部分析 N 端前 {window} aa；提交序列和原序列保留全长')
        if config['tools'][name].get('analysis_cterm'):
            window=int(config['tools'][name]['analysis_cterm'])
            item.update(analyzed_ranges=[[length-window+1,length]],
                        analysis_note=f'仅分析真实 C 端 {window} aa；原序列和提交序列保留全长')
        result[key] = item
    return result


def _failed(records, reason, unsupported=None, status='failed'):
    return {r['internal_id']: copy.deepcopy((unsupported or {}).get(r['internal_id'], dict(status=status, reason=str(reason)))) for r in records}


def _direct(job, state, records, name):
    config = state['config']
    fingerprint, provenance = _signature(config, name, state)
    cached = _cached(job, name, fingerprint) if provenance['verified'] else None
    if cached is not None:
        return cached, True
    output = job / 'raw' / name / f'attempt-{state["attempt"]:03d}'
    fasta, inputs, unsupported = _prepare(records, name, output, config)
    for item in inputs.values():
        item['map_file'] = str(output / 'input_map.json')
    if not inputs:
        return unsupported, False
    tool = config['tools'][name]
    context = _context(config, output, fasta, state['organism'])
    command = _expand(tool['command'], context)
    env = {k: str(v).format_map(context) for k, v in tool.get('env', {}).items()}
    try:
        code = _run_process(job, state, command, output, env,
                            output / tool['stdout_file'] if tool.get('stdout_file') else None)
        if code:
            raise RuntimeError(f'{name} 退出码 {code}，详见 run.log')
        path = _find_output(output, tool['output_glob'])
        results = _align(name, path, records, inputs, unsupported, config, job)
        fingerprint, provenance = _signature(config,name,state)
        atomic_json(output/'provenance.json',provenance)
        for item in results.values():
            if item['status'] == 'completed':
                item['engine'] = 'local_command'
                item['provenance'] = dict(verified=provenance['verified'],fingerprint=fingerprint,
                    manifest=str((output/'provenance.json').relative_to(job)))
        _cache(job, name, fingerprint, results, [path, output / 'input_map.json',output/'provenance.json'],
               provenance['verified'] and all(x['status'] in ('completed', 'unsupported') for x in results.values()))
        return results, False
    except Cancelled:
        raise
    except Exception as error:
        return _failed(records, error, unsupported), False


def _local_core(job, state, records, names, reason):
    state['core_execution'] = dict(engine='local_commands', reason=reason)
    results = {}
    try:
        for name in names:
            state['message'] = '运行本地 ' + name
            _save(job, state)
            results[name], _ = _direct(job, state, records, name)
    except Cancelled as error:
        error.results.update(results)
        raise
    return results, False


def _metaloc(job, state, records):
    config = state['config']
    names = ('signalp','phobius','deeplocpro') if state['organism'] in PROKARYOTES else CORE
    if any(not record['errors'] and _input_problem(record, config['tools'][name])
           for record in records for name in names):
        return _local_core(job, state, records, names, '核心模块输入范围不同，按各模型要求分别运行')
    fingerprint, provenance = _signature(config, 'core', state)
    cached = _cached(job, 'core', fingerprint) if provenance['verified'] else None
    if cached is not None:
        return cached, True
    output = job / 'raw/core'
    fasta, inputs, unsupported = _prepare(records, 'core', output, config)
    for item in inputs.values():
        item['map_file'] = str(output / 'input_map.json')
    if not inputs:
        return {name: copy.deepcopy(unsupported) for name in names}, False
    errors = _meta_problems(config,state['organism'])
    if errors:
        if config['metaloc'].get('fallback_local_commands', False):
            return _local_core(job, state, records, names, '; '.join(errors))
        return {name: _failed(records, '; '.join(errors), unsupported) for name in names}, False
    state['core_execution'] = dict(engine='metaloc',localization_tool='deeplocpro' if state['organism'] in PROKARYOTES else 'deeploc',
                                   organism_group=PROKARYOTES.get(state['organism']))
    for name in ('bin','envs'):
        target=(local_path(config['metaloc']['path'])/name).resolve()
        link=output/name
        if link.is_symlink() and link.resolve()!=target:
            link.unlink()
        if not link.exists():
            link.symlink_to(target,target_is_directory=True)
    for name in ('results', 'trace.tsv'):
        published = output / name
        if published.exists():
            previous = output / f'previous-{published.stem}-{uuid.uuid4().hex[:8]}{published.suffix}'
            published.rename(previous)
    context = _context(config, output, fasta, state['organism'])
    override = output / 'runtime.config'
    cpus, memory = max(1, int(config.get('cpus', 2))), max(1, int(config.get('memory_gb', 4)))
    override.write_text(f"process.cpus = {cpus}\nprocess.memory = '{memory} GB'\nexecutor.queueSize = 1\nconda.useMicromamba = {str(config['metaloc'].get('use_micromamba', True)).lower()}\nconda.cacheDir = '{ROOT / '.runtime/conda-envs'}'\n")
    if config['metaloc'].get('runtime_environment'):
        environment=local_path(config['metaloc']['runtime_environment']).resolve()
        with override.open('a') as stream:
            stream.write("process { withName: /.*/ { conda = '"+str(environment)+"' } }\n")
    results, files = {}, [output / 'input_map.json']
    cancelled, code = False, None
    try:
        command = _expand(config['metaloc']['command'], context)
        try:
            previous = read_json(job/'cache/core.json')
            same_environment = provenance['verified'] and previous['fingerprint'] == fingerprint
        except (OSError,ValueError,KeyError):
            same_environment = False
        if not same_environment:
            command = [arg for arg in command if arg != '-resume']
        env = {k: str(v).format_map(context) for k,v in config['metaloc'].get('env',{}).items()}
        if state['organism'] in PROKARYOTES:
            env['PROTEINLOC_DEEPLOCPRO_GROUP'] = context['deeplocpro_group']
        code = _run_process(job, state, command, output, env)
    except Cancelled:
        cancelled = True
    except Exception as error:
        return {name: _failed(records, error, unsupported) for name in names}, False
    patterns = dict(signalp='results/**/signalp_results.txt', phobius='results/**/phobius_results.txt', deeploc='results/**/deeploc_results.txt')
    patterns['deeplocpro']='results/**/deeplocpro_results.txt'
    for name in names:
        try:
            path = _find_output(output, patterns[name])
            results[name] = _align(name, path, records, inputs, unsupported, config, job,
                                   missing_status='cancelled' if cancelled else 'failed')
            files.append(path)
        except Exception as error:
            results[name] = _failed(records, error, unsupported, status='cancelled' if cancelled else 'failed')
    if code:
        state['errors'].append(f'metaLoc 退出码 {code}，已保留能够解析的独立结果')
    fingerprint, provenance = _signature(config,'core',state)
    atomic_json(output/'provenance.json',provenance)
    files.append(output/'provenance.json')
    for group in results.values():
        for item in group.values():
            if item['status']=='completed':
                item['engine']='metaloc'
                item['provenance']=dict(verified=provenance['verified'],fingerprint=fingerprint,
                    manifest=str((output/'provenance.json').relative_to(job)))
    complete = provenance['verified'] and not cancelled and not code and all(v['status'] in ('completed', 'unsupported') for group in results.values() for v in group.values())
    _cache(job, 'core', fingerprint, results, files, complete)
    if cancelled:
        raise Cancelled(results)
    return results, False


def _module_summary(results):
    counts = {status: sum(v['status'] == status for v in results.values()) for status in ('completed','failed','unsupported','cancelled')}
    if counts['completed'] == len(results):
        status = 'completed'
    elif counts['completed']:
        status = 'partial'
    elif counts['failed']:
        status = 'failed'
    elif counts['cancelled']:
        status = 'cancelled'
    else:
        status = 'unsupported'
    return dict(status=status, **counts)


def _execute_job(job, state):
    state['config'] = _upgrade_config(state['config'])
    config, organism = state['config'], state['organism']
    prok = organism in PROKARYOTES
    names = (*MODULES,'deeplocpro') if prok else MODULES
    core_names = ('signalp','phobius','deeplocpro') if prok else CORE
    records = read_json(job / 'records.json')
    predictions = {r['internal_id']: {} for r in records}
    state.update(status='running', attempt=state.get('attempt', 0)+1, started_at=now(), errors=[], modules={},
                 worker=dict(pid=os.getpid(), token=_process_identity(os.getpid())), process=None)
    _save(job, state)
    cancelled = False
    def store(name, result, cached=False):
        for key, item in result.items():
            predictions[key][name] = item
        state['modules'][name] = dict(_module_summary(result), cached=cached)
        _save(job, state)
    try:
        if organism not in EUKARYOTES and not prok:
            for name in names:
                store(name, {r['internal_id']: dict(status='unsupported', reason='未知类群或尚未适配的物种路由') for r in records})
        else:
            if prok:
                for name in ('deeploc','targetp','predgpi','er','pts','pts1_svm'):
                    store(name,{r['internal_id']:dict(status='unsupported',applicable=False,reason='此真核定位/信号模块不适用于原核蛋白') for r in records})
            if config['core_engine'] == 'metaloc':
                state['message'] = '运行 metaLoc：SignalP / Phobius / '+('DeepLocPro' if prok else 'DeepLoc')
                _save(job, state)
                core, cached = _metaloc(job, state, records)
                for name in core_names:
                    store(name, core[name], cached)
            else:
                for name in core_names:
                    state['message'] = '运行 ' + name
                    _save(job, state)
                    result, cached = _direct(job, state, records, name)
                    store(name, result, cached)
            for name in (('deeptmhmm',) if prok else ('targetp','predgpi','deeptmhmm','pts1_svm')):
                state['message'] = '运行 ' + name
                _save(job, state)
                result, cached = _direct(job, state, records, name)
                store(name, result, cached)
    except Cancelled as error:
        for name, result in error.results.items():
            store(name, result)
        cancelled = True
    for record in records:
        key = record['internal_id']
        if organism in EUKARYOTES:
            predictions[key].update(scan_signals(record, predictions[key], config.get('rules')))
        for name in names:
            predictions[key].setdefault(name, dict(status='cancelled' if cancelled else 'failed', reason='任务取消' if cancelled else '模块未返回结果'))
        record['predictions'] = predictions[key]
        record['annotations'] = annotate(predictions[key])
        statuses = [p['status'] for p in predictions[key].values() if p.get('applicable') is not False]
        record['status'] = 'invalid' if record['errors'] else 'completed' if all(s == 'completed' for s in statuses) else 'partial'
    for name in names:
        state['modules'][name] = dict(_module_summary({r['internal_id']: predictions[r['internal_id']][name] for r in records}),
                                     cached=state['modules'].get(name, {}).get('cached', False))
    state['message'] = '保存结果与导出表格'
    _save(job, state)
    state['export'] = export_results(records, job)
    complete = all(r['status'] == 'completed' for r in records) and not state['errors']
    some = any(p['status'] == 'completed' for values in predictions.values() for p in values.values())
    state.update(status='cancelled' if cancelled else 'completed' if complete else 'partial' if some else 'failed',
                 message='已取消，保留已有结果' if cancelled else '全部分析完成' if complete else '分析不完整，请查看模块状态和日志',
                 completed_proteins=sum(r['status'] == 'completed' for r in records), finished_at=now(), process=None)
    _save(job, state)


def run_worker(runs_dir):
    runs_dir = Path(runs_dir).resolve()
    lock = _lock(runs_dir)
    if lock is None:
        return
    try:
        for state in list_jobs(runs_dir):
            _recover(Path(state['path']), state)
        while True:
            queued = sorted((s for s in list_jobs(runs_dir) if s['status'] == 'queued'), key=lambda s: s['id'])
            if not queued:
                return
            state = queued[0]
            job = Path(state.pop('path'))
            execution_lock = None
            try:
                while execution_lock is None and not (job/'cancel.request').exists():
                    execution_lock = _execution_lock()
                    if execution_lock is None:
                        state['message'] = '等待其他任务目录中的分析完成'
                        _save(job,state)
                        time.sleep(.2)
                if (job / 'cancel.request').exists():
                    state.update(status='cancelled', message='任务已取消')
                    _save(job, state)
                    continue
                if file_hash(job / 'original.fasta') != state['input_sha256'] or file_hash(job / 'records.json') != state['records_sha256']:
                    raise ValueError('输入或 ID 映射完整性检查失败')
                _execute_job(job, state)
            except Exception as error:
                state.update(status='failed', process=None, message=str(error))
                state.setdefault('errors', []).append(str(error))
                _save(job, state)
                with (job / 'run.log').open('a') as log:
                    traceback.print_exc(file=log)
            finally:
                if execution_lock is not None:
                    if state['status'] in TERMINAL:
                        (ROOT/'.runtime/active.json').unlink(missing_ok=True)
                    execution_lock.close()
    finally:
        lock.close()


def main():
    parser = argparse.ArgumentParser(description='ProteinLoc 本地亚细胞定位分析')
    sub = parser.add_subparsers(dest='action', required=True)
    doctor = sub.add_parser('doctor', help='检查模型和环境')
    doctor.add_argument('--config')
    worker = sub.add_parser('worker')
    worker.add_argument('--runs-dir', required=True)
    submit = sub.add_parser('submit')
    submit.add_argument('fasta')
    submit.add_argument('--organism', default='fungi', choices=(*EUKARYOTES, *PROKARYOTES, 'unknown'))
    submit.add_argument('--config')
    submit.add_argument('--runs-dir')
    submit.add_argument('--wait', action='store_true')
    for name in ('resume', 'cancel', 'status'):
        cmd = sub.add_parser(name)
        cmd.add_argument('job')
        if name == 'resume':
            cmd.add_argument('--config')
    args = parser.parse_args()
    if args.action == 'doctor':
        report = check_environment(load_config(args.config))
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return 0 if all(x['ready'] for x in report) else 2
    if args.action == 'worker':
        run_worker(args.runs_dir)
    elif args.action == 'submit':
        config = load_config(args.config)
        job = create_job(Path(args.fasta).read_text(encoding='utf-8-sig'), args.organism, config, args.runs_dir)
        print(job, flush=True)
        ensure_worker(job.parent)
        if args.wait:
            state = wait_job(job)
            print(json.dumps(state, ensure_ascii=False, indent=2))
            return 0 if state['status'] == 'completed' else 2
    elif args.action == 'resume':
        resume_job(args.job, load_config(args.config))
        ensure_worker(Path(args.job).parent)
    elif args.action == 'cancel':
        cancel_job(args.job)
    elif args.action == 'status':
        state = read_json(Path(args.job) / 'job.json')
        _recover(Path(args.job), state)
        print(json.dumps(state, ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
