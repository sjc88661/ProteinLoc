"""Real available-model validation. Missing models remain failed, never simulated."""
import argparse
import json
import random
import resource
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from parsers import STANDARD_AA, atomic_json
from pipeline import create_job, load_config, read_json, run_worker


def references():
    chunks, metadata = [], []
    for path in sorted((ROOT/'examples/reference').glob('P*.json')):
        data = read_json(path)
        accession = data['primaryAccession']
        sequence = data['sequence']['value']
        location_comments = [x for x in data.get('comments',[]) if x.get('commentType') == 'SUBCELLULAR LOCATION']
        metadata.append(dict(accession=accession, organism=data['organism'], length=len(sequence),
                             entry_audit=data.get('entryAudit'), location_comments=location_comments,
                             source=f'https://rest.uniprot.org/uniprotkb/{accession}.json'))
        chunks.append(f'>{accession} {data.get("uniProtkbId",accession)}\n{sequence}\n')
    if not chunks:
        raise RuntimeError('Missing downloaded UniProt reference JSON')
    text = ''.join(chunks)
    (ROOT/'examples/yeast_reference.fasta').write_text(text,encoding='utf-8')
    atomic_json(ROOT/'examples/reference/metadata.json',metadata)
    return text


def execute(text,label):
    started = time.monotonic()
    job = create_job(text,'fungi',load_config(),ROOT/'validation/runs',species_name=label)
    run_worker(job.parent)
    state = read_json(job/'job.json')
    result = dict(label=label,job=str(job.relative_to(ROOT)),count=state['count'],status=state['status'],
                  elapsed_seconds=round(time.monotonic()-started,3),modules=state['modules'],
                  python_peak_rss_kb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  child_peak_rss_kb=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                  errors=state['errors'])
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--count',type=int,default=50)
    p.add_argument('--skip-reference',action='store_true',help='复用已完成的参考样本记录，只执行本次批量验收')
    args=p.parse_args()
    results=([r for r in read_json(ROOT/'validation/real_runs.json') if r['count']==3 and r['status']=='completed']
             if args.skip_reference else [execute(references(),'公开酿酒酵母参考蛋白；非独立精度测试')])
    atomic_json(ROOT/'validation/real_runs.json',results)
    rng=random.Random(20261110)
    alphabet='ACDEFGHIKLMNPQRSTVWY'
    chunks=[]
    for i in range(args.count):
        sequence='M'+''.join(rng.choices(alphabet,k=60))+('SKL' if i%10==0 else 'HDEL' if i%10==1 else 'AAA')
        chunks.append(f'>synthetic_stress_{i+1:05d}\n{sequence}\n')
    results.append(execute(''.join(chunks),f'{args.count} 条合成序列批量流程验收；真实模块执行'))
    atomic_json(ROOT/'validation/real_runs.json',results)


if __name__=='__main__':
    main()
