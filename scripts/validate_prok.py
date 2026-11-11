"""Run at most 50 synthetic inputs through the actual local prokaryotic workflow."""
import argparse
import json
import random
import resource
import sys
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from parsers import atomic_json
from pipeline import create_job,load_config,read_json,run_worker


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--count',type=int,default=50)
    args=parser.parse_args()
    if not 1 <= args.count <= 50:
        parser.error('本次验收只允许 1–50 条')
    rng=random.Random(20261110)
    chunks=[]
    for index in range(args.count):
        sequence='M'+''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY',k=60))+('SKL' if index%10==0 else 'HDEL' if index%10==1 else 'AAA')
        chunks.append(f'>synthetic_stress_{index+1:05d}\n{sequence}\n')
    started=time.monotonic()
    job=create_job(''.join(chunks),'gram_positive',load_config(),ROOT/'validation/runs',
                   species_name=f'{args.count} 条原核合成工程样本；真实本地流程')
    print('job='+str(job.relative_to(ROOT)),flush=True)
    run_worker(job.parent)
    state=read_json(job/'job.json')
    result=dict(job=str(job.relative_to(ROOT)),count=state['count'],status=state['status'],
                elapsed_seconds=round(time.monotonic()-started,3),organism=state['organism'],
                modules=state['modules'],core_execution=state.get('core_execution'),errors=state['errors'],
                python_peak_rss_kb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                child_peak_rss_kb=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    atomic_json(ROOT/'validation/deeplocpro_run.json',result)
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
    records=read_json(job/'results.json')
    assert state['status']=='completed' and len(records)==args.count
    assert all(r['status']=='completed' for r in records)
    for name in ('signalp','phobius','deeplocpro','deeptmhmm'):
        assert state['modules'][name]['completed']==args.count and not state['modules'][name]['cached']
    for record in records:
        scores=record['predictions']['deeplocpro']['scores']
        assert scores['Outer Membrane']==scores['Periplasmic']==0
        assert all(record['predictions'][name].get('applicable') is False for name in ('deeploc','targetp','predgpi','er','pts'))
    invocations=list((job/'raw/core/results').glob('**/proteinloc_invocation.json'))
    assert len(invocations)==1 and read_json(invocations[0])['group']=='positive'


if __name__=='__main__':
    main()
