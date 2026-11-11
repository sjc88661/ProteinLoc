"""Synthetic external-tool outputs for integration tests only. NOT a predictor."""
import argparse
import csv
import json
import sys
import subprocess
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from parsers import read_fasta

p = argparse.ArgumentParser()
p.add_argument('--tool', required=True)
p.add_argument('--input', required=True)
p.add_argument('--output', required=True)
p.add_argument('--delay', type=float, default=0)
p.add_argument('--fail', action='store_true')
p.add_argument('--extra', action='store_true')
p.add_argument('--omit-last', action='store_true')
p.add_argument('--spawn-stubborn', action='store_true')
p.add_argument('--escape-session', action='store_true')
p.add_argument('--model-file')
args = p.parse_args()
time.sleep(args.delay)
if args.fail:
    sys.exit(7)
rows = read_fasta(Path(args.input).read_text())
ids = [r['original_id'] for r in rows][::-1]
if args.omit_last:
    ids = ids[:-1]
if args.extra:
    ids.append('UNEXPECTED_ID')
out = Path(args.output)
out.mkdir(parents=True, exist_ok=True)
if args.spawn_stubborn:
    ready=out/'stubborn.ready'
    code=f"import os,pathlib,signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); pathlib.Path({str(ready)!r}).write_text(str(os.getpid())); time.sleep(60)"
    subprocess.Popen([sys.executable,'-c',code],start_new_session=args.escape_session)
    time.sleep(60)
if args.tool == 'signalp':
    (out / 'prediction_results.txt').write_text('# SignalP-6.0\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n' + ''.join(f'{i}\tOTHER\t0.99\t0.01\t\n' for i in ids))
elif args.tool == 'phobius':
    (out / 'phobius_results.txt').write_text('SEQENCE ID TM SP PREDICTION\n' + ''.join(f'{i} 0 0 i\n' for i in ids))
elif args.tool == 'deeploc':
    label=Path(args.model_file).read_text().strip() if args.model_file else 'Cytoplasm'
    (out / 'results_fixture.csv').write_text('Protein_ID,Localizations,Cytoplasm\n' + ''.join(f'{i},{label},0.9\n' for i in ids))
elif args.tool == 'deeplocpro':
    (out/'results_fixture.csv').write_text(',ACC,Localization,Cell wall & surface,Extracellular,Cytoplasmic,Cytoplasmic Membrane,Outer Membrane,Periplasmic\n'+''.join(f'{n},{i},Periplasmic,0.01,0.02,0.03,0.04,0.10,0.80\n' for n,i in enumerate(ids)))
elif args.tool == 'targetp':
    (out / 'targetp_summary.targetp2').write_text('# TargetP-2.0\n# ID\tPrediction\tOTHER\tSP\tmTP\tCS Position\n' + ''.join(f'{i}\tOTHER\t0.9\t0.05\t0.05\t\n' for i in ids))
elif args.tool == 'predgpi':
    (out / 'predgpi.json').write_text(json.dumps([{'accession': i, 'features': []} for i in ids]))
elif args.tool == 'pts1_svm':
    sequences={r['original_id']:r['sequence'] for r in rows}
    (out/'predictions.csv').write_text('id,prediction,is_pts1,decision_function,length,c_terminal_14\n'+''.join(f'{i},0,False,-0.12345678912345678,{len(sequences[i])},{sequences[i][-14:]}\n' for i in ids))
elif args.tool == 'deeptmhmm':
    sequences={r['original_id']:r['sequence'] for r in rows}
    (out/'prediction').mkdir(exist_ok=True)
    (out/'prediction/predicted_topologies.3line').write_text(''.join(f'>{i} | GLOB\n{sequences.get(i,"A"*50)}\n'+('I'*len(sequences.get(i,'A'*50)))+'\n' for i in ids))
