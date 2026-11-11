"""Independent numerical and internal-CV audit of the supplied PTS1 NPZ weights."""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import warnings
from pathlib import Path

os.environ.setdefault('OMP_NUM_THREADS','2')
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
ROOT=Path(__file__).resolve().parents[1]
SOURCE=(ROOT/'.runtime/pts1-upstream').resolve()
sys.path[:0]=[str(ROOT),str(SOURCE),str(ROOT/'tools/pts1')]
import numpy as np
import sklearn
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedKFold,StratifiedGroupKFold
from sklearn.svm import SVC
from threadpoolctl import threadpool_limits
from Bio import SeqIO
from pts1_prediction_tool.pts1_prediction import PTS1_Predictor
from predict_pts1_numpy import NumpyPTS1Predictor
from parsers import atomic_json


def metrics(truth,predictions):
    tn,fp,fn,tp=map(int,confusion_matrix(truth,predictions,labels=[0,1]).ravel())
    return dict(tp=tp,tn=tn,fp=fp,fn=fn,sensitivity=tp/(tp+fn),specificity=tn/(tn+fp),
                precision=tp/(tp+fp) if tp+fp else 0,accuracy=(tp+tn)/(tp+tn+fp+fn))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--baseline',action='store_true',help='Inspect the supplied predictor before compatibility fixes; no CV')
    args=parser.parse_args()
    started=time.monotonic()
    warnings.filterwarnings('ignore',message='X has feature names')
    with threadpool_limits(limits=2):
        reference=PTS1_Predictor()
        model=reference.trained_model
        portable=NumpyPTS1Predictor(ROOT/'tools/pts1/svm_weights.npz')
        positives=list(SeqIO.parse(SOURCE/'pts1_prediction_tool/datasets/pts1_filtered_peroxisomal_sequences.fasta','fasta'))
        negatives=list(SeqIO.parse(SOURCE/'pts1_prediction_tool/datasets/raw_non_peroxisomal_sequences.fasta','fasta'))
        pos,pos_y=reference._prepare_sequences_for_machine_learning(positives,1)
        neg,neg_y=reference._prepare_sequences_for_machine_learning(negatives,0)
        _,labels,features=reference._merge_and_convert(neg_y,neg,pos_y,pos)
        x,y=features.to_numpy(dtype=np.float64),np.ravel(labels).astype(int)
        samples=[str(r.seq) for r in positives if len(r.seq)>=14]+[str(r.seq) for r in negatives if len(r.seq)>=14][:500]
        expected=np.array([model.decision_function(reference._prepare_sequences(s).to_numpy())[0] for s in samples])
        observed=np.array([portable.decision_function(s) for s in samples])
        with np.load(ROOT/'tools/pts1/svm_weights.npz',allow_pickle=False) as weights:
            coefficient_delta=float(np.max(np.abs(model.dual_coef_-weights['dual_coef'])))
            intercept_delta=float(np.max(np.abs(model.intercept_-weights['intercept'])))
        assert model.support_vectors_.shape==(1308,308)
        assert np.array_equal(model.support_vectors_,portable.sv)
        assert coefficient_delta<1e-10 and intercept_delta<1e-10
        assert np.max(np.abs(expected-observed))<1e-10
        assert np.array_equal(expected>0,observed>0)
        edge_sequences=['MMMMKLSKMLLLSLSKLSKLSKLSKL','ACDEFGHIKLMNSKL','ACDEFGHIKLMBSKL',
                        'ACDEFGHIKLMZSKL','ACDEFGHIKLMJSKL','ACDEFGHIKLMOSKL','ACDEFGHIKLMUSKL','ACDEFGHIKLMXSKL']
        edges=[]
        for sequence in edge_sequences:
            matrix=reference._prepare_sequences(sequence).to_numpy()
            value=float(model.decision_function(matrix)[0])
            other=portable.decision_function(sequence)
            edges.append(dict(sequence=sequence,reference=value,numpy=other,absolute_error=abs(value-other),
                              feature_ones=int(matrix.sum()),prediction=int(value>0)))
        report=dict(upstream_commit='0c1041905d12b96b72063aec5dccfd6073f31435',
                    weight_sha256=hashlib.sha256((ROOT/'tools/pts1/svm_weights.npz').read_bytes()).hexdigest(),
                    environment=dict(python=sys.version,numpy=np.__version__,sklearn=sklearn.__version__),
                    training=dict(positives=len(pos),negatives=len(neg),features=x.shape[1],support_vectors=model.support_vectors_.shape[0],
                                  support_per_class=model.n_support_.tolist(),raw_positive_records=len(positives),raw_negative_records=len(negatives)),
                    numerical=dict(records=len(samples),max_absolute_error=float(np.max(np.abs(expected-observed))),
                                   prediction_mismatches=int(np.count_nonzero((expected>0)!=(observed>0))),
                                   coefficient_max_error=coefficient_delta,intercept_max_error=intercept_delta),
                    edge_cases=edges)
        if not args.baseline:
            assert max(row['absolute_error'] for row in edges)<1e-10
            groups=np.array([''.join(row) for row in pos+neg])
            report['unique_c_terminal_14_groups']=len(set(groups))
            report['repeated_c_terminal_14_records']=len(groups)-len(set(groups))
            for name,splitter in [('stratified_cv',StratifiedKFold(5,shuffle=True,random_state=42)),
                                  ('tail_grouped_cv',StratifiedGroupKFold(5,shuffle=True,random_state=42))]:
                predictions=np.zeros_like(y)
                folds=[]
                split = splitter.split(x,y,groups) if name=='tail_grouped_cv' else splitter.split(x,y)
                for fold,(train,test) in enumerate(split,1):
                    if name=='tail_grouped_cv':
                        assert not set(groups[train]) & set(groups[test])
                    estimator=SVC(kernel='rbf',C=2.,gamma=.092,random_state=42)
                    estimator.fit(x[train],y[train])
                    predictions[test]=estimator.predict(x[test])
                    folds.append(dict(fold=fold,n_test=len(test),**metrics(y[test],predictions[test])))
                report[name]=dict(overall=metrics(y,predictions),folds=folds)
            environment=dict(os.environ,PYTHONPATH=str(SOURCE))
            tests=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(SOURCE/'tests'),
                                  '-p','pts1_prediction_tests.py','-v'],capture_output=True,text=True,env=environment)
            (ROOT/'validation/pts1_upstream_tests.log').write_text(tests.stdout+tests.stderr)
            assert tests.returncode==0,tests.stdout+tests.stderr
            report['upstream_tests']=dict(exit_code=0,count=4)
            report['interpretation']='Internal cross-validation only; exact C-tail grouping does not exclude homology or training-selection bias.'
        report['seconds']=round(time.monotonic()-started,3)
        target=ROOT/'validation'/('pts1_baseline.json' if args.baseline else 'pts1_verification.json')
        atomic_json(target,report)
        print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
