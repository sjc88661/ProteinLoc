"""Run the unchanged DeepLocPro model with explicit local device/group settings."""
import argparse
import json
import os
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main(argv=None):
    parser=argparse.ArgumentParser(description='Local DeepLocPro 1.0')
    parser.add_argument('-f','--fasta')
    parser.add_argument('-o','--output',default='./outputs')
    parser.add_argument('-g','--group',choices=['any','negative','positive','archaea'],
                        default=os.environ.get('PROTEINLOC_DEEPLOCPRO_GROUP','any'))
    parser.add_argument('-d','--device',choices=['cpu','cuda'],default='cpu')
    parser.add_argument('-p','--plot',action='store_true')
    parser.add_argument('--version',action='store_true')
    args=parser.parse_args(argv)
    source=(ROOT/'.runtime/deeplocpro-source').resolve()
    if not (source/'DeepLocPro/__init__.py').is_file():
        parser.error('DeepLocPro 未安装；运行 bash scripts/python.sh scripts/install_deeplocpro.py')
    # Upstream embed_batch otherwise moves tokens to CUDA even for a CPU model.
    if args.device=='cpu':
        os.environ['CUDA_VISIBLE_DEVICES']='-1'
    os.environ.setdefault('TORCH_HOME',str((ROOT/'.runtime/model-cache/torch').resolve()))
    sys.path.insert(0,str(source))
    if args.version:
        import importlib.metadata as metadata
        import torch
        print(json.dumps(dict(tool='DeepLocPro',version='1.0.0',commit='32c0b37c7a31a0686ce022a697f48998fc2997b1',
                              python=sys.version,torch=torch.__version__,esm=metadata.version('fair-esm'))))
        return
    if not args.fasta:
        parser.error('-f/--fasta 必填')
    Path(args.output).mkdir(parents=True,exist_ok=True)
    invocation=dict(tool='DeepLocPro',commit='32c0b37c7a31a0686ce022a697f48998fc2997b1',
                    group=args.group,device=args.device,fasta=str(Path(args.fasta).resolve()),
                    source=str(source),cuda_visible_devices=os.environ.get('CUDA_VISIBLE_DEVICES'))
    (Path(args.output)/'proteinloc_invocation.json').write_text(json.dumps(invocation,indent=2),encoding='utf-8')
    from DeepLocPro.deeplocpro import main as predict
    predict(args)


if __name__=='__main__':
    main()
