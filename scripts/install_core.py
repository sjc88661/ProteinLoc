"""Install the two user-supplied DTU packages into an isolated, pinned Linux environment."""
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
NATIVE=(ROOT/'.runtime/venv').resolve().parent
prefix=NATIVE/'core'
environment=os.environ.copy()
environment['TORCH_HOME']=str(ROOT/'.runtime/model-cache/torch')
environment['PIP_CONSTRAINT']=str(ROOT/'tools/core-constraints.txt')
environment['HF_HOME']=str(ROOT/'.runtime/model-cache/huggingface')
if not (prefix/'bin/python').exists():
    subprocess.run([str(ROOT/'.runtime/bin/micromamba'),'create','--yes','--prefix',str(prefix),
                    '--file',str(ROOT/'tools/core-environment.yml')],env=environment,check=True)
link=ROOT/'.runtime/core'
if not link.exists():
    link.symlink_to(prefix,target_is_directory=True)
phobius=prefix/'share/phobius-1.01-5'
if not (phobius/'phobius.model').exists():
    shutil.copytree(ROOT/'.runtime/phobius/phobius',phobius,symlinks=True,dirs_exist_ok=True)
source=NATIVE/'core-src'
source.mkdir(parents=True,exist_ok=True)
for archive,folder in [('signalp-6.0i.fast.tar.gz','signalp6_fast'),('deeploc-2.1.All.tar.gz','deeploc2_package')]:
    if not (source/folder).exists():
        with tarfile.open(ROOT/'tools/assets'/archive) as tar:
            tar.extractall(source,filter='data')
python=prefix/'bin/python'
subprocess.run([str(python),'-m','pip','install','-r',str(ROOT/'tools/core-pip-requirements.txt')],
               env=environment,check=True)
for package in [source/'signalp6_fast/signalp-6-package',source/'deeploc2_package']:
    subprocess.run([str(python),'-m','pip','install','--no-build-isolation','--no-deps',str(package)],
                   env=environment,check=True)
weights=source/'signalp6_fast/signalp-6-package/models'
target=Path(subprocess.check_output([str(python),'-c','import signalp; from pathlib import Path; print(Path(signalp.__file__).parent/"model_weights")'],text=True).strip())
for item in weights.iterdir():
    if item.is_file() and not (target/item.name).exists():
        shutil.copy2(item,target/item.name)
print('Pinned local SignalP / DeepLoc / Phobius environment installed.',flush=True)
