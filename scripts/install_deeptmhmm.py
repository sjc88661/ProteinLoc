"""Prepare the unchanged academic DeepTMHMM package with isolated legacy dependencies."""
import os
import shutil
import subprocess
import zipfile
from pathlib import Path, PurePosixPath

ROOT=Path(__file__).resolve().parents[1]
native=(ROOT/'.runtime/venv').resolve().parent
source=native/'DeepTMHMM-Academic-License-v1.0'
if not (source/'.proteinloc-extracted').exists():
    with zipfile.ZipFile(ROOT/'tools/assets/DeepTMHMM-Academic-License-v1.0.zip') as archive:
        members=[]
        for member in archive.infolist():
            path=PurePosixPath(member.filename)
            if path.is_absolute() or '..' in path.parts:
                raise ValueError('Unsafe model ZIP member')
            if path.parts[0]=='DeepTMHMM-Academic-License-v1.0':
                members.append(member)
        archive.extractall(native,members=members)
    (source/'.proteinloc-extracted').touch()
prefix=native/'deeptmhmm'
if not (prefix/'bin/python').exists():
    subprocess.run([str(native/'core/bin/python'),'-m','venv','--system-site-packages',str(prefix)],check=True)
environment=os.environ.copy()
environment['PYTHONPATH']=''
environment.pop('PIP_CONSTRAINT',None)
subprocess.run([str(prefix/'bin/python'),'-m','pip','install',
                'numpy==1.23.4','biopython==1.79','fair-esm==0.4.0','h5py==3.7.0','PeptideBuilder==1.1.0'],
               env=environment,check=True)
for name,target in [('deeptmhmm',prefix),('deeptmhmm-source',source)]:
    link=ROOT/'.runtime'/name
    if not link.exists():
        link.symlink_to(target,target_is_directory=True)
print('Local DeepTMHMM academic package prepared.',flush=True)
