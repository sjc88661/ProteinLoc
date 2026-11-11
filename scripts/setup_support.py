"""Download public workflow dependencies only; never request DTU licenses or upload sequences."""
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNTIME = ROOT / '.runtime'
META_COMMIT = 'c10d8035603134a0c9916dcbcebdcf51d0559e95'
PREDGPI_COMMIT = '35ac13788771aaf06b7174852476baea9321ae36'


def download(url, target):
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + '.download')
    with urllib.request.urlopen(url, timeout=90) as response, temporary.open('wb') as output:
        shutil.copyfileobj(response, output)
    os.replace(temporary, target)


def main():
    bin_dir = RUNTIME / 'bin'
    bin_dir.mkdir(parents=True, exist_ok=True)
    native = Path(sys.prefix).resolve().parent
    fonts = native / 'fonts'
    fonts.mkdir(parents=True, exist_ok=True)
    for name in ('msyh.ttc','msyhbd.ttc','consola.ttf','segoeui.ttf'):
        source = Path('/mnt/c/Windows/Fonts') / name
        target = fonts / name
        if source.is_file() and not target.exists():
            target.symlink_to(source)
    font_cache = native / 'font-cache'
    font_cache.mkdir(parents=True, exist_ok=True)
    (RUNTIME/'fonts.conf').write_text('<?xml version="1.0"?>\n<!DOCTYPE fontconfig SYSTEM "fonts.dtd">\n'
        '<fontconfig><include ignore_missing="yes">/etc/fonts/fonts.conf</include>'
        f'<dir>{fonts}</dir><cachedir>{font_cache}</cachedir></fontconfig>',encoding='utf-8')
    tk=native/'tk'
    environment=['export FONTCONFIG_FILE='+shlex.quote(str(RUNTIME/'fonts.conf'))]
    if (tk/'usr/lib/python3.12/lib-dynload').exists():
        libraries=str(tk/'usr/lib/x86_64-linux-gnu')+':'+str(tk/'usr/lib')
        pythonpath=str(tk/'usr/lib/python3.12')+':'+str(tk/'usr/lib/python3.12/lib-dynload')
        environment.extend(['export LD_LIBRARY_PATH='+shlex.quote(libraries)+'"${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"',
            'export PYTHONPATH='+shlex.quote(pythonpath)+'"${PYTHONPATH:+:$PYTHONPATH}"',
            'export TCL_LIBRARY='+shlex.quote(str(tk/'usr/share/tcltk/tcl8.6')),
            'export TK_LIBRARY='+shlex.quote(str(tk/'usr/share/tcltk/tk8.6'))])
    (RUNTIME/'env.sh').write_text('\n'.join(environment)+'\n',encoding='utf-8')
    for name in ('mamba', 'conda-envs', 'nextflow'):
        target = native / name
        target.mkdir(parents=True, exist_ok=True)
        link = RUNTIME / name
        if not link.exists() and not link.is_symlink():
            link.symlink_to(target, target_is_directory=True)
    (native/'model-cache').mkdir(parents=True,exist_ok=True)
    if not (RUNTIME/'model-cache').exists():
        (RUNTIME/'model-cache').symlink_to(native/'model-cache',target_is_directory=True)
    targetp_archive=ROOT/'tools/assets/targetp-2.0.Linux.tar.gz'
    if targetp_archive.exists():
        expected='e3767a5b674233be6d19b7f8f67826c096ccdec16d607329a6eb758055142cfe'
        targetp_root=native/'targetp'
        if not (targetp_root/'TargetP-2.0/bin/targetp').exists():
            digest=hashlib.sha256()
            with targetp_archive.open('rb') as stream:
                for block in iter(lambda:stream.read(1024*1024),b''):
                    digest.update(block)
            if digest.hexdigest()!=expected:
                raise RuntimeError('TargetP GitHub LFS checksum mismatch')
            targetp_root.mkdir(parents=True,exist_ok=True)
            with tarfile.open(targetp_archive) as tar:
                tar.extractall(targetp_root,filter='data')
        if not (RUNTIME/'targetp').exists():
            (RUNTIME/'targetp').symlink_to(targetp_root,target_is_directory=True)
    phobius_archive=ROOT/'tools/assets/phobius101_linux.tgz'
    if phobius_archive.exists():
        phobius_root=native/'phobius'
        if not (phobius_root/'phobius/phobius.model').exists():
            if hashlib.sha256(phobius_archive.read_bytes()).hexdigest()!='caf6256d04f8977578b8f0914ee4f4a5bf7bd82af1e67af3932c7762a6627f3e':
                raise RuntimeError('Phobius author package checksum mismatch')
            phobius_root.mkdir(parents=True,exist_ok=True)
            with tarfile.open(phobius_archive) as tar:
                tar.extractall(phobius_root,filter='data')
        decoder=phobius_root/'phobius/decodeanhmm'
        if not decoder.is_symlink():
            decoder.rename(decoder.with_name('decodeanhmm.32bit'))
            decoder.symlink_to('decodeanhmm.64bit')
        (phobius_root/'phobius/decodeanhmm.64bit').chmod(0o755)
        if not (RUNTIME/'phobius').exists():
            (RUNTIME/'phobius').symlink_to(phobius_root,target_is_directory=True)
    meta = ROOT / 'upstream/metaLoc'
    if not (meta / 'main.nf').exists():
        subprocess.run(['git','clone','https://github.com/scottc-bio/metaLoc.git',str(meta)],check=True)
        subprocess.run(['git','-C',str(meta),'checkout',META_COMMIT],check=True)
    actual = subprocess.check_output(['git','-C',str(meta),'rev-parse','HEAD'],text=True).strip()
    if actual != META_COMMIT:
        raise RuntimeError('metaLoc commit changed; do not overwrite an existing checkout automatically')
    nextflow = bin_dir / 'nextflow'
    if not nextflow.exists():
        download('https://github.com/nextflow-io/nextflow/releases/download/v25.10.2/nextflow-25.10.2-dist', nextflow)
    nextflow.chmod(0o755)
    micromamba = bin_dir / 'micromamba'
    if not micromamba.exists():
        archive = RUNTIME / 'micromamba.tar.bz2'
        download('https://micro.mamba.pm/api/micromamba/linux-64/latest', archive)
        with tarfile.open(archive) as tar:
            source = tar.extractfile('bin/micromamba')
            if source is None:
                raise RuntimeError('micromamba archive is incomplete')
            with micromamba.open('wb') as output:
                shutil.copyfileobj(source,output)
    micromamba.chmod(0o755)
    predgpi = ROOT / 'upstream/predgpi-master'
    if not (predgpi / 'predgpi.py').exists():
        subprocess.run(['git','clone','https://github.com/BolognaBiocomp/predgpi.git',str(predgpi)],check=True)
        subprocess.run(['git','-C',str(predgpi),'checkout',PREDGPI_COMMIT],check=True)
    desktop = Path.home() / '.local/share/applications/proteinloc.desktop'
    desktop.parent.mkdir(parents=True, exist_ok=True)
    desktop.write_text('[Desktop Entry]\nType=Application\nName=ProteinLoc\nComment=Protein localization analysis\n'
                       f'Exec=bash "{ROOT}/run.sh"\nPath={ROOT}\nTerminal=false\nCategories=Science;Biology;\n', encoding='utf-8')
    lock_path = ROOT/'upstream.lock.json'
    lock = json.loads(lock_path.read_text(encoding='utf-8')) if lock_path.exists() else {}
    lock.update({'metaloc': {'url':'https://github.com/scottc-bio/metaLoc','commit':actual},
            'predgpi': {'url':'https://github.com/BolognaBiocomp/predgpi','distribution':'official source archive',
                        'commit':PREDGPI_COMMIT,
                        'source_sha256':hashlib.sha256((predgpi/'predgpi.py').read_bytes()).hexdigest()},
            'nextflow': {'version':'25.10.2','sha256':hashlib.sha256(nextflow.read_bytes()).hexdigest()},
            'micromamba': {'version':subprocess.check_output([str(micromamba),'--version'],text=True).strip(),
                           'sha256':hashlib.sha256(micromamba.read_bytes()).hexdigest()}})
    if targetp_archive.exists():
        lock['targetp']={'url':'https://media.githubusercontent.com/media/JJAlmagro/TargetP-2.0/master/targetp-2.0.Linux/targetp-2.0.Linux.tar.gz',
                        'sha256':'e3767a5b674233be6d19b7f8f67826c096ccdec16d607329a6eb758055142cfe','version':'2.0'}
    if phobius_archive.exists():
        lock['phobius']={'url':'https://phobius.sbc.su.se/download/phobius101_linux.tgz',
                        'sha256':hashlib.sha256(phobius_archive.read_bytes()).hexdigest(),'version':'1.01',
                        'decoder':'bundled decodeanhmm.64bit; original 32-bit decoder retained'}
    for name,filename,version in [('signalp','signalp-6.0i.fast.tar.gz','6.0 fast'),
                                  ('deeploc','deeploc-2.1.All.tar.gz','2.1 Fast'),
                                  ('deeptmhmm','DeepTMHMM-Academic-License-v1.0.zip','1.0 Academic')]:
        package=ROOT/'tools/assets'/filename
        if package.exists():
            digest=hashlib.sha256()
            with package.open('rb') as stream:
                for block in iter(lambda:stream.read(1024*1024),b''):
                    digest.update(block)
            lock[name]=dict(version=version,source='Locally supplied authorized package: '+filename,
                            bytes=package.stat().st_size,sha256=digest.hexdigest())
    lock['core_environment']=dict(file='tools/core-environment.yml',
            sha256=hashlib.sha256((ROOT/'tools/core-environment.yml').read_bytes()).hexdigest(),
            constraints_sha256=hashlib.sha256((ROOT/'tools/core-constraints.txt').read_bytes()).hexdigest())
    (ROOT/'upstream.lock.json').write_text(json.dumps(lock,indent=2),encoding='utf-8')
    print('Public workflow dependencies prepared; licensed model packages still checked by doctor.')


if __name__ == '__main__':
    main()
