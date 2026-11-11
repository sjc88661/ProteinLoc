"""Fetch public author sources, TargetP GitHub LFS, and the Phobius academic package."""
import hashlib
import json
import shutil
import sys
import urllib.request
import zipfile
from pathlib import Path, PurePosixPath

ROOT=Path(__file__).resolve().parents[1]
SOURCES={
    'predgpi-master':'https://codeload.github.com/BolognaBiocomp/predgpi/zip/refs/heads/master',
    'signalp-6.0-main':'https://codeload.github.com/fteufel/signalp-6.0/zip/refs/heads/main',
    'phobius_www-master':'https://codeload.github.com/statisticalbiotechnology/phobius_www/zip/refs/heads/master',
    'TargetP-2.0-master':'https://codeload.github.com/JJAlmagro/TargetP-2.0/zip/refs/heads/master',
}
TARGETP_URL='https://media.githubusercontent.com/media/JJAlmagro/TargetP-2.0/master/targetp-2.0.Linux/targetp-2.0.Linux.tar.gz'
TARGETP_SHA='e3767a5b674233be6d19b7f8f67826c096ccdec16d607329a6eb758055142cfe'
PHOBIUS_URL='https://phobius.sbc.su.se/download/phobius101_linux.tgz'
PHOBIUS_SHA='caf6256d04f8977578b8f0914ee4f4a5bf7bd82af1e67af3932c7762a6627f3e'


def fetch(url,path):
    temporary=path.with_suffix(path.suffix+'.download')
    path.parent.mkdir(parents=True,exist_ok=True)
    with urllib.request.urlopen(url,timeout=180) as response, temporary.open('wb') as output:
        shutil.copyfileobj(response,output)
    temporary.replace(path)


def main():
    report=[]
    for name,url in SOURCES.items():
        directory=ROOT/'upstream'/name
        if not directory.is_dir():
            archive=ROOT/'.runtime'/(name+'.zip')
            fetch(url,archive)
            with zipfile.ZipFile(archive) as zipped:
                for entry in zipped.infolist():
                    relative=PurePosixPath(entry.filename)
                    if relative.is_absolute() or '..' in relative.parts:
                        raise ValueError('Unsafe GitHub ZIP member')
                zipped.extractall(ROOT/'upstream')
        report.append({'source':url,'directory':str(directory.relative_to(ROOT)),'downloaded':directory.is_dir()})
    targetp=ROOT/'tools/assets/targetp-2.0.Linux.tar.gz'
    if not targetp.exists():
        fetch(TARGETP_URL,targetp)
    digest=hashlib.sha256()
    with targetp.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            digest.update(block)
    if digest.hexdigest()!=TARGETP_SHA:
        raise ValueError('TargetP author LFS checksum mismatch')
    report.append({'source':TARGETP_URL,'sha256':digest.hexdigest(),'bytes':targetp.stat().st_size,'complete_release':True})
    phobius=ROOT/'tools/assets/phobius101_linux.tgz'
    if not phobius.exists():
        fetch(PHOBIUS_URL,phobius)
    digest=hashlib.sha256(phobius.read_bytes()).hexdigest()
    if digest!=PHOBIUS_SHA:
        raise ValueError('Phobius author package checksum mismatch')
    report.append({'source':PHOBIUS_URL,'sha256':digest,'bytes':phobius.stat().st_size,
                   'complete_release':True,'terms':'private study, education or nonprofit research'})
    path=ROOT/'validation/github_downloads.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))
    print('Source downloads do not imply all model weights are present; see MODEL_SETUP.md.')


if __name__=='__main__':
    main()
