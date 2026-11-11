"""Prepare pinned, unmodified DeepLocPro source and the official ESM-2 cache."""
import hashlib
import json
import shutil
import sys
import zipfile
from pathlib import Path, PurePosixPath

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from parsers import atomic_json

COMMIT='32c0b37c7a31a0686ce022a697f48998fc2997b1'
ASSETS={
    'deeplocpro-32c0b37.zip':'b1b8b80ade332e2f984b34501389edd68fd9660662bace51cb29a267c95e6ab3',
    'esm2_t33_650M_UR50D.pt':'ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c',
    'esm2_t33_650M_UR50D-contact-regression.pt':'8ffe6edbd4173dc8d45c2cd5cb27d43aad77ec26b4c768200c58ae1f96693575'}


def digest(path):
    value=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):
            value.update(chunk)
    return value.hexdigest()


def main():
    native=(ROOT/'.runtime/core').resolve().parent
    manifest={}
    for filename,expected in ASSETS.items():
        path=ROOT/'tools/assets'/filename
        if digest(path)!=expected:
            raise ValueError('安装文件 SHA256 不匹配：'+filename)
        manifest[filename]=dict(bytes=path.stat().st_size,sha256=expected)
    source=native/('deeplocpro-'+COMMIT)
    if not (source/'.proteinloc-source.json').exists():
        with zipfile.ZipFile(ROOT/'tools/assets/deeplocpro-32c0b37.zip') as archive:
            for member in archive.infolist():
                path=PurePosixPath(member.filename)
                if path.is_absolute() or '..' in path.parts or path.parts[0]!='deeplocpro-'+COMMIT:
                    raise ValueError('不安全或错误的 DeepLocPro ZIP 路径')
            archive.extractall(native)
        atomic_json(source/'.proteinloc-source.json',dict(commit=COMMIT,sha256=ASSETS['deeplocpro-32c0b37.zip']))
    if len(list((source/'DeepLocPro/models/checkpoints').glob('model_*')))!=20:
        raise ValueError('DeepLocPro 分类器权重不完整')
    link=ROOT/'.runtime/deeplocpro-source'
    if link.is_symlink() and link.resolve()!=source:
        link.unlink()
    if not link.exists():
        link.symlink_to(source,target_is_directory=True)
    cache=(ROOT/'.runtime/model-cache').resolve()/'torch/hub/checkpoints'
    cache.mkdir(parents=True,exist_ok=True)
    for filename,expected in ASSETS.items():
        if not filename.endswith('.pt'):
            continue
        target=cache/filename
        if not target.exists() or digest(target)!=expected:
            staged=cache/(filename+'.installing')
            shutil.copyfile(ROOT/'tools/assets'/filename,staged)
            if digest(staged)!=expected:
                raise ValueError('复制后的基础权重校验失败：'+filename)
            staged.replace(target)
    launcher=ROOT/'.runtime/bin/deeplocpro'
    launcher.write_text('#!'+str(native/'core/bin/python')+'\nimport sys\nsys.path.insert(0,'+repr(str(ROOT))+')\nfrom scripts.deeplocpro_local import main\nmain()\n')
    launcher.chmod(0o755)
    lock=json.loads((ROOT/'upstream.lock.json').read_text())
    lock['deeplocpro']=dict(url='https://github.com/Jaimomar99/deeplocpro',commit=COMMIT,
                           version='1.0.0',license='CC BY-NC-SA 4.0',assets=manifest)
    atomic_json(ROOT/'upstream.lock.json',lock)
    atomic_json(ROOT/'validation/deeplocpro_install.json',dict(source=str(source),assets=manifest,
                model_url='https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt'))
    print('DeepLocPro source, 20 classifiers and ESM-2 cache prepared.',flush=True)


if __name__=='__main__':
    main()
