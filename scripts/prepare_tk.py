"""Ubuntu 24.04 private Tk/Xft runtime; no dpkg database or system changes."""
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
native=Path(sys.prefix).resolve().parent
cache=ROOT/'.runtime/tk-debs'
cache.mkdir(parents=True,exist_ok=True)
packages=('python3-tk','libtcl8.6','libtk8.6','blt','tk8.6-blt2.5','libxss1')
missing=[name for name in packages if not list(cache.glob(name+'_*.deb'))]
if missing:
    release=Path('/etc/os-release').read_text()
    if 'Ubuntu' not in release or '24.04' not in release:
        raise RuntimeError('私有 Tk 自动部署当前针对 Ubuntu 24.04；其他发行版请安装对应 python3-tk 并配置运行环境')
    subprocess.run(['apt-get','download',*missing],cwd=cache,check=True)
target=native/'tk'
target.mkdir(parents=True,exist_ok=True)
for package in cache.glob('*.deb'):
    subprocess.run(['dpkg-deb','-x',str(package),str(target)],check=True)
print('Private Tk/Xft ready:',target)
