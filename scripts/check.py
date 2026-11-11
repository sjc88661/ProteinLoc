"""Record fresh regression-test evidence; no scientific model is simulated here."""
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
folder=ROOT/'validation'
folder.mkdir(exist_ok=True)
command=[sys.executable,'-m','unittest','discover','-s','tests','-v']
result=subprocess.run(command,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
(folder/'tests.log').write_text(result.stdout,encoding='utf-8')
matched=re.search(r'Ran (\d+) tests? in ([\d.]+)s',result.stdout)
evidence={'time_utc':datetime.now(timezone.utc).isoformat(),'command':command,'exit_code':result.returncode,
          'tests':int(matched[1]) if matched else None,'seconds':float(matched[2]) if matched else None,
          'passed':result.returncode==0}
(folder/'tests.json').write_text(json.dumps(evidence,indent=2),encoding='utf-8')
print(result.stdout,end='')
print(json.dumps(evidence,indent=2))
raise SystemExit(result.returncode)
