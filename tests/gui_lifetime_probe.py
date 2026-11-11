"""Run in its own interpreter because a Tcl thread-lifetime bug aborts the process."""
import gc
import sys
import tempfile
import time
import tkinter as tk
import weakref
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from app import ProteinLocApp
from parsers import read_fasta

with tempfile.TemporaryDirectory() as temp:
    root=tk.Tk()
    root.withdraw()
    app=ProteinLocApp(root,runs_dir=Path(temp)/'runs')
    app.filtered=read_fasta(''.join(f'>protein_{i}\nM'+'A'*50+'\n' for i in range(2000)))
    for record in app.filtered:
        record.update(status='partial',predictions={},annotations=[])
    folder=app.export_to(Path(temp)/'export')
    messages=app.messages
    reference=weakref.ref(app)
    root.destroy()
    del app,root
    gc.collect()
    retained=reference() is not None
    deadline=time.monotonic()+20
    while not (folder/'summary.xlsx').exists() and time.monotonic()<deadline:
        time.sleep(.05)
    if retained:
        print('Background export retained Tk widgets after window destruction',flush=True)
        raise SystemExit(1)
    if not (folder/'summary.xlsx').exists():
        raise SystemExit(2)
    print('Window destroyed safely while export continued',flush=True)
