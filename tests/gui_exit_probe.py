"""Do not keep the interpreter alive manually: export must survive mainloop exit."""
import sys
import tkinter as tk
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from app import ProteinLocApp
from parsers import read_fasta

folder=Path(sys.argv[1])
root=tk.Tk()
root.withdraw()
app=ProteinLocApp(root,runs_dir=folder/'runs')
app.filtered=read_fasta(''.join(f'>protein_{i}\nM'+'A'*50+'\n' for i in range(1000)))
for record in app.filtered:
    record.update(status='partial',predictions={},annotations=[])
print(app.export_to(folder),flush=True)
root.destroy()
