"""Exercise saved results in the actual GUI without starting a prediction."""
import argparse
import json
import sys
import time
import tkinter as tk
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from app import ProteinLocApp
from parsers import atomic_json
from pipeline import read_json


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('job',type=Path)
    parser.add_argument('--report',type=Path,default=ROOT/'validation/gui_capacity.json')
    args=parser.parse_args()
    job=args.job.resolve()
    expected=read_json(job/'job.json')['count']
    root=tk.Tk()
    root.withdraw()
    app=ProteinLocApp(root,runs_dir=job.parent)
    app.load_job(job)
    started=time.monotonic()
    while len(app.records)!=expected and time.monotonic()-started<60:
        root.update()
        time.sleep(.02)
    loaded=time.monotonic()-started
    assert len(app.records)==expected
    assert len(app.tree.get_children())==min(200,expected)
    app.turn_page(1)
    next_page_first='P00000201' if expected>200 else 'P00000001'
    assert app.tree.get_children()[0]==next_page_first
    app.search.set('synthetic_stress_00001')
    app.apply_filter()
    assert len(app.filtered)==1
    export=app.export_to(ROOT/'validation/gui-export')
    while not (export/'summary.xlsx').is_file() and time.monotonic()-started<90:
        root.update()
        time.sleep(.02)
    assert len(read_json(export/'results.json'))==1
    assert (export/'summary.xlsx').stat().st_size>0
    report=dict(job=str(job.relative_to(ROOT)),records=len(app.records),load_seconds=round(loaded,3),
                page_size=200,next_page_first_id=next_page_first,filtered_records=1,
                export=str(export.relative_to(ROOT)),excel_written=True)
    root.destroy()
    atomic_json(args.report,report)
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
