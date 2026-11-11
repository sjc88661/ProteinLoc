"""Run the unmodified academic DeepTMHMM package locally in its model directory."""
import os
import sys
from pathlib import Path

source=Path(__file__).resolve().parents[1]/'.runtime/deeptmhmm-source'
os.chdir(source.resolve())
os.execv(sys.executable,[sys.executable,str(source.resolve()/'predict.py'),*sys.argv[1:]])
