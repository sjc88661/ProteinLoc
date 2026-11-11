"""Run the unchanged upstream PredGPI with NumPy's removed builtin-int alias.

np.int was exactly Python's int (not a different numerical model). This shim
restores that name only, preserving the original algorithm and model files.
"""
import os
import importlib
import sys
from pathlib import Path

import numpy

if 'int' not in numpy.__dict__:
    numpy.int = int

upstream = Path(__file__).resolve().parents[1] / 'upstream/predgpi-master'
os.environ.setdefault('PREDGPI_HOME', str(upstream))
sys.path.insert(0,str(upstream))
sys.argv[0] = str(upstream/'predgpi.py')
# Ordinary import avoids the reproducible long-batch crash of runpy.run_path
# with this runtime. Execute the unchanged official CLI and prediction code.
importlib.import_module('predgpi').main()
