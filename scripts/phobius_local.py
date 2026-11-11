"""Use the official Phobius CLI via stdin so input paths never enter its legacy shell command."""
import os
import sys
from pathlib import Path

script=Path(__file__).resolve().parents[1]/'.runtime/phobius/phobius/phobius.pl'
if sys.argv[1:]==['--version']:
    os.execvp('perl',['perl',str(script.resolve()),'-help'])
if len(sys.argv)!=2:
    raise SystemExit('usage: phobius_local.py input.fasta | --version')
with Path(sys.argv[1]).open('rb') as fasta:
    os.dup2(fasta.fileno(),0)
os.execvp('perl',['perl',str(script.resolve()),'-short'])
