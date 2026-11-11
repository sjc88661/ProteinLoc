import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from parsers import parse_phobius
from pipeline import ROOT


@unittest.skipUnless((ROOT/'.runtime/phobius/phobius/phobius.model').exists(),'需要作者提供的 Phobius 模型')
class PhobiusRealTests(unittest.TestCase):
    def test_official_opsin_fixture_uses_64bit_model_and_handles_input_path_spaces(self):
        with tempfile.TemporaryDirectory(prefix='phobius test ') as temp:
            fasta=Path(temp)/'input proteins.fasta'
            fasta.write_bytes((ROOT/'.runtime/phobius/phobius/OPSD_SHEEP').read_bytes())
            result=subprocess.run([sys.executable,str(ROOT/'scripts/phobius_local.py'),str(fasta)],
                                  capture_output=True,text=True,timeout=30)
            self.assertEqual(result.returncode,0,result.stderr)
            raw=Path(temp)/'phobius.txt'
            raw.write_text(result.stdout)
            prediction=parse_phobius(raw)['OPSD_SHEEP']
            self.assertEqual(prediction['tm_count'],7)
            self.assertEqual(prediction['tm_regions'],[[40,63],[75,95],[115,133],[153,173],[203,223],[254,277],[283,308]])
            self.assertFalse(prediction['signal_peptide'])
            self.assertEqual(prediction['c_terminal_side'],'cytoplasmic')


if __name__=='__main__':
    unittest.main()
