import importlib.util
import hashlib
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from pipeline import ROOT, create_job, load_config, read_json, run_worker


@unittest.skipUnless(importlib.util.find_spec('numpy'), '需要真实 PredGPI 的 numpy 依赖')
class PredGpiRealTests(unittest.TestCase):
    def test_fifty_real_predictions_match_official_golden_output(self):
        # Reuse the first 50 records of the previously verified official output.
        rng=random.Random(20261110)
        chunks=[]
        for i in range(50):
            seq='M'+''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY',k=60))+('SKL' if i%10==0 else 'HDEL' if i%10==1 else 'AAA')
            chunks.append(f'>P{i+1:08d}\n{seq}\n')
        with tempfile.TemporaryDirectory() as temp:
            fasta,output=Path(temp)/'input.faa',Path(temp)/'predgpi.json'
            fasta.write_text(''.join(chunks))
            result=subprocess.run([sys.executable,str(ROOT/'scripts/predgpi_compat.py'),'-f',str(fasta),'-o',str(output),'-m','json'],
                                  capture_output=True,text=True,timeout=180)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(hashlib.sha256(output.read_bytes()).hexdigest(),
                             '3407e6610bda565d712f334df153f43b0d5d09e00e1d14ceb0db2093cddfa0cc')

    def test_official_full_and_cut_fixture_run_with_current_numpy(self):
        text=(ROOT/'upstream/predgpi-master/testdata/test.fasta').read_text()
        config=load_config()
        config['core_engine']='commands'
        for name,tool in config['tools'].items():
            if name!='predgpi':
                tool.update(command=['/not-installed-for-this-model-test'],provenance_paths=[],watch_paths=[],
                            version_command=['/not-installed-for-this-model-test'])
        with tempfile.TemporaryDirectory() as temp:
            job=create_job(text,'other_euk',config,Path(temp)/'runs')
            run_worker(job.parent)
            records=read_json(job/'results.json')
            self.assertEqual([r['original_id'] for r in records],['5NTD_BOOMI','5NTD_BOOMI_cut'])
            for record in records:
                prediction=record['predictions']['predgpi']
                self.assertEqual(prediction['status'],'completed',prediction)
                self.assertIsInstance(prediction['gpi_anchored'],bool)


if __name__=='__main__':
    unittest.main()
