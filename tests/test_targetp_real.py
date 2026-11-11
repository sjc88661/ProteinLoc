import tempfile
import unittest
from pathlib import Path

from pipeline import ROOT, create_job, load_config, read_json, run_worker


@unittest.skipUnless((ROOT/'.runtime/targetp/TargetP-2.0/bin/targetp').exists(),'需要作者 GitHub 的 TargetP Linux 发行包')
class TargetPRealTests(unittest.TestCase):
    def test_reference_hsp60_local_output_is_normalized_and_traceable(self):
        config=load_config()
        config['core_engine']='commands'
        for name,tool in config['tools'].items():
            if name!='targetp':
                tool.update(command=['/not-installed-for-this-model-test'],provenance_paths=[],watch_paths=[],
                            version_command=['/not-installed-for-this-model-test'])
        with tempfile.TemporaryDirectory() as temp:
            job=create_job((ROOT/'examples/yeast_reference.fasta').read_text(),'fungi',config,Path(temp)/'runs')
            run_worker(job.parent)
            data={r['original_id']:r for r in read_json(job/'results.json')}
            for record in data.values():
                self.assertEqual(record['predictions']['targetp']['status'],'completed')
            hsp60=data['P19882']['predictions']['targetp']
            self.assertEqual(hsp60['prediction'],'mTP')
            self.assertEqual(hsp60['raw_prediction'],'MT')
            self.assertEqual(hsp60['cleavage_after'],21)
            self.assertTrue(hsp60['provenance']['verified'])


if __name__=='__main__':
    unittest.main()
