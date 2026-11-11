import importlib.util
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
import parsers
import pipeline

ROOT=Path(__file__).resolve().parents[1]


class PTS1Tests(unittest.TestCase):
    def test_unknown_category_is_zero_one_hot_not_x(self):
        spec=importlib.util.spec_from_file_location('portable_pts1',ROOT/'tools/pts1/predict_pts1_numpy.py')
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'toy.npz'
            np.savez(path,support_vector_codes=np.zeros((1,14),dtype=np.int8),
                     aa_alphabet=np.array(list('ACDEFGHIKLMNPQRSTVWYXU')),
                     dual_coef=np.array([[1.]]),intercept=np.array([0.]),gamma=np.array(.1),
                     c_terminal_length=np.array(14),classes=np.array([0,1]),C=np.array(2.),
                     kernel=np.array('rbf'),n_support=np.array([0,1]))
            model=module.NumpyPTS1Predictor(path)
            # One missing-category block contributes squared distance 1; X is a distinct one-hot block.
            self.assertAlmostEqual(model.decision_function('B'+'A'*13),0.9048374180359595,places=14)
            self.assertAlmostEqual(model.decision_function('X'+'A'*13),0.8187307530779818,places=14)

    def test_signed_margin_and_candidate_are_not_converted_to_probability(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'prediction.csv'
            path.write_text('id,prediction,is_pts1,decision_function,length,c_terminal_14\nP00000001,1,True,1.234567891234567,26,SLSKLSKLSKLSKL\n')
            self.assertTrue(hasattr(parsers,'parse_pts1_svm'),'缺少 PTS1 SVM 独立输出解析')
            result=parsers.parse_pts1_svm(path)['P00000001']
            self.assertTrue(result['pts1_candidate'])
            self.assertEqual(result['scores']['decision_function'],1.234567891234567)
            self.assertEqual(result['output_tail'],'SLSKLSKLSKLSKL')
            self.assertIn('概率',result['interpretation'])

    def test_tail_model_preserves_u_and_only_rejects_unmodeled_c_terminal_residues(self):
        config=pipeline.load_config()
        self.assertIn('pts1_svm',config['tools'])
        records=parsers.read_fasta('>selenocysteine\nM'+'A'*15+'U'+'A'*10+'SKL\n>ambiguous_tail\nM'+'A'*15+'B'+'A'*10+'SKL\n>ambiguous_nterm\nB'+'A'*20+'SKL\n>short\nMASLK\n')
        with tempfile.TemporaryDirectory() as temp:
            _,inputs,unsupported=pipeline._prepare(records,'pts1_svm',Path(temp),config)
        self.assertEqual(inputs['P00000001']['submitted_sequence'],records[0]['sequence'])
        self.assertEqual(inputs['P00000001']['changed_residues'],[])
        self.assertIn('P00000002',unsupported)
        self.assertIn('P00000003',inputs)
        self.assertIn('P00000004',unsupported)


if __name__=='__main__':
    unittest.main()
