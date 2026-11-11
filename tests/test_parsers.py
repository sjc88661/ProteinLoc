import importlib
import json
import tempfile
import unittest
from pathlib import Path


class DataTests(unittest.TestCase):
    def setUp(self):
        try:
            self.p = importlib.import_module('parsers')
        except ImportError:
            self.fail('FASTA/结果解析和基序注释尚未实现')

    def file(self, text):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        path = Path(folder.name) / 'result.txt'
        path.write_text(text, encoding='utf-8')
        return path

    def test_duplicate_headers_keep_distinct_ids_and_invalid_records(self):
        rows = self.p.read_fasta('>dup first\nACD\n>dup second\nMZZ\n>bad\nMA1\n>empty\n')
        self.assertEqual(len(rows), 4)
        self.assertEqual([r['internal_id'] for r in rows], ['P00000001', 'P00000002', 'P00000003', 'P00000004'])
        self.assertEqual(rows[0]['header'], 'dup first')
        self.assertEqual(rows[1]['sequence'], 'MZZ')
        self.assertTrue(rows[0]['warnings'])
        self.assertTrue(rows[1]['warnings'])
        self.assertTrue(rows[2]['errors'])
        self.assertTrue(rows[3]['errors'])

    def test_signalp_keeps_scores_and_uses_ids_not_order(self):
        path = self.file('# SignalP-6.0\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\nP00000002\tOTHER\t0.98\t0.02\t\nP00000001\tSP(Sec/SPI)\t0.01\t0.99\tCS pos: 20-21. Pr: 0.9\n')
        rows = self.p.parse_signalp(path)
        self.assertTrue(rows['P00000001']['signal_peptide'])
        self.assertEqual(rows['P00000001']['cleavage_after'], 20)
        self.assertFalse(rows['P00000002']['signal_peptide'])
        self.assertEqual(rows['P00000001']['scores']['SP(Sec/SPI)'], 0.99)

    def test_duplicate_model_ids_are_rejected(self):
        path = self.file('Protein_ID,Localizations\nP1,Cytoplasm\nP1,Nucleus\n')
        with self.assertRaises(ValueError):
            self.p.parse_deeploc(path)

    def test_phobius_retains_exact_tm_regions_and_tail(self):
        path = self.file('SEQENCE ID TM SP PREDICTION\nP1 2 Y n2-8c13/14o30-50i70-90o\nP2 0 0 i\n')
        rows = self.p.parse_phobius(path)
        self.assertEqual(rows['P1']['tm_regions'], [[30, 50], [70, 90]])
        self.assertEqual(rows['P1']['c_terminal_side'], 'non_cytoplasmic')
        self.assertFalse(rows['P2']['signal_peptide'])

    def test_phobius_incomplete_topology_is_not_success(self):
        with self.assertRaises(ValueError):
            self.p.parse_phobius(self.file('ID TM SP PREDICTION\nP1 2 0 o30-50i\n'))

    def test_deeptmhmm_keeps_alpha_beta_signal_and_residue_coordinates(self):
        rows=self.p.parse_deeptmhmm(self.file('>P2 | BETA\nAAAAAAAAAAAA\nOOBBBPPPBBBO\n>P1 | SP+TM\nMAAAAAAAAA\nSSSOOMMMII\n'))
        self.assertEqual(rows['P2']['beta_regions'],[[3,5],[9,11]])
        self.assertEqual(rows['P2']['tm_count'],0)
        self.assertEqual(rows['P1']['tm_regions'],[[6,8]])
        self.assertEqual(rows['P1']['signal_regions'],[[1,3]])
        self.assertTrue(rows['P1']['signal_peptide'])
        self.assertEqual(rows['P1']['c_terminal_side'],'cytoplasmic')

    def test_deeptmhmm_incomplete_or_duplicate_output_is_not_negative(self):
        for text in ['>P1 | GLOB\nMAAAA\nIII\n','>P1 | GLOB\nMA\nII\n>P1 | GLOB\nMA\nII\n']:
            with self.assertRaises(ValueError):
                self.p.parse_deeptmhmm(self.file(text))

    def test_deeploc_membrane_types_from_actual_21_header_are_retained(self):
        rows=self.p.parse_deeploc(self.file('Protein_ID,Localizations,Signals,Membrane types,Peripheral,Soluble\nP1,Plastid,Chloroplast transit peptide,Peripheral|Soluble,0.6678,0.6309\n'))
        self.assertEqual(rows['P1']['membrane_types'],['Peripheral','Soluble'])
        self.assertEqual(rows['P1']['scores']['Soluble'],0.6309)

    def test_deeploc_multilabel_and_targetp_cleavage(self):
        rows = self.p.parse_deeploc(self.file('Protein_ID,Localizations,Cytoplasm,Nucleus\nP2,"Cytoplasm;Nucleus",0.8,0.7\n'))
        self.assertEqual(rows['P2']['localizations'], ['Cytoplasm', 'Nucleus'])
        rows = self.p.parse_targetp(self.file('# TargetP-2.0\n# ID\tPrediction\tOTHER\tSP\tmTP\tCS Position\nP2\tmTP\t0.1\t0.1\t0.8\tCS pos: 31-32. Pr: 0.8\n'))
        self.assertEqual(rows['P2']['cleavage_after'], 31)

    def test_targetp_local_mt_label_is_normalized_without_losing_raw_label(self):
        rows=self.p.parse_targetp(self.file('# TargetP-2.0\n# ID\tPrediction\tOTHER\tSP\tMT\tCS Position\nP2\tMT\t0.1\t0.1\t0.8\tCS pos: 21-22. Pr: 0.8\n'))
        self.assertEqual(rows['P2']['prediction'],'mTP')
        self.assertEqual(rows['P2']['raw_prediction'],'MT')

    def test_predgpi_negative_is_explicit_and_omega_is_not_offset(self):
        data = [{'accession': 'P2', 'sequence': {'length': 100}, 'features': []},
                {'accession': 'P1', 'features': [{'description': 'GPI-anchor', 'begin': 83, 'end': 83, 'score': 0.7}]}]
        rows = self.p.parse_predgpi(self.file(json.dumps(data)))
        self.assertFalse(rows['P2']['gpi_anchored'])
        self.assertEqual(rows['P1']['omega_site'], 83)
        self.assertEqual(rows['P1']['scores']['category_score'], 0.7)

    def test_er_terminal_motifs_require_background_for_interpretation(self):
        record = self.p.read_fasta('>p\nM' + 'A' * 40 + 'HDEL')[0]
        result = self.p.scan_signals(record, {})
        self.assertEqual(result['er']['matches'][0]['sequence'], 'HDEL')
        self.assertFalse(result['er']['matches'][0]['background_supported'])
        supported = self.p.scan_signals(record, {'signalp': {'status': 'completed', 'signal_peptide': True}, 'phobius': {'status': 'completed', 'tm_count': 0}})
        self.assertTrue(supported['er']['matches'][0]['background_supported'])
        internal = self.p.read_fasta('>p\nMHDELA')[0]
        self.assertEqual(self.p.scan_signals(internal, {})['er']['matches'], [])

    def test_pts2_window_and_pts1_terminal_position(self):
        record = self.p.read_fasta('>p\nMRLAAAAAHL' + 'A' * 25 + 'SKL')[0]
        pts = self.p.scan_signals(record, {})['pts']['matches']
        self.assertEqual([(x['kind'], x['start'], x['end']) for x in pts], [('PTS1_classic', 36, 38), ('PTS2_candidate', 2, 10)])
        internal = self.p.read_fasta('>p\nM' + 'A' * 30 + 'RLAAAAAHLSKLA')[0]
        self.assertEqual(self.p.scan_signals(internal, {})['pts']['matches'], [])

    def test_exports_keep_every_record_and_protect_excel_formulas(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        records = self.p.read_fasta('>=formula\nMAAA\n>bad\nMA1\n')
        for r in records:
            r['predictions'] = {}
            r['annotations'] = []
            r['status'] = 'partial'
        self.p.export_results(records, Path(folder.name))
        import csv
        with open(Path(folder.name) / 'summary.tsv', encoding='utf-8') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['Protein ID'], "'=formula")
        data = json.loads((Path(folder.name) / 'results.json').read_text())
        self.assertEqual(data[0]['original_id'], '=formula')


if __name__ == '__main__':
    unittest.main()
