import importlib
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch


@unittest.skipUnless(os.environ.get('DISPLAY'), '需要 WSLg 或 X11 显示环境')
class GuiTests(unittest.TestCase):
    def test_length_and_tm_columns_sort_by_numeric_value(self):
        import app
        import tkinter as tk
        from parsers import read_fasta, summary_row
        with tempfile.TemporaryDirectory() as temp:
            root = tk.Tk()
            root.withdraw()
            window = app.ProteinLocApp(root, runs_dir=Path(temp))
            try:
                records = read_fasta('>short\n' + 'A'*10 + '\n>long\n' + 'A'*100 + '\n>middle\n' + 'A'*22)
                for record, count in zip(records, (2, 10, 1)):
                    record.update(status='completed', predictions={'phobius': dict(status='completed', tm_count=count)})
                window.records = records
                window.rows = {record['internal_id']: summary_row(record) for record in records}
                for column, ascending, descending in (
                    ('Length', ['short', 'middle', 'long'], ['long', 'middle', 'short']),
                    ('TM', ['middle', 'short', 'long'], ['long', 'short', 'middle'])):
                    window.sort_by(column)
                    self.assertEqual([r['original_id'] for r in window.filtered], ascending)
                    window.sort_by(column)
                    self.assertEqual([r['original_id'] for r in window.filtered], descending)
            finally:
                root.destroy()

    def test_switching_jobs_clears_previous_protein_details(self):
        import app
        import tkinter as tk
        from parsers import read_fasta, summary_row
        with tempfile.TemporaryDirectory() as temp:
            root = tk.Tk()
            root.withdraw()
            window = app.ProteinLocApp(root, runs_dir=Path(temp))
            try:
                record = read_fasta('>old\nMAAA')[0]
                window.records = window.filtered = [record]
                window.rows = {record['internal_id']: summary_row(record)}
                window.render_page()
                window.tree.selection_set(record['internal_id'])
                window.show_detail()
                window.set_text(window.log_text, 'previous job log')
                self.assertIn('old', window.detail_text.get('1.0', 'end'))
                window.load_job(Path(temp) / 'new-job')
                self.assertNotIn('"original_id": "old"', window.detail_text.get('1.0', 'end'))
                self.assertNotIn('MAAA', window.sequence_text.get('1.0', 'end'))
                self.assertEqual(window.log_text.get('1.0', 'end').strip(), '')
                self.assertEqual(window.feature_canvas.find_all(), ())
            finally:
                root.destroy()

    def test_switching_output_directory_clears_previous_results(self):
        import app
        import tkinter as tk
        from parsers import read_fasta, summary_row
        with tempfile.TemporaryDirectory() as temp:
            root = tk.Tk()
            root.withdraw()
            window = app.ProteinLocApp(root, runs_dir=Path(temp) / 'first')
            try:
                record = read_fasta('>old\nMAAA')[0]
                window.records = window.filtered = [record]
                window.rows = {record['internal_id']: summary_row(record)}
                window.selected_job = Path(temp) / 'first' / 'old-job'
                window.page = 2
                window.render_page()
                with patch.object(app.filedialog, 'askdirectory', return_value=str(Path(temp) / 'empty')):
                    window.choose_output()
                self.assertEqual(window.records, [])
                self.assertEqual(window.filtered, [])
                self.assertEqual(window.tree.get_children(), ())
                self.assertIsNone(window.selected_job)
                self.assertEqual(window.page, 0)
            finally:
                root.destroy()

    def test_prokaryotic_groups_and_localization_filter_are_available(self):
        import app
        import tkinter as tk
        from parsers import read_fasta,summary_row
        for name in ('gram_negative','gram_positive','archaea','prokaryote'):
            self.assertIn(name,app.GROUPS.values())
        with tempfile.TemporaryDirectory() as temp:
            root=tk.Tk()
            root.withdraw()
            window=app.ProteinLocApp(root,runs_dir=Path(temp))
            try:
                record=read_fasta('>outer\nMAAAAA\n')[0]
                record.update(status='completed',annotations=[],predictions={'deeplocpro':dict(status='completed',localizations=['Outer Membrane'])})
                window.records=[record]
                window.rows={record['internal_id']:summary_row(record)}
                window.location_filter.set('Outer Membrane')
                window.apply_filter()
                self.assertEqual(len(window.filtered),1)
                self.assertIn('DeepLocPro',window.tree.item(record['internal_id'],'values'))
                cytoplasm=read_fasta('>cytoplasm\nMAAAAA\n')[0]
                membrane=read_fasta('>membrane\nMAAAAA\n')[0]
                membrane['internal_id']='P00000002'
                for item,label in ((cytoplasm,'Cytoplasmic'),(membrane,'Cytoplasmic Membrane')):
                    item.update(status='completed',annotations=[],predictions={'deeplocpro':dict(status='completed',localizations=[label])})
                window.records=[cytoplasm,membrane]
                window.rows={item['internal_id']:summary_row(item) for item in window.records}
                window.location_filter.set('Cytoplasmic')
                window.apply_filter()
                self.assertEqual([item['original_id'] for item in window.filtered],['cytoplasm'])
            finally:
                root.destroy()

    def test_sorting_keeps_the_selected_proteins_details(self):
        import app
        import tkinter as tk
        from parsers import read_fasta, summary_row
        with tempfile.TemporaryDirectory() as temp:
            root = tk.Tk()
            root.withdraw()
            window = app.ProteinLocApp(root, runs_dir=Path(temp))
            try:
                window.records = read_fasta('>first\nMAAAA\n>second\nMKK\n')
                window.rows = {r['internal_id']: summary_row(r) for r in window.records}
                window.apply_filter()
                window.tree.selection_set('P00000002')
                root.update()
                window.sort_by('Length')
                window.sort_by('Length')
                root.update()
                self.assertEqual(window.tree.selection(), ('P00000002',))
                self.assertIn('"original_id": "second"', window.detail_text.get('1.0', 'end'))
                self.assertIn('MKK', window.sequence_text.get('1.0', 'end'))
            finally:
                root.destroy()

    def test_paging_and_empty_filter_update_the_details(self):
        import app
        import tkinter as tk
        from parsers import read_fasta, summary_row
        with tempfile.TemporaryDirectory() as temp:
            root = tk.Tk()
            root.withdraw()
            window = app.ProteinLocApp(root, runs_dir=Path(temp))
            try:
                window.records = read_fasta('>first\nMAAAA\n>second\nMKK\n')
                window.rows = {r['internal_id']: summary_row(r) for r in window.records}
                window.page_size = 1
                window.apply_filter()
                window.turn_page(1)
                self.assertEqual(window.tree.selection(), ('P00000002',))
                self.assertIn('MKK', window.sequence_text.get('1.0', 'end'))
                window.search.set('no-matching-protein')
                window.apply_filter()
                self.assertEqual(window.tree.selection(), ())
                self.assertNotIn('"original_id": "second"', window.detail_text.get('1.0', 'end'))
                self.assertNotIn('MKK', window.sequence_text.get('1.0', 'end'))
                self.assertTrue(window.detail_text.get('1.0', 'end').strip())
                self.assertEqual(window.feature_canvas.find_all(), ())
            finally:
                root.destroy()

    def test_export_finishes_when_application_interpreter_exits(self):
        with tempfile.TemporaryDirectory() as temp:
            result=subprocess.run([sys.executable,str(Path(__file__).with_name('gui_exit_probe.py')),temp],capture_output=True,text=True,timeout=30)
            self.assertEqual(result.returncode,0,result.stderr)
            exports=list(Path(temp).glob('export-*/summary.xlsx'))
            self.assertEqual(len(exports),1,'关闭应用中止了已经请求的导出')

    def test_closing_window_during_export_does_not_leave_tk_objects_in_worker_thread(self):
        probe=Path(__file__).with_name('gui_lifetime_probe.py')
        result=subprocess.run([sys.executable,str(probe)],capture_output=True,text=True,timeout=30)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_history_filter_details_and_filtered_export_use_real_saved_records(self):
        try:
            app_module = importlib.import_module('app')
        except ImportError:
            self.fail('Tkinter 应用尚未实现')
        import tkinter as tk
        from parsers import atomic_json, read_fasta, scan_signals
        from pipeline import create_job, load_config, read_json
        with tempfile.TemporaryDirectory() as temp:
            runs = Path(temp) / 'runs'
            text = '>ER protein\nM' + 'A'*45 + 'HDEL\n>PTS protein\nM' + 'A'*45 + 'SKL\n'
            job = create_job(text, 'fungi', load_config(), runs)
            records = read_fasta(text)
            for record in records:
                record['predictions'] = scan_signals(record, {})
                record['status'], record['annotations'] = 'partial', []
            atomic_json(job / 'results.json', records)
            state = read_json(job / 'job.json')
            state['status'] = 'partial'
            atomic_json(job / 'job.json', state)
            root = tk.Tk()
            root.withdraw()
            app = app_module.ProteinLocApp(root, runs_dir=runs)
            try:
                app.load_job(job)
                deadline = time.monotonic() + 15
                while len(app.records) != 2 and time.monotonic() < deadline:
                    root.update()
                    time.sleep(.02)
                self.assertEqual(len(app.records), 2)
                self.assertEqual(len(app.tree.get_children()), 2)
                self.assertEqual(app.tree.selection(), ('P00000001',))
                self.assertIn('ER_luminal_retrieval', app.detail_text.get('1.0', 'end'))
                self.assertIn('HDEL', app.sequence_text.get('1.0', 'end'))
                self.assertTrue(app.feature_canvas.find_all())
                app.search.set('SKL')
                app.apply_filter()
                self.assertEqual(len(app.filtered), 1)
                self.assertEqual(app.filtered[0]['original_id'], 'PTS')
                self.assertEqual(app.tree.selection(), ('P00000002',))
                self.assertIn('PTS1_classic', app.detail_text.get('1.0', 'end'))
                exported = app.export_to(Path(temp) / 'exports')
                deadline = time.monotonic() + 15
                while not (exported / 'summary.xlsx').exists() and time.monotonic() < deadline:
                    root.update()
                    time.sleep(.02)
                self.assertTrue((exported / 'summary.tsv').exists())
                self.assertTrue((exported / 'summary.xlsx').exists())
                self.assertEqual(len(read_json(exported / 'results.json')), 1)
            finally:
                root.destroy()


if __name__ == '__main__':
    unittest.main()
