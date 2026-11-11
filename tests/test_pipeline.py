import copy
import importlib
import json
import os
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class PipelineTests(unittest.TestCase):
    def setUp(self):
        try:
            self.p = importlib.import_module('pipeline')
        except ImportError:
            self.fail('持久化任务和本地执行流程尚未实现')
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.runs = Path(self.temp.name) / 'runs'
        self.cfg = self.p.load_config()
        self.cfg['core_engine'] = 'commands'
        self.cfg['metaloc']['provenance_paths']=[str(ROOT/'tests/fixture_tool.py')]
        self.cfg['metaloc']['version_command']=[sys.executable,'--version']
        self.cfg['metaloc']['assets']={key:str(ROOT/'tests/fixture_tool.py') for key in self.cfg['metaloc']['assets']}
        for name, tool in self.cfg['tools'].items():
            tool['command'] = [sys.executable, str(ROOT / 'tests/fixture_tool.py'), '--tool', name, '--input', '{input}', '--output', '{output}']
            tool.pop('stdout_file', None)
            tool.pop('env', None)
            tool.pop('watch_paths', None)
            tool.pop('required_files', None)
            tool['provenance_paths'] = [str(ROOT / 'tests/fixture_tool.py')]
            tool['version_command'] = [sys.executable, '--version']
        self.text = '>same first\nM' + 'A'*45 + 'HDEL\n>same second\nM' + 'A'*45 + 'SKL\n'

    def job(self, cfg=None):
        return self.p.create_job(self.text, 'fungi', cfg or self.cfg, self.runs)

    def state(self, job):
        return self.p.read_json(job / 'job.json')

    def test_complete_pipeline_preserves_duplicate_ids_and_reordered_predictions(self):
        job = self.job()
        self.p.run_worker(self.runs)
        self.assertEqual(self.state(job)['status'], 'completed')
        data = self.p.read_json(job / 'results.json')
        self.assertEqual(len(data), 2)
        self.assertEqual(data[0]['header'], 'same first')
        self.assertEqual(data[1]['predictions']['pts']['matches'][0]['kind'], 'PTS1_classic')
        self.assertEqual(set(data[0]['predictions']), {'signalp','phobius','deeploc','targetp','predgpi','deeptmhmm','er','pts','pts1_svm'})
        self.assertTrue((job / 'summary.tsv').exists())

    def test_prokaryotic_job_completes_required_modules_and_keeps_eukaryote_functions_inapplicable(self):
        config=copy.deepcopy(self.cfg)
        job=self.p.create_job(self.text,'gram_negative',config,self.runs)
        self.p.run_worker(self.runs)
        self.assertEqual(self.state(job)['status'],'completed')
        records=self.p.read_json(job/'results.json')
        self.assertEqual(records[0]['predictions']['deeplocpro']['localizations'],['Periplasmic'])
        for name in ('deeploc','targetp','predgpi','er','pts','pts1_svm'):
            self.assertFalse(records[0]['predictions'][name]['applicable'])
            self.assertEqual(records[0]['predictions'][name]['status'],'unsupported')
            self.assertFalse((job/'raw'/name).exists())
        self.assertEqual(self.state(job)['completed_proteins'],2)

    def test_missing_prokaryotic_localization_row_remains_partial(self):
        config=copy.deepcopy(self.cfg)
        config['tools']['deeplocpro']['command'].append('--omit-last')
        job=self.p.create_job(self.text,'archaea',config,self.runs)
        self.p.run_worker(self.runs)
        self.assertEqual(self.state(job)['status'],'partial')
        records=self.p.read_json(job/'results.json')
        prediction=records[0]['predictions'].get('deeplocpro',{})
        self.assertEqual(prediction.get('status'),'failed')
        self.assertNotIn('localizations',prediction)

    def test_metaloc_prokaryotic_branch_receives_specific_group_and_other_signalp_mode(self):
        config=copy.deepcopy(self.cfg)
        config['core_engine']='metaloc'
        config['metaloc']['env'].pop('PROTEINLOC_DEEPLOCPRO_GROUP',None)
        job=self.p.create_job(self.text,'gram_positive',config,self.runs)
        state=self.state(job)
        state['attempt']=1
        records=self.p.read_json(job/'records.json')
        calls=[]
        def invoke(job,state,command,cwd,env):
            calls.append((command,env))
            return 0
        with patch.object(self.p,'_meta_problems',return_value=[]), patch.object(self.p,'_run_process',side_effect=invoke):
            results,_=self.p._metaloc(job,state,records)
        command,environment=calls[0]
        self.assertEqual(command[command.index('--organism')+1],'other')
        self.assertEqual(environment.get('PROTEINLOC_DEEPLOCPRO_GROUP'),'positive')
        self.assertIn('deeplocpro',results)
        self.assertNotIn('deeploc',results)

    def test_long_prokaryotic_input_retains_other_tools_without_sending_it_to_deeplocpro(self):
        config=copy.deepcopy(self.cfg)
        config['core_engine']='metaloc'
        config['metaloc']['command']=[sys.executable,'-c','raise SystemExit(99)']
        config['tools']['deeplocpro']['max_length']=60
        text='>long_bacterium\nM'+'A'*60+'\n'
        job=self.p.create_job(text,'gram_negative',config,self.runs)
        with patch.object(self.p,'_meta_problems',return_value=[]):
            self.p.run_worker(self.runs)
        record=self.p.read_json(job/'results.json')[0]
        self.assertEqual(record['predictions']['deeplocpro']['status'],'unsupported')
        self.assertIn('61 aa',record['predictions']['deeplocpro']['reason'])
        self.assertIn('60 aa 上限',record['predictions']['deeplocpro']['reason'])
        self.assertEqual(record['predictions']['signalp']['status'],'completed')
        self.assertEqual(record['predictions']['phobius']['status'],'completed')
        self.assertEqual(self.state(job)['status'],'partial')

    def test_identical_prepared_fasta_keeps_mtime_for_nextflow_resume(self):
        records = self.p.read_fasta(self.text)
        folder = Path(self.temp.name) / 'prepared'
        fasta, _, _ = self.p._prepare(records, 'core', folder, self.cfg)
        original_mtime = 1_600_000_000_000_000_000
        os.utime(fasta, ns=(original_mtime, original_mtime))
        self.p._prepare(records, 'core', folder, self.cfg)
        self.assertEqual(fasta.stat().st_mtime_ns, original_mtime)
        records[0]['sequence'] += 'V'
        records[0]['sequence_length'] += 1
        self.p._prepare(records, 'core', folder, self.cfg)
        self.assertIn('HDELV\n', fasta.read_text())
        self.assertNotEqual(fasta.stat().st_mtime_ns, original_mtime)

    def test_signalp_internal_window_is_distinct_from_full_submitted_sequence(self):
        self.text='>long\nM'+'A'*99+'\n'
        job=self.job()
        self.p.run_worker(self.runs)
        record=self.p.read_json(job/'results.json')[0]
        prediction=record['predictions']['signalp']
        self.assertEqual(prediction['submitted_length'],100)
        self.assertEqual(prediction['analyzed_ranges'],[[1,70]])
        self.assertEqual(len(record['sequence']),100)

    def test_failure_is_partial_and_resume_reuses_unchanged_modules(self):
        self.cfg['tools']['targetp']['command'].append('--fail')
        job = self.job()
        self.p.run_worker(self.runs)
        self.assertEqual(self.state(job)['status'], 'partial')
        data = self.p.read_json(job / 'results.json')
        self.assertEqual(data[0]['predictions']['targetp']['status'], 'failed')
        self.assertFalse('prediction' in data[0]['predictions']['targetp'])
        self.assertEqual(data[0]['predictions']['predgpi']['status'], 'completed')
        signalp_attempts = len(list((job / 'raw/signalp').iterdir()))
        self.cfg['tools']['targetp']['command'].remove('--fail')
        self.p.resume_job(job, self.cfg)
        self.p.run_worker(self.runs)
        self.assertEqual(self.state(job)['status'], 'completed')
        self.assertEqual(len(list((job / 'raw/signalp').iterdir())), signalp_attempts)

    def test_unknown_output_id_fails_module_instead_of_misassigning(self):
        self.cfg['tools']['deeploc']['command'].append('--extra')
        job = self.job()
        self.p.run_worker(self.runs)
        data = self.p.read_json(job / 'results.json')
        self.assertEqual(data[0]['predictions']['deeploc']['status'], 'failed')
        self.assertEqual(self.state(job)['status'], 'partial')

    def test_missing_output_row_is_not_negative(self):
        self.cfg['tools']['predgpi']['command'].append('--omit-last')
        job = self.job()
        self.p.run_worker(self.runs)
        data = self.p.read_json(job / 'results.json')
        self.assertEqual(data[0]['predictions']['predgpi']['status'], 'failed')
        self.assertNotIn('gpi_anchored', data[0]['predictions']['predgpi'])

    def test_resume_refuses_tampered_original_input(self):
        job = self.job()
        self.p.run_worker(self.runs)
        (job / 'original.fasta').write_text('>changed\nMAAAA\n')
        with self.assertRaises(ValueError):
            self.p.resume_job(job, self.cfg)

    def test_missing_real_tools_preserve_local_annotation_and_failure_records(self):
        cfg = self.p.load_config()
        cfg['metaloc']['assets'] = {k: '/definitely/missing/'+k for k in cfg['metaloc']['assets']}
        cfg['metaloc']['provenance_paths']=[str(ROOT/'tests/fixture_tool.py')]
        cfg['metaloc']['version_command']=[sys.executable,'--version']
        for name in ('signalp','deeploc','targetp','deeptmhmm'):
            cfg['tools'][name]['command'][0]='/definitely/missing/'+name
            cfg['tools'][name]['provenance_paths']=[]
            cfg['tools'][name]['version_command']=['/definitely/missing/'+name]
        job = self.job(cfg)
        self.p.run_worker(self.runs)
        state = self.state(job)
        self.assertEqual(state['status'], 'partial')
        data = self.p.read_json(job / 'results.json')
        self.assertEqual(data[0]['predictions']['er']['status'], 'completed')
        self.assertEqual(data[0]['predictions']['signalp']['status'], 'failed')

    def test_unready_metaloc_uses_available_local_core_tools(self):
        config=copy.deepcopy(self.cfg)
        config['core_engine']='metaloc'
        config['metaloc']['fallback_local_commands']=True
        config['metaloc']['assets']={k:'/definitely/missing/'+k for k in config['metaloc']['assets']}
        job=self.job(config)
        state=self.state(job)
        records=self.p.read_json(job/'records.json')
        def direct(job,state,records,name):
            return {r['internal_id']:dict(status='completed' if name=='phobius' else 'failed') for r in records},False
        with patch.object(self.p,'_direct',side_effect=direct) as run:
            results,cached=self.p._metaloc(job,state,records)
        self.assertEqual(results['phobius'][records[0]['internal_id']]['status'],'completed')
        self.assertEqual([call.args[-1] for call in run.call_args_list],list(self.p.CORE))
        self.assertEqual(state['core_execution']['engine'],'local_commands')
        self.assertFalse(cached)

    def test_cancelled_fallback_keeps_completed_core_predictions(self):
        config = copy.deepcopy(self.cfg)
        config['core_engine'] = 'metaloc'
        job = self.job(config)
        run_process = self.p._run_process

        def cancel_second_tool(job, state, command, cwd, *args, **kwargs):
            if 'phobius' in command:
                raise self.p.Cancelled()
            return run_process(job, state, command, cwd, *args, **kwargs)

        with patch.object(self.p, '_meta_problems', return_value=['missing metaLoc dependency']), \
             patch.object(self.p, '_run_process', side_effect=cancel_second_tool):
            self.p.run_worker(self.runs)
        records = self.p.read_json(job / 'results.json')
        self.assertEqual(self.state(job)['status'], 'cancelled')
        self.assertEqual(records[0]['predictions']['signalp']['status'], 'completed')
        self.assertEqual(records[0]['predictions']['phobius']['status'], 'cancelled')
        self.assertEqual(self.state(job)['modules']['signalp']['completed'], 2)

    def test_cancelled_module_summary_keeps_cancelled_status(self):
        result = {'P00000001': dict(status='cancelled', reason='任务取消')}
        self.assertEqual(self.p._module_summary(result)['status'], 'cancelled')

    def test_cancelled_metaloc_keeps_published_independent_outputs(self):
        config = copy.deepcopy(self.cfg)
        config['core_engine'] = 'metaloc'
        job = self.job(config)

        def publish_then_cancel(job, state, command, cwd, env):
            folder = cwd / 'results' / 'signalp'
            folder.mkdir(parents=True)
            (folder / 'signalp_results.txt').write_text(
                '# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n'
                'P00000001\tOTHER\t0.99\t0.01\t\n'
                'P00000002\tOTHER\t0.99\t0.01\t\n')
            raise self.p.Cancelled()

        with patch.object(self.p, '_meta_problems', return_value=[]), \
             patch.object(self.p, '_run_process', side_effect=publish_then_cancel):
            self.p.run_worker(self.runs)
        records = self.p.read_json(job / 'results.json')
        self.assertEqual(self.state(job)['status'], 'cancelled')
        self.assertEqual(records[0]['predictions']['signalp']['status'], 'completed')
        self.assertEqual(records[0]['predictions']['signalp']['engine'], 'metaloc')
        self.assertEqual(records[0]['predictions']['phobius']['status'], 'cancelled')
        self.assertFalse(self.p.read_json(job / 'cache/core.json')['complete'])

    def test_short_eukaryotic_input_uses_each_core_tools_own_range(self):
        config = copy.deepcopy(self.cfg)
        config['core_engine'] = 'metaloc'
        config['metaloc']['command'] = [sys.executable, '-c', 'raise SystemExit(99)']
        text = '>short\nMAAAA\n>normal\nM' + 'A'*49 + '\n'
        job = self.p.create_job(text, 'fungi', config, self.runs)
        with patch.object(self.p, '_meta_problems', return_value=[]):
            self.p.run_worker(self.runs)
        short, normal = self.p.read_json(job / 'results.json')
        self.assertEqual(short['predictions']['phobius']['status'], 'completed')
        self.assertEqual(normal['predictions']['phobius']['status'], 'completed')
        for name in ('signalp', 'deeploc'):
            self.assertEqual(short['predictions'][name]['status'], 'unsupported')
            self.assertEqual(normal['predictions'][name]['status'], 'completed')
        self.assertEqual(normal['status'], 'completed')

    def test_doctor_reports_local_fallback_readiness_per_tool(self):
        config=copy.deepcopy(self.cfg)
        config['core_engine']='metaloc'
        config['metaloc']['fallback_local_commands']=True
        config['metaloc']['assets']={k:'/definitely/missing/'+k for k in config['metaloc']['assets']}
        for name in ('signalp','deeploc'):
            config['tools'][name]['command'][0]='/definitely/missing/'+name
        report={r['module']:r for r in self.p.check_environment(config)}
        self.assertTrue(report['phobius']['ready'])
        self.assertEqual(report['phobius']['engine'],'local_command')
        self.assertTrue(report['phobius']['notes'])
        self.assertFalse(report['signalp']['ready'])

    def test_metaloc_uses_absolute_environment_in_its_independent_working_directory(self):
        config=copy.deepcopy(self.cfg)
        config['core_engine']='metaloc'
        environment=Path(self.temp.name)/'model-environment'
        config['metaloc']['runtime_environment']=str(environment)
        job=self.job(config)
        state=self.state(job)
        records=self.p.read_json(job/'records.json')
        with patch.object(self.p,'_signature',return_value=('test',dict(verified=False))), \
             patch.object(self.p,'_meta_problems',return_value=[]), \
             patch.object(self.p,'_run_process',return_value=1):
            self.p._metaloc(job,state,records)
        override=(job/'raw/core/runtime.config').read_text()
        self.assertIn("conda = '"+str(environment)+"'",override)

    def test_metaloc_retry_archives_the_previous_execution_trace(self):
        config = copy.deepcopy(self.cfg)
        config['core_engine'] = 'metaloc'
        job = self.job(config)
        output = job / 'raw/core'
        output.mkdir(parents=True)
        (output / 'trace.tsv').write_text('previous execution trace')

        def write_current_trace(job, state, command, cwd, env):
            (cwd / 'trace.tsv').write_text('current execution trace')
            return 1

        with patch.object(self.p, '_signature', return_value=('fixture', dict(verified=False))), \
             patch.object(self.p, '_meta_problems', return_value=[]), \
             patch.object(self.p, '_run_process', side_effect=write_current_trace):
            self.p._metaloc(job, self.state(job), self.p.read_json(job / 'records.json'))
        previous = list(output.glob('previous-trace-*.tsv'))
        self.assertEqual(len(previous), 1)
        self.assertEqual(previous[0].read_text(), 'previous execution trace')
        self.assertEqual((output / 'trace.tsv').read_text(), 'current execution trace')

    def test_metaloc_can_stage_auxiliary_scripts_from_independent_task_directory(self):
        config=copy.deepcopy(self.cfg)
        config['core_engine']='metaloc'
        job=self.job(config)
        state=self.state(job)
        records=self.p.read_json(job/'records.json')
        with patch.object(self.p,'_signature',return_value=('test',dict(verified=False))), \
             patch.object(self.p,'_meta_problems',return_value=[]), \
             patch.object(self.p,'_run_process',return_value=1):
            self.p._metaloc(job,state,records)
        self.assertTrue((job/'raw/core/bin/merge_results.py').is_file())

    def test_cancel_terminates_external_process_and_keeps_recoverable_state(self):
        self.cfg['tools']['signalp']['command'].extend(['--delay', '30'])
        job = self.job()
        worker = subprocess.Popen([sys.executable, str(ROOT / 'pipeline.py'), 'worker', '--runs-dir', str(self.runs)])
        try:
            deadline = time.monotonic() + 10
            while time.monotonic() < deadline:
                if self.state(job).get('process'):
                    break
                time.sleep(.1)
            self.assertTrue(self.state(job).get('process'))
            self.p.cancel_job(job)
            worker.wait(timeout=10)
            self.assertEqual(self.state(job)['status'], 'cancelled')
            self.assertFalse(self.state(job).get('process'))
        finally:
            if worker.poll() is None:
                worker.terminate()
                worker.wait(timeout=5)

    def test_two_workers_do_not_execute_same_jobs_twice(self):
        one, two = self.job(), self.job()
        workers = [subprocess.Popen([sys.executable, str(ROOT / 'pipeline.py'), 'worker', '--runs-dir', str(self.runs)]) for _ in range(2)]
        try:
            for worker in workers:
                worker.wait(timeout=30)
        finally:
            for worker in workers:
                if worker.poll() is None:
                    worker.terminate()
                    worker.wait(timeout=5)
        self.assertEqual([self.state(x)['status'] for x in (one, two)], ['completed', 'completed'])
        self.assertEqual(len(list((one / 'raw/signalp').iterdir())), 1)
        self.assertEqual(len(list((two / 'raw/signalp').iterdir())), 1)

    def stubborn_job(self):
        self.cfg['tools']['signalp']['command'].append('--spawn-stubborn')
        job=self.job()
        worker=subprocess.Popen([sys.executable,str(ROOT/'pipeline.py'),'worker','--runs-dir',str(self.runs)])
        ready=job/'raw/signalp/attempt-001/stubborn.ready'
        deadline=time.monotonic()+30
        while not ready.exists() and time.monotonic()<deadline:
            time.sleep(.1)
        self.assertTrue(ready.exists())
        child=int(ready.read_text())
        child_identity=self.p._process_identity(child)
        info=self.state(job)['process']
        def cleanup():
            try:
                os.killpg(info['pid'],signal.SIGKILL)
            except ProcessLookupError:
                pass
            if self.p._process_identity(child)==child_identity:
                try:
                    os.kill(child,signal.SIGKILL)
                except ProcessLookupError:
                    pass
            if worker.poll() is None:
                worker.kill()
            worker.wait(timeout=5)
        self.addCleanup(cleanup)
        return job,worker,child

    def test_cancellation_kills_child_even_after_wrapper_exits_on_term(self):
        job,worker,child=self.stubborn_job()
        self.p.cancel_job(job)
        worker.wait(timeout=15)
        self.assertIsNone(self.p._process_identity(child),'工具的子进程仍在运行')
        self.assertEqual(self.state(job)['status'],'cancelled')

    def test_cancellation_follows_owned_child_into_new_session(self):
        self.cfg['tools']['signalp']['command'].append('--escape-session')
        job,worker,child=self.stubborn_job()
        self.p.cancel_job(job)
        worker.wait(timeout=15)
        self.assertIsNone(self.p._process_identity(child),'新 session 内的模型进程没有停止')

    def test_resume_stops_orphan_group_before_requeue(self):
        job,worker,child=self.stubborn_job()
        worker.kill()
        worker.wait(timeout=5)
        self.p.resume_job(job,self.cfg)
        self.assertIsNone(self.p._process_identity(child),'续跑遗留旧工具进程')
        self.assertEqual(self.state(job)['status'],'queued')

    def test_wait_detects_dead_worker_instead_of_hanging(self):
        self.assertTrue(hasattr(self.p,'wait_job'),'缺少能恢复失联 worker 的等待入口')
        job,worker,child=self.stubborn_job()
        worker.kill()
        worker.wait(timeout=5)
        state=self.p.wait_job(job)
        self.assertEqual(state['status'],'interrupted')
        self.assertIsNone(self.p._process_identity(child))

    def test_new_task_directory_recovers_previous_global_orphan(self):
        job,worker,child=self.stubborn_job()
        worker.kill()
        worker.wait(timeout=5)
        config=copy.deepcopy(self.cfg)
        config['tools']['signalp']['command'].remove('--spawn-stubborn')
        other=Path(self.temp.name)/'new-root'
        next_job=self.p.create_job(self.text,'fungi',config,other)
        self.p.run_worker(other)
        self.assertIsNone(self.p._process_identity(child),'跨目录启动留下了旧模型计算')
        self.assertEqual(self.state(job)['status'],'interrupted')
        self.assertEqual(self.state(next_job)['status'],'completed')

    def test_different_task_directories_share_execution_limit(self):
        self.cfg['tools']['signalp']['command'].extend(['--delay','10'])
        one=self.job()
        other=Path(self.temp.name)/'other-runs'
        two=self.p.create_job(self.text,'fungi',self.cfg,other)
        workers=[]
        try:
            workers.append(subprocess.Popen([sys.executable,str(ROOT/'pipeline.py'),'worker','--runs-dir',str(self.runs)]))
            deadline=time.monotonic()+20
            while not self.state(one).get('process') and time.monotonic()<deadline:
                time.sleep(.1)
            self.assertTrue(self.state(one).get('process'))
            workers.append(subprocess.Popen([sys.executable,str(ROOT/'pipeline.py'),'worker','--runs-dir',str(other)]))
            deadline=time.monotonic()+3
            while time.monotonic()<deadline:
                self.assertNotEqual(self.state(two)['status'],'running','两个目录同时执行了模型')
                time.sleep(.1)
        finally:
            self.p.cancel_job(one)
            self.p.cancel_job(two)
            for worker in workers:
                try:
                    worker.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    worker.kill()
                    worker.wait(timeout=5)

    def test_model_bytes_changed_with_same_size_and_mtime_invalidate_cache(self):
        model=Path(self.temp.name)/'model.bin'
        model.write_text('Nucleus')
        stat=model.stat()
        self.cfg['tools']['deeploc']['command'].extend(['--model-file',str(model)])
        self.cfg['tools']['deeploc']['provenance_paths'].append(str(model))
        job=self.job()
        self.p.run_worker(self.runs)
        self.assertEqual(self.p.read_json(job/'results.json')[0]['predictions']['deeploc']['localizations'],['Nucleus'])
        model.write_text('Plastid')
        os.utime(model,ns=(stat.st_atime_ns,stat.st_mtime_ns))
        self.p.resume_job(job,self.cfg)
        self.p.run_worker(self.runs)
        self.assertEqual(self.p.read_json(job/'results.json')[0]['predictions']['deeploc']['localizations'],['Plastid'])

    def test_nonexecuted_flag_names_on_path_do_not_change_model_cache_identity(self):
        folder=Path(self.temp.name)/'bin'
        folder.mkdir()
        unrelated=folder/'--output'
        unrelated.write_text('first unrelated executable')
        unrelated.chmod(0o755)
        job=self.job()
        state=self.state(job)
        with patch.dict(os.environ,{'PATH':str(folder)+os.pathsep+os.environ.get('PATH','')}):
            first,_=self.p._signature(self.cfg,'signalp',state)
            unrelated.write_text('changed unrelated executable')
            second,_=self.p._signature(self.cfg,'signalp',state)
        self.assertEqual(first,second,'未执行的参数不应当成 PATH 中的依赖程序')

    def test_invalid_version_probe_disables_cache_without_aborting_predictions(self):
        self.cfg['tools']['signalp']['version_command']=[sys.executable,'-c','print({unknown_placeholder})']
        provenance=self.p._provenance(self.cfg,'signalp')
        self.assertFalse(provenance['verified'])
        self.assertIn('version_error',provenance)

    def test_provenance_scans_canonical_model_directory(self):
        native=Path(self.temp.name)/'native-models'
        native.mkdir()
        (native/'weights.bin').write_bytes(b'actual model')
        alias=Path(self.temp.name)/'mounted-alias'
        alias.symlink_to(native,target_is_directory=True)
        self.cfg['tools']['signalp']['provenance_paths']=[str(alias)]
        provenance=self.p._provenance(self.cfg,'signalp')
        self.assertIn(str(native/'weights.bin'),provenance['files'])

    def test_executable_resolution_preserves_the_virtual_environment(self):
        native=Path(self.temp.name)/'native-env'
        (native/'bin').mkdir(parents=True)
        (native/'bin/python').symlink_to(sys.executable)
        alias=Path(self.temp.name)/'mounted-env'
        alias.symlink_to(native,target_is_directory=True)
        command=self.p._expand([str(alias/'bin/python'),'-c','print(1)'],{})
        self.assertEqual(command[0],str(native/'bin/python'))
        self.assertNotEqual(command[0],str(Path(sys.executable).resolve()))


if __name__ == '__main__':
    unittest.main()
