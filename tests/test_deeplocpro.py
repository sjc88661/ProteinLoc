"""Contract checks for the real DeepLocPro CSV and prokaryote routing."""
import csv
import os
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import parsers
import pipeline

HEADER=',ACC,Localization,Cell wall & surface,Extracellular,Cytoplasmic,Cytoplasmic Membrane,Outer Membrane,Periplasmic\n'


class DeepLocProTests(unittest.TestCase):
    def test_actual_old_config_migrates_builtin_organism_modes_for_prokaryotes(self):
        root=Path(__file__).resolve().parents[1]
        old=pipeline.read_json(root/'tests/fixtures/legacy_config.json')
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'old.json'
            path.write_text(json.dumps(old))
            migrated=pipeline.load_config(path)
        context=pipeline._context(migrated,organism='gram_positive')
        context.update(work='work',override='runtime.config')
        core=pipeline._expand(migrated['metaloc']['command'],context)
        signalp=pipeline._expand(migrated['tools']['signalp']['command'],context)
        self.assertEqual(core[core.index('--organism')+1],'other')
        self.assertEqual(signalp[signalp.index('--organism')+1],'other')

    def test_resume_without_replacement_config_upgrades_persisted_old_task(self):
        root=Path(__file__).resolve().parents[1]
        old=pipeline.read_json(root/'tests/fixtures/legacy_config.json')
        with tempfile.TemporaryDirectory() as temp:
            job=pipeline.create_job('>old\nM'+'A'*49+'\n','fungi',pipeline.load_config(),Path(temp))
            state=pipeline.read_json(job/'job.json')
            state.update(config=old,status='interrupted')
            parsers.atomic_json(job/'job.json',state)
            pipeline.resume_job(job)
            migrated=pipeline.read_json(job/'job.json')['config']
        self.assertIn('pts1_svm',migrated['tools'])
        self.assertIn('deeplocpro',migrated['tools'])
        context=pipeline._context(migrated,organism='gram_negative')
        command=pipeline._expand(migrated['tools']['signalp']['command'],context)
        self.assertEqual(command[command.index('--organism')+1],'other')

    def test_custom_program_organism_arguments_are_not_rewritten(self):
        config=pipeline.load_config()
        config['metaloc']['command']=['/custom/runner','--organism','euk']
        config['tools']['signalp']['command']=['/custom/secretory-model','--organism','eukarya']
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'custom.json'
            path.write_text(json.dumps(config))
            loaded=pipeline.load_config(path)
        self.assertEqual(loaded['metaloc']['command'],['/custom/runner','--organism','euk'])
        self.assertEqual(loaded['tools']['signalp']['command'],['/custom/secretory-model','--organism','eukarya'])

    def test_existing_eukaryotic_config_gets_only_the_new_optional_tool_default(self):
        config=pipeline.load_config()
        config['tools'].pop('deeplocpro')
        config['tools']['deeploc']['command']=['/custom/deeploc','--fasta','{input}']
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'existing-config.json'
            path.write_text(json.dumps(config))
            loaded=pipeline.load_config(path)
        self.assertIn('deeplocpro',loaded['tools'])
        self.assertEqual(loaded['tools']['deeploc']['command'],['/custom/deeploc','--fasta','{input}'])

    def test_cpu_entry_hides_cuda_before_upstream_import_and_receives_metaloc_group(self):
        with tempfile.TemporaryDirectory() as temp:
            home=Path(temp)
            package=home/'.runtime/deeplocpro-source/DeepLocPro'
            package.mkdir(parents=True)
            (package/'__init__.py').write_text('')
            (package/'deeplocpro.py').write_text("import json,os\nfrom pathlib import Path\nvisible_at_import=os.environ.get('CUDA_VISIBLE_DEVICES')\ndef main(args):\n    (Path(args.output)/'observed.json').write_text(json.dumps(dict(device=args.device,group=args.group,visible_at_import=visible_at_import)))\n")
            output=home/'output/observed.json'
            code=f"from pathlib import Path; from scripts import deeplocpro_local as adapter; adapter.ROOT=Path({str(home)!r}); adapter.main(['-f','unused.faa','-o',{str(output.parent)!r},'-d','cpu'])"
            environment=dict(os.environ,CUDA_VISIBLE_DEVICES='0',PROTEINLOC_DEEPLOCPRO_GROUP='positive')
            process=subprocess.run([sys.executable,'-c',code],capture_output=True,text=True,env=environment)
            self.assertEqual(process.returncode,0,process.stdout+process.stderr)
            observed=json.loads(output.read_text())
            self.assertEqual(observed,dict(device='cpu',group='positive',visible_at_import='-1'))
            invocation=json.loads((output.parent/'proteinloc_invocation.json').read_text())
            self.assertEqual(invocation['group'],'positive')
            self.assertEqual(invocation['device'],'cpu')

    def parse(self, rows):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'results.csv'
            path.write_text(HEADER+rows)
            self.assertTrue(hasattr(parsers,'parse_deeplocpro'),'缺少 DeepLocPro 独立解析器')
            return parsers.parse_deeplocpro(path)

    def test_indexed_csv_preserves_actual_ids_six_scores_and_multiclass_label(self):
        result=self.parse('0,P00000002,Periplasmic,0.01,0.02,0.03,0.04,0.10,0.80\n1,P00000001,Cytoplasmic,0,0,1,0,0,0\n')
        self.assertEqual(result['P00000002']['localizations'],['Periplasmic'])
        self.assertEqual(result['P00000002']['scores']['Periplasmic'],0.80)
        self.assertEqual(len(result['P00000002']['scores']),6)
        self.assertEqual(result['P00000001']['raw']['ACC'],'P00000001')
        self.assertEqual(result['P00000001']['tool'],'DeepLocPro')

    def test_duplicate_ids_and_invalid_probabilities_are_not_successful_predictions(self):
        valid='0,P00000001,Cytoplasmic,0,0,1,0,0,0\n'
        for rows in (valid+valid, valid.replace(',1,0,0,0',',nan,0,0,0'),
                     valid.replace(',1,0,0,0',',0.2,0,0,0'),valid.replace('Cytoplasmic,','invented,')):
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                self.parse(rows)

    def test_cli_group_and_signalp_organism_are_selected_from_task_source(self):
        config=pipeline.load_config()
        for organism,group in [('gram_negative','negative'),('gram_positive','positive'),('archaea','archaea'),('prokaryote','any')]:
            with self.subTest(organism=organism):
                context=pipeline._context(config,organism=organism)
                self.assertEqual(context.get('deeplocpro_group'),group)
                self.assertEqual(context.get('metaloc_organism'),'other')
                self.assertEqual(context.get('signalp_organism'),'other')
        context=pipeline._context(config,organism='fungi')
        self.assertEqual(context.get('metaloc_organism'),'euk')
        self.assertEqual(context.get('signalp_organism'),'eukarya')

    def test_prokaryotic_localization_and_model_are_used_in_the_shared_summary(self):
        record=parsers.read_fasta('>bacterium\nMAAAAA\n')[0]
        record['predictions']={'deeploc':{'status':'unsupported','applicable':False},
            'deeplocpro':{'status':'completed','tool':'DeepLocPro','localizations':['Outer Membrane']},
            'er':{'status':'unsupported','applicable':False}}
        row=parsers.summary_row(record)
        self.assertEqual(row['Localization'],'Outer Membrane')
        self.assertEqual(row['Localization model'],'DeepLocPro')
        self.assertEqual(row['ER signal'],'不适用')


if __name__=='__main__':
    unittest.main()
