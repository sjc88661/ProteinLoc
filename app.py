"""简洁 Linux/Tkinter 界面；计算由 pipeline.py 的独立进程执行。"""
import argparse
import json
import os
import queue
import shutil
import subprocess
import threading
import time
import tkinter as tk
from datetime import datetime
from pathlib import Path
from tkinter import filedialog, messagebox, scrolledtext, ttk

from parsers import MODULES, atomic_json, export_results, summary_row
from pipeline import (ROOT, TERMINAL, PROKARYOTES, _recover, cancel_job, check_environment, create_job,
                      ensure_worker, list_jobs, load_config, local_path, read_json, resume_job)

if (ROOT / '.runtime/fonts.conf').exists():
    os.environ.setdefault('FONTCONFIG_FILE', str(ROOT / '.runtime/fonts.conf'))

GROUPS = {'真菌 / 酵母': 'fungi', '植物 / 含质体真核': 'plant', '动物': 'animal',
          '其他真核': 'other_euk', '革兰氏阴性菌':'gram_negative', '革兰氏阳性菌':'gram_positive',
          '古菌':'archaea', '原核（分组未知）':'prokaryote', '类群未知': 'unknown'}
STATUS = {'queued':'排队', 'running':'分析中', 'completed':'全部完成', 'partial':'部分结果',
          'failed':'失败', 'cancelled':'已取消', 'interrupted':'中断待续跑', 'invalid':'输入无效'}
VISIBLE = ('Protein ID','Length','Status','Localization','Localization model','SP','TM','DeepTMHMM','TargetP','GPI','ER signal','PTS1 SVM','PTS candidate')


class ProteinLocApp:
    def __init__(self, root, config_path=None, runs_dir=None):
        self.root = root
        self.config_path = Path(config_path) if config_path else ROOT / 'config.json'
        config = load_config(self.config_path)
        self.remember_session = runs_dir is None
        saved = {}
        if self.remember_session:
            try:
                saved = read_json(ROOT/'.runtime/gui-state.json')
            except (OSError,ValueError):
                pass
        if runs_dir is None and saved.get('last_runs_dir'):
            runs_dir = saved['last_runs_dir']
        self.runs_dir = Path(runs_dir) if runs_dir else local_path(config.get('runs_dir', 'runs'))
        self.known_runs_dirs = {self.runs_dir.resolve()}
        self.known_runs_dirs.update(Path(p) for p in saved.get('known_runs_dirs',[]) if Path(p).is_dir())
        self.records, self.filtered, self.jobs = [], [], []
        self.rows, self.selected_job = {}, None
        self.page, self.page_size, self.loaded_stamp = 0, 200, None
        self.messages = queue.Queue()
        self.loading, self.sort_key, self.sort_reverse = False, None, False
        self.root.title('ProteinLoc · 蛋白质亚细胞定位分析')
        self.root.geometry('1180x840')
        self.root.minsize(940, 700)
        self.root.protocol('WM_DELETE_WINDOW', self.close)
        self.build_gui()
        self.render_page()
        self.poll()

    def build_gui(self):
        style = ttk.Style(self.root)
        style.theme_use('clam')
        style.configure('TButton', padding=(9, 5))
        style.configure('Treeview', rowheight=26)
        style.configure('Horizontal.TProgressbar', background='#5b8da8')
        style.configure('Treeview.Heading', font=('sans-serif', 10, 'bold'))
        base = ttk.Frame(self.root, padding=12)
        base.pack(fill='both', expand=True)
        ttk.Label(base, text='ProteinLoc', font=('sans-serif', 17, 'bold')).pack(anchor='w')
        ttk.Label(base, text='本地批量分析 · 定位预测与信号证据', foreground='#526578').pack(anchor='w', pady=(0,8))
        input_frame = ttk.LabelFrame(base, text='输入蛋白 FASTA', padding=8)
        input_frame.pack(fill='x')
        self.input_text = tk.Text(input_frame, height=4, wrap='none', font=('monospace',10))
        self.input_text.pack(side='left', fill='x', expand=True)
        buttons = ttk.Frame(input_frame)
        buttons.pack(side='right', padx=(8,0))
        ttk.Button(buttons, text='导入 FASTA', command=self.import_fasta).pack(fill='x')
        ttk.Button(buttons, text='清空输入', command=lambda: self.input_text.delete('1.0','end')).pack(fill='x', pady=(5,0))
        settings = ttk.Frame(base)
        settings.pack(fill='x', pady=8)
        ttk.Label(settings, text='类群').pack(side='left')
        self.group = tk.StringVar(value='真菌 / 酵母')
        ttk.Combobox(settings, textvariable=self.group, values=list(GROUPS), width=19, state='readonly').pack(side='left', padx=5)
        ttk.Label(settings, text='物种名称（可选）').pack(side='left', padx=(6,0))
        self.species = tk.StringVar()
        ttk.Entry(settings, textvariable=self.species, width=22).pack(side='left', padx=5)
        ttk.Button(settings, text='环境检查', command=self.doctor).pack(side='right')
        ttk.Button(settings, text='高级设置', command=self.settings).pack(side='right', padx=5)
        output = ttk.Frame(base)
        output.pack(fill='x', pady=(0,8))
        ttk.Label(output, text='任务目录').pack(side='left')
        self.output = tk.StringVar(value=str(self.runs_dir))
        ttk.Entry(output, textvariable=self.output).pack(side='left', fill='x', expand=True, padx=5)
        ttk.Button(output, text='选择目录', command=self.choose_output).pack(side='left')
        controls = ttk.Frame(base)
        controls.pack(fill='x')
        ttk.Button(controls, text='开始分析 / 加入队列', command=self.submit).pack(side='left')
        ttk.Button(controls, text='取消所选任务', command=self.cancel).pack(side='left', padx=6)
        ttk.Button(controls, text='续跑 / 重试', command=self.resume).pack(side='left')
        ttk.Button(controls, text='打开结果目录', command=self.open_folder).pack(side='right')
        history = ttk.Frame(base)
        history.pack(fill='x', pady=8)
        ttk.Label(history, text='任务记录').pack(side='left')
        self.history = ttk.Combobox(history, state='readonly')
        self.history.pack(side='left', fill='x', expand=True, padx=5)
        self.history.bind('<<ComboboxSelected>>', self.select_history)
        self.progress = ttk.Progressbar(base, maximum=len(MODULES))
        self.progress.pack(fill='x')
        self.status = tk.StringVar(value='准备就绪')
        ttk.Label(base, textvariable=self.status, wraplength=1100).pack(fill='x', pady=5)
        filters = ttk.Frame(base)
        filters.pack(fill='x', pady=(3,5))
        ttk.Label(filters, text='筛选').pack(side='left')
        self.search = tk.StringVar()
        search = ttk.Entry(filters, textvariable=self.search, width=25)
        search.pack(side='left', fill='x', expand=True, padx=5)
        search.bind('<Return>', lambda e: self.apply_filter())
        self.state_filter = tk.StringVar(value='全部状态')
        self.location_filter = tk.StringVar(value='全部定位')
        self.sp_filter = tk.StringVar(value='SP 全部')
        self.gpi_filter = tk.StringVar(value='GPI 全部')
        for variable, values, width in (
            (self.state_filter, ['全部状态','全部完成','部分结果','输入无效'], 11),
            (self.location_filter, ['全部定位','Cytoplasm','Nucleus','Extracellular','Mitochondrion','Cell membrane','Endoplasmic reticulum','Golgi apparatus','Lysosome/Vacuole','Peroxisome','Chloroplast','Plastid','Cytoplasmic','Cytoplasmic Membrane','Outer Membrane','Periplasmic','Cell wall & surface'], 18),
            (self.sp_filter, ['SP 全部','SP Yes','SP No'], 10),
            (self.gpi_filter, ['GPI 全部','GPI Yes','GPI No'], 11)):
            box = ttk.Combobox(filters, textvariable=variable, values=values, width=width, state='readonly')
            box.pack(side='left', padx=2)
            box.bind('<<ComboboxSelected>>', lambda e: self.apply_filter())
        ttk.Button(filters, text='应用', command=self.apply_filter).pack(side='left', padx=4)
        table = ttk.Frame(base)
        table.pack(fill='both', expand=True)
        self.tree = ttk.Treeview(table, columns=VISIBLE, show='headings', height=8)
        for key in VISIBLE:
            self.tree.heading(key, text={'Localization model':'Model','PTS candidate':'PTS motifs'}.get(key,key), command=lambda k=key: self.sort_by(k))
            self.tree.column(key, width=140 if key in ('Protein ID','Localization') else 120 if key=='Localization model' else 110 if key=='DeepTMHMM' else 80, minwidth=55,
                             stretch=key in ('ER signal','PTS candidate','Localization'))
        scroll_y, scroll_x = ttk.Scrollbar(table, orient='vertical', command=self.tree.yview), ttk.Scrollbar(table, orient='horizontal', command=self.tree.xview)
        self.tree.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
        self.tree.grid(row=0,column=0,sticky='nsew')
        scroll_y.grid(row=0,column=1,sticky='ns')
        scroll_x.grid(row=1,column=0,sticky='ew')
        table.rowconfigure(0,weight=1)
        table.columnconfigure(0,weight=1)
        self.tree.bind('<<TreeviewSelect>>', self.show_detail)
        paging = ttk.Frame(base)
        paging.pack(fill='x', pady=5)
        self.count_label = ttk.Label(paging, text='尚无结果')
        self.count_label.pack(side='left')
        ttk.Button(paging, text='上一页', command=lambda: self.turn_page(-1)).pack(side='left', padx=(12,4))
        ttk.Button(paging, text='下一页', command=lambda: self.turn_page(1)).pack(side='left')
        ttk.Button(paging, text='导出筛选结果', command=self.export_dialog).pack(side='right')
        notebook = ttk.Notebook(base)
        notebook.pack(fill='both', pady=(4,0))
        evidence = ttk.Frame(notebook)
        sequence = ttk.Frame(notebook)
        logs = ttk.Frame(notebook)
        notebook.add(evidence, text='蛋白证据与分数')
        notebook.add(sequence, text='序列与特征')
        notebook.add(logs, text='任务日志')
        self.detail_text = scrolledtext.ScrolledText(evidence, height=8, wrap='word', font=('monospace',10))
        self.detail_text.pack(fill='both', expand=True)
        self.feature_canvas = tk.Canvas(sequence, height=100, bg='white', highlightthickness=0)
        self.feature_canvas.pack(fill='x')
        self.sequence_text = scrolledtext.ScrolledText(sequence, height=4, font=('monospace',10))
        self.sequence_text.pack(fill='both', expand=True)
        self.log_text = scrolledtext.ScrolledText(logs, height=8, wrap='word', font=('monospace',9))
        self.log_text.pack(fill='both', expand=True)

    @staticmethod
    def set_text(widget, content):
        widget.configure(state='normal')
        widget.delete('1.0','end')
        widget.insert('1.0',content)
        widget.configure(state='disabled')

    def import_fasta(self):
        path = filedialog.askopenfilename(filetypes=[('FASTA','*.fasta *.faa *.fa *.txt'),('全部文件','*')])
        if path:
            try:
                text = Path(path).read_text(encoding='utf-8-sig')
                self.input_text.delete('1.0','end')
                self.input_text.insert('1.0',text)
            except (OSError, UnicodeError) as error:
                messagebox.showerror('读取失败', str(error))

    def choose_output(self):
        path = filedialog.askdirectory(initialdir=self.runs_dir)
        if path:
            self.output.set(path)
            self.runs_dir = Path(path)
            self.known_runs_dirs.add(self.runs_dir.resolve())
            self.save_session()
            self.load_job(None)

    def submit(self):
        try:
            self.runs_dir = Path(self.output.get()).expanduser().resolve()
            self.known_runs_dirs.add(self.runs_dir)
            self.save_session()
            self.runs_dir.mkdir(parents=True, exist_ok=True)
            job = create_job(self.input_text.get('1.0','end').strip(), GROUPS[self.group.get()],
                             load_config(self.config_path), self.runs_dir, self.species.get())
            ensure_worker(self.runs_dir)
            self.load_job(job)
            self.status.set('任务已加入队列：' + job.name)
        except (OSError, ValueError, KeyError) as error:
            messagebox.showerror('无法开始分析', str(error))

    def save_session(self):
        if self.remember_session:
            atomic_json(ROOT/'.runtime/gui-state.json',dict(last_runs_dir=str(self.runs_dir),
                known_runs_dirs=sorted(str(p) for p in self.known_runs_dirs)))

    def load_job(self, job):
        self.selected_job = Path(job) if job is not None else None
        self.records, self.filtered, self.rows = [], [], {}
        self.loaded_stamp = None
        self.page = 0
        self.history.set('')
        self.progress['value'] = 0
        self.status.set('正在读取任务…' if self.selected_job else '准备就绪')
        self.set_text(self.log_text, '')
        self.render_page()
        self.load_results()

    def select_history(self, event=None):
        index = self.history.current()
        if 0 <= index < len(self.jobs):
            self.load_job(self.jobs[index]['path'])

    def load_results(self):
        if not self.selected_job or self.loading:
            return
        path = self.selected_job / 'results.json'
        if not path.is_file():
            return
        stamp = (str(path), path.stat().st_mtime_ns)
        if stamp == self.loaded_stamp:
            return
        self.loading = True
        messages = self.messages
        def load():
            try:
                messages.put(('results', stamp, read_json(path)))
            except (OSError, ValueError) as error:
                messages.put(('error', str(error)))
            finally:
                messages.put(('loaded',))
        threading.Thread(target=load, daemon=True).start()

    def poll(self):
        try:
            while True:
                message = self.messages.get_nowait()
                if message[0] == 'results':
                    _, stamp, records = message
                    if self.selected_job and stamp[0] == str(self.selected_job / 'results.json'):
                        self.records, self.loaded_stamp = records, stamp
                        self.rows = {r['internal_id']: summary_row(r) for r in records}
                        self.apply_filter()
                elif message[0] == 'loaded':
                    self.loading = False
                elif message[0] == 'environment':
                    self.show_environment(message[1])
                elif message[0] == 'exported':
                    self.status.set('导出完成：' + str(message[1]))
                elif message[0] == 'error':
                    self.status.set(message[1])
        except queue.Empty:
            pass
        try:
            self.jobs = list_jobs(self.runs_dir)
            for state in self.jobs:
                if state['status'] == 'running':
                    _recover(Path(state['path']), state)
            self.history['values'] = [f"{s['id']} · {s['species_name'] or s['organism']} · {s['count']} 条 · {STATUS.get(s['status'],s['status'])}" for s in self.jobs]
            for i, state in enumerate(self.jobs):
                if self.selected_job and Path(state['path']) == self.selected_job:
                    self.history.current(i)
                    self.status.set(f"{STATUS.get(state['status'],state['status'])} · {state.get('message','')} · {state.get('completed_proteins',0)}/{state['count']} 条完整")
                    self.progress['value'] = len(state.get('modules',{}))
                    self.progress['maximum'] = len(state.get('modules',{})) if state['status'] in TERMINAL else len(MODULES)+(state['organism'] in PROKARYOTES)
                    self.load_results()
                    log = self.selected_job / 'run.log'
                    if log.exists():
                        with log.open('rb') as stream:
                            stream.seek(max(0, log.stat().st_size-12000))
                            self.set_text(self.log_text, stream.read().decode('utf-8', errors='replace'))
                    break
            if self.selected_job is None and self.jobs:
                self.load_job(self.jobs[0]['path'])
            if any(s['status'] == 'queued' for s in self.jobs):
                ensure_worker(self.runs_dir)
        except (OSError, ValueError, KeyError) as error:
            self.status.set('读取任务状态失败：' + str(error))
        self.root.after(700, self.poll)

    def apply_filter(self):
        query = self.search.get().strip().casefold()
        self.filtered = []
        for record in self.records:
            row = self.rows.get(record['internal_id']) or summary_row(record)
            if query and query not in ' '.join(str(x) for x in row.values()).casefold():
                continue
            if self.state_filter.get() != '全部状态' and STATUS.get(record['status']) != self.state_filter.get():
                continue
            if self.location_filter.get() != '全部定位' and self.location_filter.get() not in str(row['Localization']).split('; '):
                continue
            if self.sp_filter.get() != 'SP 全部' and row['SP'] != self.sp_filter.get().split()[-1]:
                continue
            if self.gpi_filter.get() != 'GPI 全部' and row['GPI'] != self.gpi_filter.get().split()[-1]:
                continue
            self.filtered.append(record)
        if self.sort_key:
            def sort_value(record):
                value = self.rows[record['internal_id']].get(self.sort_key, '')
                if isinstance(value, (int, float)):
                    return (0, value)
                return (1, str(value).casefold())
            self.filtered.sort(key=sort_value, reverse=self.sort_reverse)
        self.page = 0
        self.render_page()

    def sort_by(self, key):
        self.sort_reverse = self.sort_key == key and not self.sort_reverse
        self.sort_key = key
        self.apply_filter()

    def turn_page(self, delta):
        pages = max(1, (len(self.filtered)+self.page_size-1)//self.page_size)
        self.page = max(0, min(pages-1, self.page+delta))
        self.render_page()

    def render_page(self):
        previous_selection = self.tree.selection()
        self.tree.delete(*self.tree.get_children())
        for record in self.filtered[self.page*self.page_size:(self.page+1)*self.page_size]:
            row = self.rows[record['internal_id']]
            self.tree.insert('', 'end', iid=record['internal_id'], values=[row[k] for k in VISIBLE])
        visible = self.tree.get_children()
        if visible:
            selected = next((key for key in previous_selection if key in visible), visible[0])
            self.tree.selection_set(selected)
            self.tree.focus(selected)
            self.tree.see(selected)
        pages = max(1, (len(self.filtered)+self.page_size-1)//self.page_size)
        self.count_label.configure(text=f'筛选 {len(self.filtered)} / 共 {len(self.records)} 条 · 第 {self.page+1}/{pages} 页')
        self.show_detail()

    def show_detail(self, event=None):
        selected = self.tree.selection()
        if not selected:
            message = '请在上方表格中选择一条蛋白。' if self.filtered else '当前没有可显示的蛋白结果。'
            for widget in (self.detail_text, self.sequence_text):
                self.set_text(widget, message)
            self.feature_canvas.delete('all')
            return
        record = next((r for r in self.records if r['internal_id'] == selected[0]), None)
        if record is None:
            return
        self.set_text(self.detail_text, json.dumps(record, ensure_ascii=False, indent=2))
        sequence = record['sequence']
        self.set_text(self.sequence_text, '\n'.join(f'{i+1:>7}  {sequence[i:i+80]}' for i in range(0,len(sequence),80)))
        canvas = self.feature_canvas
        canvas.configure(height=135)
        canvas.delete('all')
        width, length = max(600,canvas.winfo_width()), max(1,len(sequence))
        left, right = 90, width-30
        scale = (right-left)/length
        canvas.create_line(left,25,right,25,fill='#59748a',width=2)
        canvas.create_text(left,10,text='1',anchor='w')
        canvas.create_text(right,10,text=str(len(sequence)),anchor='e')
        predictions = record.get('predictions',{})
        features = []
        sp = predictions.get('signalp',{})
        if sp.get('status') == 'completed' and sp.get('signal_peptide') and sp.get('cleavage_after'):
            features.append(('SP',1,sp['cleavage_after'],'#59aa7c',40))
        for start,end in predictions.get('phobius',{}).get('tm_regions',[]):
            features.append(('TM',start,end,'#df9a45',55))
        for key,label,color,y in [('tm_regions','DTM α','#ce7e1a',70),('beta_regions','DTM β','#b85c79',85)]:
            for start,end in predictions.get('deeptmhmm',{}).get(key,[]):
                features.append((label,start,end,color,y))
        omega = predictions.get('predgpi',{}).get('omega_site')
        if omega:
            features.append(('GPI',omega,omega,'#a26ab8',100))
        for name in ('er','pts'):
            for item in predictions.get(name,{}).get('matches',[]):
                features.append((name.upper(),item['start'],item['end'],'#568bd1',115))
        if not features:
            canvas.create_text(left, 65, text='当前结果暂无可绘制的特征区段。', anchor='w', fill='#526578')
        for label,start,end,color,y in features:
            x1,x2=left+(start-1)*scale,left+end*scale
            canvas.create_rectangle(x1,y-5,max(x1+3,x2),y+4,fill=color,outline='')
            canvas.create_text(45,y,text=label,fill=color)

    def cancel(self):
        if self.selected_job:
            cancel_job(self.selected_job)

    def resume(self):
        if self.selected_job:
            try:
                resume_job(self.selected_job,load_config(self.config_path))
                ensure_worker(self.runs_dir)
                self.loaded_stamp=None
            except (OSError,ValueError) as error:
                messagebox.showerror('无法续跑',str(error))

    def export_to(self, directory):
        folder = Path(directory) / ('export-' + datetime.now().strftime('%Y%m%d-%H%M%S-%f'))
        records = list(self.filtered)
        messages = self.messages
        def export():
            try:
                export_results(records, folder)
                messages.put(('exported',folder))
            except Exception as error:
                messages.put(('error','导出失败：'+str(error)))
        threading.Thread(target=export,daemon=False).start()
        return folder

    def export_dialog(self):
        if not self.filtered:
            messagebox.showinfo('导出','当前没有筛选结果')
            return
        folder=filedialog.askdirectory(title='选择导出位置，将自动建立独立子目录')
        if folder:
            self.status.set('正在导出筛选结果…')
            self.export_to(folder)

    def open_folder(self):
        folder=self.selected_job or self.runs_dir
        try:
            if shutil.which('explorer.exe'):
                path=subprocess.check_output(['wslpath','-w',str(folder)],text=True).strip()
                subprocess.Popen(['explorer.exe',path],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            else:
                subprocess.Popen(['xdg-open',str(folder)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        except (OSError,subprocess.SubprocessError) as error:
            messagebox.showinfo('结果目录',str(folder)+'\n'+str(error))

    def doctor(self):
        messages, config_path = self.messages, self.config_path
        def check():
            try:
                messages.put(('environment',check_environment(load_config(config_path))))
            except Exception as error:
                messages.put(('error',str(error)))
        self.status.set('正在检查本地工具和模型文件…')
        threading.Thread(target=check,daemon=True).start()

    def show_environment(self,report):
        window=tk.Toplevel(self.root)
        window.title('环境检查')
        text=scrolledtext.ScrolledText(window,width=100,height=24,wrap='word')
        text.pack(fill='both',expand=True,padx=10,pady=10)
        lines=['检查程序、安装包和本地规则；真实预测是否成功由任务验收确认。\n']
        for row in report:
            lines.append(f"{row['module']} · {row['version']} · {'文件/命令已就绪' if row['ready'] else '需要配置'}")
            lines.extend('  '+issue for issue in row['issues'])
        lines.extend(['','授权安装包及配置方法：'+str(ROOT/'README.md'),'缺少工具时任务保留明确失败状态及其他可用结果。'])
        self.set_text(text,'\n'.join(lines))

    def settings(self):
        window=tk.Toplevel(self.root)
        window.title('高级设置 · 本地工具路径和参数')
        text=scrolledtext.ScrolledText(window,width=105,height=28,font=('monospace',10))
        text.pack(fill='both',expand=True,padx=10,pady=10)
        text.insert('1.0',json.dumps(load_config(self.config_path),ensure_ascii=False,indent=2))
        def save():
            try:
                config=json.loads(text.get('1.0','end'))
                temp=self.config_path.with_suffix('.draft.json')
                atomic_json(temp,config)
                try:
                    load_config(temp)
                finally:
                    temp.unlink(missing_ok=True)
                atomic_json(self.config_path,config)
                window.destroy()
            except (OSError,ValueError,KeyError,TypeError) as error:
                messagebox.showerror('设置无效',str(error),parent=window)
        ttk.Button(window,text='保存设置（续跑时重新检查缓存）',command=save).pack(pady=(0,10))

    def close(self):
        active=[s for folder in self.known_runs_dirs for s in list_jobs(folder) if s['status'] in ('queued','running')]
        if not active:
            self.root.destroy()
            return
        window=tk.Toplevel(self.root)
        window.title('分析任务仍在运行')
        window.transient(self.root)
        window.grab_set()
        ttk.Label(window,text=f'当前有 {len(active)} 个运行或排队任务。\n后台继续需要保持 WSL 运行。',padding=18).pack()
        def leave(cancel=False):
            if cancel:
                for state in active:
                    cancel_job(state['path'])
            else:
                for folder in self.known_runs_dirs:
                    ensure_worker(folder)
            self.root.destroy()
        buttons=ttk.Frame(window,padding=10)
        buttons.pack()
        ttk.Button(buttons,text='后台继续并关闭窗口',command=leave).pack(side='left',padx=5)
        ttk.Button(buttons,text='取消任务并关闭',command=lambda:leave(True)).pack(side='left',padx=5)
        ttk.Button(buttons,text='返回',command=window.destroy).pack(side='left',padx=5)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--config')
    parser.add_argument('--runs-dir')
    args=parser.parse_args()
    root=tk.Tk()
    ProteinLocApp(root,args.config,args.runs_dir)
    root.mainloop()


if __name__=='__main__':
    main()
