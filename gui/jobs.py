"""Persisted background jobs survive browser reruns without duplicate inference."""
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import uuid

ROOT=Path(__file__).resolve().parent.parent
JOBS=Path(os.environ.get('EVOPHI_GUI_JOBS',ROOT/'.cache/gui/jobs')).resolve()


def write(path,value):
    path=Path(path);tmp=path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)


def folder(identifier):
    if not re.fullmatch(r'[a-f0-9]{32}',str(identifier)):raise ValueError('Invalid job identifier.')
    return JOBS/identifier


def read(identifier):
    path=folder(identifier)
    state=json.loads((path/'job.json').read_text())
    pid_file=path/'worker.pid'
    if state['status'] in ['queued','running'] and pid_file.exists():
        alive=True
        try:
            pid=int(pid_file.read_text());os.kill(pid,0)
            proc=Path('/proc')/str(pid)
            if proc.exists():
                alive=(proc/'stat').read_text().split()[2]!='Z' and identifier in (proc/'cmdline').read_bytes().decode(errors='replace')
        except (ProcessLookupError,FileNotFoundError):alive=False
        if not alive:
            state.update(status='failed',error='The background calculation was interrupted. Retry to reuse cached work.')
            write(path/'job.json',state)
    progress=path/'progress.json'
    if progress.exists():
        try:state['progress']=json.loads(progress.read_text())
        except json.JSONDecodeError:pass
    return state


def create(frame,device='auto',msa=None,msa_filtered=False):
    JOBS.mkdir(parents=True,exist_ok=True)
    identifier=uuid.uuid4().hex;path=folder(identifier);path.mkdir()
    frame.to_csv(path/'variants.csv',index=False)
    if msa is not None:(path/'alignment.a3m').write_bytes(msa)
    state=dict(id=identifier,status='queued',created=time.time(),variants=len(frame),proteins=int(frame.Sequence.nunique()),
               device=device,msa=msa is not None,msa_filtered=msa_filtered,error=None)
    write(path/'job.json',state)
    with (path/'worker.log').open('a') as log:
        process=subprocess.Popen([sys.executable,str(ROOT/'gui/worker.py'),identifier],cwd=ROOT,
            stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    state['worker_pid']=process.pid
    # The worker also updates this file; persist its PID separately to avoid a race.
    (path/'worker.pid').write_text(str(process.pid))
    return identifier


def recent(limit=20):
    if not JOBS.exists():return []
    items=[]
    for path in JOBS.glob('*/job.json'):
        try:items.append(read(path.parent.name))
        except (OSError,ValueError):continue
    return sorted(items,key=lambda item:item['created'],reverse=True)[:limit]


def cancel(identifier):
    (folder(identifier)/'cancel.request').touch()


def retry(identifier):
    import pandas as pd
    path=folder(identifier);state=read(identifier)
    alignment=path/'alignment.a3m'
    return create(pd.read_csv(path/'variants.csv',keep_default_na=False),state['device'],
        alignment.read_bytes() if alignment.exists() else None,state['msa_filtered'])


def log_tail(identifier,limit=8000):
    path=folder(identifier)/'pipeline.log'
    if not path.exists():return ''
    with path.open('rb') as stream:
        stream.seek(max(0,path.stat().st_size-limit))
        return stream.read().decode(errors='replace')
