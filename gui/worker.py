"""Run one validated pipeline job, serializing access to shared model/MSA caches."""
import fcntl
import json
import os
import signal
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
from gui import jobs


def stop(process):
    try:os.killpg(process.pid,signal.SIGTERM)
    except ProcessLookupError:pass
    try:process.wait(timeout=15)
    except subprocess.TimeoutExpired:
        try:os.killpg(process.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        process.wait()


def run(identifier):
    path=jobs.folder(identifier);state=jobs.read(identifier)
    def update(status,**extra):
        state.update(status=status,updated=time.time(),**extra);jobs.write(path/'job.json',state)
    process=None
    try:
        with (jobs.JOBS.parent/'pipeline.lock').open('a') as lock:
            while True:
                if (path/'cancel.request').exists():update('cancelled');return
                try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);break
                except BlockingIOError:time.sleep(.5)
            if (path/'cancel.request').exists():update('cancelled');return
            update('running',started=time.time())
            command=[sys.executable,'-u',str(ROOT/'inference/run.py'),'--input',str(path/'variants.csv'),
                '--output',str(path/'predictions.csv'),'--device',state['device'],'--progress-file',str(path/'progress.json')]
            if state['msa']:command+=['--msa',str(path/'alignment.a3m')]
            if state['msa_filtered']:command+=['--msa-filtered']
            with (path/'pipeline.log').open('w') as log:
                process=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                while process.poll() is None:
                    if (path/'cancel.request').exists():
                        stop(process)
                        update('cancelled');return
                    time.sleep(.5)
            if process.returncode==0 and (path/'predictions.csv').exists():
                update('complete',finished=time.time())
            else:
                progress=path/'progress.json'
                error=json.loads(progress.read_text()).get('error') if progress.exists() else None
                update('failed',error=error or 'The calculation stopped unexpectedly. Open the run log for details.')
    except Exception as error:
        if process is not None and process.poll() is None:stop(process)
        update('failed',error=str(error))


if __name__=='__main__':run(sys.argv[1])
