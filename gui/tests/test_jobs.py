from pathlib import Path
import json
import tempfile
import time
import unittest
from unittest.mock import patch
from gui import jobs
from gui.worker import run


class JobTests(unittest.TestCase):
    def test_interrupted_worker_can_be_recognized(self):
        with tempfile.TemporaryDirectory() as tmp,patch.object(jobs,'JOBS',Path(tmp)):
            identifier='a'*32;path=jobs.folder(identifier);path.mkdir()
            jobs.write(path/'job.json',dict(id=identifier,status='running',created=time.time()))
            (path/'worker.pid').write_text('999999999')
            state=jobs.read(identifier)
            self.assertEqual(state['status'],'failed')
            self.assertIn('interrupted',state['error'])

    def test_cancel_queued_job_never_launches_inference(self):
        with tempfile.TemporaryDirectory() as tmp,patch.object(jobs,'JOBS',Path(tmp)/'jobs'):
            identifier='b'*32;path=jobs.folder(identifier);path.mkdir(parents=True)
            jobs.write(path/'job.json',dict(id=identifier,status='queued',created=time.time()))
            jobs.cancel(identifier)
            with patch('gui.worker.subprocess.Popen') as launch:run(identifier)
            launch.assert_not_called()
            self.assertEqual(jobs.read(identifier)['status'],'cancelled')


if __name__=='__main__':unittest.main()
