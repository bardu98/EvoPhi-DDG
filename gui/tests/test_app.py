from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch

from streamlit.testing.v1 import AppTest
from gui import jobs

ROOT=Path(__file__).resolve().parents[2]


def button(app,label):
    return next(item for item in app.button if item.label==label)


class InterfaceTests(unittest.TestCase):
    def test_private_download_error_shows_browser_solution(self):
        with tempfile.TemporaryDirectory() as tmp,patch.object(jobs,'JOBS',Path(tmp)):
            identifier='a'*32;path=jobs.folder(identifier);path.mkdir()
            jobs.write(path/'job.json',dict(id=identifier,status='failed',created=time.time(),variants=1,proteins=1,
                error="Command '['gh', 'api']' returned non-zero exit status 4."))
            app=AppTest.from_file(str(ROOT/'gui/app.py'),default_timeout=30)
            app.query_params['job']=identifier
            app.run()
            self.assertEqual(len(app.exception),0)
            self.assertTrue(any('GitHub authentication is needed' in item.value for item in app.error))
            self.assertFalse(any('exit status 4' in item.value for item in app.error))
            self.assertTrue(any('Verify and save checkpoint' in item.value for item in app.markdown))

    def test_example_and_validation_before_launch(self):
        with tempfile.TemporaryDirectory() as tmp,patch.object(jobs,'JOBS',Path(tmp)),patch.object(jobs,'create') as create:
            app=AppTest.from_file(str(ROOT/'gui/app.py'),default_timeout=30).run()
            self.assertEqual(len(app.exception),0)
            button(app,'Load a real example').click().run()
            button(app,'Check input').click().run()
            self.assertTrue(any('Input valid: 1' in item.value for item in app.success))
            app.text_area(key='mutation_text').set_value('A13V').run()
            button(app,'Calculate predictions').click().run()
            self.assertTrue(any('sequence residue 13 is K' in item.value for item in app.error))
            create.assert_not_called()
            self.assertEqual(len(app.exception),0)


if __name__=='__main__':unittest.main()
