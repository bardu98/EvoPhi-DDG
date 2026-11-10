import subprocess
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
from inference.weights import ensure_weights


class PrivateDownloadTests(unittest.TestCase):
    def test_missing_github_login_has_actionable_message_and_removes_partial(self):
        failure=subprocess.CalledProcessError(4,['gh','api'],stderr=b'Authentication required')
        with (tempfile.TemporaryDirectory() as tmp,
              patch('inference.weights.urllib.request.urlopen',side_effect=OSError('private asset')),
              patch('shutil.which',return_value='/bin/gh'),
              patch('inference.weights.subprocess.run',side_effect=failure)):
            with self.assertRaisesRegex(RuntimeError,'GitHub CLI is not signed in') as result:
                ensure_weights(tmp)
            self.assertIn('Model files',str(result.exception))
            self.assertIn('gh auth login',str(result.exception))
            self.assertFalse(list(Path(tmp).glob('*.part')))
            self.assertFalse(list(Path(tmp).glob('*.pt')))

    def test_other_github_failure_explains_access_without_asserting_missing_login(self):
        failure=subprocess.CalledProcessError(1,['gh','api'],stderr=b'HTTP 404')
        with (tempfile.TemporaryDirectory() as tmp,
              patch('inference.weights.urllib.request.urlopen',side_effect=OSError('private asset')),
              patch('shutil.which',return_value='/bin/gh'),
              patch('inference.weights.subprocess.run',side_effect=failure)):
            with self.assertRaisesRegex(RuntimeError,'Check that the signed-in account has access'):
                ensure_weights(tmp)
            self.assertFalse(list(Path(tmp).glob('*.part')))


if __name__=='__main__':unittest.main()
