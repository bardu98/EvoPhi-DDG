from pathlib import Path
import os
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]


class LauncherTests(unittest.TestCase):
    def test_first_setup_and_second_launch_reuses_environment(self):
        with tempfile.TemporaryDirectory(prefix='evophi setup ') as tmp:
            path=Path(tmp);conda=path/'conda';log=path/'calls.log';prefix=path/'app environment'
            conda.write_text('''#!/usr/bin/env bash
set -e
echo setup >> "$EVOPHI_LAUNCHER_LOG"
while [[ $# -gt 0 ]]; do
  if [[ "$1" == --prefix ]]; then prefix="$2";break;fi
  shift
done
mkdir -p "$prefix/bin"
printf '#!/usr/bin/env bash\nprintf "launch\\n" >> "$EVOPHI_LAUNCHER_LOG"\n' > "$prefix/bin/python"
chmod +x "$prefix/bin/python"
''')
            conda.chmod(0o755)
            env=dict(os.environ,CONDA_EXE=str(conda),EVOPHI_APP_PREFIX=str(prefix),EVOPHI_LAUNCHER_LOG=str(log))
            for _ in range(2):
                subprocess.run(['bash',str(ROOT/'gui/start.sh'),'--no-browser'],env=env,check=True,capture_output=True,text=True)
            self.assertEqual(log.read_text().splitlines(),['setup','launch','launch'])


if __name__=='__main__':unittest.main()
