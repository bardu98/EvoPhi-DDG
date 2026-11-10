"""Start the local interface on loopback and open its browser page."""
from pathlib import Path
import argparse
import socket
import subprocess
import sys
import time
import urllib.request
import webbrowser
import shutil

ROOT=Path(__file__).resolve().parent.parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port',type=int,default=8501)
    parser.add_argument('--no-browser',action='store_true')
    args=parser.parse_args()
    with socket.socket() as sock:
        try:sock.bind(('127.0.0.1',args.port));port=args.port
        except OSError:sock.bind(('127.0.0.1',0));port=sock.getsockname()[1]
    command=[sys.executable,'-m','streamlit','run',str(ROOT/'gui/app.py'),
        '--server.address=127.0.0.1',f'--server.port={port}','--server.headless=true',
        '--browser.gatherUsageStats=false','--server.maxUploadSize=600','--client.toolbarMode=minimal',
        '--theme.base=light','--theme.primaryColor=#207e71','--theme.backgroundColor=#ffffff']
    server=subprocess.Popen(command,cwd=ROOT)
    url=f'http://localhost:{port}'
    try:
        ready=False
        for _ in range(120):
            if server.poll() is not None:return server.returncode
            try:
                with urllib.request.urlopen(f'http://127.0.0.1:{port}/_stcore/health',timeout=1) as response:
                    ready=response.status==200
            except OSError:pass
            if ready:break
            time.sleep(.25)
        if not ready:raise RuntimeError('The interface did not become ready. Check the startup messages above.')
        print(f'\nEvoPhi-DDG is ready: {url}\nKeep this window open while using the interface.\n',flush=True)
        if not args.no_browser:
            proc=Path('/proc/version')
            is_wsl=proc.exists() and 'microsoft' in proc.read_text().lower()
            try:
                if is_wsl and shutil.which('cmd.exe'):
                    subprocess.Popen(['cmd.exe','/c','start','',url],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
                else:webbrowser.open(url)
            except OSError:print('Open the address above in your browser.',flush=True)
        return server.wait()
    except KeyboardInterrupt:return 0
    finally:
        if server.poll() is None:
            server.terminate()
            try:server.wait(timeout=10)
            except subprocess.TimeoutExpired:server.kill();server.wait()


if __name__=='__main__':raise SystemExit(main())
