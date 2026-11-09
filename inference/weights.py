"""Fetch and verify the ready-to-use checkpoint assets once, with atomic writes."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import urllib.request

HERE=Path(__file__).resolve().parent


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def ensure_weights(folder):
    manifest=json.loads((HERE/'weights.json').read_text())
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    paths={}
    for name,item in manifest['assets'].items():
        target=folder/item['filename'];paths[name]=target
        if target.exists():
            if sha256(target)!=item['sha256']:
                raise ValueError(f'Checkpoint checksum mismatch: {target}. Remove it and rerun to download again.')
            print(f'Checkpoint ready: {name}',flush=True);continue
        print(f'Downloading {name} checkpoint ({item["bytes"]/1024**2:.0f} MiB)…',flush=True)
        tmp=target.with_suffix(target.suffix+'.part')
        # gh uses the user's own account when the repository/release is private.
        # Once the repository is public, HTTPS downloads need no GitHub login.
        try:
            request=urllib.request.Request(item['url'],headers={'User-Agent':'EvoPhi-DDG/1.0'})
            with urllib.request.urlopen(request,timeout=60) as source,tmp.open('wb') as out:
                while block:=source.read(4*1024*1024):out.write(block)
        except Exception as download_error:
            import shutil
            if not shutil.which('gh'):
                raise RuntimeError('Cannot download checkpoint. While this repository is private, download the release assets '
                    f'from {manifest["release_url"]} into {folder}, or authenticate GitHub CLI with gh auth login. '
                    'For offline use, copy both checkpoint files into that folder.') from download_error
            asset_id=item['asset_id']
            with tmp.open('wb') as out:
                subprocess.run(['gh','api',f'repos/bardu98/EvoPhi-DDG/releases/assets/{asset_id}',
                                '-H','Accept: application/octet-stream'],stdout=out,check=True)
        if sha256(tmp)!=item['sha256']:
            tmp.unlink(missing_ok=True)
            raise ValueError(f'Downloaded {name} checkpoint failed its SHA256 check.')
        os.replace(tmp,target)
    return paths,manifest
