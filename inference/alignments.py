"""ColabFold UniRef-only searches, resumable tickets and the fixed HHfilter protocol."""
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import time

HERE=Path(__file__).resolve().parent
FILTER={'id':90,'diff':1000,'cov':0,'qid':0,'qsc':-20,'maxseq':65535}


def read_msa(path):
    import re
    rows=[];name=None;parts=[]
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith('#'):continue
        if line.startswith('>'):
            if name is not None:rows.append((name,re.sub('[a-z.*]','',''.join(parts))))
            name=line[1:];parts=[]
        else:parts.append(line.strip())
    if name is not None:rows.append((name,re.sub('[a-z.*]','',''.join(parts))))
    return rows


def validate(rows,sequence):
    if not rows or rows[0][1]!=sequence:
        raise ValueError('The first alignment row must match the complete WT sequence exactly.')
    if any(len(s)!=len(sequence) or set(s)-set('ABCDEFGHIJKLMNOPQRSTUVWXYZ-') for _,s in rows):
        raise ValueError('Alignment rows have inconsistent lengths or unsupported characters.')


def request(session,method,url,**kwargs):
    import requests
    for attempt in range(5):
        try:
            response=session.request(method,url,timeout=(15,90),**kwargs)
            if response.status_code in [429,500,502,503,504]:
                response.raise_for_status()
            response.raise_for_status();return response
        except requests.RequestException:
            if attempt==4:raise
            time.sleep(min(60,5*2**attempt))


def search(sequence,folder,server,max_wait,poll_seconds):
    import requests
    archive=folder/'result.tar.gz';ticket_file=folder/'ticket.json'
    raw=folder/'raw.a3m'
    if raw.exists():
        validate(read_msa(raw),sequence);return raw
    session=requests.Session();session.headers['User-Agent']='EvoPhi-DDG/1.0 (+https://github.com/bardu98/EvoPhi-DDG)'
    if not archive.exists():
        if ticket_file.exists():ticket=json.loads(ticket_file.read_text())
        else:
            # An accepted ticket is persisted and reused on restart. Never resubmit accepted jobs.
            for attempt in range(12):
                ticket=request(session,'POST',server+'/ticket/msa',data={'q':'>101\n'+sequence+'\n','mode':'all'}).json()
                if ticket.get('id') and ticket.get('status') not in ['UNKNOWN','RATELIMIT']:break
                if ticket.get('status') in ['ERROR','MAINTENANCE']:raise RuntimeError(f'MSA service: {ticket}')
                time.sleep(min(60,10+attempt*5))
            if not ticket.get('id') or ticket.get('status') in ['UNKNOWN','RATELIMIT']:
                raise RuntimeError('The MSA server did not accept the job; try again later.')
            ticket_file.write_text(json.dumps(ticket,indent=2)+'\n')
        start=time.monotonic();previous=None
        while True:
            status=ticket.get('status')
            if status!=previous:print(f'MSA server: {status} (ticket {ticket["id"]})',flush=True);previous=status
            if status=='COMPLETE':break
            if status in ['ERROR','MAINTENANCE']:raise RuntimeError(f'MSA service: {ticket}')
            if time.monotonic()-start>max_wait:
                raise TimeoutError('MSA search is still pending. Rerun the same command to resume the saved ticket.')
            time.sleep(poll_seconds)
            ticket=request(session,'GET',server+'/ticket/'+ticket['id']).json()
            ticket_file.write_text(json.dumps(ticket,indent=2)+'\n')
        data=request(session,'GET',server+'/result/download/'+ticket['id']).content
        tmp=archive.with_suffix('.part');tmp.write_bytes(data);tmp.replace(archive)
    # Read only the expected member; never extract arbitrary archive paths.
    with tarfile.open(archive,'r:gz') as tar:
        members=[m for m in tar.getmembers() if Path(m.name).name=='uniref.a3m']
        if len(members)!=1:raise ValueError('MSA response does not contain one uniref.a3m file.')
        text=tar.extractfile(members[0]).read().decode()
    blocks=[s.strip()+'\n' for s in text.split('\x00') if s.strip()]
    candidates=[]
    for block in blocks:
        candidate=folder/'candidate.a3m';candidate.write_text(block)
        if read_msa(candidate) and read_msa(candidate)[0][1]==sequence:candidates.append(block)
    (folder/'candidate.a3m').unlink(missing_ok=True)
    if len(candidates)!=1:raise ValueError('Could not identify the requested query in the MSA response.')
    raw.write_text(candidates[0]);validate(read_msa(raw),sequence)
    return raw


def get_alignment(sequence,cache,hhfilter,server,max_wait,poll_seconds,provided=None,already_filtered=False):
    key=hashlib.sha256((sequence+'\n'+server+'\n'+json.dumps(FILTER,sort_keys=True)).encode()).hexdigest()
    # Include the input bytes when users supply an alignment, avoiding stale-context reuse.
    if provided:key+='_'+hashlib.sha256(Path(provided).read_bytes()).hexdigest()[:16]
    key+='_ready' if already_filtered else '_hhfilter'
    folder=Path(cache)/'msa'/key;folder.mkdir(parents=True,exist_ok=True)
    target=folder/'alignment.fasta'
    if target.exists():
        rows=read_msa(target);validate(rows,sequence)
        meta=folder/'provenance.json'
        metadata=json.loads(meta.read_text()) if meta.exists() else {}
        if metadata.get('alignment_sha256') and hashlib.sha256(target.read_bytes()).hexdigest()!=metadata['alignment_sha256']:
            raise ValueError(f'Cached alignment checksum mismatch: {target}')
        return rows,dict(metadata,cache=str(target),reused=True)
    raw=Path(provided).resolve() if provided else search(sequence,folder,server,max_wait,poll_seconds)
    validate(read_msa(raw),sequence)
    temporary=folder/'alignment.tmp.fasta'
    if already_filtered:
        rows=read_msa(raw)
        temporary.write_text(''.join(f'>{n}\n{s}\n' for n,s in rows))
        version='user-provided, already filtered'
    else:
        binary=shutil.which(hhfilter)
        if binary is None:raise FileNotFoundError('HHfilter not found. Activate the environment from inference/environment.yml.')
        filtered=folder/'filtered.a3m'
        command=[binary,'-i',str(raw),'-o',str(filtered)]
        for k,v in FILTER.items():command += ['-'+k,str(v)]
        with (folder/'filter.log').open('w') as log:
            subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT)
            subprocess.run(['perl',str(HERE/'vendor/reformat.pl'),'-v','1','-r','a3m','fas',str(filtered),str(temporary)],
                           check=True,stdout=log,stderr=subprocess.STDOUT)
        version=subprocess.run([binary,'-h'],capture_output=True,text=True).stdout.splitlines()[0]
        rows=read_msa(temporary)
    validate(rows,sequence)
    temporary.replace(target)
    if len(rows)<2:print('Note: no homolog rows remain; the model will use the query alone.',flush=True)
    metadata=dict(query_sha256=hashlib.sha256(sequence.encode()).hexdigest(),n_rows=len(rows),query_length=len(sequence),
        hhfilter=version,filter_parameters=FILTER,mode='all',result='uniref.a3m only',server=server,
        alignment_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),already_filtered=already_filtered,
        database_release='Not attested by individual service jobs; see documentation for database provenance.')
    (folder/'provenance.json').write_text(json.dumps(metadata,indent=2)+'\n')
    return rows,dict(cache=str(target),reused=False,**metadata)
