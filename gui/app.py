"""Browser interface for EvoPhi-DDG, JanusDDG and their fixed ensemble."""
from pathlib import Path
import hashlib
import importlib.util
import json
import shutil
import sys
import time
from datetime import datetime

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
import pandas as pd
import streamlit as st
from gui import inputs,jobs

st.set_page_config(page_title='EvoPhi-DDG',page_icon='E',layout='wide',initial_sidebar_state='expanded')
st.markdown('''<style>
.block-container {max-width:1180px; padding-top:3.8rem;}
h1 {letter-spacing:-0.045em; font-weight:650 !important;}
h2,h3 {letter-spacing:-0.025em;}
[data-testid="stSidebar"] {background:#f3f6f5;}
[data-testid="stSidebar"] .block-container {padding-top:2rem;}
.eyebrow {font-size:.76rem; letter-spacing:.14em; text-transform:uppercase; color:#28736b; font-weight:650; margin-bottom:.6rem;}
.lead {color:#596a68; font-size:1.05rem; max-width:780px; margin-bottom:1.5rem;}
.step-line {padding:.45rem 0; border-bottom:1px solid #edf0ef; color:#596a68;}
.step-line.active {color:#176f64; font-weight:650;}
</style>''',unsafe_allow_html=True)


def decode(upload):
    try:return upload.getvalue().decode('utf-8-sig')
    except UnicodeDecodeError:raise ValueError('Upload a plain-text FASTA or UTF-8 CSV file.')


@st.cache_data(ttl=60)
def readiness():
    missing=[]
    for module,name in [('torch','PyTorch'),('esm','fair-esm')]:
        if importlib.util.find_spec(module) is None:missing.append(name)
    for command,name in [('hhfilter','HHfilter'),('perl','Perl')]:
        if shutil.which(command) is None:missing.append(name)
    gpu=False
    if 'PyTorch' not in missing:
        import torch
        gpu=torch.cuda.is_available()
    return missing,gpu


def save_checkpoint(upload,manifest):
    expected=next((item for item in manifest['assets'].values() if item['filename']==upload.name),None)
    if expected is None:raise ValueError('Select evophi-ddg.pt or janusddg.pt from the checkpoint release.')
    folder=ROOT/'inference/weights';folder.mkdir(exist_ok=True)
    target=folder/expected['filename'];tmp=target.with_suffix('.upload')
    h=hashlib.sha256();upload.seek(0)
    with tmp.open('wb') as stream:
        while block:=upload.read(4*1024*1024):h.update(block);stream.write(block)
    if h.hexdigest()!=expected['sha256']:
        tmp.unlink(missing_ok=True);raise ValueError('This file does not match the released checkpoint. Download the file again.')
    tmp.replace(target)


if 'active_job' not in st.session_state:
    candidate=st.query_params.get('job')
    try:st.session_state.active_job=candidate if candidate and jobs.folder(candidate).exists() else None
    except ValueError:st.session_state.active_job=None

with st.sidebar:
    st.markdown('### EvoPhi-DDG')
    st.caption('Protein stability · local workspace')
    st.divider()
    missing,gpu=readiness()
    if missing:st.warning('Missing: '+', '.join(missing)+'. Start the app with the guided launcher.')
    else:st.success('Prediction environment ready')
    st.caption('NVIDIA GPU available' if gpu else 'CPU inference available')
    device_label=st.selectbox('Compute device',['Automatic','CPU','NVIDIA GPU'] if gpu else ['Automatic','CPU'],key='device')
    device={'Automatic':'auto','CPU':'cpu','NVIDIA GPU':'cuda'}[device_label]
    st.caption('Automatic uses the GPU when available. CPU inference takes longer.')
    history=jobs.recent()
    if history:
        st.divider()
        labels={item['id']:datetime.fromtimestamp(item['created']).strftime('%d %b, %H:%M')+
            f' · {item["variants"]} variant(s) · {item["status"]}' for item in history}
        selected=st.selectbox('Previous runs',[None]+list(labels),format_func=lambda key:'Select a run' if key is None else labels[key])
        if st.button('Open selected run',disabled=selected is None):
            st.session_state.active_job=selected;st.query_params['job']=selected;st.rerun()
    st.divider()
    with st.expander('Model files'):
        manifest=json.loads((ROOT/'inference/weights.json').read_text())
        for name,item in manifest['assets'].items():
            present=(ROOT/'inference/weights'/item['filename']).exists()
            st.write(('Available · ' if present else 'Download on first run · ')+('EvoPhi-DDG' if name=='evophi' else 'JanusDDG'))
        st.caption('While the repository is private, you can download the two checkpoints in your browser and add them here.')
        st.link_button('Open checkpoint downloads',manifest['release_url'])
        checkpoint=st.file_uploader('Add a downloaded checkpoint',type=['pt'],key='checkpoint_upload')
        if st.button('Verify and save checkpoint',disabled=checkpoint is None):
            try:save_checkpoint(checkpoint,manifest);st.success('Checkpoint verified and saved.')
            except ValueError as error:st.error(str(error))

st.markdown('<div class="eyebrow">Sequence to stability</div>',unsafe_allow_html=True)
st.title('Protein stability predictions')
st.markdown('<div class="lead">Upload a protein or a batch of variants. Get EvoPhi-DDG, JanusDDG and their ensemble in one downloadable table.</div>',unsafe_allow_html=True)

mode=st.radio('Input format',['One protein','Batch CSV','FASTA + variant table'],horizontal=True,key='input_mode')
active=None
if st.session_state.active_job:
    try:active=jobs.read(st.session_state.active_job)
    except (OSError,ValueError):st.session_state.active_job=None
busy=active is not None and active['status'] in ['queued','running']

uploaded_msa=None;msa_filtered=False
if mode=='One protein':
    left,right=st.columns([1.6,1],gap='large')
    with left:
        if st.button('Load a real example',disabled=busy):
            st.session_state.sequence_text=(ROOT/'inference/examples/protein.fasta').read_text()
            st.session_state.mutation_text='K13A'
            st.session_state.fasta_revision=st.session_state.get('fasta_revision',0)+1
        sequence_text=st.text_area('WT sequence or FASTA',height=165,key='sequence_text',placeholder='Paste the complete wild-type amino-acid sequence')
        fasta=st.file_uploader('Or upload one FASTA file',type=['fasta','fa','faa','fas'],key='single_fasta_'+str(st.session_state.get('fasta_revision',0)))
        mutation_text=st.text_area('Variants — one per line',height=105,key='mutation_text',placeholder='K13A\nL22V\nK13A_L22V')
        st.caption('Positions start at 1 in the exact sequence supplied above. Join sites in one multiple mutant with underscores.')
    with right:
        with st.container(border=True):
            st.markdown('#### Your prediction')
            st.write('Each input variant receives three ΔΔG estimates in kcal/mol.')
            st.write('**Positive:** predicted stabilization')
            st.write('**Negative:** predicted destabilization')
            st.caption('For example, K13A changes lysine at position 13 to alanine. K13A_L22V is one double mutant.')
        with st.expander('Use an existing alignment'):
            uploaded_msa=st.file_uploader('WT alignment',type=['a3m','fasta','fas','fa'],key='alignment_upload')
            msa_filtered=st.checkbox('Already filtered with the intended HHfilter protocol',value=False)
            st.caption('The first row must match the complete WT sequence. An uploaded alignment avoids a new online homology search.')
elif mode=='Batch CSV':
    csv_file=st.file_uploader('Variant CSV',type=['csv','tsv'],key='batch_csv')
    st.caption('Use Sequence and MTS columns; ID is optional. One row represents one variant. Several rows can share the same WT sequence.')
    st.download_button('Download an example CSV',(ROOT/'inference/examples/variants.csv').read_bytes(),file_name='variants_example.csv',mime='text/csv')
else:
    first,second=st.columns(2)
    with first:batch_fasta=st.file_uploader('Protein FASTA',type=['fasta','fa','faa','fas'],key='batch_fasta')
    with second:variant_file=st.file_uploader('Variant table',type=['csv','tsv'],key='fasta_variants')
    st.caption('The table needs ID and MTS columns. Each ID must match the first word after > in the FASTA header.')
    st.download_button('Download a variant-table template','ID,MTS\nprotein_A,K13A\nprotein_A,K13A_L22V\nprotein_B,A1V\n',
                       file_name='variant_table_template.csv',mime='text/csv')


def current_input():
    if mode=='One protein':return inputs.single_input(decode(fasta) if fasta else sequence_text,mutation_text)
    if mode=='Batch CSV':
        if csv_file is None:raise ValueError('Upload the variant CSV first.')
        return inputs.batch_input(decode(csv_file))
    if batch_fasta is None or variant_file is None:raise ValueError('Upload both the protein FASTA and the variant table.')
    return inputs.batch_input(decode(variant_file),decode(batch_fasta))

st.caption('New homology searches send WT sequences to the public ColabFold service. First use also downloads the model weights; later runs reuse cached files.')
check,run,_=st.columns([1,1.35,3])
validate_clicked=check.button('Check input',disabled=busy,width='stretch')
run_clicked=run.button('Calculate predictions',type='primary',disabled=busy,width='stretch')
if validate_clicked or run_clicked:
    try:
        frame=current_input()
        if uploaded_msa:
            from inference.alignments import validate as validate_msa,read_msa
            import tempfile
            with tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp)/'alignment.a3m';path.write_bytes(uploaded_msa.getvalue())
                validate_msa(read_msa(path),frame.iloc[0].Sequence)
        if run_clicked:
            if missing:raise ValueError('The prediction environment is incomplete. Close this app and use the guided launcher.')
            identifier=jobs.create(frame,device,uploaded_msa.getvalue() if uploaded_msa else None,msa_filtered)
            st.session_state.active_job=identifier;st.query_params['job']=identifier;st.rerun()
        else:
            st.success(f'Input valid: {len(frame)} variant(s), {frame.Sequence.nunique()} distinct WT sequence(s).')
            preview=frame[['ID','MTS']].copy();preview['Sequence length']=frame.Sequence.str.len()
            st.dataframe(preview,hide_index=True,width='stretch')
    except ValueError as error:st.error(str(error))


@st.fragment(run_every=2)
def progress_panel():
    identifier=st.session_state.get('active_job')
    if not identifier:return
    state=jobs.read(identifier);path=jobs.folder(identifier)
    st.divider()
    st.subheader('Run progress')
    st.caption(f'{state["variants"]} variant(s) · {state["proteins"]} WT sequence(s) · '+datetime.fromtimestamp(state['created']).strftime('%d %b %Y, %H:%M'))
    status=state['status'];progress=state.get('progress',{});stage=progress.get('stage','checking')
    if status in ['complete','cancelled','failed'] and st.session_state.get('last_terminal_job')!=(identifier,status):
        st.session_state.last_terminal_job=(identifier,status)
        st.rerun()
    names={'checking':'Checking inputs and environment','models':'Preparing model weights','alignments':'Finding and filtering homologs',
           'evophi':'Predicting with EvoPhi-DDG','janus':'Predicting with JanusDDG','complete':'Predictions complete'}
    if status in ['queued','running']:
        label='Waiting for the other local prediction job' if status=='queued' else names.get(stage,'Calculation in progress')
        completed,total=progress.get('completed',0),progress.get('total',0)
        with st.status(label,state='running',expanded=True):
            if status=='running' and total and stage in ['alignments','evophi','janus']:
                st.write(f'{completed} of {total} '+('protein alignments ready' if stage=='alignments' else 'variants evaluated'))
            st.caption('You may refresh or close this browser tab. Keep the computer running; the calculation continues in the background.')
            if stage=='alignments':st.caption('The remote server queue determines search time. No time estimate is assumed.')
            if stage=='janus':st.caption('On first use, this stage also downloads the ESM2 backbone.')
        if st.button('Stop this run',key='cancel_run'):
            jobs.cancel(identifier);st.info('Stopping the calculation… Cached alignments will be retained.')
    elif status=='complete':
        st.success('Predictions ready')
        frame=pd.read_csv(path/'predictions.csv')
        columns=['ID','MTS','DDG_EvoPhi_DDG','DDG_JanusDDG','DDG_Ensemble']
        st.dataframe(frame[columns],hide_index=True,width='stretch',column_config={
            'DDG_EvoPhi_DDG':st.column_config.NumberColumn('EvoPhi-DDG',format='%.2f'),
            'DDG_JanusDDG':st.column_config.NumberColumn('JanusDDG',format='%.2f'),
            'DDG_Ensemble':st.column_config.NumberColumn('Ensemble',format='%.2f')})
        st.caption('ΔΔG in kcal/mol. Positive values indicate stabilization. The CSV retains full numerical precision.')
        a,b=st.columns([1,2])
        a.download_button('Download predictions CSV',(path/'predictions.csv').read_bytes(),file_name='predictions.csv',mime='text/csv',type='primary')
        provenance=path/'predictions.provenance.json'
        if provenance.exists():b.download_button('Download run details',provenance.read_bytes(),file_name='prediction_details.json',mime='application/json')
    elif status=='cancelled':
        st.info('This run was stopped. Its completed alignment searches remain cached.')
    else:
        st.error(state.get('error') or 'The calculation stopped before completion.')
        if 'download' in (state.get('error') or '').lower():st.caption('For this private preview, you can add the released files under Model files in the sidebar.')
    if status in ['failed','cancelled']:
        if st.button('Retry using cached work',key='retry_run'):
            new=jobs.retry(identifier);st.session_state.active_job=new;st.query_params['job']=new;st.rerun()
    with st.expander('Run log'):
        st.code(jobs.log_tail(identifier) or 'The calculation is starting.',language=None)


progress_panel()
