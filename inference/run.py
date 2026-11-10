"""One command: WT sequence and substitutions -> MSAs -> both models -> ensemble."""
from pathlib import Path
import argparse
import json
import os
import sys
import time

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
from inference.inputs import sequence,read_fasta,variant,window


def arguments():
    parser=argparse.ArgumentParser(description=__doc__)
    source=parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--sequence',help='Complete WT amino-acid sequence')
    source.add_argument('--fasta',type=Path,help='FASTA with one complete WT sequence')
    source.add_argument('--input',type=Path,help='CSV with ID, Sequence and MTS columns')
    parser.add_argument('--mutations',nargs='+',help='Separate variants by spaces; join sites of one variant with underscores')
    parser.add_argument('--output',type=Path,default=Path('predictions.csv'))
    parser.add_argument('--cache',type=Path,default=ROOT/'.cache/inference')
    parser.add_argument('--weights',type=Path,default=Path(__file__).resolve().parent/'weights')
    parser.add_argument('--device',choices=['auto','cpu','cuda'],default='auto')
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--hhfilter',default='hhfilter',help='HHfilter executable or full path')
    parser.add_argument('--msa-server',default='https://api.colabfold.com')
    parser.add_argument('--max-msa-wait',type=int,default=7200,help='Seconds before pausing; rerun resumes the saved ticket')
    parser.add_argument('--poll-seconds',type=int,default=10)
    parser.add_argument('--msa',type=Path,help='Optional existing MSA; valid for inputs with one WT sequence')
    parser.add_argument('--msa-filtered',action='store_true',help='Explicitly use supplied MSAs as already filtered; skip HHfilter')
    parser.add_argument('--esm2-checkpoint',type=Path,help='Optional local ESM2 checkpoint for offline use')
    parser.add_argument('--progress-file',type=Path,help='Optional machine-readable progress file for the local interface')
    return parser.parse_args()


def load_inputs(args):
    import pandas as pd
    if args.input:
        if args.mutations:raise ValueError('Put mutations in the CSV MTS column when using --input.')
        frame=pd.read_csv(args.input,keep_default_na=False)
        if not {'Sequence','MTS'}.issubset(frame):raise ValueError('Input CSV needs Sequence and MTS columns; ID is optional.')
        if 'ID' not in frame:frame['ID']=[f'variant_{i+1}' for i in range(len(frame))]
    else:
        if not args.mutations:raise ValueError('Supply --mutations, for example --mutations K13A K13A_L22V.')
        wt=read_fasta(args.fasta) if args.fasta else sequence(args.sequence)
        frame=pd.DataFrame([dict(ID=f'variant_{i+1}',Sequence=wt,MTS=m) for i,m in enumerate(args.mutations)])
    if frame.empty:raise ValueError('No variants supplied.')
    if args.msa and frame.Sequence.nunique()!=1:raise ValueError('--msa requires one WT sequence; use an MSA column for multiple proteins.')
    rows=[]
    for index,row in frame.iterrows():
        wt=sequence(row.Sequence);mutant,mts,changes=variant(wt,row.MTS)
        left,right=window(len(wt),changes)
        msa=args.msa
        if not msa and 'MSA' in frame and row.MSA:
            msa=Path(row.MSA)
            if not msa.is_absolute():msa=args.input.parent/msa
        if args.msa_filtered and msa is None:raise ValueError('--msa-filtered requires a supplied alignment for every record.')
        rows.append(dict(ID=str(row.ID),Sequence=wt,MTS=mts,mutant=mutant,left=left,right=right,msa=msa))
    return rows


def main():
    args=arguments()
    def progress(stage,completed=0,total=0,error=None):
        if args.progress_file:
            args.progress_file.parent.mkdir(parents=True,exist_ok=True)
            target=args.progress_file;tmp=target.with_suffix('.tmp')
            tmp.write_text(json.dumps(dict(stage=stage,completed=completed,total=total,error=error,updated=time.time()))+'\n')
            tmp.replace(target)
    try:
        records=load_inputs(args)  # Reject invalid coordinates before network activity.
        progress('checking',total=len(records))
        import shutil
        if not args.msa_filtered and not shutil.which(args.hhfilter):
            raise FileNotFoundError('HHfilter is missing. Create/activate the environment from inference/environment.yml.')
        import torch,pandas as pd
        from inference.weights import ensure_weights,sha256
        from inference.alignments import get_alignment
        from inference.models import EvoPhi,Janus
        torch.set_num_threads(args.threads)
        torch.manual_seed(42)
        torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
        device='cuda' if args.device=='auto' and torch.cuda.is_available() else ('cpu' if args.device=='auto' else args.device)
        if device=='cuda' and not torch.cuda.is_available():raise ValueError('CUDA unavailable; use --device cpu or install a CUDA-enabled PyTorch build.')
        print(f'Predicting {len(records)} variant(s) on {device}. Positive DDG means stabilizing; units are kcal/mol.',flush=True)
        args.cache.mkdir(parents=True,exist_ok=True)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        progress('models',total=2)
        weights,manifest=ensure_weights(args.weights)
        alignments={};provenance=[]
        n_contexts=len({(row['Sequence'],str(row['msa'])) for row in records})
        progress('alignments',total=n_contexts)
        for row in records:
            key=(row['Sequence'],str(row['msa']))
            if key not in alignments:
                print(f'Preparing MSA for {row["ID"]} ({len(row["Sequence"])} residues)…',flush=True)
                alignments[key],meta=get_alignment(row['Sequence'],args.cache,args.hhfilter,args.msa_server.rstrip('/'),
                    args.max_msa_wait,args.poll_seconds,row['msa'],args.msa_filtered)
                provenance.append(meta)
                progress('alignments',completed=len(alignments),total=n_contexts)
        print('Running EvoPhi-DDG…',flush=True)
        progress('evophi',total=len(records))
        model=EvoPhi(weights['evophi'],device)
        output=[]
        for row in records:
            value,depth=model.predict(row['Sequence'],row['mutant'],alignments[(row['Sequence'],str(row['msa']))],row['left'],row['right'])
            output.append(dict(ID=row['ID'],Sequence=row['Sequence'],MTS=row['MTS'],DDG_EvoPhi_DDG=value,
                               msa_depth=depth,msa_window_start=row['left']+1,msa_window_end=row['right']))
            print(f'{row["ID"]} {row["MTS"]}: EvoPhi-DDG {value:.2f}',flush=True)
            progress('evophi',completed=len(output),total=len(records))
        del model
        if device=='cuda':torch.cuda.empty_cache()
        print('Running JanusDDG (first use also downloads ESM2)…',flush=True)
        progress('janus',total=len(records))
        model=Janus(weights['janus'],device,args.cache,args.esm2_checkpoint)
        partial=args.output.with_name(args.output.stem+'.partial.csv')
        for index,(row,result) in enumerate(zip(records,output)):
            value=model.predict(row['Sequence'],row['mutant'])
            result['DDG_JanusDDG']=value
            result['DDG_Ensemble']=(result['DDG_EvoPhi_DDG']+value)/2
            pd.DataFrame(output[:index+1]).to_csv(partial,index=False)
            print(f'{row["ID"]} {row["MTS"]}: JanusDDG {value:.2f}, ensemble {result["DDG_Ensemble"]:.2f}',flush=True)
            progress('janus',completed=index+1,total=len(records))
        # Preserve full numerical precision in CSV; console display uses two decimals.
        partial.replace(args.output)
        report=dict(models='Published frozen JanusDDG and zero-hidden-bias EvoPhi-DDG; no training or calibration',
            units='kcal/mol',positive='stabilizing',ensemble='equal-weight arithmetic mean',device=device,
            torch=torch.__version__,seed=42,encoder_dtype='FP32',potential_dtype='FP64',max_tokens=16000,
            query_pooling='Published [1:-1], omits final biological residue',
            checkpoints={name:sha256(path) for name,path in weights.items()},release=manifest['release_url'],alignments=provenance)
        args.output.with_suffix('.provenance.json').write_text(json.dumps(report,indent=2)+'\n')
        print(f'Predictions saved to {args.output}',flush=True)
        progress('complete',completed=len(records),total=len(records))
    except Exception as error:
        progress('failed',error=str(error))
        print(f'Error: {error}',file=sys.stderr,flush=True)
        if os.environ.get('EVOPHI_DEBUG'):raise
        return 1
    return 0


if __name__=='__main__':raise SystemExit(main())
