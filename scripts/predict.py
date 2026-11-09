"""Score single or multiple substitutions using prefiltered, fixed-context MSAs."""
from pathlib import Path
import argparse
import io
import re
import sys
import zipfile
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
from evophi import potential_weights, single, multiple


def read_alignment(path):
    records=[]
    name=None
    seq=[]
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        if line.startswith('>'):
            if name is not None:
                records.append((name,re.sub(r'[a-z.*]','',''.join(seq))))
            name=line[1:];seq=[]
        else:
            if name is None:
                raise ValueError('Alignment must be FASTA or A3M')
            seq.append(line.strip())
    if name is not None:
        records.append((name,re.sub(r'[a-z.*]','',''.join(seq))))
    if not records or len(set(len(s) for _,s in records))!=1:
        raise ValueError('MSA rows must have equal length after removing insertions')
    return records


def load_checkpoint(path,torch):
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            members=[n for n in archive.namelist() if Path(n).name=='MSA_Finetuning.pt']
            if members:
                return torch.load(io.BytesIO(archive.read(members[0])),map_location='cpu',weights_only=True)['MODEL_STATE']
    return torch.load(path,map_location='cpu',weights_only=True)['MODEL_STATE']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path,help='CSV: code, wt_seq, mut_seq, msa; optional janus_score')
    parser.add_argument('--checkpoint',required=True,type=Path,help='Published MSA_Finetuning.pt or enclosing ZIP')
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--device',default='cuda')
    parser.add_argument('--base-encoder',type=Path,help='Optional local ESM MSA Transformer base checkpoint; otherwise fair-esm downloads it')
    args=parser.parse_args()
    import torch
    import esm
    torch.set_num_threads(4)
    torch.backends.cuda.matmul.allow_tf32=False
    torch.backends.cudnn.allow_tf32=False
    frame=pd.read_csv(args.input)
    required={'code','wt_seq','mut_seq','msa'}
    if not required.issubset(frame.columns) or not frame.code.is_unique:
        raise ValueError('Require unique codes and columns: '+', '.join(sorted(required)))
    state=load_checkpoint(args.checkpoint,torch)
    if args.base_encoder:
        encoder,alphabet=esm.pretrained.load_model_and_alphabet_local(str(args.base_encoder))
    else:
        encoder,alphabet=esm.pretrained.esm_msa1b_t12_100M_UR50S()
    prefix='msa_transformer.'
    encoder.load_state_dict({k[len(prefix):]:v for k,v in state.items() if k.startswith(prefix)},strict=True)
    encoder.eval().to(args.device).requires_grad_(False)
    converter=alphabet.get_batch_converter()
    q=potential_weights(state['fc1.weight'].numpy(),state['fc2.weight'].numpy())
    output=[]
    for row in frame.itertuples():
        a,b=row.wt_seq,row.mut_seq
        if len(a)!=len(b) or not (2<=len(a)<=1022) or any(x not in 'ACDEFGHIKLMNPQRSTVWY' for x in a+b):
            raise ValueError(f'{row.code}: provide equal-length canonical-AA query windows of 2–1022 residues')
        sites=[i for i,(x,y) in enumerate(zip(a,b)) if x!=y]
        if not sites:
            raise ValueError(f'{row.code}: no substitutions')
        path=Path(row.msa)
        if not path.is_absolute():
            path=args.input.parent/path
        alignment=read_alignment(path)
        if alignment[0][1]!=a:
            raise ValueError(f'{row.code}: first MSA row must be the exact WT query')
        depth=min(len(alignment),16000//(len(a)+1))
        homologs=alignment[1:depth]
        features={}

        def phi(sequence,position):
            if sequence not in features:
                _,_,tokens=converter([[('query',sequence)]+homologs])
                with torch.inference_mode():
                    reps=encoder(tokens.to(args.device),repr_layers=[12])['representations'][12][0,0]
                # Match published pooling exactly: BOS present, EOS absent.
                features[sequence]=(reps.detach().float().cpu().numpy(),
                                    reps[1:-1].mean(0).detach().float().cpu().numpy())
            reps,mean=features[sequence]
            return np.concatenate((reps[position+1],mean))

        forward=[];reverse=[]
        for position in sites:
            ai=a[:position]+b[position]+a[position+1:]
            bi=b[:position]+a[position]+b[position+1:]
            forward.append(single(phi(a,position),phi(ai,position),q))
            reverse.append(single(phi(b,position),phi(bi,position),q))
        prediction=multiple(forward,reverse)
        result=dict(code=row.code,n_substitutions=len(sites),depth=depth,length=len(a),evophi_score=prediction)
        if hasattr(row,'janus_score'):
            result.update(janus_score=float(row.janus_score),ensemble_score=(prediction+float(row.janus_score))/2)
        output.append(result)
        print(f'{row.code}: {prediction:.2f} kcal/mol',flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(output).to_csv(args.output,index=False)


if __name__=='__main__':
    main()
