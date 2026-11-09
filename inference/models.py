"""Frozen EvoPhi-DDG and native JanusDDG inference with portable state dictionaries."""
from argparse import Namespace
import hashlib
from pathlib import Path
import numpy as np


class EvoPhi:
    def __init__(self,checkpoint,device):
        import torch,esm
        saved=torch.load(checkpoint,map_location='cpu',weights_only=True)
        if saved.get('model')!='EvoPhi-DDG':raise ValueError('Not an EvoPhi-DDG checkpoint.')
        self.device=device
        self.alphabet=esm.Alphabet.from_architecture('msa_transformer')
        self.encoder=esm.MSATransformer(Namespace(**saved['msa_args']),self.alphabet)
        self.encoder.load_state_dict(saved['encoder_state'],strict=True)
        self.encoder.eval().requires_grad_(False).to(device)
        self.converter=self.alphabet.get_batch_converter()
        self.q=saved['potential_weights'].double().numpy()
        self.context=None;self.features={}

    def predict(self,wt,mutant,alignment,left,right):
        import torch
        from evophi import single,multiple
        a,b=wt[left:right],mutant[left:right]
        sites=[p for p,(x,y) in enumerate(zip(a,b)) if x!=y]
        depth=min(len(alignment),16000//(len(a)+1))
        homologs=[(name,seq[left:right]) for name,seq in alignment[1:depth]]
        key=(a,left,right,tuple(homologs))
        if key!=self.context:self.context=key;self.features={}

        def phi(seq,p):
            if seq not in self.features:
                _,_,tokens=self.converter([[('query',seq)]+homologs])
                with torch.inference_mode():
                    reps=self.encoder(tokens.to(self.device),repr_layers=[12])['representations'][12][0,0]
                self.features[seq]=(reps.detach().cpu().float().numpy(),reps[1:-1].mean(0).detach().cpu().float().numpy())
            reps,mean=self.features[seq]
            return np.concatenate((reps[p+1],mean))
        forward=[];reverse=[]
        for p in sites:
            ai=a[:p]+b[p]+a[p+1:];bi=b[:p]+a[p]+b[p+1:]
            forward.append(single(phi(a,p),phi(ai,p),self.q))
            reverse.append(single(phi(b,p),phi(bi,p),self.q))
        return multiple(forward,reverse),depth


class Janus:
    def __init__(self,checkpoint,device,cache,esm2_checkpoint=None):
        import torch,esm
        from inference.vendor.janus_model import Cross_Attention_DDG,TransformerRegression
        saved=torch.load(checkpoint,map_location='cpu',weights_only=True)
        if saved.get('model')!='JanusDDG':raise ValueError('Not a JanusDDG checkpoint.')
        self.head=Cross_Attention_DDG(TransformerRegression)
        self.head.load_state_dict(saved['state_dict'],strict=True)
        self.head.eval().requires_grad_(False).to(device)
        self.device=device
        self.cache=Path(cache)/'esm2_embeddings';self.cache.mkdir(parents=True,exist_ok=True)
        if esm2_checkpoint:
            self.encoder,alphabet=esm.pretrained.load_model_and_alphabet_local(str(esm2_checkpoint))
        else:
            torch.hub.set_dir(str(Path(cache)/'esm_models'))
            self.encoder,alphabet=esm.pretrained.esm2_t33_650M_UR50D()
        self.encoder.eval().requires_grad_(False).to(device)
        self.converter=alphabet.get_batch_converter()

    def embedding(self,sequence):
        import torch
        file=self.cache/(hashlib.sha256(sequence.encode()).hexdigest()+'.npy')
        if file.exists():
            feature=np.load(file,allow_pickle=False)
            if feature.shape!=(len(sequence),1280):raise ValueError(f'Invalid ESM2 cache: {file}')
            return torch.from_numpy(feature).unsqueeze(0).to(self.device)
        _,_,tokens=self.converter([('query',sequence)])
        with torch.inference_mode():
            reps=self.encoder(tokens.to(self.device),repr_layers=[33])['representations'][33][0,1:-1]
        feature=reps.detach().cpu().float().numpy()
        tmp=file.with_suffix('.part.npy');np.save(tmp,feature);tmp.replace(file)
        return torch.from_numpy(feature).unsqueeze(0).to(self.device)

    def predict(self,wt,mutant):
        import torch
        a=self.embedding(wt);b=self.embedding(mutant)
        with torch.inference_mode():
            result=self.head(a,b,torch.tensor([len(wt)],device=self.device,dtype=torch.float32))
        return float(result.item())
