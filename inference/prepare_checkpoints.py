"""Maintainer utility: export portable weights from the verified published checkpoints."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
import zipfile

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--msa-checkpoint',required=True,type=Path)
    parser.add_argument('--janus-checkpoint',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    import torch
    import numpy as np
    from inference.weights import sha256
    from inference.vendor import janus_model
    import __main__
    torch.set_num_threads(4)
    original_janus='989f454dc53f890923f5364147c3c6d54612630c85d30d60445921d5b5940810'
    original_msa='30c5534f2985e7026b3a714e0089879bcbd52f14d6e023a2fdded57af97476f2'
    assert sha256(args.janus_checkpoint)==original_janus
    with zipfile.ZipFile(args.msa_checkpoint) as archive:
        payload=archive.read('MSA_Finetuning.pt')
    assert hashlib.sha256(payload).hexdigest()==original_msa
    state=torch.load(io.BytesIO(payload),map_location='cpu',weights_only=True)['MODEL_STATE']
    encoder={key.removeprefix('msa_transformer.'):value for key,value in state.items() if key.startswith('msa_transformer.')}
    q=state['fc1.weight'].double().numpy().T@state['fc2.weight'][0].double().numpy()/2
    config=dict(embed_dim=768,ffn_embed_dim=3072,attention_heads=12,layers=12,dropout=.1,
                attention_dropout=.1,activation_dropout=0.,max_tokens=16384,max_tokens_per_msa=16384,
                max_positions=1024,embed_positions_msa=True,
                embed_positions_msa_dim=encoder['msa_position_embedding'].shape[-1])
    args.output.mkdir(parents=True,exist_ok=True)
    torch.save(dict(format_version=1,encoder_state=encoder,msa_args=config,potential_weights=torch.from_numpy(q),
                    source_sha256=original_msa,model='EvoPhi-DDG',license='CC-BY-4.0'),args.output/'evophi-ddg.pt')
    # Original release pickles classes under __main__; only this SHA-verified conversion uses pickle.
    for name in ['Cross_Attention_DDG','TransformerRegression','SinusoidalPositionalEncoding']:
        setattr(__main__,name,getattr(janus_model,name))
    original=torch.load(args.janus_checkpoint,map_location='cpu',weights_only=False).eval()
    converted=janus_model.Cross_Attention_DDG(janus_model.TransformerRegression)
    converted.load_state_dict(original.state_dict(),strict=True)
    converted.eval();original.base_ddg.speach_att_type=False
    with torch.inference_mode():
        rng=torch.Generator().manual_seed(42)
        a=torch.randn(1,79,1280,generator=rng);b=torch.randn(1,79,1280,generator=rng);length=torch.tensor([79.])
        torch.testing.assert_close(original(a,b,length),converted(a,b,length),rtol=0,atol=0)
    torch.save(dict(format_version=1,state_dict=original.state_dict(),source_sha256=original_janus,
                    model='JanusDDG',license='MIT'),args.output/'janusddg.pt')
    assets={}
    for name,filename in [('evophi','evophi-ddg.pt'),('janus','janusddg.pt')]:
        path=args.output/filename
        assets[name]=dict(filename=filename,bytes=path.stat().st_size,sha256=sha256(path))
    (args.output/'checkpoint_manifest.json').write_text(json.dumps(dict(format_version=1,assets=assets,
        sources=dict(msa_sha256=original_msa,janus_sha256=original_janus)),indent=2)+'\n')
    print(json.dumps(assets,indent=2))


if __name__=='__main__':main()
