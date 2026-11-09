"""Validate sequence-coordinate substitutions before any downloads or MSA searches."""
import re
from pathlib import Path

AMINO_ACIDS=set('ACDEFGHIKLMNPQRSTVWY')


def sequence(value):
    value=''.join(str(value).split()).upper()
    if not 20<=len(value)<=3719:
        raise ValueError('Provide a sequence of 20–3719 amino acids (JanusDDG architecture limits).')
    if set(value)-AMINO_ACIDS:
        raise ValueError('Sequences must contain only the 20 standard amino acids; no gaps or stop symbols.')
    return value


def read_fasta(path):
    names=[];parts=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):names.append(line[1:].strip())
        elif line.strip():parts.append(line.strip())
    if len(names)!=1:
        raise ValueError('Use a FASTA containing exactly one WT sequence, or use a CSV for multiple proteins.')
    return sequence(''.join(parts))


def variant(wt,mutations):
    wt=sequence(wt);changes=[];seen=set()
    for text in re.split(r'[_ ,;]+',str(mutations).strip().upper()):
        match=re.fullmatch(r'([ACDEFGHIKLMNPQRSTVWY])([1-9][0-9]*)([ACDEFGHIKLMNPQRSTVWY])',text)
        if not match:raise ValueError(f'Invalid mutation {text!r}; use K13A or K13A_L22V.')
        old,number,new=match.groups();position=int(number)-1
        if position>=len(wt):raise ValueError(f'{text}: position exceeds sequence length {len(wt)}.')
        if wt[position]!=old:raise ValueError(f'{text}: sequence residue {number} is {wt[position]}, not {old}.')
        if old==new:raise ValueError(f'{text}: source and destination residues are identical.')
        if position in seen:raise ValueError(f'{text}: a variant cannot change the same site twice.')
        seen.add(position);changes.append((position,new))
    changes.sort()
    if changes[-1][0]-changes[0][0]>=1022:
        raise ValueError('All substitutions in a variant must fit in one 1022-residue MSA window.')
    mutant=list(wt)
    for pos,aa in changes:mutant[pos]=aa
    return ''.join(mutant),'_'.join(f'{wt[p]}{p+1}{aa}' for p,aa in changes),changes


def window(length,changes):
    if length<=1022:return 0,length
    midpoint=(changes[0][0]+changes[-1][0])//2
    left=max(0,min(length-1022,midpoint-511))
    if not all(left<=p<left+1022 for p,_ in changes):
        left=max(0,min(length-1022,changes[0][0]))
    return left,left+1022
