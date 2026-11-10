"""Parse biological inputs and explain row errors before starting a job."""
import io
import re
import pandas as pd
from inference.inputs import sequence,variant


def fasta_records(text):
    result={};identifier=None;parts=[]
    def finish():
        if identifier in result:raise ValueError(f'The FASTA identifier {identifier!r} occurs more than once.')
        result[identifier]=sequence(''.join(parts))
    for line in text.splitlines():
        line=line.strip()
        if not line:continue
        if line.startswith('>'):
            if identifier is not None:finish()
            header=line[1:].strip()
            if not header:raise ValueError('Each FASTA record needs a name after >.')
            identifier=header.split()[0];parts=[]
        else:
            if identifier is None:raise ValueError('A FASTA file must start with > followed by a protein identifier.')
            parts.append(line)
    if identifier is not None:finish()
    if not result:raise ValueError('The FASTA file contains no proteins.')
    return result


def table(text):
    try:
        frame=pd.read_csv(io.StringIO(text),sep=None,engine='python',keep_default_na=False)
    except Exception as error:raise ValueError('The variant table could not be read. Use a CSV with a header row.') from error
    frame.columns=[str(col).strip() for col in frame.columns]
    if not frame.columns.is_unique:raise ValueError('The variant table contains duplicate column names.')
    if frame.empty:raise ValueError('The variant table contains no variants.')
    return frame


def single_input(text,mutations):
    if text.lstrip().startswith('>'):
        records=fasta_records(text)
        if len(records)!=1:raise ValueError('For several FASTA records, choose FASTA + variant table.')
        name,wt=next(iter(records.items()))
    else:name,wt='protein',sequence(text)
    variants=[item.strip() for item in re.split(r'[\n;]+',mutations) if item.strip()]
    if not variants:raise ValueError('Enter at least one variant, for example K13A.')
    return validate(pd.DataFrame([dict(ID=name,Sequence=wt,MTS=item) for item in variants]))


def batch_input(csv_text,fasta_text=None):
    frame=table(csv_text)
    if fasta_text is not None:
        if not {'ID','MTS'}.issubset(frame):raise ValueError('With a FASTA file, the table needs ID and MTS columns.')
        records=fasta_records(fasta_text)
        missing=set(frame.ID.astype(str))-set(records)
        if missing:raise ValueError('These table IDs are absent from the FASTA: '+', '.join(sorted(missing)[:8]))
        frame['Sequence']=[records[str(identifier)] for identifier in frame.ID]
    if not {'Sequence','MTS'}.issubset(frame):raise ValueError('The CSV needs Sequence and MTS columns. ID is optional.')
    return validate(frame)


def validate(frame):
    rows=[]
    for index,row in frame.reset_index(drop=True).iterrows():
        name=str(row.get('ID',f'protein_{index+1}')).strip() or f'protein_{index+1}'
        try:
            wt=sequence(row.Sequence);_,mts,_=variant(wt,row.MTS)
        except ValueError as error:raise ValueError(f'Row {index+1} ({name}): {error}') from error
        rows.append(dict(ID=name,Sequence=wt,MTS=mts))
    if not rows:raise ValueError('No variants were supplied.')
    return pd.DataFrame(rows)
