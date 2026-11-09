"""Consistent manuscript emphasis; metric direction and displayed ties matter."""
import numpy as np


def bold(text):
    return r'\textbf{'+(r'\boldmath ' if '$' in str(text) else '')+str(text)+'}'


def formatted(value,fmt='.2f'):
    if fmt=='d':return str(int(value))
    if fmt=='corr' and 1e-14<1-value<.005:
        mantissa,exponent=f'{1-value:.2e}'.split('e')
        return rf'$1-{mantissa}\times10^{{{int(exponent)}}}$'
    if 0<abs(value)<.005:
        mantissa,exponent=f'{value:.2e}'.split('e')
        if fmt.startswith('+') and value>0:mantissa='+'+mantissa
        return rf'${mantissa}\times10^{{{int(exponent)}}}$'
    return format(value,'+.2f' if fmt.startswith('+') else '.2f')


def displayed_number(value,fmt):
    if fmt=='d':return int(value)
    if fmt=='corr' and 1e-14<1-value<.005:return 1-float(f'{1-value:.2e}')
    return float(f'{value:.2e}') if 0<abs(value)<.005 else float(f'{value:.2f}')


def best_cells(values,formats,directions):
    """Return formatted cells and row winners; None marks descriptive columns.

    Ties at the precision actually printed in the table receive equal emphasis.
    """
    values=np.asarray(values,dtype=float)
    cells=[[formatted(x,fmt)
            for x,fmt in zip(row,formats)] for row in values]
    winning=np.zeros(values.shape,dtype=bool)
    for col,direction in enumerate(directions):
        if direction is None:continue
        rounded=np.array([displayed_number(row[col],formats[col]) for row in values])
        target=rounded.max() if direction=='max' else rounded.min()
        winning[:,col]=rounded==target
    return [[bold(cell) if winning[i,j] else cell for j,cell in enumerate(row)]
            for i,row in enumerate(cells)],winning.any(axis=1)


def performance_rows(labels,values,formats,directions):
    cells,winning=best_cells(values,formats,directions)
    return [' & '.join([bold(label) if win else label]+row)+r' \\'
            for label,row,win in zip(labels,cells,winning)]
