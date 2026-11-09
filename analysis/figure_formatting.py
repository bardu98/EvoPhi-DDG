"""Display-only model names and two-decimal manuscript figure formatting."""
import math

from matplotlib.ticker import FuncFormatter, MaxNLocator

MSA_ORIGINAL='MSAesm_ddG (orig.)'
MSA_ANTI='MSAesm_ddG (AS)'
LABELS={
    'MSAesm original':MSA_ORIGINAL,'MSA original':MSA_ORIGINAL,
    'MSAesm anti':MSA_ANTI,'MSA anti':MSA_ANTI,'Antisymmetrized':MSA_ANTI,
    'JanusDDG':'JanusDDG','Ensemble':'Ensemble',
    'EvoPhi-DDG':'EvoPhi-DDG','MSAesm_anti':MSA_ANTI,
    'Mean EvoPhi':'Ensemble (EvoPhi-DDG)',
    'Mean original MSA':'Mean (MSAesm_ddG orig.)',
    'Mean anti':'Mean (MSAesm_ddG AS)',
}


def correlation(value):
    if 1e-14<1-value<.005:
        mantissa,exponent=f'{1-value:.2e}'.split('e')
        return rf'$1-{mantissa}\times10^{{{int(exponent)}}}$'
    return f'{value:.2f}'


def label(name):
    return LABELS.get(name,name)


def two_decimals(value,_position=None):
    if abs(value)<.005:value=0.
    return f'{value:.2f}'


def scientific_two_decimals(value,_position=None):
    if value<=0:return ''
    exponent=math.floor(math.log10(value)+1e-12)
    mantissa=value/10**exponent
    return rf'${mantissa:.2f}\times10^{{{exponent}}}$'


def format_axes(fig,skip=()):
    """Format continuous axes; explicitly preserve count/categorical axes.

    Log-axis labels retain small errors with two-decimal scientific mantissas.
    Numerical source values, plotted positions and statistical results are unchanged.
    """
    skip=set(skip)
    for ax in fig.axes:
        for name,axis in [('x',ax.xaxis),('y',ax.yaxis)]:
            if (ax,name) in skip:continue
            if axis.get_scale()=='symlog':
                axis.set_major_formatter(FuncFormatter(lambda value,_: scientific_two_decimals(value) if 0<value<.01 else f'{value:.2f}'))
            elif axis.get_scale()=='log':
                axis.set_major_formatter(FuncFormatter(scientific_two_decimals))
                axis.set_minor_formatter(FuncFormatter(lambda *_:''))
            else:
                axis.set_major_formatter(FuncFormatter(two_decimals))
                # Longer two-decimal ticks need fewer columns in small panels.
                if name=='x':axis.set_major_locator(MaxNLocator(nbins=3))
