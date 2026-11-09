"""Generate benchmark tables using the frozen JanusDDG experimental labels."""
import json

import numpy as np
import pandas as pd

from analysis import HERE, dataset_frames, score
from table_formatting import performance_rows


def metric_values(frame, prediction):
    result = score(frame.y, prediction)
    return [result[key] for key in ('pearson', 'spearman', 'mae', 'rmse')]


def table_block(name, frame, models):
    rows = [rf'\multicolumn{{5}}{{l}}{{\textbf{{{name}}} ($n={len(frame)}$)}} \\']
    labels = [label for label, _ in models]
    values = [metric_values(frame, prediction) for _, prediction in models]
    rows.extend(performance_rows(labels, values, ['.2f'] * 4,
                                 ['max', 'max', 'min', 'min']))
    return '\n'.join(rows)


def main():
    frames = dataset_frames(pd.read_csv(HERE / 'data/s669.csv'),
                            pd.read_csv(HERE / 'data/ptmul_nr.csv'))
    reference = json.loads((HERE / 'generated/statistics.json').read_text())
    summary = []
    details = []
    for name, frame in frames:
        models = [('JanusDDG', frame.J),
                  ('EvoPhi-DDG', frame.M),
                  ('Equal-weight ensemble', (frame.J + frame.M) / 2)]
        for key, (_, prediction) in zip(('J', 'M', 'E'), models):
            result = score(frame.y, prediction)
            for metric in ('pearson', 'spearman', 'mae', 'rmse'):
                if not np.isclose(result[metric], reference[name][key][metric],
                                  rtol=0, atol=1e-10):
                    raise ValueError(f'Metric mismatch: {name}, {key}, {metric}')
        summary.append(table_block(name, frame, models))
        if name != 'PTmul-NR':
            models.insert(1, (r'MSAesm\_ddG original', frame.M_original))
            models.insert(2, (r'MSAesm\_ddG antisymmetrized', frame.M_anti))
            body = table_block(name, frame, models)
        else:
            # Preserve the published-bias control and construction sensitivities.
            models.insert(1, (r'MSAesm\_ddG antisymmetrized', frame.M_anti))
            models.extend([
                (r'EvoPhi-DDG: sequential, both directions', frame.M_sequence),
                (r'EvoPhi-DDG: all mutation orders', frame.M_orders)])
            values = []
            for _, prediction in models:
                full = score(frame.y, prediction)
                subset = score(frame.loc[frame.retained, 'y'], prediction[frame.retained])
                values.append([full['pearson'], full['mae'], full['rmse'], subset['pearson']])
            body = '\n'.join(performance_rows(
                [label for label, _ in models], values, ['.2f'] * 4,
                ['max', 'min', 'min', 'max']))
        header = (r'Predictor & Pearson & Spearman & MAE & RMSE'
                  if name != 'PTmul-NR' else
                  r'Predictor & Pearson$_{82}$ & MAE$_{82}$ & RMSE$_{82}$ & Pearson$_{67}$')
        label = {'S669': 's669', 'S461': 's461', 'S96': 's96', 'PTmul-NR': 'ptmul'}[name]
        caption = (rf'\textbf{{Detailed evaluation on {name}.}} '
                   r'All predictors use the experimental labels distributed with JanusDDG. '
                   r'MAE and RMSE are in kcal/mol. Bold marks the best values and associated methods, including ties at two decimals.')
        if name == 'PTmul-NR':
            caption += (' The full set contains 82 variants; the training-restricted subset contains 67. '
                        'The primary EvoPhi-DDG prediction uses zero hidden biases and component-and-whole antisymmetrization; '
                        'the other constructions are sensitivity analyses.')
        details.append('\n'.join([
            r'\noindent\begin{minipage}{\linewidth}',
            r'\section*{Supplementary Tables}' if name == 'S669' else '',
            r'\centering\small',
            r'\captionof{table}{' + caption + '}', rf'\label{{tab:{label}}}',
            r'\begin{tabular}{lrrrr}\toprule', header + r'\\\midrule',
            body, r'\bottomrule\end{tabular}', r'\end{minipage}\par']))
    header = (r'\begin{tabular}{lrrrr}\toprule'
              '\n' + r'Predictor & Pearson & Spearman & MAE & RMSE\\\midrule')
    (HERE / 'generated/benchmark_summary_rows.tex').write_text(
        header + '\n' + '\n\\midrule\n'.join(summary) +
        '\n' + r'\bottomrule\end{tabular}' + '\n')
    # Keep each caption and table together, below the supplementary heading.
    # Two tables per page avoid float reordering and stretched float spacing.
    (HERE / 'generated/benchmark_details.tex').write_text(
        '\n\\bigskip\n'.join(details[:2]) + '\n\\clearpage\n' +
        '\n\\bigskip\n'.join(details[2:]) + '\n')
    print('Generated summary and detailed tables; all primary metrics match statistics.json.')


if __name__ == '__main__':
    main()
