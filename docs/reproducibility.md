# Reproduction guide

## What the reproduction command does

`python scripts/reproduce.py` runs the original paper analyses adapted to repository-relative paths. It needs NumPy, pandas, SciPy and Matplotlib, with the tested versions pinned in `requirements.txt`. It does not need PyTorch, a GPU, an MSA server or the original workspace.

The input scores and pathway predictions are frozen. This reproduces the statistics and plotted results from those scores, rather than claiming a new end-to-end execution of alignment searches and neural-network inference. The optional inference workflow is documented separately.

`python scripts/verify.py` checks SHA256 hashes for immutable inputs and compares regenerated outputs with `reference/`. Floating-point results are compared with absolute tolerance `1e-10` and relative tolerance `1e-9`. PDF byte identity is not required: plotting libraries can change metadata and rendering without changing numerical results.

## Figure sources

Names refer to files rather than manuscript figure numbers, so the links remain valid if the paper is reordered.

| Figure | Inputs | Generator |
| --- | --- | --- |
| `accuracy` | `data/s669.csv`, `s461.csv`, `s96.csv`, `ptmul_nr.csv` | `analysis/analysis.py` |
| `physics` | Benchmark scores and `data/transitivity.csv` | `analysis/analysis.py` |
| `multisite_transitivity` | `data/s669_multisite_msa_paths.csv.gz`, `s669_multisite_janus_paths.csv.gz` | `analysis/s669_multisite_paper.py` |
| `classification` | All four benchmark tables, applying the experimental exclusion band | `analysis/classification.py` |
| `complementarity` | Benchmark predictions and squared-error decomposition | `analysis/analysis.py` |
| `ensemble_ablation` | `data/ensemble_ablation_predictions.csv` and same-site edge/path arrays | `analysis/ensemble_ablation.py` |
| `S2450_bias_comparison` | `data/model_controls/S2450_predictions.csv` | `analysis/bias_comparison.py` |
| `accuracy_s96`, `classification_s96` | S96 inputs only; supplemental individual exports of the combined panels | `analysis/refresh_figures.py` |

`msaesm.png`, `MSA_JanusDDG.png` and `ensemble_schematic.tex` are static model diagrams. They are supplied rather than regenerated from numerical data.

## Table sources

| Table content | LaTeX output in `generated/` | Numerical source |
| --- | --- | --- |
| Benchmark overview | `benchmark_summary_rows.tex` | `statistics.json`, `results/regression.csv` |
| Detailed single-substitution benchmarks | `s669_rows.tex`, `s461_rows.tex`, `s96_rows.tex`, `benchmark_details.tex` | Benchmark CSVs and `statistics.json` |
| Multiple substitutions and alternative constructions | `ptmul_rows.tex`, `benchmark_details.tex` | `data/ptmul_nr.csv` |
| Classification | `classification_rows.tex` | `classification_statistics.json`, `results/classification.csv` |
| Same-site transitivity | `transitivity_rows.tex` | `data/physics/`, `data/transitivity_repeats.csv` |
| Multisite transitivity | `multisite_rows.tex` | Compressed S669 path-level files |
| Original-head and antisymmetry controls | `ensemble_ablation_rows.tex` | `ensemble_ablation_statistics.json` |
| Paired Pearson differences | `bootstrap_rows.tex` | `statistics.json`, `results/bootstrap.csv` |
| Literature context and recalculated S669 scores | `literature_rows.tex` | `data/literature_comparisons.csv`, `data/s669.csv` |

The standalone manuscript text and its revision markup are not included. Historical literature values are transcribed separately in `data/literature_comparisons.csv` with citation keys. The generated literature table retains those values and recalculates only the rows from this study; it does not present the historical predictors as matched reruns.

## Statistical settings

- Regression uses Pearson, Spearman, mean absolute error and root mean squared error on all direct variants: S669 669, S461 461, S96 96 and PTmul-NR 82.
- Bootstrap groups are complete WT sequences: 87, 46, 14 and 13 groups, respectively. Entire groups are resampled together, with repeated experimental measurements retained.
- Regression comparisons use 5,000 paired protein-bootstrap draws, NumPy `default_rng(42)`, with percentile intervals at 2.5% and 97.5%. Differences are computed before rounding.
- Classification retains only experimental `abs(y) > 0.50` kcal/mol. Stabilizing labels have `y > 0.50`; predicted stabilizing means score `> 0`. No predictions are discarded according to magnitude. Classification uncertainty uses 1,000 paired protein-bootstrap draws with seed 42.
- The classification retained counts are recorded in `generated/classification_statistics.json`; individual-model and ensemble scores use exactly the same retained records within each benchmark.
- S669 same-site paths use 1, 3, 5, 7 or 9 intermediates and ten repetitions. The frozen draws are in `data/physics/paths.npz`; their reuse avoids any dependency on random-number consumption order.
- S669 multisite paths use 2, 4, 6, 8 or 10 intermediates, ten repetitions and seed 42. Each target mutation is surrounded by temporary substitutions at distinct other sites and their reversions. The protocol contains 92 query/window contexts, 196,196 states, 229,558 directed edges and 33,450 paths.
- Multisite uncertainty uses 5,000 paired bootstrap draws of complete WT groups, preserving all ten repetitions. The reported closure RMSE is the mean of the ten repetition-wise RMSEs, not the RMSE pooled over all repetitions.
- Physical consistency uses predictions only; no experimental intermediate-state labels enter closure metrics.

These intervals describe the supplied datasets and procedures. They do not account for the full history of exploratory benchmark inspection.

## Display precision

Calculations use unrounded source values. Standard displayed metrics have two decimals. Small nonzero errors are shown in scientific notation with a two-decimal mantissa. Near-unit physical correlations are written as `1 − epsilon` when rounding would hide a nonzero deficit. Counts remain integers. A floating-point defect close to machine precision is distinguished from a mathematical guarantee.

## Automated checks

```bash
python -m unittest discover -s tests
python scripts/reproduce.py
python scripts/verify.py
```

The tests check the equivalence between the zero-bias antisymmetric ReLU head and its linear potential, same-site closure, reversal of the multiple-substitution construction and equal-weight combination. GitHub Actions runs these checks and uploads the regenerated artifacts.
