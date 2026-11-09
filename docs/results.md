# Benchmark results

All methods use the frozen experimental labels distributed with JanusDDG. The ensemble is the equal-weight mean of JanusDDG and EvoPhi-DDG. Errors are in kcal/mol. Values below use two decimals; CSV and JSON files retain full precision.

| Benchmark | Variants | Model | Pearson | Spearman | MAE | RMSE |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| S669 | 669 | JanusDDG | 0.55 | 0.55 | 0.97 | 1.39 |
| S669 | 669 | EvoPhi-DDG | 0.53 | 0.55 | 0.99 | 1.39 |
| S669 | 669 | Ensemble | 0.59 | 0.61 | 0.91 | 1.33 |
| S461 | 461 | JanusDDG | 0.69 | 0.66 | 0.71 | 0.96 |
| S461 | 461 | EvoPhi-DDG | 0.66 | 0.65 | 0.73 | 0.97 |
| S461 | 461 | Ensemble | 0.75 | 0.72 | 0.65 | 0.88 |
| S96 | 96 | JanusDDG | 0.52 | 0.50 | 1.50 | 2.10 |
| S96 | 96 | EvoPhi-DDG | 0.48 | 0.47 | 1.53 | 2.14 |
| S96 | 96 | Ensemble | 0.55 | 0.52 | 1.45 | 2.04 |
| PTmul-NR | 82 | JanusDDG | 0.61 | 0.57 | 1.57 | 2.06 |
| PTmul-NR | 82 | EvoPhi-DDG | 0.66 | 0.58 | 1.45 | 1.95 |
| PTmul-NR | 82 | Ensemble | 0.66 | 0.61 | 1.47 | 1.95 |

S461 is a curation sensitivity analysis of S669, rather than an independent test. S96 also shares biological variants with S669; see the overlap audit in `generated/`. The ensemble and EvoPhi-DDG both round to Pearson 0.66 on PTmul-NR; the ensemble does not exceed EvoPhi-DDG in unrounded Pearson on the full set.

## Numerical source files

- `results/regression.csv`: individual models, ensemble and multiple-substitution sensitivity analyses.
- `results/classification.csv`: individual classifiers and ensemble after the experimental exclusion band.
- `results/bootstrap.csv`: paired differences in Pearson and protein-bootstrap intervals.
- `results/same_site_transitivity.csv`: same-site closure summaries.
- `results/multisite_transitivity.csv`: S669 multisite closure summaries.
- `results/ensemble_controls.csv`: original-head, projected-head and zero-bias controls.
- `generated/`: LaTeX tables and complete statistical outputs, including additional intervals.
