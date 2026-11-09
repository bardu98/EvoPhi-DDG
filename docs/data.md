# Data files

## Experimental labels and exact inputs

`data/experimental_labels.csv` is the frozen experimental reference. It contains `dataset`, unique internal `code`, experimental `y` in kcal/mol, and the original JanusDDG source file, row and variant identifier. Repeated S669 measurements are retained. The same mapping is applied to every predictor through `analysis/label_policy.py`.

`data/benchmark_inputs.csv` contains the 1,308 benchmark records across S669, S461, S96 and PTmul-NR, with `protein`, `mutations`, `wt_seq`, `mut_seq` and `y`. Positions in `mutations` are one-based positions in the supplied sequence, rather than PDB residue numbering. The WT residue and experimental value were checked against each original source row when this file was assembled. S461 is represented with its own curated labels and contexts; do not replace these by the labels for overlapping S669 entries.

## Prediction tables

`data/s669.csv`, `s461.csv`, `s96.csv` and `ptmul_nr.csv` share these fields:

| Column | Meaning |
| --- | --- |
| `code` | Unique variant key used for alignment and label checks |
| `group` | Complete WT-sequence group for paired protein bootstrap |
| `y` | Experimental value from JanusDDG |
| `J` | Published-checkpoint JanusDDG score |
| `M` | Primary EvoPhi-DDG score |
| `M_anti` | Antisymmetrized MSAesm_ddG with published nonzero hidden biases |
| `M_original` | Original directional MSAesm_ddG score, when present |
| `M_reverse` | Original reverse-direction score, when present |
| `retained` | Membership in the specified training-similarity-filtered sensitivity subset |

`y_msa_reference`, where present, is an audit column only and is never used as an accuracy target. PTmul-NR also contains the sequential and all-order constructions and their published-hidden-bias controls. `M_previous` denotes the forward endpoint-star construction used in earlier controls; it is not the primary multiple-substitution score.

The ensemble is calculated from `J` and `M` rather than relying on rounded values. The MSA-derived primary scores are independently checked against `data/model_controls/` before analysis.

## Pathway data

`data/physics/` contains same-site directed-edge arrays, the site-to-benchmark mapping and frozen random intermediate residues. Array dimensions are site context, source amino acid and destination amino acid, using the order `ACDEFGHIKLMNPQRSTVWY`. The EvoPhi array is the potential difference computed from the frozen embeddings and head weights. These are edge-level scores, not the original high-dimensional encoder embeddings.

`data/s669_multisite_msa_paths.csv.gz` and `s669_multisite_janus_paths.csv.gz` contain direct scores, pathway sums and signed closure errors for the S669 multisite protocol. The ensemble path values are checked against the mean of component values. The signed closure error is `path_sum − direct`.

The older `multisite_*` PTmul graph files describe the exploratory three-edge detour analysis. They are retained for provenance and are distinct from the primary S669 multisite results shown in the paper. They are not substituted for the S669 files by the reproduction command.

## Provenance and checksums

- `data/checksums.json`: SHA256 hashes for immutable reproduction inputs.
- `generated/label_policy_provenance.json`: original JanusDDG file hashes and mapping checks.
- `generated/model_policy_provenance.json`: zero-hidden-bias prediction provenance.
- `generated/s461_provenance.json`, `s461_mapping.csv`: context and curation checks.
- `generated/s96_provenance.json`, `s96_s669_overlap.csv`: S96 alignment depth and benchmark overlap audit.
- `data/s669_multisite_protocol.json`, `s669_multisite_janus_verification.json`: multisite protocol and checkpoint verification.

Derived outputs, such as classification-selected records and repeat-level summaries, are rewritten when analyses are rerun. They are compared numerically to the reference results rather than treated as immutable inputs.
