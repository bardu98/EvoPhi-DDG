# Sources and licensing

## Source models

| Component | Source | Revision used in the original workspace |
| --- | --- | --- |
| MSAesm_ddG | [RitAreaSciencePark/PLM4Muts](https://github.com/RitAreaSciencePark/PLM4Muts) | `6bbbdeda629d3b62e2f26d7fd3231979a3a1784a` |
| JanusDDG | [compbiomed-unito/JanusDDG](https://github.com/compbiomed-unito/JanusDDG) | `7ef53fde72c9b9b715b0a9f2f80910ebf927bd59` |

The primary predictions use their published checkpoint weights. Checkpoint hashes and acquisition details are in [inference.md](inference.md). Ready-to-use converted checkpoints are supplied through repository releases for the [unified inference pipeline](../inference/README.md). The full upstream repositories, training datasets, virtual environments and temporary inference caches are not bundled here.

The MSAesm_ddG source uses AGPL-3.0-or-later notices. The analysis and EvoPhi-DDG code in this repository is released under the same license, supplied as `LICENSE`. JanusDDG is distributed under MIT; its upstream license is preserved in `docs/JanusDDG-LICENSE`. Third-party software retains its own license.

Benchmark measurements and distributed prediction tables retain their original provenance and any original data-license terms. The code license does not relicense third-party data or model weights. `generated/label_policy_provenance.json` identifies the JanusDDG distribution files used for labels, and `references.bib` provides the literature references. The files under `data/physics/` and the compressed S669 multisite files are derived prediction artifacts from this work.

## Benchmark provenance

- S669: Pancotti et al. (2022), evaluated with the experimental labels and sequence mapping distributed with JanusDDG.
- S461: curated S669 subset described by Hernández et al. (2023), evaluated using the curated JanusDDG-distributed labels and contexts.
- S96: DDGun web-server benchmark, Montanucci et al. (2022), using the JanusDDG-distributed labels and sequence mapping.
- PTmul-NR: multiple-substitution benchmark distributed with JanusDDG; the source bibliography and upstream papers document its construction.

Literature comparator values are distinct from predictions rerun in this project. Neither historical literature values nor alternative S669 labels are substituted into the matched empirical comparisons.

## Scope

The repository is an analysis and inference companion to the EvoPhi-DDG manuscript in preparation. It covers the four benchmark evaluations, component controls, classification, complementarity and physical-consistency tests. The supplied code and numerical outputs document these constructions and comparisons.
