# Predict stability changes from your sequence

This workflow starts from a WT protein sequence and one or more variants. It finds homologs, filters the alignment, evaluates EvoPhi-DDG and JanusDDG, and writes their predictions and equal-weight ensemble to a CSV. Neither model is trained or calibrated on your inputs.

The ready-to-use checkpoints are supplied as [repository release assets](https://github.com/bardu98/EvoPhi-DDG/releases/tag/checkpoints-v1), rather than large files in the Git history. The program downloads them automatically on first use and verifies their SHA256 hashes. The EvoPhi-DDG file already contains the fine-tuned MSA encoder and potential readout. The original ESM2 backbone needed by JanusDDG is downloaded automatically from fair-esm.

## Install once

Use Linux or macOS with Conda/Miniforge installed. On Windows, use WSL2. From the repository directory:

```bash
conda env create -f inference/environment.yml
conda activate evophi-ddg
```

This installs the Python dependencies, HHfilter and Perl. A CUDA-capable NVIDIA GPU is recommended; CPU inference is supported but slower. The first prediction downloads approximately 0.49 GB of model checkpoints plus approximately 2.6 GB for ESM2. Allow extra disk space for the Python environment and its GPU libraries. Subsequent predictions reuse these downloads.

While this repository is private, downloading release assets requires access to it. If needed, authenticate GitHub CLI once with `gh auth login`. Once the repository is public, checkpoint downloads do not require a GitHub account. For offline use, put the two release files into `inference/weights/` and use an existing alignment and local ESM2 checkpoint as described below.

## Predict a single substitution

Save your complete WT sequence in a FASTA file, for example `my_protein.fasta`. Mutation positions start at **1 in that exact sequence**. `K13A` means replacing the lysine at position 13 with alanine.

```bash
python inference/run.py --fasta my_protein.fasta --mutations K13A --output predictions.csv
```

To try the included real example:

```bash
python inference/run.py --fasta inference/examples/protein.fasta --mutations K13A --output predictions.csv
```

The program automatically chooses a GPU if available. Add `--device cpu` to request CPU inference explicitly.

## Predict multiple substitutions and multiple variants

Within one variant, join substitutions with underscores. Separate distinct variants with spaces:

```bash
python inference/run.py --fasta my_protein.fasta --mutations K13A L22V K13A_L22V --output predictions.csv
```

This example requests three predictions: each single substitution and the double substitution. Replace the example mutations with changes that match your own sequence. The program checks WT residues and positions before starting any alignment search.

EvoPhi-DDG uses the primary endpoint-star and whole-variant antisymmetric construction for multiple substitutions. JanusDDG evaluates the complete WT and fully mutated sequence directly. The program does not replace the primary EvoPhi-DDG construction by an arbitrary sequential mutation order.

## Predict variants from a CSV

For multiple proteins or a longer variant list, create a CSV with these columns:

| Column | Content |
| --- | --- |
| `ID` | Your protein or variant identifier; optional |
| `Sequence` | Complete WT amino-acid sequence |
| `MTS` | Substitution, such as `K13A`, or multiple substitutions, such as `K13A_L22V` |

```bash
python inference/run.py --input my_variants.csv --output predictions.csv
```

`inference/examples/variants.csv` provides a working single- and multiple-substitution example. No experimental labels are needed.

## Read the results

The output includes:

| Column | Meaning |
| --- | --- |
| `ID`, `Sequence`, `MTS` | The validated input |
| `DDG_EvoPhi_DDG` | EvoPhi-DDG prediction |
| `DDG_JanusDDG` | JanusDDG prediction |
| `DDG_Ensemble` | Equal-weight mean of the two predictions |
| `msa_depth` | Alignment rows used by EvoPhi-DDG, including the query |
| `msa_window_start`, `msa_window_end` | One-based inclusive MSA query window in the complete sequence |

All stability changes are in **kcal/mol**, with **positive values indicating stabilization**. The console shows two decimals; the CSV retains full precision. A companion `predictions.provenance.json` records model hashes, alignment provenance, parameters and device.

## Alignment search, filtering and caching

The default uses the public ColabFold MMseqs2 service, requests `mode=all`, and uses only its unpaired UniRef result. Sequences are sent to that external service when a new search is needed. You can supply your own alignment instead when you want to keep a sequence local.

Filtering uses HHfilter 3.3.0 with `-id 90 -diff 1000 -cov 0 -qid 0 -qsc -20 -maxseq 65535`, followed by the original `reformat.pl -r` conversion. The model keeps the first post-filtering rows within 16,000 tokens. Query pooling preserves the published `[1:-1]` convention. Only query sequences change between mutant states; homologs, their order, retained depth and window remain fixed.

Alignments and accepted search tickets are cached under `.cache/inference/`. If a server search times out, rerun the same command: it reuses the accepted ticket instead of submitting a duplicate job. Server queues determine alignment-search time. Changing the cache location with `--cache` creates a separate cache. Searches and filtering are done on the complete WT sequence; the fixed MSA window is selected afterwards for model inference.

Supported inputs contain 20–3,719 standard amino acids. EvoPhi-DDG uses at most 1,022 residues, while JanusDDG retains the complete sequence. For long proteins, the MSA window is centered on the altered site or mutation span; all changed sites must fit in the same window. Wider variants are rejected rather than silently dropping substitutions. Positions are sequence positions, not PDB residue numbers. Insertions, deletions and ambiguous amino acids are not supported.

New online searches can differ from the paper's frozen alignments as databases evolve. The exact database release is not attested by individual job archives; see [database provenance](../docs/inference.md).

## Use existing alignments or run offline

To supply a raw A3M/FASTA alignment and still apply the standard HHfilter protocol:

```bash
python inference/run.py --fasta my_protein.fasta --mutations K13A --msa my_alignment.a3m --output predictions.csv
```

The first row must match the complete WT sequence. If the alignment has already been filtered with the intended protocol, add `--msa-filtered`. For CSV inputs, an optional `MSA` column can specify a separate alignment path per protein, relative to the CSV directory.

For a complete offline example, with downloaded weights and the standard ESM2 checkpoint available locally:

```bash
python inference/run.py --fasta inference/examples/protein.fasta --mutations K13A \
  --msa examples/alignment.fasta --msa-filtered \
  --weights inference/weights \
  --esm2-checkpoint /path/to/esm2_t33_650M_UR50D.pt \
  --output predictions.csv
```

Fair-esm also expects the matching ESM2 contact-regression file beside the local backbone checkpoint. There is no separate base MSA Transformer download: the EvoPhi-DDG checkpoint already includes its fine-tuned encoder.

## Verification and source attribution

The exported JanusDDG state dictionary was checked against the SHA-verified original pickled model. The unified inference was also checked on real S96 and PTmul-NR inputs against the frozen benchmark predictions. Small FP32 differences across hardware are expected. The supplied filter reproduced the archived S96 alignment exactly, and the complete online search-to-prediction route was exercised on a real query.

```bash
python -m unittest discover -s inference/tests
```

`prepare_checkpoints.py` documents the checkpoint conversion; ordinary users do not need to run it. `weights.json` records checkpoint hashes, source hashes and release download details. Source code from JanusDDG retains its MIT license. The original MSA weights and the derived EvoPhi-DDG checkpoint are CC BY 4.0, with attribution to Cuturello et al. (2024); the JanusDDG source and checkpoint are attributed to Barducci et al. (2026). The original GPL reformat script and the relevant license texts are preserved in `vendor/`. See the main [source notes](../docs/sources.md) and bibliography.
