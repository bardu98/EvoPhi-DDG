# Inference and alignment preparation

## Published weights

The reproduction workflow uses frozen scores, so it does not need model weights. To run inference on new inputs, obtain the original MSAesm_ddG `MSA_Finetuning.pt` checkpoint from [Zenodo record 14026821](https://zenodo.org/records/14026821). The original checkpoint used for this work has SHA256:

```text
30c5534f2985e7026b3a714e0089879bcbd52f14d6e023a2fdded57af97476f2
```

This hash refers to the `MSA_Finetuning.pt` archive member, not the ZIP bytes. The published `MSA_Finetuning.zip` has MD5 `f0d56f503029ef7a1dcea58df4923fa6`. The local archive filename is immaterial.

The JanusDDG published final checkpoint is `models/JanusDDG_fine_tuned.pth` in its upstream distribution. The SHA256 used here is:

```text
989f454dc53f890923f5364147c3c6d54612630c85d30d60445921d5b5940810
```

No locally retrained or fitted checkpoint is used in the primary benchmark results.

## Alignment procedure

The benchmark MSAs were obtained through the public ColabFold MMseqs2 MSA service using `/ticket/msa`, `mode=all`, and the unpaired UniRef result. Environmental hits were not added to the MSA used here.

The ColabFold server configuration and database history identify this backend as UniRef30_2302. However, the saved job archives do not contain the runtime database path or database hash. Thus the exact backend release is supported by server documentation, rather than independently attested in each archived job. See the [server configuration](https://github.com/sokrypton/ColabFold/blob/main/MsaServer/config.json) and [database history](https://github.com/sokrypton/ColabFold/wiki/MSA-Server-Database-History). Do not interpret UniRef30 as a minimum query identity threshold.

We filter the resulting MSA with HHfilter from HH-suite 3.3 using:

```bash
hhfilter -i query.a3m -o query.filtered.a3m \
  -id 90 -diff 1000 -cov 0 -qid 0 -qsc -20 -maxseq 65535
```

The original pipeline then used `reformat.pl -r` to remove insertion columns. The inference script can also remove A3M lowercase insertions directly, preserving aligned row order. HHfilter reduces redundancy; it does not perform the homology search. `-diff 1000` is not a hard cap of 1,000 rows, and `-id 90` is an upper pairwise redundancy threshold, not a lower query identity threshold.

After filtering, retain the first rows in their existing order, including the query. With query length `L`, use:

```text
depth = min(number_of_filtered_rows, floor(16000 / (L + 1)))
```

The query must be first. The model uses at most 1,022 biological residues per window. For long benchmark proteins, the original mutation-centered windows and coordinates are recorded in `data/physics/variants.csv`; JanusDDG uses its full-sequence context. Matching the original context matters for reproducing scores. A new online search can return different homologs as server databases change.

## Run EvoPhi-DDG on prepared inputs

Install the optional dependencies, selecting an appropriate PyTorch CUDA build for your machine if required:

```bash
pip install -r requirements-inference.txt
```

Prepare a CSV containing `code`, `wt_seq`, `mut_seq` and `msa`. MSA paths can be absolute or relative to the CSV. The first MSA row must exactly match `wt_seq`. The script accepts equal-length query windows of 2–1,022 canonical amino acids and detects every changed site. It supports both single and multiple substitutions.

```bash
python scripts/predict.py variants.csv \
  --checkpoint weights/MSA_Finetuning.pt \
  --output predictions.csv \
  --device cuda
```

The checkpoint can also be the ZIP containing `MSA_Finetuning.pt`. Fair-esm normally downloads the base MSA Transformer before loading the complete fine-tuned encoder state. For offline use, supply the base checkpoint with `--base-encoder weights/esm_msa1b_t12_100M_UR50S.pt`; its matching contact-regression checkpoint must be in the same directory as expected by fair-esm.

The encoder runs in evaluation mode with frozen weights, layer 12, FP32 and TF32 disabled. Potential weights and readout arithmetic use float64. The mean representation preserves the published `[1:-1]` pooling convention, including omission of the last biological residue. Homolog rows remain fixed; only the query changes between states. This script is an independent implementation of the EvoPhi readout and endpoint projection, not a copy of the original training program.

If the input CSV also contains `janus_score`, the output includes the equal-weight ensemble. Those scores must refer to the same variants, sign convention and intended sequence contexts. Use the upstream JanusDDG inference program and published checkpoint to obtain them; this repository does not replace its ESM2 encoder and architecture.

This command performs inference from supplied MSAs. It does not submit new alignment jobs or automatically crop proteins. The paper's cached-score reproduction remains the checked route for rebuilding all published numerical artifacts.

## Supplied inference example

`examples/variants.csv` and `examples/alignment.fasta` contain one real S96 input with its filtered MSA. Run the command above with `examples/variants.csv` as the input. `examples/expected_predictions.csv` records a CPU check against the published checkpoint. The score is approximately −0.75 kcal/mol. Its difference from the frozen benchmark score was `1.62e-06` kcal/mol, consistent with FP32 encoder arithmetic across devices; do not require byte-identical inference results across hardware.

`examples/multiple_variants.csv` and `examples/multiple.a3m` provide a PTmul-NR example with four substitutions. Its CPU score is approximately 0.96 kcal/mol, differing from the frozen score by `9.29e-06` kcal/mol. The expected output is in `examples/expected_multiple_predictions.csv`, and both checks are recorded in `generated/inference_checks.json`.
