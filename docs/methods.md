# Model construction

## From MSA representations to a state vector

The frozen, published fine-tuned MSA Transformer receives the query and its homologous aligned sequences. We use layer 12 and the representations of the query row only. Each residue has 768 features.

At mutation position `p`, the state vector is the concatenation of the query residue representation at `p` and the mean query representation, giving 1,536 features. The published implementation pools the token slice `[1:-1]`. Its alphabet has a beginning-of-sequence token and no end-of-sequence token, so this slice omits the final biological residue. The reproduction preserves this convention; changing it would change the model inputs.

For a fixed site and homolog context `C`, write the vector as `phi_p,C(A)`. Only the query changes between sequence states. The homolog rows, their order, retained depth and window are fixed.

## Original head and antisymmetric control

The original head acts on the difference between two state vectors:

```text
x = phi(A) - phi(B)
M(A,B) = vᵀ ReLU(Wx + b) + c
M_AS(A,B) = [M(A,B) - M(B,A)] / 2
```

`b` is a vector of hidden-unit biases; `c` is the scalar output bias. The projection cancels `c` and enforces antisymmetry, but nonzero `b` generally prevents exact same-site transitivity. This published-hidden-bias control is retained in the supplied data as `M_anti`.

## EvoPhi-DDG

Setting the hidden biases `b` to zero yields

```text
EvoPhi(A,B) = vᵀ [ReLU(Wx) - ReLU(-Wx)] / 2
            = vᵀ W [phi(A) - phi(B)] / 2
            = U(A) - U(B).
```

The identity `ReLU(z) - ReLU(-z) = z` makes this a potential difference. No parameters are fitted. `W`, `v`, encoder weights and alignment contexts are retained from the published model. Antisymmetry and same-site transitivity follow algebraically for a fixed site and context. Numerical residuals can remain because of floating-point arithmetic.

The potential depends on the site and context. It is not an identified absolute folding free energy and does not establish global transitivity across different mutation sites.

## Multiple substitutions

Let `A` be the complete WT and `B` the sequence containing all `k` substitutions. Let `Ai` contain only substitution `i` introduced into `A`, and let `Bminus_i` restore only site `i` in `B`.

```text
F(A,B) = sum_i EvoPhi(A,Ai)
F(B,A) = sum_i EvoPhi(B,Bminus_i)
P(A,B) = [F(A,B) - F(B,A)] / 2
```

All terms in the first sum start from the same WT background; all terms in the second start from the same fully mutated background. `P(B,A) = -P(A,B)` exactly in exact arithmetic. For one substitution, the construction reduces to the single-substitution score. The same homolog rows are retained for every query state within a variant's context.

The primary multiple-substitution prediction uses this endpoint construction. Two sensitivity analyses instead sum sequential single-site changes in each direction: one uses a fixed site order; the other averages over all `k!` orders. Both then apply whole-variant projection. These alternatives are retained as `M_sequence` and `M_orders` and are not substituted for the primary score.

The ensemble-control table also compares the directional WT star of the original raw model, whole-variant projection of the raw endpoint stars, component-and-whole projection retaining the published hidden biases, and the zero-hidden-bias EvoPhi-DDG construction. The first comparison changes endpoint-background treatment as well as projection; it is not a pure single-head ablation.

## JanusDDG and the ensemble

JanusDDG uses the published final checkpoint and its native whole-variant inference procedure. For the benchmark analyses, the released and verified predictions are provided directly; no new fitting is performed.

```text
Ensemble(A,B) = [JanusDDG(A,B) + EvoPhi_or_P(A,B)] / 2
```

The ensemble inherits exact antisymmetry. Its signed closure error on an identical path is the mean of the two signed component closure errors. On same-site paths, EvoPhi-DDG has zero closure error in exact arithmetic, so the ensemble closure error is half JanusDDG's. On multisite paths, EvoPhi-DDG can have nonzero closure error and no analogous exact guarantee applies.
