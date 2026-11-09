"""EvoPhi-DDG readout and endpoint construction, without encoder dependencies."""
import numpy as np


def potential_weights(fc1_weight, fc2_weight):
    """Return W.T @ v / 2 from the published head; neither bias enters it."""
    w = np.asarray(fc1_weight, dtype=np.float64)
    v = np.asarray(fc2_weight, dtype=np.float64).reshape(-1)
    if w.ndim != 2 or w.shape[0] != len(v):
        raise ValueError("Incompatible head weight dimensions")
    return w.T @ v / 2


def single(phi_a, phi_b, weights):
    """Score A -> B using fixed site and homolog context."""
    a, b = np.asarray(phi_a), np.asarray(phi_b)
    if a.shape != b.shape or a.shape[-1] != len(weights):
        raise ValueError("State vectors and potential weights must have matching dimensions")
    return (a.astype(np.float64) - b.astype(np.float64)) @ weights


def multiple(forward_components, reverse_components):
    """Project endpoint stars: A->Ai and B->B-minus-i, each from a fixed background."""
    f, r = np.asarray(forward_components), np.asarray(reverse_components)
    if f.ndim != 1 or f.shape != r.shape or not len(f):
        raise ValueError("Supply one forward and reverse component for each changed site")
    return float((f.sum() - r.sum()) / 2)


def ensemble(janus, evophi):
    """Fixed equal-weight mean; values must describe the same variants and direction."""
    j, e = np.asarray(janus), np.asarray(evophi)
    if j.shape != e.shape:
        raise ValueError("JanusDDG and EvoPhi-DDG scores must have the same shape")
    return (j + e) / 2
