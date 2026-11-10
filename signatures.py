"""Mutational signatures from scratch, using only numpy.

Functions here follow the code shown in the five-part post
"5 cancer genomics projects you can do at home with Claude"
(@thebiosciencedesk). Names are kept close to the slides so you can
compare them line by line.

Matrix file format (tab-separated):
    header:  Type <TAB> sample1 <TAB> sample2 ...
    rows:    channel label <TAB> counts per sample
"""

import numpy as np


def tsv(path):
    """Read a tab-separated matrix. Returns (column names, row labels, values)."""
    with open(path) as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh]
    cols = rows[0][1:]
    labels = [r[0] for r in rows[1:]]
    values = np.array([[float(v) for v in r[1:]] for r in rows[1:]])
    return cols, labels, values


def cos(a, b):
    """Cosine similarity between two vectors."""
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))


def nmf(V, k, seed, it=4000):
    """Non-negative matrix factorisation with Lee & Seung multiplicative updates.

    Returns W (channels x k, signatures) and H (k x samples, exposures).
    """
    rs = np.random.default_rng(seed)
    W = rs.random((V.shape[0], k)) + 0.1
    H = rs.random((k, V.shape[1])) + 0.1
    for _ in range(it):
        H *= (W.T @ V) / (W.T @ W @ H + 1e-9)
        W *= (V @ H.T) / (W @ H @ H.T + 1e-9)
    return W, H


def pair(A, B):
    """Greedy one-to-one matching of signatures in A to signatures in B.

    Returns the mean cosine over matched pairs (1.0 = identical sets).
    """
    k, used, tot = A.shape[1], set(), 0.0
    S = [[cos(A[:, i], B[:, j]) for j in range(k)] for i in range(k)]
    for i in sorted(range(k), key=lambda i: -max(S[i])):
        j = max((j for j in range(k) if j not in used), key=lambda j: S[i][j])
        used.add(j)
        tot += S[i][j]
    return tot / k


def channel_fraction_tcw(TYPE, M):
    """Step 2: the C>T / C>G channels with a T on the left and A/T on the right.

    Returns (n_tcw, n_cg, observed share, expected share).
    """
    CG = [i for i, t in enumerate(TYPE) if t[2:5] in ("C>T", "C>G")]
    TCW = [i for i in CG if TYPE[i][0] == "T" and TYPE[i][6] in "AT"]
    OBS = float(M[TCW].sum() / M[CG].sum())
    EXP = len(TCW) / float(len(CG))
    return len(TCW), len(CG), OBS, EXP
