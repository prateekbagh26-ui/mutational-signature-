"""Run the five steps from the post on the 21 breast cancer genomes.

    python run_analysis.py
"""

import numpy as np

from signatures import channel_fraction_tcw, cos, nmf, pair, tsv

BREAST = "data/breast21_sbs96.txt"
COSMIC = "data/cosmic_v3.4_sbs_grch37.txt"


def step1_open(M, SAMP):
    print("STEP 1 - open 21 genomes")
    N = M.sum(0)
    print(f"  channels x samples: {M.shape}")
    print(f"  total mutations:    {int(M.sum()):,}")
    print(f"  largest tumour:     {SAMP[int(N.argmax())]} with {int(N.max()):,} "
          f"({N.max() / M.sum():.0%} of all)")
    print(f"  ratio to quietest:  {N.max() / N.min():.1f}x")


def step2_see_it(TYPE, M):
    print("STEP 2 - see APOBEC before fitting")
    n_tcw, n_cg, obs, exp = channel_fraction_tcw(TYPE, M)
    print(f"  C>T / C>G channels: {n_cg}, of which TCA/TCT: {n_tcw}")
    print(f"  observed share in TCA/TCT: {obs:.1%}  (expected by chance: {exp:.1%})")
    print(f"  fold enrichment: {obs / exp:.1f}x")


def step3_factorise(M):
    print("STEP 3 - factorise into K=3 signatures (seed 0)")
    W, H = nmf(M, 3, 0)
    fit = cos(M.ravel(), (W @ H).ravel())
    print(f"  reconstruction cosine: {fit:.5f}")
    return W, H


def step4_ask_cosmic(TYPE_COS, SIG, COS, W):
    print("STEP 4 - compare with COSMIC")
    WN = W / W.sum(0)
    best = []
    for j in range(WN.shape[1]):
        s = np.array([cos(WN[:, j], COS[:, i]) for i in range(COS.shape[1])])
        o = np.argsort(-s)
        best.append([SIG[o[0]], round(float(s[o[0]]), 3),
                     SIG[o[1]], round(float(s[o[1]]), 3)])
        print(f"  signature {j + 1}: best {SIG[o[0]]} ({s[o[0]]:.3f}), "
              f"runner-up {SIG[o[1]]} ({s[o[1]]:.3f})")
    return best


def step5_stability(M):
    print("STEP 5 - how many signatures? fit vs. repeatability")
    FIT, STAB = [], []
    for k in (2, 3, 4, 5, 6):
        R = [nmf(M, k, s) for s in range(8)]
        FIT.append(max(cos(M.ravel(), (w @ h).ravel()) for w, h in R))
        S = [w / w.sum(0) for w, _ in R]
        STAB.append(min(pair(S[a], S[b]) for a in range(8) for b in range(a + 1, 8)))
    for k, f, s in zip((2, 3, 4, 5, 6), FIT, STAB):
        print(f"  K={k}: fit {f:.5f}   repeatability {s:.3f}")
    return FIT, STAB


def main():
    SAMP, TYPE, M = tsv(BREAST)
    SIG, T2, COS = tsv(COSMIC)
    assert TYPE == T2, "channel order differs between breast data and COSMIC"
    step1_open(M, SAMP)
    step2_see_it(TYPE, M)
    W, _H = step3_factorise(M)
    step4_ask_cosmic(TYPE, SIG, COS, W)
    step5_stability(M)


if __name__ == "__main__":
    main()
