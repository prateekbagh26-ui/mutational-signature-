"""Convert the 21 breast-cancer genome matrix from MATLAB format to tab-separated text.

Source: 21_breast_WGS_substitutions.mat, shipped inside the SigProfilerAssignment
and SigProfilerExtractor pip packages. Run this once:

    pip install scipy
    python prepare_data.py /path/to/21_breast_WGS_substitutions.mat data/breast21_sbs96.txt
"""

import sys

import numpy as np
import scipy.io as sio


def main(mat_path, out_path):
    d = sio.loadmat(mat_path)
    counts = d["originalGenomes"]  # 96 channels x 21 samples
    names = [str(x[0]) for x in d["sampleNames"].ravel()]
    ctx = [str(x[0]).strip("[]'") for x in d["subtypes"]]  # trinucleotide context
    mut = [str(x[0]).strip("[]'") for x in d["types"]]  # mutation class, e.g. C>A

    labels = []
    for c, m in zip(ctx, mut):
        ref, _alt = m.split(">")
        assert c[1] == ref, f"context {c} does not match mutation {m}"
        labels.append(f"{c[0]}[{m}]{c[2]}")

    with open(out_path, "w") as fh:
        fh.write("Type\t" + "\t".join(names) + "\n")
        for lab, row in zip(labels, counts):
            fh.write(lab + "\t" + "\t".join(str(int(v)) for v in row) + "\n")
    print(f"wrote {out_path}: {counts.shape[0]} channels x {counts.shape[1]} samples, "
          f"{int(counts.sum())} mutations")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
