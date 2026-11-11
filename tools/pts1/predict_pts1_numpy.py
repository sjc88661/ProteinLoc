#!/usr/bin/env python3
"""Standalone PTS1 predictor that needs ONLY numpy + svm_weights.npz.

The upstream package (https://github.com/oKoch/pts1-prediction-tool) retrains an
sklearn SVC from FASTA data on every predictor construction. This script loads the
exported weights (support-vector codes, dual coefficients, intercept, gamma) and
reproduces the identical RBF-SVM decision function in plain numpy, so no
scikit-learn / pandas / biopython installation is required.

Examples
--------
    python3 predict_pts1_numpy.py --sequence MMMMKLSKMLLLSLSKLSKLSKLSKL
    python3 predict_pts1_numpy.py --fasta proteins.fasta > predictions.csv
"""

import argparse
import csv
import hashlib
import json
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_NPZ = os.path.join(HERE, "svm_weights.npz")


class NumpyPTS1Predictor:
    def __init__(self, npz_path=DEFAULT_NPZ):
        w = np.load(npz_path)
        self.alphabet = [str(a) for a in w["aa_alphabet"]]
        self.index = {a: i for i, a in enumerate(self.alphabet)}
        self.k = int(w["c_terminal_length"])
        codes = w["support_vector_codes"].astype(np.int64)
        n_sv = codes.shape[0]
        self.sv = np.zeros((n_sv, len(self.alphabet) * self.k), dtype=np.float64)
        for j in range(self.k):
            self.sv[np.arange(n_sv), j * len(self.alphabet) + codes[:, j]] = 1.0
        self.dual_coef = w["dual_coef"].astype(np.float64).reshape(-1)
        self.intercept = float(np.ravel(w["intercept"])[0])
        self.gamma = float(w["gamma"])
        self.C = float(w['C'])
        self._sv_sqnorm = float(self.k)  # every support vector has k ones
        w.close()

    def decision_function(self, aa_seq):
        seq = str(aa_seq).upper().replace("*", "")
        if len(seq) < self.k:
            raise ValueError(f"sequence shorter than {self.k} residues: {aa_seq!r}")
        tail = seq[-self.k:]
        row = np.zeros(len(self.alphabet) * self.k, dtype=np.float64)
        unknown = 0
        for j, aa in enumerate(tail):
            idx = self.index.get(aa)
            if idx is None:  # Upstream pandas categorical unknowns have an all-zero block.
                unknown += 1
                continue
            row[j * len(self.alphabet) + idx] = 1.0
        with warnings.catch_warnings(), np.errstate(all="ignore"):
            warnings.simplefilter("ignore")
            dots = self.sv @ row
            sq_dist = self._sv_sqnorm + float(row.sum()) - 2.0 * dots
            kern = np.exp(-self.gamma * sq_dist)
            return float(self.dual_coef @ kern + self.intercept)

    def predict(self, aa_seq):
        return 1 if self.decision_function(aa_seq) > 0 else 0

    def check_for_pts1(self, aa_seq):
        return self.predict(aa_seq) == 1


def read_fasta(path):
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def main(argv=None):
    ap = argparse.ArgumentParser(description="PTS1 predictor using exported SVM weights (numpy only)")
    ap.add_argument("--npz", default=DEFAULT_NPZ, help="path to svm_weights.npz")
    ap.add_argument("--sequence", help="single amino-acid sequence (C-terminus last)")
    ap.add_argument("--fasta", help="FASTA file with one or more sequences")
    ap.add_argument("--out", help="CSV output path (default: stdout)")
    ap.add_argument('--version',action='store_true')
    args = ap.parse_args(argv)

    model = NumpyPTS1Predictor(args.npz)
    if args.version:
        with open(args.npz,'rb') as stream:
            digest=hashlib.sha256(stream.read()).hexdigest()
        print(json.dumps(dict(tool='oKoch PTS1 SVM',upstream_version='1.0.2',weights='user-reproduced',
                              numpy=np.__version__,n_support=model.sv.shape[0],c_terminal_length=model.k,
                              C=model.C,gamma=model.gamma,sha256=digest)))
        return
    records = []
    if args.sequence:
        records.append(("input", args.sequence))
    if args.fasta:
        records.extend(read_fasta(args.fasta))
    if not records:
        ap.error("provide --sequence and/or --fasta")

    out = open(args.out, "w", newline="") if args.out else sys.stdout
    writer = csv.writer(out)
    writer.writerow(["id", "prediction", "is_pts1", "decision_function", "length", "c_terminal_14"])
    for name, seq in records:
        seq = str(seq).upper().replace("*", "")
        if len(seq) < model.k:
            writer.writerow([name, "error: sequence too short", "", "", len(seq), ""])
            continue
        dec = model.decision_function(seq)
        pred = int(dec > 0)
        writer.writerow([name, pred, bool(pred), f"{dec:.17g}", len(seq), seq[-model.k:]])
    if args.out:
        out.close()
        print(f"wrote {args.out} ({len(records)} sequences)")


if __name__ == "__main__":
    main()
