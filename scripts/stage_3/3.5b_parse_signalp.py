#!/usr/bin/env python3
"""
3.5b_parse_signalp.py
Convert native SignalP6 prediction output into a simple accession /
YES|NO table for downstream PCuAC functional filtering.

Used for 3.5_signalp_pcuac.sh after signalp6 finishes.

Inputs:
- SignalP6 prediction_results.txt (native --format txt output,
  found in --output_dir from the signalp6 run)

Outputs:
- out/signalp_results.tsv  (accession, signal_n_terminal_peptide)

This is a pure format conversion: no biological logic is applied.

Important note!
SignalP is only used here for N-terminal signal peptide validation of
PCuAC candidates, not for any feature calculation. Stage 4 feature
extraction (MW, pI, GRAVY, composition) is entirely separate.
"""

import argparse
import csv
import os
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--signalp-dir", default="out/signalp_pcuac",
                    help="SignalP6 --output_dir from 3.5_signalp_pcuac.sh")
    ap.add_argument("--input", default=None,
                    help="Override: path to the native SignalP results file "
                         "(default: <signalp-dir>/prediction_results.txt)")
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()

    in_path = args.input or os.path.join(args.signalp_dir, "prediction_results.txt")

    if not os.path.exists(in_path):
        sys.exit(f"ERROR: SignalP results file not found: {in_path}")

    os.makedirs(args.outdir, exist_ok=True)
    out_path = os.path.join(args.outdir, "signalp_results.tsv")

    n_yes = 0
    n_no = 0
    rows_written = 0
    header_fields = None

    with open(in_path) as fh, open(out_path, "w", newline="") as ofh:
        writer = csv.writer(ofh, delimiter="\t")
        writer.writerow(["accession", "signal_n_terminal_peptide"])

        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if header_fields is None and "Prediction" in line:
                    header_fields = line.lstrip("#").strip().split("\t")
                continue

            fields = line.split("\t")
            if len(fields) < 2:
                continue

            accession = fields[0].strip()
            prediction = fields[1].strip()

            has_signal = prediction.upper().startswith("SP")
            if has_signal:
                n_yes += 1
            else:
                n_no += 1

            writer.writerow([accession, "YES" if has_signal else "NO"])
            rows_written += 1

    if rows_written == 0:
        sys.exit(f"ERROR: no prediction rows parsed from {in_path}. "
                  f"Check the SignalP6 output format.")

    if header_fields and len(header_fields) > 1 and "prediction" not in header_fields[1].lower():
        print(f"WARNING: expected column 2 to be 'Prediction', found "
              f"'{header_fields[1]}' - verify SignalP6 output format.", file=sys.stderr)

    print("-" * 50)
    print("SignalP output parsed")
    print(f"Input: {in_path}")
    print(f"Rows written: {rows_written}  (signal peptide: {n_yes}, none: {n_no})")
    print(f"wrote {out_path}")
    print("-" * 50)


if __name__ == "__main__":
    main()