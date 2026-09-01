#!/usr/bin/env python3
"""
3.4_identify_nirk.py
Identify confirmed NirK sequences using the residue checks from Script 3.3.

3.3 identifies conserved NirK features:
- Type 1 copper site (His-Cys HC motif conservation)
- Type 2 copper site (conserved His residue signature)
- Catalytic residues

This script applies the final NirK filter:
NirK = T1_complete AND T2_complete AND length_ok AND (catalytic_ok unless --no-catalytic)

It then reads the 3.3 residue-status table > selects sequences that pass
the filter > extracts them from the 3.3 FASTA file.

Inputs:
- out/nirk_site_status.csv
- out/nirk_bonafide.fasta

Outputs:
- out/nirk.fasta
- out/nirk_accessions.txt

No new alignments are performed; this script only filters existing results.
"""

import argparse
import csv
import os
import sys


def read_fasta(path):
    hid, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if hid is not None:
                    yield hid, "".join(buf)
                hid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if hid is not None:
            yield hid, "".join(buf)

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--status", default="out/nirk_site_status.csv",
                    help="3.3 Nir site-status CSV")
    ap.add_argument("--fasta", default="out/nirk_bonafide.fasta",
                    help="3.3 Nir bona fide FASTA (sequences to subset)")
    ap.add_argument("--outdir", default="out")

    # Argument to apply if you want output for NirK with no catalytic dyad
    # nir_cat defined in 3.3
    ap.add_argument("--no-catalytic", action="store_true",
                    help="relax: require T1 AND T2 only (no catalytic dyad)")
    args = ap.parse_args()

    if not os.path.exists(args.status):
        sys.exit(f"ERROR: status CSV not found: {args.status}")

    require_cat = not args.no_catalytic

    fasta_lengths = {}
    for hid, seq in read_fasta(args.fasta):
        fasta_lengths[hid] = len(seq.replace("-", ""))

    nirk_accs = set()
    n_rows = 0
    n_both = 0  # T1 AND T2 (3.3 bonafide)
    n_cat_of_both = 0
    length_pass = 0
    length_fail = 0

    with open(args.status, newline="") as fh:
        r = csv.DictReader(fh)
        for row in r:
            n_rows += 1

            t1 = row.get("T1_complete") == "1"
            t2 = row.get("T2_complete") == "1"
            cat = row.get("catalytic_ok") == "1"

            acc = row["accession"]

            seq_len = fasta_lengths.get(acc, 0)
            length_ok = 300 <= seq_len <= 996

            if length_ok:
                length_pass += 1
            else:
                length_fail += 1

            if t1 and t2:
                n_both += 1
                if cat:
                    n_cat_of_both += 1
            is_nirk = (
                    t1
                    and t2
                    and length_ok
                    and (cat or not require_cat)
            )
            if is_nirk:
                nirk_accs.add(row["accession"])

    # Subset the FASTA (only sequences whose accession is a NirK)
    out_fasta = os.path.join(args.outdir, "nirk.fasta")
    out_list = os.path.join(args.outdir, "nirk_accessions.txt")
    written = 0
    with open(out_fasta, "w") as ffh, open(out_list, "w") as lfh:
        for hid, seq in read_fasta(args.fasta):
            if hid in nirk_accs:
                ffh.write(f">{hid}\n{seq}\n")
                lfh.write(hid + "\n")
                written += 1

    print("-" * 50)
    print("NirK identification - final checks")
    print("-" * 50)
    print(f"Definition: T1 AND T2" + (" AND catalytic" if require_cat else " (catalytic NOT required)"))
    print(f"Status rows examined: {n_rows}")
    print(f"Both copper sites: {n_both} (3.3 bonafide)")
    print(f"of which, {n_cat_of_both} are catalytic-positive")
    print(f"Length-pass: {length_pass}")
    print(f"Length-fail: {length_fail}")
    print()
    print(f"NirK (this filter): {len(nirk_accs)}")
    print(f"Sequences written: {written} to {out_fasta}")
    if written != len(nirk_accs):
        print(f"NOTE: {len(nirk_accs) - written} NirK accession(s) not found in the FASTA."
              f"(expected 0 if --fasta is the matching 3.3 bonafide set)")
    print(f"accession list: {out_list}")
    print("-" * 50)


if __name__ == "__main__":
    main()