#!/usr/bin/env python3
"""
3.5_identify_functional_pcu.py
Final PCuAC classification: combine the Cu(I) motif result, the length
filter, and the SignalP N-terminal signal peptide prediction into the
final functional PCuAC definition.

Functional PCuAC = Cu(I) motif AND length requirement AND SignalP signal peptide

Cu(I) motif:
    His69, Met80, His103, Met105

Length:
    100-350 aa

The length check is re-applied here defensively even though it was
already applied upstream in 3.5_pcuac_candidates.py.

The Cu(I) motif is re-read from --status (not re-derived from --fasta)
since it was already computed once in 3.3_tag_sites.py; no biological
logic is duplicated, only the existing status column is combined with
the SignalP result.

No new alignments are performed; this script only filters existing results.

Inputs:
- out/pcuac_site_status.csv
- out/pcuac_signalp_candidates.fasta
- SignalP prediction (out/signalp_results.tsv)

Outputs:
- out/pcuac_functional.fasta
- out/pcuac_functional_accessions.txt
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
    ap.add_argument("--status", default="out/pcuac_site_status.csv",
                    help="3.3 Pcu site-status CSV")
    ap.add_argument("--fasta", default="out/pcuac_signalp_candidates.fasta",
                    help="3.5_pcuac_candidates.py candidate FASTA "
                         "(motif + length already applied)")
    ap.add_argument("--outdir", default="out")
    ap.add_argument("--signalp", required=True,
                    help="SignalP prediction")
    args = ap.parse_args()

    if not os.path.exists(args.status):
        sys.exit(f"ERROR: status CSV not found: {args.status}")
    if not os.path.exists(args.signalp):
        sys.exit(f"ERROR: signalP prediction not found: {args.signalp}")

    signal_status = {}

    with open(args.signalp, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            signal_status[row["accession"]] = (
            row["signal_n_terminal_peptide"] == "YES")

    fasta_sequences = {}
    for hid, seq in read_fasta(args.fasta):
        fasta_sequences[hid] = seq

    func_accs = set()

    n_rows = 0
    n_pcuac = 0

    length_pass = 0
    length_fail = 0

    signal_pass = 0
    signal_fail = 0

    with open(args.status, newline="") as fh:
        r = csv.DictReader(fh)
        for row in r:
            acc = row["accession"]
            if acc not in fasta_sequences:
                # not in the candidate set passed in (already failed
                # motif or length upstream) - not part of this run
                continue

            n_rows += 1
            motif = row.get("pcuac_motif_present") == "1"
            seq = fasta_sequences.get(acc, "")
            seq_len = len(seq.replace("-", ""))
            length_ok = 100 <= seq_len <= 350
            signal_ok = signal_status.get(acc, False)

            if length_ok:
                length_pass += 1
            else:
                length_fail += 1
            if signal_ok:
                signal_pass += 1
            else:
                signal_fail += 1
            is_pcuac = (
                    motif and length_ok and signal_ok
            )
            if is_pcuac:
                n_pcuac += 1
                func_accs.add(acc)

    out_fasta = os.path.join(args.outdir, "pcuac_functional.fasta")
    out_list = os.path.join(args.outdir, "pcuac_functional_accessions.txt")
    written = 0
    with open(out_fasta, "w") as ffh, open(out_list, "w") as lfh:
        for hid, seq in read_fasta(args.fasta):
            if hid in func_accs:
                ffh.write(f">{hid}\n{seq}\n")
                lfh.write(hid + "\n")
                written += 1

    print("-" * 50)
    print("PCuAC identification (Cu(I) motif)")
    print()
    print(f"Definition: His69, Met80, His103, Met105 ")
    print(f"Candidates examined (post motif+length filter): {n_rows}")
    print(f"PCuAC passing motif + length + signal filter: {n_pcuac}")
    print(f"Length-pass (defensive re-check): {length_pass}")
    print(f"Length-fail (defensive re-check): {length_fail}")
    if length_fail > 0:
        print(f"WARNING: {length_fail} candidate(s) failed the length re-check - "
              f"verify --status/--fasta correspond to the same run.")
    print(f"Signal-pass: {signal_pass}")
    print(f"Signal-fail: {signal_fail}")
    print(f"Sequences written: {written}  -> {out_fasta}")
    if written != len(func_accs):
        print(f"NOTE: {len(func_accs) - written} accession(s) not found in the FASTA "
              f"(expected 0 if --fasta is the matching candidate set)")
    print(f"accession list: {out_list}")
    print("-" * 50)


if __name__ == "__main__":
    main()