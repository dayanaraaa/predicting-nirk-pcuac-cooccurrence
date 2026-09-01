#!/usr/bin/env python3
"""
3.5_pcuac_candidates.py
Build the PCuAC SignalP candidate set: apply the existing PCuAC Cu(I)
motif requirement and the existing 100-350 aa length requirement to the
3.3 bona fide PCuAC set, before SignalP is run.

This script performs no SignalP-related classification and does not
apply the final PCuAC functional definition.

Inputs:
- out/pcuac_site_status.csv
- out/pcuac_bonafide.fasta

Outputs:
- out/pcuac_signalp_candidates.fasta
- out/pcuac_signalp_candidates_accessions.txt
- out/pipeline_metadata/3.5_pcuac_candidates_metadata.txt

Dependencies:
- Python 3
"""

import argparse
import csv
import os
import sys
import platform
import datetime

PIPELINE_VERSION = "stage3.v2"


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


def write_pipeline_metadata(outdir, script_name, input_files, output_files, notes=""):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{script_name}_metadata.txt")
    with open(path, "w") as fh:
        fh.write("Stage: 3\n")
        fh.write(f"Script: {script_name}\n")
        fh.write(f"Pipeline version: {PIPELINE_VERSION}\n\n")
        fh.write(f"Date generated: {datetime.datetime.now().isoformat(timespec='seconds')}\n\n")
        fh.write(f"Python: {platform.python_version()}\n\n")
        fh.write("Input files:\n")
        for f in input_files:
            fh.write(f"- {f}\n")
        fh.write("\nOutput files:\n")
        for f in output_files:
            fh.write(f"- {f}\n")
        if notes:
            fh.write("\n" + notes.strip() + "\n")
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--status", default="out/pcuac_site_status.csv",
                    help="3.3 PCuAC site-status CSV")
    ap.add_argument("--fasta", default="out/pcuac_bonafide.fasta",
                    help="3.3 PCuAC bona fide FASTA (sequences to subset)")
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()

    if not os.path.exists(args.status):
        sys.exit(f"ERROR: status CSV not found: {args.status}")

    fasta_lengths = {}
    for hid, seq in read_fasta(args.fasta):
        fasta_lengths[hid] = len(seq.replace("-", ""))

    n_input = 0
    n_motif_pass = 0
    n_length_pass = 0
    candidate_accs = set()

    with open(args.status, newline="") as fh:
        r = csv.DictReader(fh)
        for row in r:
            n_input += 1
            acc = row["accession"]

            motif = row.get("pcuac_motif_present") == "1"
            seq_len = fasta_lengths.get(acc, 0)
            length_ok = 100 <= seq_len <= 350

            if motif:
                n_motif_pass += 1
            if length_ok:
                n_length_pass += 1

            if motif and length_ok:
                candidate_accs.add(acc)

    out_fasta = os.path.join(args.outdir, "pcuac_signalp_candidates.fasta")
    out_list = os.path.join(args.outdir, "pcuac_signalp_candidates_accessions.txt")
    written = 0
    with open(out_fasta, "w") as ffh, open(out_list, "w") as lfh:
        for hid, seq in read_fasta(args.fasta):
            if hid in candidate_accs:
                ffh.write(f">{hid}\n{seq}\n")
                lfh.write(hid + "\n")
                written += 1

    n_removed = n_input - written

    print("-" * 50)
    print("PCuAC SignalP candidate filtering (motif + length only)")
    print("-" * 50)
    print(f"Input sequences (3.3 bona fide): {n_input}")
    print(f"Passing PCuAC motif filter: {n_motif_pass}")
    print(f"Passing length filter (100-350 aa): {n_length_pass}")
    print(f"SignalP candidates written (motif AND length): {written}")
    print(f"Removed before SignalP: {n_removed}")
    print(f"Output FASTA: {out_fasta}")
    if written != len(candidate_accs):
        print(f"NOTE: {len(candidate_accs) - written} candidate accession(s) not found in the FASTA.")
    print(f"accession list: {out_list}")
    print("-" * 50)

    write_pipeline_metadata(
        outdir="out/pipeline_metadata",
        script_name="3.5_pcuac_candidates",
        input_files=[args.status, args.fasta],
        output_files=[out_fasta, out_list],
        notes=(
            "Filters applied: PCuAC Cu(I) motif present AND length 100-350 aa.\n"
            "No new biological criteria - identical filters to those used in\n"
            "3.5_identify_functional_pcu.py, applied here only to reduce SignalP input."
        ),
    )


if __name__ == "__main__":
    main()