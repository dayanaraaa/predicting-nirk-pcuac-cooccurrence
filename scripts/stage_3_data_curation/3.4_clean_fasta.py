#!/usr/bin/env python3
"""
Objective:
Clean the bona fide FASTA files produced in Stage 3.3. This step removes
sequences that are too short, too long, contain too many gaps or ambiguous
amino acids, or are duplicates. The output FASTAs are used for downstream
annotation and feature extraction.

Inputs:
- out/nir_bonafide.fasta
- out/pcu_bonafide.fasta

Outputs:
- out/<name>_clean.fasta = cleaned sequences (degapped)
- out/<name>_clean_ids.txt = one accession per line
- out/<name>_rejected.tsv = rejected sequences + reason

Parameters:
Literature-supported sequence length ranges:
- NirK: 300–996 aa
- PCuAC: 100–350 aa

Additional quality-control thresholds:
--max-gap = maximum allowed gap fraction
--max-ambig = maximum ambiguous AA count
"""

import argparse
import os
import sys

AMBIG_AA = set("BJOUZ")


def read_fasta(path):
    """
    Yield (accession, sequence) from a FASTA file.
    """
    if not os.path.exists(path):
        sys.exit(f"File not found: {path}")
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq)
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if header:
        yield header, "".join(seq)


def clean_fasta(
    input_fasta,
    out_fasta,
    out_ids,
    out_rejected,
    *,
    min_len,
    max_len,
    max_gap,
    max_ambig,
    label,
):
    """
    Apply length, gap, ambiguity, and duplicate filters.
    """
    seen_ids = set()
    seen_seqs = set()

    counts = {
        "label": label,
        "n_in": 0,
        "n_kept": 0,
        "n_dup_id": 0,
        "n_dup_seq": 0,
        "n_short": 0,
        "n_long": 0,
        "n_gap": 0,
        "n_ambig": 0,
    }

    with open(out_fasta, "w") as ff, \
         open(out_ids, "w") as fi, \
         open(out_rejected, "w") as fr:

        fr.write("accession\treason\tlength\tgap_frac\tambig_count\n")

        for acc, seq in read_fasta(input_fasta):
            counts["n_in"] += 1

            seq_u = seq.upper()
            degapped = seq_u.replace("-", "").replace(".", "")
            length = len(degapped)
            gap_frac = (len(seq_u) - length) / max(len(seq_u), 1)
            ambig_n = sum(aa in AMBIG_AA for aa in degapped)

            def reject(reason):
                fr.write(f"{acc}\t{reason}\t{length}\t{gap_frac:.3f}\t{ambig_n}\n")

            if acc in seen_ids:
                counts["n_dup_id"] += 1
                reject("duplicate_id")
                continue

            if length < min_len:
                counts["n_short"] += 1
                reject("too_short")
                continue

            if length > max_len:
                counts["n_long"] += 1
                reject("too_long")
                continue

            if gap_frac > max_gap:
                counts["n_gap"] += 1
                reject("high_gap")
                continue

            if ambig_n > max_ambig:
                counts["n_ambig"] += 1
                reject("high_ambig")
                continue

            if degapped in seen_seqs:
                counts["n_dup_seq"] += 1
                reject("duplicate_seq")
                continue

            seen_ids.add(acc)
            seen_seqs.add(degapped)
            counts["n_kept"] += 1

            ff.write(f">{acc}\n{degapped}\n")
            fi.write(acc + "\n")

    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", choices=["nir", "pcu", "both"], default="both")
    ap.add_argument("--nir-fasta", default="out/nir_bonafide.fasta")
    ap.add_argument("--pcu-fasta", default="out/pcu_bonafide.fasta")
    ap.add_argument("--outdir", default="out")

    ap.add_argument("--max-gap", type=float, default=0.10)
    ap.add_argument("--max-ambig", type=int, default=2)

    # Literature-supported sequence length ranges
    ap.add_argument("--nir-min", type=int, default=300)
    ap.add_argument("--nir-max", type=int, default=996)

    ap.add_argument("--pcu-min", type=int, default=100)
    ap.add_argument("--pcu-max", type=int, default=350)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    targets = []
    if args.name in ("nir", "both"):
        targets.append(("nir", args.nir_fasta, args.nir_min, args.nir_max))
    if args.name in ("pcu", "both"):
        targets.append(("pcu", args.pcu_fasta, args.pcu_min, args.pcu_max))

    for name, fasta, min_len, max_len in targets:
        out_fasta = os.path.join(args.outdir, f"{name}_clean.fasta")
        out_ids = os.path.join(args.outdir, f"{name}_clean_ids.txt")
        out_rej = os.path.join(args.outdir, f"{name}_rejected.tsv")

        summary = clean_fasta(
            fasta,
            out_fasta,
            out_ids,
            out_rej,
            min_len=min_len,
            max_len=max_len,
            max_gap=args.max_gap,
            max_ambig=args.max_ambig,
            label=name.upper(),
        )

        kept_pct = (
            f"{100 * summary['n_kept'] / summary['n_in']:.1f}%"
            if summary["n_in"] else "N/A"
        )
        print(f"{name.upper()}: kept {summary['n_kept']} / {summary['n_in']} ({kept_pct})")

        print(f"  cleaned FASTA: {out_fasta}")
        print(f"  cleaned IDs:   {out_ids}")
        print(f"  rejected:      {out_rej}")
        print()

if __name__ == "__main__":
    main()
