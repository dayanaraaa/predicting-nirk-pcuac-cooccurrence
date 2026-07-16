#!/usr/bin/env python3
"""
Objective:
Clean and standardise BLAST hits.

Inputs:
- out/nir_seed_refseq_hits.tsv
- out/pcu_seed_refseq_hits.tsv

Outputs:
- out/nir_blast_metadata_clean.tsv
- out/pcu_blast_metadata_clean.tsv

Purpose:
Read raw BLAST TSVs, clean accessions, and keep only the best hit per accession (highest bitscore).
Then export clean metadata for script 3.2.
"""
import csv
import os
import re

# q = query sequence ID
# s = subject taxID
BLAST_COLS = [
    "qseqid","sseqid","pident","length","mismatch","gapopen",
    "qstart","qend","sstart","send","evalue","bitscore",
    "staxids","sscinames","stitle"
]

def clean_acc(acc):
    """
    Standardise accession IDs.
    """
    acc = acc.strip()
    if "|" in acc:
        parts = acc.split("|")
        return parts[1] if len(parts) >= 2 else parts[-1]
    return acc

def to_float(x, default=0.0):
    try:
        return float(x)
    except:
        return default

def read_blast(path):
    """
    Read BLAST TSVs and keep only the best hit per accession.
    Best = highest bitscore.
    """
    best = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(BLAST_COLS):
                continue

            d = dict(zip(BLAST_COLS, parts))
            acc = clean_acc(d["sseqid"])
            score = to_float(d["bitscore"])

            if acc not in best or score > to_float(best[acc]["bitscore"]):
                d["accession"] = acc
                best[acc] = d
    return best

def write_clean_table(hits, outpath):
    """
    Write clean BLAST metadata table.
    """
    cols = [
        "accession",
        "qseqid",
        "pident",
        "length",
        "bitscore",
        "evalue",
        "staxids",
        "sscinames",
        "stitle"
    ]
    with open (outpath, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for acc, d in hits.items():
            w.writerow([
                d.get("accession", ""),
                d.get("qseqid", ""),
                d.get("pident", ""),
                d.get("length", ""),
                d.get("bitscore", ""),
                d.get("evalue", ""),
                d.get("staxids", ""),
                d.get("sscinames", ""),
                d.get("stitle", ""),
            ])

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--nir", required=True, help="nir_seed_refseq_hits.tsv")
    ap.add_argument("--pcu", required=True, help="pcu_seed_refseq_hits.tsv")
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Read + Clean BLAST hits
    nir_hits = read_blast(args.nir)
    pcu_hits = read_blast(args.pcu)

    # Print summary
    print(f"NirK accessions retained: {len(nir_hits):,}")
    print(f"PCuAC accessions retained: {len(pcu_hits):,}")

    # Write the cleaned metadata tables
    write_clean_table(
        nir_hits,
        os.path.join(args.outdir, "nir_blast_metadata_clean.tsv")
    )
    write_clean_table(
        pcu_hits,
        os.path.join(args.outdir, "pcu_blast_metadata_clean.tsv")
    )

if __name__ == "__main__":
    main()