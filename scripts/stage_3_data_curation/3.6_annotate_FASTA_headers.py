#!/usr/bin/env python3
"""
Objective:
Add organism names and NirK/PCuAC class labels to cleaned FASTA sequences
using the accession-level labels generated in Stage 3.5.

Inputs:
- out/nir_clean.fasta
- out/pcu_clean.fasta
- out/sequence_label_map.csv
- out/species_master_table.csv

Outputs:
- out/nir_annotated.fasta
- out/pcu_annotated.fasta
- out/3.6_annotation_report.tsv

The sequence_label_map.csv file is the primary source for accession-level
labels because it contains all cleaned sequences, not only representative
accessions.
"""

import argparse
import csv
import os
import sys

BLAST_COLS = [
    "accession",
    "qseqid",
    "bitscore",
    "evalue",
    "staxids",
    "sscinames",
    "stitle",
]


CLASS_LABEL = {
    "both": "NirK+ PCuAC+",
    "nirK_only": "NirK+ PCuAC-",
    "pcuac_only": "NirK- PCuAC+",
}


def read_fasta(path):
    if not os.path.exists(path):
        sys.exit(f"Missing FASTA: {path}")
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

def strip_version(acc):
    return acc.split(".")[0]


def load_species_table(path):
    if not os.path.exists(path):
        sys.exit(f"Missing species table: {path}")
    acc_to_org = {}
    acc_to_class = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            sp = row["species"].strip()
            cls = row["class"].strip()
            for key in ("nirK_rep_accession", "pcuac_rep_accession"):
                acc = row.get(key, "").strip()
                if acc and acc != "N/A":
                    base = strip_version(acc)
                    acc_to_org[base] = sp
                    acc_to_class[base] = cls
    return acc_to_org, acc_to_class


def load_sequence_label_map(path):
    """
    Load Stage 3.5's accession-level label map.
    """
    if not os.path.exists(path):
        print(f"WARNING: {path} not found. Falling back to representative-only "
              f"non-representative sequences will get guessed labels.")
        return {}
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            base = strip_version(row["accession"].strip())
            out[base] = {
                "species": row["species"].strip(),
                "class": row["class"].strip(),
                "protein": row["protein"].strip(),
                "is_representative": row.get("is_representative", "0").strip(),
            }
    return out


def load_blast(path):
    """
    Fallback organism names from BLAST.
    """
    if not os.path.exists(path):
        return {}
    out = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < len(BLAST_COLS):
                continue
            d = dict(zip(BLAST_COLS, parts))
            acc = strip_version(d["accession"])
            org = d["sscinames"].split(";")[0].strip()
            if acc and org:
                out[acc] = org
    return out


def annotate(in_fasta, out_fasta, source, label_map, acc_to_org, acc_to_class, fallback, report):
    written = 0
    guessed = 0
    with open(out_fasta, "w") as fh:
        for acc, seq in read_fasta(in_fasta):
            base = strip_version(acc)
            mapped = label_map.get(base)

            if mapped is not None:
                # true species-derived label
                # for this exact accession, whether or not it's the
                # species' chosen representative.
                organism = mapped["species"]
                raw_cls = mapped["class"]
                is_rep = mapped["is_representative"]
                provenance = "label_map"
            else:
                organism = (
                    acc_to_org.get(base)
                    or fallback.get(base)
                    or "Unknown organism"
                )
                raw_cls = acc_to_class.get(base)
                if not raw_cls:
                    raw_cls = "nirK_only" if source == "nir" else "pcuac_only"
                is_rep = "1" if base in acc_to_class else "0"
                provenance = "guessed"
                guessed += 1

            label = CLASS_LABEL.get(raw_cls, raw_cls)
            fh.write(f">{acc} | {organism} | {label}\n{seq}\n")
            written += 1
            report.append({
                "accession": acc,
                "organism": organism,
                "class": label,
                "source": source,
                "is_representative": is_rep,
                "label_source": provenance,
            })
    if guessed:
        print(f"  {source.upper()}: {guessed}/{written} sequences fell back to a "
              f"guessed label (not found in sequence_label_map.csv)")
    return written


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nir-fasta", default="out/nir_clean.fasta")
    ap.add_argument("--pcu-fasta", default="out/pcu_clean.fasta")

    ap.add_argument("--species-table", default="out/species_master_table.csv")
    ap.add_argument("--sequence-label-map", default="out/sequence_label_map.csv")

    ap.add_argument("--nir-blast", default="out/nir_blast_metadata_clean.tsv")
    ap.add_argument("--pcu-blast", default="out/pcu_blast_metadata_clean.tsv")
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    label_map = load_sequence_label_map(args.sequence_label_map)
    acc_to_org, acc_to_class = load_species_table(args.species_table)

    fallback = {}
    fallback.update(load_blast(args.nir_blast))
    fallback.update(load_blast(args.pcu_blast))

    report = []

    targets = [
        ("nir", args.nir_fasta, os.path.join(args.outdir, "nir_annotated.fasta")),
        ("pcu", args.pcu_fasta, os.path.join(args.outdir, "pcu_annotated.fasta")),
    ]

    for source, in_fa, out_fa in targets:
        if os.path.exists(in_fa):
            n = annotate(in_fa, out_fa, source, label_map, acc_to_org, acc_to_class, fallback, report)
            print(f"{source.upper()}: {n} sequences annotated: {out_fa}")

    report_path = os.path.join(args.outdir, "3.6_annotation_report.tsv")
    with open(report_path, "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["accession", "organism", "class", "source", "is_representative", "label_source"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(report)

    print(f"Annotation report written: {report_path}")

if __name__ == "__main__":
    main()
