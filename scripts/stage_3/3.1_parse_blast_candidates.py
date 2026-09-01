#!/usr/bin/env python3
"""
3.1_parse_blast_candidates.py
Parse NirK and PCuAC BLAST candidate hits and generate
clean accession lists and metadata tables for sequence retrieval.

Inputs:
- NirK BLAST output
- PCuAC BLAST output

Outputs:
- nirk_accessions_clean.txt
- pcuac_accessions_clean.txt
- nirk_candidate_metadata.tsv
- pcuac_candidate_metadata.tsv
"""

import argparse
import csv
import os
import re
import sys


TSV_COLS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "staxids",
    "sscinames",
    "stitle",
]


def normalize_header(value):
    return re.sub(r"[^a-z0-9]", "", value.lower())


def strip_accession(raw):
    """
    Remove BLAST accession wrappers.
    Example:
    ref|WP_012345678.1| -> WP_012345678.1
    """
    return re.sub(r"^[a-z]+\|", "", raw).rstrip("|")

def split_taxids(raw):
    raw = (raw or "").strip()
    if not raw or raw.upper() == "N/A":
        return set()
    return {taxid for taxid in re.split(r"[;,]\s*", raw) if taxid}

def to_float(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return 0.0

def read_hits(path, query_family):
    if not os.path.exists(path):
        sys.exit(f"ERROR: file not found: {path}")

    with open(path, newline="") as fh:
        first_line = fh.readline()

    is_web = "," in first_line and (
        "taxid" in normalize_header(first_line)
        or "scientificname" in normalize_header(first_line)
    )

    rows = []

    if is_web:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            header_map = {normalize_header(k): k for k in reader.fieldnames or []}

            def get_column(*names):
                for name in names:
                    if name in header_map:
                        return header_map[name]
                return None

            tax_col = get_column("taxid", "taxids")
            org_col = get_column("scientificname", "organism")
            acc_col = get_column("accession", "acc")
            pid_col = get_column("perident", "peridentity", "identity")
            evalue_col = get_column("evalue", "expectvalue")
            score_col = get_column("maxscore", "bitscore", "score")
            length_col = get_column("acclen", "accessionlength", "length")

            for row in reader:
                if not any(row.values()):
                    continue

                rows.append({
                    "accession": (row.get(acc_col) or "").strip(),
                    "query_accession": "",
                    "query_family": query_family,
                    "taxids": split_taxids(row.get(tax_col)),
                    "organism": (row.get(org_col) or "N/A").strip(),
                    "pident": row.get(pid_col, ""),
                    "evalue": row.get(evalue_col, ""),
                    "bitscore": row.get(score_col, ""),
                    "length": row.get(length_col, ""),
                })

    else:
        with open(path, newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")

            for line in reader:
                if len(line) < len(TSV_COLS):
                    continue

                row = {name: line[i] for i, name in enumerate(TSV_COLS)}

                rows.append({
                    "accession": strip_accession(row["sseqid"]),
                    "query_accession": row["qseqid"],
                    "query_family": query_family,
                    "taxids": split_taxids(row["staxids"]),
                    "organism": row["sscinames"] if row["sscinames"] != "N/A" else row["stitle"],
                    "pident": row["pident"],
                    "evalue": row["evalue"],
                    "bitscore": row["bitscore"],
                    "length": row["length"],
                })
    return rows

def collapse_accessions(rows):
    best = {}

    for row in rows:
        accession = row["accession"]
        if not accession:
            continue
        if accession not in best or to_float(row["bitscore"]) > to_float(best[accession]["bitscore"]):
            best[accession] = row
    return [best[accession] for accession in sorted(best)]

def write_outputs(rows, prefix, outdir):
    accessions_path = os.path.join(outdir, f"{prefix}_accessions_clean.txt")
    metadata_path = os.path.join(outdir, f"{prefix}_candidate_metadata.tsv")

    with open(accessions_path, "w") as fh:
        for row in rows:
            fh.write(row["accession"] + "\n")

    with open(metadata_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")

        writer.writerow([
            "accession",
            "query_accession",
            "query_family",
            "taxids",
            "organism",
            "pident",
            "evalue",
            "bitscore",
            "length",
        ])

        for row in rows:
            writer.writerow([
                row["accession"],
                row["query_accession"],
                row["query_family"],
                ";".join(sorted(row["taxids"])),
                row["organism"],
                row["pident"],
                row["evalue"],
                row["bitscore"],
                row["length"],
            ])
    return accessions_path, metadata_path

def main():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--nirk", required=True, help="NirK BLAST hit file")
    parser.add_argument("--pcuac", required=True, help="PCuAC BLAST hit file")
    parser.add_argument("--outdir", default="out", help="Output directory")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    nirk_hits = collapse_accessions(read_hits(args.nirk, "nirk"))
    pcuac_hits = collapse_accessions(read_hits(args.pcuac, "pcuac"))

    outputs = []
    outputs.extend(write_outputs(nirk_hits, "nirk", args.outdir))
    outputs.extend(write_outputs(pcuac_hits, "pcuac", args.outdir))

    print("-" * 50)
    print("BLAST candidates parsed")
    print(f"NirK candidates: {len(nirk_hits)}")
    print(f"PCuAC candidates: {len(pcuac_hits)}")
    print("-" * 50)

    for output in outputs:
        print(f"wrote {output}")


if __name__ == "__main__":
    main()