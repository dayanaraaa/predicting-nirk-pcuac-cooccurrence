#!/usr/bin/env python3
"""
Objective:
Evaluate conserved copper-binding residues in NirK and PCuAC hits aligned
to their seed proteins. Sequences passing ligand criteria are marked as
bona fide homologues candidates.

Inputs:
MAFFT --keeplength alignment containing:
- seed sequence
- all BLAST hits for NirK or PCuAC

Outputs:
- out/<name>_site_status.csv
    Residue table with ligand positions, boolean flags, and category labels.

- out/<name>_bonafide.fasta
    Ungapped sequences passing the bona fide criteria.

Parameters:
NirK:
- min length 300 aa
- Type‑1 ligands: His134, Cys175, His183, Met188
- Type‑2 ligands: His139, His174, His329
- catalytic dyad: Asp137, His280
- bona fide = Type‑1 AND Type‑2

PCuAC:
- min length 100 aa
- Cu(I) motif: His69, Met80, His103, Met105
- Cu(II) proxy: His count ≥ 3 in C‑terminal(last 30 aa)
- bona fide = Cu(I) motif present
"""

import argparse
import csv
import os
import sys

NIR_SITES = {
    "T1_134": (134, "H"), "T1_175": (175, "C"),
    "T1_183": (183, "H"), "T1_188": (188, "M"),
    "T2_139": (139, "H"), "T2_174": (174, "H"),
    "T2_329": (329, "H"),
    "cat_137": (137, "D"), "cat_280": (280, "H"),
}
NIR_T1 = ["T1_134", "T1_175", "T1_183", "T1_188"]
NIR_T2 = ["T2_139", "T2_174", "T2_329"]
NIR_CAT = ["cat_137", "cat_280"]

PCU_SITES = {
    "His69": (69, "H"), "Met80": (80, "M"),
    "His103": (103, "H"), "Met105": (105, "M"),
}
PCU_MOTIF = ["His69", "Met80", "His103", "Met105"]

SEED_IDS = {
    "Q02219", "sp|Q02219|ANIA_NEIGO",
    "A0AAQ1E0N0", "tr|A0AAQ1E0N0|A0AAQ1E0N0_NEIGO",
}

MIN_LEN = {"nir": 300, "pcu": 100}


def read_fasta(path):
    if not os.path.exists(path):
        sys.exit(f"Alignment not found: {path}")
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


def aa_at(seq, pos1):
    """
    Residue positions are defined relative to MAFFT seed-aligned reference sequence 
    """
    idx = pos1 - 1
    return seq[idx] if 0 <= idx < len(seq) else "-"

def process_nir(records, status_fh, bona_fh):
    writer = csv.writer(status_fh)
    writer.writerow(
        ["accession"] + list(NIR_SITES.keys()) +
        ["T1_complete", "T2_complete", "catalytic_ok", "category", "bonafide"]
    )

    examined = passed = 0
    counts = {}

    for hid, seq in records:
        if hid in SEED_IDS:
            continue

        ungapped = seq.replace("-", "")
        if len(ungapped) < MIN_LEN["nir"]:
            continue
        if "*" in ungapped or "X" in ungapped:
            continue

        examined += 1

        res = {k: aa_at(seq, pos) for k, (pos, _) in NIR_SITES.items()}
        t1 = all(res[k] == NIR_SITES[k][1] for k in NIR_T1)
        t2 = all(res[k] == NIR_SITES[k][1] for k in NIR_T2)
        cat_ok = all(res[k] == NIR_SITES[k][1] for k in NIR_CAT)
        bonafide = t1 and t2

        if t1 and t2:
            category = "both"
        elif t1:
            category = "t1_only"
        elif t2:
            category = "t2_only"
        else:
            category = "neither"

        counts[category] = counts.get(category, 0) + 1

        writer.writerow(
            [hid] + [res[k] for k in NIR_SITES] +
            [int(t1), int(t2), int(cat_ok), category, int(bonafide)]
        )

        if bonafide:
            passed += 1
            bona_fh.write(f">{hid}\n{ungapped}\n")

    return examined, passed, counts


def process_pcu(records, status_fh, bona_fh):
    writer = csv.writer(status_fh)
    writer.writerow(
        ["accession"] + list(PCU_SITES.keys()) +
        ["motif_ok", "his_tail_count", "category", "bonafide"]
    )

    examined = passed = 0
    counts = {}

    for hid, seq in records:
        if hid in SEED_IDS:
            continue

        ungapped = seq.replace("-", "")
        if len(ungapped) < MIN_LEN["pcu"]:
            continue
        if "*" in ungapped or "X" in ungapped:
            continue

        examined += 1

        res = {k: aa_at(seq, pos) for k, (pos, _) in PCU_SITES.items()}
        motif_ok = all(res[k] == PCU_SITES[k][1] for k in PCU_MOTIF)
        his_tail = ungapped[-30:].count("H")
        site2 = his_tail >= 3

        if motif_ok and site2:
            category = "both"
        elif motif_ok:
            category = "site1_only"
        elif site2:
            category = "site2_only"
        else:
            category = "neither"

        bonafide = motif_ok
        counts[category] = counts.get(category, 0) + 1

        writer.writerow(
            [hid] + [res[k] for k in PCU_SITES] +
            [int(motif_ok), his_tail, category, int(bonafide)]
        )

        if bonafide:
            passed += 1
            bona_fh.write(f">{hid}\n{ungapped}\n")

    return examined, passed, counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True, choices=["nir", "pcu"])
    ap.add_argument("--aln", required=True)
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    records = list(read_fasta(args.aln))
    if not records:
        sys.exit("Alignment contains no sequences.")

    status_path = os.path.join(args.outdir, f"{args.name}_site_status.csv")
    bona_path = os.path.join(args.outdir, f"{args.name}_bonafide.fasta")

    with open(status_path, "w", newline="") as sfh, open(bona_path, "w") as ffh:
        if args.name == "nir":
            examined, passed, counts = process_nir(records, sfh, ffh)
        else:
            examined, passed, counts = process_pcu(records, sfh, ffh)
    pct = f"{100 * passed / examined:.1f}%" if examined else "N/A"

    print(f"{args.name.upper()} residue filter:")
    print(f"examined: {examined}")
    print(f"bona fide: {passed} ({pct})")
    print(f"status table: {status_path}")
    print(f" bona fide FASTA: {bona_path}")


if __name__ == "__main__":
    main()

