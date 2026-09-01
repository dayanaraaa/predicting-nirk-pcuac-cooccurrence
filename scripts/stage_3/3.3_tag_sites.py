#!/usr/bin/env python3
"""
3.3_tag_sites.py
Checks conserved copper-binding residues in MAFFT --keeplength alignments
and labels NirK/PCuAC candidate sequences based on functional site presence.

Called by 3.3_residue_filter.sh after candidate sequences are aligned onto
the reference seed coordinate frame.

For NirK, AniA (seed: Q02219):
- Type 1 copper site: His, Cys (HC motif), His, Met
- Type 2 copper site: three conserved His residues (sequence signature only;
  this does not confirm the structural, partly inter-subunit, T2 site)
- Catalytic Asp/His residues are also recorded as an additional confidence check

Bona fide NirK = if both Type 1 and Type 2 copper sites are present

For PCuAC, AccA (seed: A0AAQ1E0N0):
- Checks conservation of the Cu(I)-binding motif:
    His69, Met80, His103, Met105
- This script validates PCuAC copper-binding residues only.
- Signal peptide validation is performed separately.

The reference seed record is identified (robust to version suffixes and
BLAST db|ACC|NAME wrappers), always retained in the outputs, and flagged
with is_reference_seed=1. It is excluded from candidate hit statistics
(n_total/n_pass/category counts), since it is a positive control, not an
unknown hit.

Outputs:
- out/nirk_site_status.csv
- out/pcuac_site_status.csv

- out/nirk_bonafide.fasta
- out/pcuac_bonafide.fasta
"""
import argparse
import csv
import os
import sys

# 1-based seed positions for each site
NIR_SITES = {
    "T1_134": (134, "H"), "T1_175": (175, "C"), "T1_183": (183, "H"), "T1_188": (188, "M"),
    "T2_139": (139, "H"), "T2_174": (174, "H"), "T2_329": (329, "H"),
    "cat_137": (137, "D"), "cat_280": (280, "H"),
}
NIR_T1 = ["T1_134", "T1_175", "T1_183", "T1_188"]
NIR_T2 = ["T2_139", "T2_174", "T2_329"]
NIR_CAT = ["cat_137", "cat_280"]

PCU_SITES = {
    "His69": (69, "H"), "Met80": (80, "M"), "His103": (103, "H"), "Met105": (105, "M"),
}
PCU_MOTIF = ["His69", "Met80", "His103", "Met105"]

# Reference seed accessions (base accession, version stripped at match time)
NIR_SEED_ACC = "Q02219"       # AniA
PCU_SEED_ACC = "A0AAQ1E0N0"   # AccA


def read_fasta(path):
    """
    Gives (header_id, sequence).
    header_id = first token after '>'.
    """
    if not os.path.exists(path):
        sys.exit(f"ERROR: alignment not found: {path}")
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

def aa_at(seq, pos1):
    """
    Residue at 1-based position; '-' or '?' if out of range.
    """
    i = pos1 - 1
    return seq[i] if 0 <= i < len(seq) else "-"

def matches_seed(hid, seed_acc):
    """
    True if hid identifies the reference seed accession.
    Robust to version suffixes (ACC.1) and BLAST db|ACC|NAME wrappers
    (e.g. sp|Q02219|ANIA_NEIGO), since header formatting can vary
    between the manually curated seed FASTA and the BLAST-derived hits.
    """
    for token in hid.split("|"):
        if token.split(".")[0] == seed_acc:
            return True
    return False

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True, choices=["nirk", "pcuac"])
    ap.add_argument("--aln", required=True, help="MAFFT --keeplength alignment (seed + hits)")
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()

    records = list(read_fasta(args.aln))
    if not records:
        sys.exit("ERROR: no records in alignment")

    seed_acc = NIR_SEED_ACC if args.name == "nirk" else PCU_SEED_ACC
    seed_found = False
    seed_bonafide = None

    status_path = os.path.join(args.outdir, f"{args.name}_site_status.csv")
    fasta_path = os.path.join(args.outdir, f"{args.name}_bonafide.fasta")

    n_total = n_pass = 0
    cat_counts = {}

    with open(status_path, "w", newline="") as sfh, open(fasta_path, "w") as ffh:
        w = csv.writer(sfh)

        if args.name == "nirk":
            cols = (["accession"] + list(NIR_SITES.keys())
                    + ["T1_complete", "T2_complete", "catalytic_ok", "category",
                       "bonafide", "is_reference_seed"])
            w.writerow(cols)
            for hid, seq in records:
                is_seed = matches_seed(hid, seed_acc)

                res = {k: aa_at(seq, p) for k, (p, _) in NIR_SITES.items()}
                t1 = all(res[k] == NIR_SITES[k][1] for k in NIR_T1)
                t2 = all(res[k] == NIR_SITES[k][1] for k in NIR_T2)
                cat = all(res[k] == NIR_SITES[k][1] for k in NIR_CAT)
                # Category: both/t1_only/t2_only/neither
                if t1 and t2:
                    category = "both"
                elif t1:
                    category = "t1_only"
                elif t2:
                    category = "t2_only"
                else:
                    category = "neither"
                bonafide = t1 and t2

                if is_seed:
                    seed_found = True
                    seed_bonafide = bonafide
                else:
                    n_total += 1
                    cat_counts[category] = cat_counts.get(category, 0) + 1
                    if bonafide:
                        n_pass += 1

                w.writerow([hid] + [res[k] for k in NIR_SITES]
                           + [int(t1), int(t2), int(cat), category, int(bonafide), int(is_seed)])
                if bonafide or is_seed:
                    ffh.write(f">{hid}\n{seq.replace('-', '')}\n")

        else:  # PCuAC
            cols = (["accession"] + list(PCU_SITES.keys())
                    + ["pcuac_motif_present", "category", "bonafide", "is_reference_seed"])
            w.writerow(cols)
            for hid, seq in records:
                is_seed = matches_seed(hid, seed_acc)

                res = {k: aa_at(seq, p) for k, (p, _) in PCU_SITES.items()}
                pcuac_motif_present = all(res[k] == PCU_SITES[k][1] for k in PCU_MOTIF)
                bonafide = pcuac_motif_present
                category = "pcuac" if bonafide else "not_pcuac"

                if is_seed:
                    seed_found = True
                    seed_bonafide = bonafide
                else:
                    n_total += 1
                    cat_counts[category] = cat_counts.get(category, 0) + 1
                    if bonafide:
                        n_pass += 1

                w.writerow([hid] + [res[k] for k in PCU_SITES]
                           + [int(pcuac_motif_present), category, int(bonafide), int(is_seed)])
                if bonafide or is_seed:
                    ffh.write(f">{hid}\n{seq.replace('-', '')}\n")

    if not seed_found:
        sys.exit(f"ERROR: expected reference seed ({seed_acc}) not found in {args.aln}."
                  f"Check the seed FASTA header / query/{args.name}_seed.fasta.")

    print("-" * 50)
    print(f"Residue filter -- {args.name.upper()}")
    print("-" * 50)
    print()
    print(f"Reference seed detected: {seed_acc}")
    if seed_bonafide is False:
        print(f"WARNING: reference seed {seed_acc} does NOT satisfy its own site "
              f"definition in this alignment - check residue positions in this script.")
    print(f"Hits examined: {n_total}")
    print(f"Bona fide (kept): {n_pass} ({100*n_pass/n_total:.1f}%)" if n_total else "no hits")
    print(f"Category breakdown: {cat_counts}")
    print(f"Wrote {status_path}")
    print(f"Wrote {fasta_path}")
    print("-" * 50)

if __name__ == "__main__":
    main()