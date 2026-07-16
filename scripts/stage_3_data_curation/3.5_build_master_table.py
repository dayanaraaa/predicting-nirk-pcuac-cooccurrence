#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

BLAST_COLS = [
    "accession",
    "qseqid",
    "pident",
    "length",
    "bitscore",
    "evalue",
    "staxids",
    "sscinames",
    "stitle",
]

BINOMIAL_RE = re.compile(r"^([A-Z][a-z]+)\s+([a-z]+)")

def extract_species(name):
    """
    Standardise BLAST organism names to Genus species.
    """
    if not name:
        return "Unknown"
    name = name.split(";")[0].strip()
    name = re.sub(r"^(MULTISPECIES:\s*)", "", name)
    m = BINOMIAL_RE.search(name)
    if m:
        return f"{m.group(1)} {m.group(2)}"
    parts = name.split()
    return " ".join(parts[:2]) if len(parts) >= 2 else "Unknown"


def read_clean_ids(path):
    if not path or not os.path.exists(path):
        return None
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def read_status_csv(path):
    if not os.path.exists(path):
        sys.exit(f"Missing status CSV: {path}")
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            if row.get("bonafide") == "1":
                out[row["accession"]] = row
    return out


def read_blast(path):
    """
    Load BLAST TSV and keep best hit per accession.
    """
    if not os.path.exists(path):
        sys.exit(f"Missing BLAST TSV: {path}")
    best = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < len(BLAST_COLS):
                continue
            d = dict(zip(BLAST_COLS, parts))
            acc = d["accession"]
            if acc not in best or float(d["bitscore"]) > float(best[acc]["bitscore"]):
                d["accession"] = acc
                best[acc] = d
    return best


def sort_key(hit):
    ev = float(hit.get("evalue", 1.0))
    bs = float(hit.get("bitscore", 0))
    return (ev, -bs)


def build_taxid_map(blast_hits, bona_fide):
    """
    Maps taxid to the best bona fide hit.
    """
    out = {}
    for acc, hit in blast_hits.items():
        if acc not in bona_fide:
            continue
        taxids = re.split(r"[;,]", hit.get("staxids", ""))
        for tid in taxids:
            tid = tid.strip()
            if not tid:
                continue
            if tid not in out or sort_key(hit) < sort_key(out[tid]):
                hit["taxid"] = tid
                out[tid] = hit
    return out

def collapse_species(taxid_map):
    """
    Collapse taxids to species-level representatives.
    """
    groups = defaultdict(list)
    for tid, hit in taxid_map.items():
        sp = extract_species(hit.get("sscinames", ""))
        groups[sp].append((tid, hit))
    species = {}
    for sp, entries in groups.items():
        rep_tid, rep_hit = min(entries, key=lambda x: sort_key(x[1]))
        species[sp] = {"rep_taxid": rep_tid, "hit": rep_hit}
    return species

def build_sequence_label_rows(protein, accs, blast_hits, species_class, rep_accession_by_species):
    """
    Build accession-level rows for sequence_label_map.csv.
    """
    rows = []
    unresolved = 0
    for acc in sorted(accs):
        hit = blast_hits.get(acc)
        if hit is None:
            unresolved += 1
            continue
        sp = extract_species(hit.get("sscinames", ""))
        cls = species_class.get(sp, "unresolved")
        if sp not in species_class:
            unresolved += 1
        is_rep = int(acc == rep_accession_by_species.get(sp))
        rows.append({
            "accession": acc,
            "species": sp,
            "class": cls,
            "protein": protein,
            "is_representative": is_rep,
        })
    if unresolved:
        print(f"  [{protein}] WARNING: {unresolved}/{len(accs)} cleaned accessions "
              f"could not be resolved to a known species/class in sequence_label_map.csv")
    return rows

def write_sequence_label_map(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["accession", "species", "class", "protein", "is_representative"]
        )
        w.writeheader()
        w.writerows(rows)


def write_species_table(path, nir_sp, pcu_sp):
    """
    Write species-level ML table.
    """
    species = sorted(set(nir_sp) | set(pcu_sp))
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "species", "rep_taxid",
            "nirK_present", "pcuac_present",
            "nirK_rep_accession", "pcuac_rep_accession",
            "nirK_bitscore", "pcuac_bitscore",
            "class",
        ])
        for sp in species:
            nr = nir_sp.get(sp, {})
            pr = pcu_sp.get(sp, {})
            n_pres = int(bool(nr))
            p_pres = int(bool(pr))
            cls = "both" if n_pres and p_pres else ("nirK_only" if n_pres else "pcuac_only")
            rep = nr.get("rep_taxid") or pr.get("rep_taxid") or "N/A"
            n_hit = nr.get("hit", {})
            p_hit = pr.get("hit", {})
            w.writerow([
                sp, rep,
                n_pres, p_pres,
                n_hit.get("accession", "N/A"),
                p_hit.get("accession", "N/A"),
                n_hit.get("bitscore", "N/A"),
                p_hit.get("bitscore", "N/A"),
                cls,
            ])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nir-status", default="out/nir_site_status.csv")
    ap.add_argument("--pcu-status", default="out/pcu_site_status.csv")
    ap.add_argument("--nir-blast", default="out/nir_blast_metadata_clean.tsv")
    ap.add_argument("--pcu-blast", default="out/pcu_blast_metadata_clean.tsv")
    ap.add_argument("--nir-clean-ids", default="out/nir_clean_ids.txt")
    ap.add_argument("--pcu-clean-ids", default="out/pcu_clean_ids.txt")
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    nir_status = read_status_csv(args.nir_status)
    pcu_status = read_status_csv(args.pcu_status)

    nir_blast = read_blast(args.nir_blast)
    pcu_blast = read_blast(args.pcu_blast)

    nir_accs = {acc for acc, row in nir_status.items() if row.get("bonafide") == "1"}
    pcu_accs = {acc for acc, row in pcu_status.items() if row.get("bonafide") == "1"}

    nir_clean_ids = read_clean_ids(args.nir_clean_ids)
    pcu_clean_ids = read_clean_ids(args.pcu_clean_ids)

    if nir_clean_ids is None:
        print(f"WARNING: {args.nir_clean_ids} not found")
    else:
        before = len(nir_accs)
        nir_accs &= nir_clean_ids
        print(f"[nir] bona fide accessions restricted to cleaned set: {before} -> {len(nir_accs)}")

    if pcu_clean_ids is None:
        print(f"WARNING: {args.pcu_clean_ids} not found")
    else:
        before = len(pcu_accs)
        pcu_accs &= pcu_clean_ids
        print(f"[pcu] bona fide accessions restricted to cleaned set: {before} -> {len(pcu_accs)}")

    nir_taxid = build_taxid_map(nir_blast, nir_accs)
    pcu_taxid = build_taxid_map(pcu_blast, pcu_accs)

    nir_species = collapse_species(nir_taxid)
    pcu_species = collapse_species(pcu_taxid)

    nir_species_names = set(nir_species.keys())
    pcu_species_names = set(pcu_species.keys())
    overlap = nir_species_names & pcu_species_names

    print("\nDIAGNOSTIC species consistency check")
    print("Unique NirK species:", len(nir_species_names))
    print("Unique PCuAC species:", len(pcu_species_names))
    print("Overlapping species:", len(overlap))

    nir_only = list(nir_species_names - pcu_species_names)[:10]
    pcu_only = list(pcu_species_names - nir_species_names)[:10]

    print("Example NirK-only species:", nir_only)
    print("Example PCuAC-only species:", pcu_only)

    species_path = os.path.join(args.outdir, "species_master_table.csv")
    write_species_table(species_path, nir_species, pcu_species)
    print(f"Species table written: {species_path}")

    all_species = set(nir_species) | set(pcu_species)
    species_class = {
        sp: ("both" if sp in nir_species and sp in pcu_species
             else ("nirK_only" if sp in nir_species else "pcuac_only"))
        for sp in all_species
    }

    nir_rep_by_species = {sp: g["hit"].get("accession") for sp, g in nir_species.items()}
    pcu_rep_by_species = {sp: g["hit"].get("accession") for sp, g in pcu_species.items()}

    label_rows = []
    label_rows += build_sequence_label_rows("nir", nir_accs, nir_blast, species_class, nir_rep_by_species)
    label_rows += build_sequence_label_rows("pcu", pcu_accs, pcu_blast, species_class, pcu_rep_by_species)

    label_map_path = os.path.join(args.outdir, "sequence_label_map.csv")
    write_sequence_label_map(label_map_path, label_rows)
    print(f"Sequence label map written: {label_map_path} ({len(label_rows)} accessions)")

    nir_only = 0
    pcu_only = 0
    both = 0
    neither = 0

    all_species = set(nir_species) | set(pcu_species)

    for sp in all_species:
        n = sp in nir_species
        p = sp in pcu_species
        if n and p:
            both += 1
        elif n:
            nir_only += 1
        elif p:
            pcu_only += 1
        else:
            neither += 1

    print("\nSpecies Summary")
    print(f"Total species: {len(all_species)}")
    print(f"NirK only: {nir_only}")
    print(f"PCuAC only: {pcu_only}")
    print(f"Both: {both}")
    print(f"Neither (in BLAST but failed QC): {neither}")


if __name__ == "__main__":
    main()
