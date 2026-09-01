#!/usr/bin/env python3
"""
4.0_build_sequence_master_table.py
Build sequence-level and species-level NirK-PCuAC tables from validated protein accessions.

Uses cleaned accession metadata from 3.1_parse_blast_candidates.py as the single source of
truth for accession-to-organism mapping.

Species-level tables:
- Resolved species are identified from valid Genus species names.
- Unresolved annotations are retained for QC but excluded from
  downstream ML datasets by 4.1.

Representative sequence selection:
- Longest validated sequence per species and protein family.
- Ties resolved by accession order for reproducibility.
- Representative sequences are flagged for traceability.

Inputs:
- out/nirk_accessions.txt
- out/pcuac_functional_accessions.txt
- out/nirk_candidate_metadata.tsv
- out/pcuac_candidate_metadata.tsv
- out/nirk.fasta
- out/pcuac_functional.fasta

Outputs:
- out/sequence_master_table.csv (one row per validated sequence, QC/audit)
- out/species_master_table.csv (one row per species, incl. unresolved, QC)
- out/pipeline_metadata/4.0_build_sequence_master_table_metadata.txt
"""

import csv
import os
import sys
import platform
import datetime
from collections import defaultdict

PIPELINE_VERSION = "stage4.v2"


def read_accession_list(path):
    accs = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                accs.add(line)
    return accs


def read_candidate_metadata(path):
    """
    Reads the cleaned candidate metadata table from 3.1_parse_blast_candidates.py.
    Returns accession -> taxids (string), accession -> organism (string).
    """
    if not os.path.exists(path):
        sys.exit(f"ERROR: candidate metadata not found: {path}")

    acc_to_taxids = {}
    acc_to_org = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"accession", "taxids", "organism"}
        if not required.issubset(set(reader.fieldnames or [])):
            sys.exit(f"ERROR: {path} missing required column(s), "
                     f"expected {required}, found {reader.fieldnames}")
        for row in reader:
            acc = row["accession"].strip()
            if not acc:
                continue
            acc_to_taxids[acc] = row["taxids"].strip()
            acc_to_org[acc] = row["organism"].strip()
    return acc_to_taxids, acc_to_org


def read_fasta_lengths(path):
    lengths = {}
    current = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    lengths[current] = len("".join(seq))
                current = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if current:
            lengths[current] = len("".join(seq))
    return lengths

def extract_species(org, accession):
    """
    Parse a clean 'genus species' key from organism metadata.

    Returns (species_key, species_resolved).

    !!! species_resolved = False when the organism is missing or only
    identified to genus level ("Genus sp."), since these cannot be
    safely collapsed into a shared species-level observation without
    risking merging unrelated organisms together.
    """
    if not org or org == "N/A":
        return f"unresolved:{accession}", False

    fallback = None
    for name in org.split(";"):
        name = name.strip()
        parts = name.split()

        if len(parts) >= 2:
            genus, species = parts[0], parts[1]
            if species.lower() not in ("sp.", "spp."):
                return f"{genus} {species}".lower(), True
            elif fallback is None:
                fallback = name.lower()

    if fallback is not None:
        return fallback, False
    return org.lower(), False

def select_representative(rows):
    """
    Representative sequence per species/protein family: longest
    validated sequence, ties broken by accession alphabetical order.
    """
    return sorted(rows, key=lambda r: (-r["length"], r["accession"]))[0]


def write_pipeline_metadata(outdir, script_name, input_files, output_files, notes=""):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{script_name}_metadata.txt")
    with open(path, "w") as fh:
        fh.write("Stage: 4\n")
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
    nirK_list = "out/nirk_accessions.txt"
    pcuac_list = "out/pcuac_functional_accessions.txt"

    nirK_metadata = "out/nirk_candidate_metadata.tsv"
    pcuac_metadata = "out/pcuac_candidate_metadata.tsv"

    nirK_fasta = "out/nirk.fasta"
    pcuac_fasta = "out/pcuac_functional.fasta"

    if not os.path.exists(nirK_list):
        sys.exit("ERROR: NirK accession list missing.")
    if not os.path.exists(pcuac_list):
        sys.exit("ERROR: PCuAC accession list missing.")

    nirK_lengths = read_fasta_lengths(nirK_fasta)
    pcuac_lengths = read_fasta_lengths(pcuac_fasta)

    nirK_accs = sorted(read_accession_list(nirK_list))
    pcuac_accs = sorted(read_accession_list(pcuac_list))

    nirK_taxids, nirK_org = read_candidate_metadata(nirK_metadata)
    pcuac_taxids, pcuac_org = read_candidate_metadata(pcuac_metadata)

    print(f"NirK accessions: {len(nirK_accs)}")
    print(f"NirK metadata matches: {len(set(nirK_accs) & set(nirK_taxids))}")
    missing_nirK = sorted(set(nirK_accs) - set(nirK_taxids))
    if missing_nirK:
        print(f"WARNING: {len(missing_nirK)} NirK accession(s) have no 3.1 metadata "
              f"entry and will be dropped: {missing_nirK[:5]}{' ...' if len(missing_nirK) > 5 else ''}")

    print(f"PCuAC accessions: {len(pcuac_accs)}")
    print(f"PCuAC metadata matches: {len(set(pcuac_accs) & set(pcuac_taxids))}")
    missing_pcuac = sorted(set(pcuac_accs) - set(pcuac_taxids))
    if missing_pcuac:
        print(f"WARNING: {len(missing_pcuac)} PCuAC accession(s) have no 3.1 metadata "
              f"entry and will be dropped: {missing_pcuac[:5]}{' ...' if len(missing_pcuac) > 5 else ''}")

    sequence_rows = []

    for acc in nirK_accs:
        if acc not in nirK_taxids:
            continue
        species, resolved = extract_species(nirK_org.get(acc, "N/A"), acc)
        sequence_rows.append({
            "accession": acc,
            "protein_type": "NirK",
            "length": nirK_lengths.get(acc, 0),
            "taxid": nirK_taxids[acc],
            "organism_name": nirK_org.get(acc, "N/A"),
            "species": species,
            "species_resolved": int(resolved),
            "is_representative": False,
            "representative_selection_reason": "",
        })

    for acc in pcuac_accs:
        if acc not in pcuac_taxids:
            continue
        species, resolved = extract_species(pcuac_org.get(acc, "N/A"), acc)
        sequence_rows.append({
            "accession": acc,
            "protein_type": "PCuAC",
            "length": pcuac_lengths.get(acc, 0),
            "taxid": pcuac_taxids[acc],
            "organism_name": pcuac_org.get(acc, "N/A"),
            "species": species,
            "species_resolved": int(resolved),
            "is_representative": False,
            "representative_selection_reason": "",
        })

    # Collapse validated sequences into species-level table.
    # Unresolved species keys are unique per accession/annotation
    # (see extract_species), so this never merges unrelated organisms
    species_data = defaultdict(lambda: {
        "taxid": set(),
        "organism_name": None,
        "species_resolved": None,
        "NirK": [],
        "PCuAC": []
    })

    for row in sequence_rows:
        sp = row["species"]
        species_data[sp]["taxid"].add(row["taxid"])
        species_data[sp]["organism_name"] = row["organism_name"]
        species_data[sp]["species_resolved"] = row["species_resolved"]

        if row["protein_type"] == "NirK":
            species_data[sp]["NirK"].append(row)
        elif row["protein_type"] == "PCuAC":
            species_data[sp]["PCuAC"].append(row)

    species_rows = []

    for species, data in species_data.items():
        nirk_present = len(data["NirK"]) > 0
        pcuac_present = len(data["PCuAC"]) > 0

        if nirk_present and pcuac_present:
            class_label = "both"
        elif nirk_present:
            class_label = "nirK_only"
        elif pcuac_present:
            class_label = "pcuac_only"
        else:
            continue

        nirK_rep = select_representative(data["NirK"]) if nirk_present else None
        pcuac_rep = select_representative(data["PCuAC"]) if pcuac_present else None

        if nirK_rep is not None:
            nirK_rep["is_representative"] = True
            nirK_rep["representative_selection_reason"] = "longest_validated_sequence"
        if pcuac_rep is not None:
            pcuac_rep["is_representative"] = True
            pcuac_rep["representative_selection_reason"] = "longest_validated_sequence"

        species_rows.append({
            "species": species,
            "taxid": ";".join(sorted(data["taxid"])),
            "organism_name": data["organism_name"],
            "nirK_present": int(nirk_present),
            "pcuac_present": int(pcuac_present),
            "class": class_label,
            "nirK_rep_accession": nirK_rep["accession"] if nirK_rep else "NA",
            "pcuac_rep_accession": pcuac_rep["accession"] if pcuac_rep else "NA",
            "species_resolved": data["species_resolved"],
        })

    # Write sequence-level audit table
    sequence_out = "out/sequence_master_table.csv"
    with open(sequence_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "accession", "protein_type", "length", "taxid", "organism_name",
            "species", "species_resolved", "is_representative",
            "representative_selection_reason",
        ])
        writer.writeheader()
        for row in sequence_rows:
            out_row = dict(row)
            out_row["is_representative"] = "TRUE" if row["is_representative"] else "FALSE"
            writer.writerow(out_row)

    # Write species-level QC table (includes unresolved species)
    species_out = "out/species_master_table.csv"
    with open(species_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "species", "taxid", "organism_name", "nirK_present", "pcuac_present",
            "class", "nirK_rep_accession", "pcuac_rep_accession", "species_resolved",
        ])
        writer.writeheader()
        writer.writerows(species_rows)

    n_resolved = sum(r["species_resolved"] == 1 for r in species_rows)

    print("-" * 50)
    print("Sequence master table")
    print(f"NirK sequences: {sum(x['protein_type']=='NirK' for x in sequence_rows)}")
    print(f"PCuAC sequences: {sum(x['protein_type']=='PCuAC' for x in sequence_rows)}")
    print(f"Rows written: {len(sequence_rows)}")
    print(sequence_out)
    print()
    print("Species master table (QC, includes unresolved)")
    print(f"Species rows: {len(species_rows)}")
    print(f"Resolved: {n_resolved}  Unresolved: {len(species_rows) - n_resolved}")
    print(f"Both: {sum(x['class']=='both' for x in species_rows)}")
    print(f"NirK only: {sum(x['class']=='nirK_only' for x in species_rows)}")
    print(f"PCuAC only: {sum(x['class']=='pcuac_only' for x in species_rows)}")
    print(species_out)
    print("-" * 50)

    write_pipeline_metadata(
        outdir="out/pipeline_metadata",
        script_name="4.0_build_sequence_master_table",
        input_files=[nirK_list, pcuac_list, nirK_metadata, pcuac_metadata, nirK_fasta, pcuac_fasta],
        output_files=[sequence_out, species_out],
        notes=(
            "Representative selection rule: longest validated sequence per species,\n"
            "if tied, then selected by accession alphabetical order.\n"
            "species_resolved = 0 for missing organism metadata or unresolved 'Genus sp.' annotations."
        ),
    )


if __name__ == "__main__":
    main()