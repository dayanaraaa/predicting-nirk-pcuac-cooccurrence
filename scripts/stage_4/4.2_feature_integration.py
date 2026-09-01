#!/usr/bin/env python3
"""
4.2_feature_integration.py
Integrate NirK and PCuAC sequence features with the species-level co-occurrence table.

Protein features are calculated using Biopython's ProteinAnalysis module, including molecular weight,
isoelectric point, GRAVY, aromaticity, and amino acid composition.

Sequences are cleaned to retain standard amino acids before feature calculation. Any removed residues
are recorded by accession.
Sequence length features use the original FASTA sequence length and are not affected by cleaning.

Inputs:
- out/species_dataset.csv
- out/nirk.fasta
- out/pcuac_functional.fasta

Output:
- out/final_feature_table.csv
- out/pipeline_metadata/4.2_feature_integration_metadata.txt

One row = one species-level observation.
Representative accessions are taken from the species table.
"""
import sys
print(sys.executable)

import Bio
print(Bio.__version__)

import argparse
import csv
import os
import sys
import platform
import datetime
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import Bio

PIPELINE_VERSION = "stage4.v2"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path):
    sequences = {}
    acc = None
    seq = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if acc:
                    sequences[acc] = "".join(seq)
                acc = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if acc:
            sequences[acc] = "".join(seq)
    return sequences


def sanitize_sequence(seq, accession):
    """
    Strip non-standard residues (X, B, *, gaps, etc.) before
    handing the sequence to ProteinAnalysis
    """
    clean = "".join(a for a in seq.upper() if a in VALID_AA)
    n_removed = len(seq) - len(clean)
    if n_removed > 0:
        print(f"WARNING: {accession}: removed {n_removed} non-standard "
              f"residue(s) before feature calculation")
    return clean

def sequence_features(seq, prefix, accession):
    raw_len = len(seq)
    clean_seq = sanitize_sequence(seq, accession)
    analysis = ProteinAnalysis(clean_seq)

    features = {}
    features[prefix + "_length"] = raw_len
    features[prefix + "_molecular_weight"] = round(analysis.molecular_weight(), 2)
    features[prefix + "_isoelectric_point"] = round(analysis.isoelectric_point(), 3)
    features[prefix + "_gravy"] = round(analysis.gravy(), 4)
    features[prefix + "_aromaticity"] = round(analysis.aromaticity(), 5)

    # NOTE: get_amino_acids_percent() was removed in newer Biopython versions.
    # Replaced it with an equivalent manual calculation: count / length of the
    # same (cleaned) sequence ProteinAnalysis was built from. This preserves
    # the original feature definition exactly.
    aa_counts = analysis.count_amino_acids()
    aa_percent = {
        aa: count / len(clean_seq)
        for aa, count in aa_counts.items()
    }

    for aa in "ACDEFGHIKLMNPQRSTVWY":
        features[f"{prefix}_aa_frac_{aa}"] = round(aa_percent.get(aa, 0.0), 5)

    return features

def write_pipeline_metadata(outdir, script_name, input_files, output_files, notes=""):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{script_name}_metadata.txt")
    with open(path, "w") as fh:
        fh.write("Stage: 4\n")
        fh.write(f"Script: {script_name}\n")
        fh.write(f"Pipeline version: {PIPELINE_VERSION}\n\n")
        fh.write(f"Date generated: {datetime.datetime.now().isoformat(timespec='seconds')}\n\n")
        fh.write(f"Python: {platform.python_version()}\n")
        fh.write(f"Biopython: {Bio.__version__}\n\n")
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="out/species_dataset.csv")
    ap.add_argument("--nirk", default="out/nirk.fasta")
    ap.add_argument("--pcuac", default="out/pcuac_functional.fasta")
    ap.add_argument("--out", default="out/final_feature_table.csv")
    args = ap.parse_args()

    print("Loading FASTA files")
    nirK_sequences = read_fasta(args.nirk)
    pcuac_sequences = read_fasta(args.pcuac)

    print("NirK sequences:", len(nirK_sequences))
    print("PCuAC sequences:", len(pcuac_sequences))

    output_rows = []
    with open(args.species) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            features = {}

            nirK_acc = row["nirK_rep_accession"].strip()
            pcuac_acc = row["pcuac_rep_accession"].strip()

            if nirK_acc.upper() in {"NA", "N/A", ""}:
                nirK_acc = None

            if pcuac_acc.upper() in {"NA", "N/A", ""}:
                pcuac_acc = None

            if nirK_acc:
                if nirK_acc in nirK_sequences:
                    features.update(
                        sequence_features(nirK_sequences[nirK_acc], "nirK", nirK_acc)
                    )
                else:
                    print(f"WARNING: NirK accession missing from FASTA: {nirK_acc}")

            if pcuac_acc:
                if pcuac_acc in pcuac_sequences:
                    features.update(
                        sequence_features(pcuac_sequences[pcuac_acc], "pcuac", pcuac_acc)
                    )
                else:
                    print(f"WARNING: PCuAC accession missing from FASTA: {pcuac_acc}")

            row.update(features)
            output_rows.append(row)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    all_fields = []
    for row in output_rows:
        for key in row.keys():
            if key not in all_fields:
                all_fields.append(key)
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=all_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_rows)

    print()
    print("Finished")
    print("Rows written:", len(output_rows))
    print("Output:", args.out)

    write_pipeline_metadata(
        outdir="out/pipeline_metadata",
        script_name="4.2_feature_integration",
        input_files=[args.species, args.nirk, args.pcuac],
        output_files=[args.out],
        notes="Feature calculation via Bio.SeqUtils.ProtParam.ProteinAnalysis.",
    )

if __name__ == "__main__":
    main()
