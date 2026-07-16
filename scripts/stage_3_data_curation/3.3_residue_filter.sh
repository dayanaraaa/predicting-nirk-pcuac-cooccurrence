#!/bin/bash
#SBATCH -p shared # SLURM partition/queue
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=tmp:32G
#SBATCH -t 08:00:00
#SBATCH --job-name=3.3_residue_filter # SLURM job identifier
#SBATCH --output=log/3.3_residue_filter_%A_%a.out
#SBATCH --error=log/3.3_residue_filter_%A_%a.err

# SLURM email notifications are sent to the address specified below.
#SBATCH --mail-type=ALL # Enable SLURM email notifications
#SBATCH --mail-user=your_email@example.com

# Objective:
# Align NirK/PCuAC BLAST hits to their seed proteins and evaluate conserved
# copper-binding residues to identify bona fide homologues.

# Inputs:
#   query/nir_seed.fasta
#   query/pcu_seed.fasta
#   out/nir_blast_hits.fasta
#   out/pcu_blast_hits.fasta

# Outputs:
#   out/nir_site_status.csv
#   out/pcu_site_status.csv
#   out/nir_bonafide.fasta
#   out/pcu_bonafide.fasta

# Parameters:
#   TASK=1 → NirK
#   TASK=2 → PCuAC
#   MAFFT --keeplength alignment
#   Biological thresholds:
#   NirK: min length 300 aa, Type‑1/Type‑2 ligand positions, catalytic dyad
#   PCuAC: min length 100 aa, Cu(I) motif positions, His-tail ≥ 3

set -euo pipefail

if [[ -z "${TASK:-}" ]]; then
    echo "TASK not set. Use:"
    echo "  TASK=1 sbatch 3.3_residue_filter.sh   # NirK"
    echo "  TASK=2 sbatch 3.3_residue_filter.sh   # PCuAC"
    exit 1
fi

if [[ "$TASK" == 1 ]]; then
    name="nir"
elif [[ "$TASK" == 2 ]]; then
    name="pcu"
else
    echo "TASK must be 1 or 2"
    exit 1
fi

module load bioinformatics || true
module load mafft
module load python || module load python/3.10 || {
    echo "Python module not available."
    exit 1
}

root="$(pwd)"

seed="${root}/query/${name}_seed.fasta"
hits="${root}/out/${name}_blast_hits.fasta"
aln="${root}/out/aln/${name}_to_seed.aln"

mkdir -p "${root}/out/aln"

if [[ ! -f "$seed" || ! -f "$hits" ]]; then
    echo "Missing required input: $seed or $hits"
    exit 1
fi

n_hits=$(grep -c '^>' "$hits" || echo 0)
if [[ "$n_hits" -eq 0 ]]; then
    echo "No hits found for ${name}; nothing to process."
    exit 0
fi

echo "Aligning ${name} hits to seed..."
mafft \
    --thread "$SLURM_CPUS_PER_TASK" \
    --keeplength \
    --add "$hits" \
    "$seed" \
> "$aln"

echo "Tagging copper-site residues..."
python3 src/3.3_tag_sites.py \
    --name "$name" \
    --aln "$aln" \
    --outdir out

echo "Residue filtering complete for ${name}."
