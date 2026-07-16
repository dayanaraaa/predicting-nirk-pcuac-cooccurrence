#!/bin/bash
#SBATCH -p shared # SLURM partition/queue
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --gres=tmp:12G
#SBATCH -t 02:00:00
#SBATCH --job-name=3.2_fetch_sequences # SLURM job identifier
#SBATCH --output=log/3.2_fetch_sequences%A_%a.out
#SBATCH --error=log/3.2_fetch_sequences%A_%a.err
#SBATCH --array=1-2

# SLURM email notifications are sent to the address specified below.
#SBATCH --mail-type=ALL # Enable SLURM email notifications
#SBATCH --mail-user=your_email@example.com

# Objective:
# Fetch FASTA sequences for NirK and PCuAC BLAST hits
# Requires local RefSeq protein BLAST database generated in Stage 2

module load bioinformatics
module load blast/2.17.0

cd /nobackup/$USER/pcu-nir || {
  echo "Failed to enter project directory."
  exit 1
}

export BLASTDB=$PWD/db

# Array task 1 = nir, task 2 = pcu
names=("" "nir" "pcu")
name=${names[$SLURM_ARRAY_TASK_ID]}

metadata=out/${name}_blast_metadata_clean.tsv
acc=out/${name}_accessions_clean.txt
fasta=out/${name}_blast_hits.fasta

echo "[3.2] Fetching ${name} sequences from ${metadata}"

# Extract cleaned accession IDs from script 3.1 metadata
cut -f1 "$metadata" \
  | tail -n +2 \
  | sort -u \
  | grep -v '^$' > "$acc"

echo "[3.2] Unique accessions: $(wc -l < $acc)"

# Prevent empty retrievals running
if [ ! -s "$acc" ]; then
    echo "ERROR: No accession IDs found for ${name}"
    exit 1
fi

# Retrieve sequences from local RefSeq database
blastdbcmd \
  -db db/refseq_protein \
  -entry_batch "$acc" \
  -target_only \
  -out "$fasta" 2> >(tee out/${name}_blastdbcmd.err >&2)

echo "[3.2] Sequences written: $(grep -c '^>' $fasta)"
echo "[3.2] Missing IDs listed in out/${name}_blastdbcmd.err"
echo "[3.2] Done at $(date)"
