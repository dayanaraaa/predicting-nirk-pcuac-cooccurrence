#!/bin/bash
#SBATCH -p shared # SLURM partition/queue
#SBATCH -c 2                # number of CPU cores to allocate, one per thread, up to 128.
#SBATCH --mem=8G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:12G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 02:00:00     # time limit dd-hh:mm:ss
#SBATCH --job-name=3.2_fetch_sequences # SLURM job identifier
#SBATCH --output=../log/3.2_fetch_sequences%A_%a.out
#SBATCH --error=../log/3.2_fetch_sequences%A_%a.err
#SBATCH --array=1-2

# SLURM email notifications are sent to the address specified below.
#SBATCH --mail-type=ALL # Enable SLURM email notifications
#SBATCH --mail-user= removed email due to public repo.

# Fetch the actual amino-acid sequences for every BLAST hit, from the local
# refseq_protein database, using blastdbcmd. These FASTAs feed the residue
# filter (3.3) and downstream, the ML featurisation.

# Array task 1 = NirK hits
# Array task 2 = PCuAC hits

module load bioinformatics
module load blast/2.17.0

PROJECT_ROOT=/nobackup/$USER/pcu-nir
cd "$PROJECT_ROOT" || {
  echo "Failed to enter project directory: $PROJECT_ROOT"
  exit 1
}
export BLASTDB=$PROJECT_ROOT/db

names=("" "nirk" "pcuac")
name=${names[$SLURM_ARRAY_TASK_ID]}

acc=out/${name}_accessions_clean.txt
fasta=out/${name}_hits.fasta

echo "Fetching $name sequences using cleaned accession list: $acc"

if [ ! -f "$acc" ]; then
    echo "ERROR: accession list not found: $acc"
    exit 1
fi

echo "Unique accessions: $(wc -l < "$acc")"

# Pull the sequences
# -target_only keeps one FASTA record per accession even
# when that record maps to many taxids.
# Missing IDs are reported
blastdbcmd \
  -db db/refseq_protein \
  -entry_batch "$acc" \
  -target_only \
  -out "$fasta" 2> >(tee out/${name}_blastdbcmd.err >&2)

echo "Done at $(date)"
echo "Sequences written: $(grep -c '^>' $fasta)  ->  $fasta"
echo "Any IDs not found are listed in out/${name}_blastdbcmd.err"
#runtime ~minutes (local DB lookup)