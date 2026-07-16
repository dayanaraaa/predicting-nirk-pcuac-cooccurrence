#!/bin/bash
#SBATCH -p shared # You select the queue(cluster) here
#SBATCH -c 16                # number of CPU cores to allocate, one per thread, up to 128.
#SBATCH --mem=64G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:96G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 1-00:00:00   # time limit dd-hh:mm:ss (refseq is large; allow at least a day)
#SBATCH --job-name=2.2_blastp_search # This name will let you follow your job
#SBATCH --output=log/2.2_blastp_search%A_%a.out
#SBATCH --error=log/2.2_blastp_search%A_%a.err
#SBATCH --array=1-2

# SLURM email notifications are sent to the address specified below.
#SBATCH --mail-type=ALL # Enable SLURM email notifications
#SBATCH --mail-user=your_email@example.com

# This is the run that yields a
# real dataset — hundreds-to-thousands of bacterial Nir and Pcu homologues.
# Array task 1 = Nir (inputs) ; task 2 = Pcu (labels).

module load bioinformatics
module load blast/2.17.0

mkdir -p out
mkdir -p log
export BLASTDB=$PWD/db   # so -taxids expands to all bacteria + sscinames resolves

query_names=("" "nir_seed" "pcu_seed")  # index 0 empty so tasks start at 1
query=${query_names[$SLURM_ARRAY_TASK_ID]}
echo "BLASTp query = $query  (task $SLURM_ARRAY_TASK_ID) against refseq_protein"

blastp \
  -query query/${query}.fasta \
  -db db/refseq_protein \
  -taxids 2 \
  -evalue 1e-5 \
  -max_target_seqs 100000 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle" \
  -num_threads 16 \
  -out out/${query}_refseq_hits.tsv

echo "Done at $(date)"
echo "Hits for ${query}: $(wc -l < out/${query}_refseq_hits.tsv)"
#runtime  minutes-to-hours depending on hit count