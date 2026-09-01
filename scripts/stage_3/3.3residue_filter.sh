#!/bin/bash
#SBATCH -p shared # SLURM partition/queue
#SBATCH -c 8                # number of CPU cores to allocate, one per thread, up to 128.
#SBATCH --mem=16G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:24G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 04:00:00     # time limit dd-hh:mm:ss
#SBATCH --job-name=3.3residue_filter # SLURM job identifier
#SBATCH --output=../log/3.3residue_filter%A_%a.out
#SBATCH --error=../log/3.3residue_filter%A_%a.err

#SBATCH --array=1-2
# Array task 1 = NirK
# Array task 2 = PCuAC

# SLURM email notifications are sent to the address specified below.
#SBATCH --mail-type=ALL # Enable SLURM email notifications
#SBATCH --mail-user= removed email due to public repo.

# Residue filter --> pairwise-to-seed)
# For each hit, align it to the seed with
# MAFFT --add --keeplength

module load bioinformatics
module load mafft

cd ~/pcu-nir
mkdir -p out/aln

names=("" "nirk" "pcuac")
name=${names[$SLURM_ARRAY_TASK_ID]}
seed=query/${name}_seed.fasta
hits=out/${name}_hits.fasta
aln=out/aln/${name}_to_seed.aln
echo "Aligning $hits onto $seed with MAFFT --add --keeplength"

# --keeplength needs the seed as the existing alignment and the hits added to it.
mafft --thread 8 --keeplength --add "$hits" "$seed" > "$aln"
echo "Alignment done: $(grep -c '^>' $aln) records -> $aln"

# Sequence-count check: hits going in should equal (aligned records - seed
# records) coming out. Does not alter the alignment, just flags a mismatch.
n_hits=$(grep -c '^>' "$hits")
n_seed=$(grep -c '^>' "$seed")
n_aln=$(grep -c '^>' "$aln")
n_aln_hits=$((n_aln - n_seed))
echo "Record count check: hits_in=$n_hits  seed=$n_seed  aligned_excl_seed=$n_aln_hits"
if [ "$n_hits" -ne "$n_aln_hits" ]; then
  echo "WARNING: MAFFT record count mismatch for $name (expected $n_hits, got $n_aln_hits excluding seed)" >&2
fi

# Hand off to the Python 3.3 script tagger which reads the copper-site columns.
python3 scripts/3.3_tag_sites.py --name "$name" --aln "$aln" --outdir out
echo "Done at $(date)"
#runtime  MAFFT --add of tens of thousands of seqs: max.. ~1-3 h