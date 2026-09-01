#!/bin/bash
#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --gres=tmp:24G
#SBATCH -t 06:00:00
#SBATCH --job-name=3.5_signalp_pcuac
#SBATCH --output=../log/3.5_signalp_pcuac.out
#SBATCH --error=../log/3.5_signalp_pcuac.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user= removed email due to public repo.

module load bioinformatics
module load signalp6/6.0h

cd ~/pcu-nir

mkdir -p out/signalp_pcuac

signalp6 \
    --write_procs 8 \
    --fastafile out/pcuac_signalp_candidates.fasta \
    --organism other \
    --output_dir out/signalp_pcuac \
    --format txt \
    --mode fast

echo "SignalP finished: $(date)"

# Convert native SignalP6 output into the accession / YES|NO format
# consumed by 3.5_identify_functional_pcu.py
python3 scripts/3.5b_parse_signalp.py \
    --signalp-dir out/signalp_pcuac \
    --outdir out

echo "SignalP parsing finished: $(date)"