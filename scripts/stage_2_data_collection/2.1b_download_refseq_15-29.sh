#!/bin/bash
#SBATCH -p shared # SLURM partition/queue
#SBATCH -c 1                # number of CPU cores to allocate, one per thread, up to 128.
#SBATCH --mem=4G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:8G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 1-00:00:00   # time limit dd-hh:mm:ss
#SBATCH --job-name=2.1b_download_refseq_15-29 # SLURM job identifier
#SBATCH --output=../log/2.1b_download_refseq_15-29%A_%a.out
#SBATCH --error=../log/2.1b_download_refseq_15-29%A_%a.err

# SLURM email notifications are sent to the address specified below.
#SBATCH --mail-type=ALL # Enable SLURM email notifications
#SBATCH --mail-user=your_email@example.com

# Downloads RefSeq protein database volumes 15-29 into the configured database directory.
# Each tarball is downloaded, MD5-VERIFIED, extracted, then deleted.
# A ".done" marker per volume makes the script resumable: re-running it only
# re-fetches volumes that are missing or failed verification.

# Database availability checked 2026-05-31; RefSeq protein volumes 00-59 were available.

# This script is designed for execution on a SLURM-managed HPC environment.
# Database files are stored in the user's allocated scratch/nobackup directory.

mkdir -p /nobackup/$USER/pcu-nir/db

cd /nobackup/$USER/pcu-nir/db || {
  echo "Failed to enter database directory."
  exit 1
}

BASE="https://ftp.ncbi.nlm.nih.gov/blast/db"

# RefSeq protein database (links) volumes downloaded by this script
URLS=(
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.15.tar.gz"
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.16.tar.gz"
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.17.tar.gz"
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.18.tar.gz"
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.19.tar.gz"
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.20.tar.gz"
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.21.tar.gz"
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.22.tar.gz"
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.23.tar.gz"
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.24.tar.gz"
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.25.tar.gz"
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.26.tar.gz"
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.27.tar.gz"
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.28.tar.gz"
  "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.29.tar.gz"
)

ok=0; fail=0
for url in "${URLS[@]}"; do
  tb=$(basename "$url")                      # refseq_protein.NN.tar.gz
  vol=$(echo "$tb" | sed -E 's/refseq_protein\.([0-9]+)\.tar\.gz/\1/')
  flag="refseq_protein.${vol}.done"

  if [ -f "$flag" ]; then
    echo "[$(date +%H:%M:%S)] volume ${vol}: already done, skipping"
    ok=$((ok+1)); continue
  fi

  success=0
  for attempt in 1 2 3; do
    rm -f "$tb" "$tb.md5"
    echo "[$(date +%H:%M:%S)] volume ${vol}: download attempt ${attempt}"
    wget -c -q --tries=10 --waitretry=10 --retry-connrefused "$url"
    wget -q -O "$tb.md5" "${url}.md5"
    if md5sum -c "$tb.md5" >/dev/null 2>&1; then
      tar -xzf "$tb" && rm -f "$tb" "$tb.md5" && touch "$flag"
      echo "[$(date +%H:%M:%S)] volume ${vol}: OK"
      success=1; break
    fi
    echo "[$(date +%H:%M:%S)] volume ${vol}: md5 failed on attempt ${attempt}"
  done

  if [ "$success" = "1" ]; then ok=$((ok+1)); else fail=$((fail+1)); echo "volume ${vol}: GAVE UP after 3 attempts"; fi
done

echo "--------------------------------------------------"
echo "This script: ${ok} volume(s) OK, ${fail} failed."
echo "Total complete volumes on disk: $(ls refseq_protein.*.done 2>/dev/null | wc -l) / 60"
echo "--------------------------------------------------"