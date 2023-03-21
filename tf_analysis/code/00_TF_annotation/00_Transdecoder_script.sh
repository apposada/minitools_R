#!/bin/bash

# setup
set -e
x=$1

# clear the params
set --

# activation
source /mnt/sda/alberto/programs/miniconda3/bin/activate
conda activate transdecoder_venv

# running the stuff
TransDecoder.LongOrfs -t $x
blastp -query ${x}.transdecoder_dir/longest_orfs.pep -db ~/DATA/static/databases/uniprot_swissprot/uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > blastp.outfmt6
hmmscan --cpu 12  --domtblout pfam.domtblout ~/DATA/static/databases/pfam/Pfam-A.hmm ${x}.transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t $x --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only
conda deactivate

cp ${x}.transdecoder.pep ${x%.fa}.pep

echo "DONE"
