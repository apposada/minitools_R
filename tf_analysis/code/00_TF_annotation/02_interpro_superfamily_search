#!/bin/bash
# $1 is the protein sequence fasta file you want to scan
# $2 is the path to the interproscan script and db, ie wherever you have installed it 

# setup
set -e
x=$1
interproscan_path=$2

# clear the params
set --

interproscan.sh -i $x -d 20220705_${x}_interproscan_SFAM --cp
u 12 -goterms --appl SUPERFAMILY