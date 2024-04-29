#!/usr/bin/bash
####script runs PRS-CS to estimate posterior effect sizes for PRS 

############################################################################
#Multi-threading
N_THREADS=8
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

#add environment to path
export PATH="../envs/prscs/bin:$PATH"
#path to PRScs
PRScs=../software/PRScs/PRScs.py
########################PARAMETERS##########################################
#from command line
prefix=$1
CHROM=$2

#algorithm inputs
bim=./inputs/$prefix #dont include .bim 
#sumstats, N
N=$(cat ./inputs/$prefix.N.txt | awk '{print int($1+0.5)}') #prs-cs breaks if N is a floating point
gwas=./inputs/$prefix.prscs.parsed.txt
#LD reference, assuming European ancestry. Change _eur to relevant 1KGP superpop
ref=./ld_ref/ldblk_1kg_eur
#out directory

#specify output 
out_prefix=${prefix}_prscs
mkdir -p post_eff 
out_dir=./post_eff/${out_prefix}

###########################################################################
#note, uses "auto" option for determing the phi parameter. change with --phi [param] 
python $PRScs \
--ref_dir=$ref \
--bim_prefix=$bim \
--sst_file=$gwas \
--n_gwas=$N \
--out_dir=$out_dir \
--chrom=$CHROM \
--seed=44 \
--a=1 \
--b=0.5