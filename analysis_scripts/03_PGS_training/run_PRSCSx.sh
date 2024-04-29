#!/usr/bin/bash
# script runs PRS-CSx, a multi-ancestry version of PRS-CS 
############################################################################
#Multi-threading
N_THREADS=8
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

########################PARAMETERS##########################################
#from command line
trait=$1
CHROM=$2

#algorithm inputs
bim=./data/ukb.hm3.postQC

#out directory
out_dir=./data/nonUKB_sumstats/post_eff/
prefix=$trait.V1

#add environment to path
export PATH=".../envs/prscs/bin:$PATH"
#path to PRScsx
PRScsx=.../software/PRScsx/PRScsx.py
###########################################################################
EAS=./data/sumstats/parsed/$trait.EAS.prscs.parsed.txt
EUR=./data/sumstats/parsed/$trait.EUR.prscs.parsed.txt
AFR=./data/sumstats/parsed/$trait.AFR.prscs.parsed.txt
SAS=./data/sumstats/parsed/$trait.SAS.prscs.parsed.txt
AMR=./data/sumstats/parsed/$trait.AMR.prscs.parsed.txt

#sample sizes
#Note: had some issues with floating points since I provided an average sample size. Rounded to fix.
Neas=$(cat ./data/sumstats/parsed/$trait.EAS.N.txt | awk '{print int($1+0.5)}')
Neur=$(cat ./data/sumstats/parsed/$trait.EUR.N.txt | awk '{print int($1+0.5)}')
Nsas=$(cat ./data/sumstats/parsed/$trait.SAS.N.txt | awk '{print int($1+0.5)}')
Nafr=$(cat ./data/sumstats/parsed/$trait.AFR.N.txt | awk '{print int($1+0.5)}')
Namr=$(cat ./data/sumstats/parsed/$trait.AMR.N.txt | awk '{print int($1+0.5)}')

#out directory
out_dir=./data/sumstats/post_eff/
prefix=$trait

###########################################################################
python $PRScsx \
--ref_dir=./ld_ref \
--bim_prefix=$bim \
--sst_file=$EAS,$EUR,$SAS,$AFR,$AMR \
--n_gwas=$Neas,$Neur,$Nsas,$Nafr,$Namr \
--pop=EAS,EUR \
--out_dir=$out_dir \
--out_name=$prefix \
--chrom=$CHROM \
--meta=True \
--seed=44 \
--a=1 \
--b=0.5