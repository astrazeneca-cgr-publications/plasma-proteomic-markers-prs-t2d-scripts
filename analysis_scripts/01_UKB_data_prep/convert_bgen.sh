#!/usr/bin/bash
#script for preparing genotype data for PGS estimation
chr=$1

####EXTRACT UKB-PPP subjects ########
##################parameters####################
plink=../software/plink
plink2=../software/plink2
compute="--memory 7000 --threads 8"
#add R environment to path
export PATH="../envs/R/bin:$PATH"
###################ukb imputed bgen files ##########################
bgen=../datasets/ukb/Imputed/ukb_c${chr}_b0_v3.bgen
sample=../datasets/ukb/Imputed/ukb_c${chr}_b0_v3.sample
#############################################

#get ukb-ppp samples
linker=../datasets/olink_sample_map_3k_ukbsamples_forconsort_v2_linked.tsv
cut -f23 $linker | sort -u | awk '{print $1, $1}' > ukbppp_samples

##covert to pgen, keeping only specified variants##
#snps consist of hm3 variants, GWAS-signif variants, and pqtls#
snps=./data/sumstats/hm3_plus_gwas_sgnif_plus_pqtls_snp_list
out=./data/genotypes/UKBPPP.chr$chr
$plink2 --bgen $bgen 'ref-first' --sample $sample \
--keep ukbppp_samples \
--extract $snps \
--silent \
--make-pgen --out $out $compute

#filter, removing CG/AT sites and duplicates
Rscript ./scripts/01_data_prep/filter_snps.R $out.pvar 

#extract filtered snp list
$plink2 --pfile $out --extract $out.info_filter_snp_list --make-pgen --out $out.info_filtered $compute --silent

#extract pqtls since they will be used separately 
pqtls=./data/sumstats/pqtl_snp_list
$plink2 --pfile $out --extract $pqtls --make-pgen --out $out.pqtls $compute --silent

#remove unfiltered file
rm $out.pvar $out.psam $out.pgen 
