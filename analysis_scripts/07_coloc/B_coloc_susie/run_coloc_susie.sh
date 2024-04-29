# script for running coloc_susie using the coloc package #
i=$1

####parse index file #########
index=./data/coloc_data/regenie_inputs/index_cis_snplist_for_susie.txt
chr=$(awk -v i=$i 'FNR == i {print $3}' $index)
pop=$(awk -v i=$i 'FNR == i {print $2}' $index)
trait=$(awk -v i=$i 'FNR == i {print $1}' $index)
snp=$(awk -v i=$i 'FNR == i {print $4}' $index)
#####specify snps and effect alleles from gwas file #####
if [ $trait == "T2D" ]; then
    trait=T2D_BMIadj
    gwas_file=./data/coloc_data/regenie_outputs/UKB_${pop}_BMIadj_step2_BT_${trait}_${snp}_${trait}.regenie.gz
elif [ $trait == "BMI" ]; then
    gwas_file=./data/coloc_data/regenie_outputs/UKB_${pop}_BMIonly_step2_QT_${trait}_${snp}_${trait}.regenie.gz 
else
    gwas_file=./data/coloc_data/regenie_outputs/UKB_${pop}_step2_BT_${trait}_${snp}_${trait}.regenie.gz
fi

#prepare temporary file#
if [ -f $gwas_file ]; then
    zcat $gwas_file | awk -v OFS='\t' 'NR >1 {print $3,$5}' > ./data/ld/inputs/chr${chr}_${pop}_${trait}_${snp}.window
    snps=./data/ld/inputs/chr${chr}_${pop}_${trait}_${snp}.window
else
    exit 1
fi

##################plink parameters####################
plink=../software/plink
plink2=../software/plink2
compute="--memory 7000 --threads 8"
###################ukb imputed bgen files ##########################
#bgen=../datasets/ukb/Imputed/ukb_imp_chr${chr}_v3.bgen
#sample=../datasets/ukb/Imputed/ukb_imp_chr${chr}_v3.sample
bgen=../datasets/ukb/Imputed/ukb_c${chr}_b0_v3.bgen
sample=../datasets/ukb/Imputed/ukb_c${chr}_b0_v3.sample

###########prepared pgen file for ld#############
pgen=./data/ld/inputs/ukb_chr${chr}
#############################################
ukbppp_samples=./data/ukbppp_samples
disc_samples=./data/ukbppp_disc_samples

####extract snps and ukb-ppp cases from bgen file and make bed file, only needed once per chrom #####
if [ ! -f $pgen.pgen ]; then
    out=./data/ld/inputs/ukb_chr${chr}
    chr_snps=./data/ld/inputs/chr${chr}_snps_for_ld
    $plink2 --bgen $bgen 'ref-first' --sample $sample \
    --keep $ukbppp_samples \
    --extract $chr_snps \
    --make-pgen 'fill-missing-from-dosage' \
    --out $out $compute \
    --silent 
fi 


####extract locus and convert to bed file ######
bfile=./data/ld/inputs/ukb_chr${chr}_${pop}_${trait}_${snp}
$plink2 --pfile $pgen \
--keep $ukbppp_samples \
--rm-dup 'exclude-all' \
--extract $snps \
--make-bed \
--out $bfile $compute \
--silent


###calculate ld matrix for all ukb-ppp subjects ####
#subset to mac of 50 as variants rarer than that do not have summary stats
out=./data/ld/outputs/ukb_chr${chr}_${pop}_${trait}_${snp}_combined
$plink --bfile $bfile \
--a1-allele $snps 2 1 \
--mac 50 \
--r square gz \
--make-just-bim \
--out $out $compute \
--silent
###calculate ld matrix for just disc subjects #####
out=./data/ld/outputs/ukb_chr${chr}_${pop}_${trait}_${snp}_disc
$plink --bfile $bfile \
--keep $disc_samples \
--a1-allele $snps 2 1 \
--mac 50 \
--r 'square' 'gz' \
--make-just-bim \
--out $out $compute \
--silent

#clean up by removing binary fileset and temporary file with snps
rm $bfile.* $snps 

##now run coloc ####
export PATH="../envs/MR/bin:$PATH"

Rscript ./scripts/06_coloc_MR/run_coloc_susie.R $i

###done ####

