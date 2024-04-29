#script for clumping snps in order to idenify instruments for reverse MR
chr=$1
##################parameters####################
plink=../software/plink
plink2=../software/plink2
compute="--memory 7000 --threads 8"
###################ukb imputed bgen files ##########################
#bgen=../datasets/ukb/Imputed/ukb_imp_chr${chr}_v3.bgen
#sample=../datasets/ukb/Imputed/ukb_imp_chr${chr}_v3.sample
bgen=../datasets/ukb/Imputed/ukb_c${chr}_b0_v3.bgen
sample=../datasets/ukb/Imputed/ukb_c${chr}_b0_v3.sample
#############################################

#get list of snps by chromosome
snps=./data/reverse_mr_data/sumstats/GWAS_signif_snp_list_intersection
awk -v c=$chr '($1 == c) {print $2}' $snps > temp_chr${chr}_snps

#create temporary plink file set
disc=./data/ukbppp_disc_samples

#use UKB-PPP discovery samples as represantative of LD structure for clumping
out=UKB_for_LD.chr$chr 
$plink2 --bgen $bgen 'ref-first' --sample $sample \
--keep $disc \
--extract temp_chr${chr}_snps \
--silent \
--make-pgen 'fill-missing-from-dosage' \
--out $out $compute

#convert to binary plink1 file set
$plink2 --pfile $out --make-bed --out $out.p1 \
--rm-dup 'exclude-all' \
--silent $compute 


#clump for each trait
for trait in BMI T2D T2D_BMIadj CAD CKD NASH_AST_ALT; do

    sumstats=./data/reverse_mr_data/sumstats/${trait}_GWAS_signif_results.tsv.gz

    $plink --bfile $out.p1 --clump $sumstats \
    --clump-r2 0.001 --clump-kb 500 \
    --out ./data/reverse_mr_data/clumped/${trait}_chr$chr \
    --silent $compute  

done

#done, clean up
rm $out.* temp_chr${chr}_snps 
