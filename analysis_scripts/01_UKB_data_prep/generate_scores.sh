###script for generating genetic scores for UKB-PPP participants#####
cd ./data/pgs

#plink2 parameters
plink2=../software/plink2
compute="--memory 30000 --threads 8"

#genotype file
pfile=../genotypes/UKBPPP.hm3_plus.info_filtered
hwe_file=../genotypes/hwe_fail_snp_list 


###scores from PRS-CS####
scores=$(ls ../sumstats/post_eff/)

for score in $scores; do
out="UKBPPP_$(echo $score | sed 's/_pst_eff_a1_b0.5_phiauto.txt.gz/_prscs/g' | sed 's/_prscs_prscs/_prscs/g')"

zcat ../sumstats/post_eff/$score | cut -f2,4,6 > temp_score 
$plink2 --pfile $pfile --exclude $hwe_file --score temp_score --out $out $compute --silent 

done

###genome-wide scores ####
#NASH
score=../sumstats/scores/anstee_nash_grs.txt
out=UKBPPP_NASH_GRS
$plink2 --pfile $pfile --exclude $hwe_file --score $score 1 3 4 --out $out $compute --silent

#NASH-biomarker defined
score=../sumstats/scores/vujkovic_alt_grs.txt
out=UKBPPP_NASH_biomarker_GRS
$plink2 --pfile $pfile --exclude $hwe_file --score $score 1 4 6 --out $out $compute --silent

#BMI
score=../sumstats/scores/BMI.gwas_signif.PGS000841.txt.gz
zgrep -v "#" $score > temp_score 
out=UKBPPP_BMI_GWAS_signif
$plink2 --pfile $pfile --exclude $hwe_file --score temp_score --out $out $compute --silent

#CAD
score=../sumstats/scores/CAD.gwas_signif.PGS000019.txt.gz
zgrep -v "#" $score |  cut -f1,4,5 > temp_score 
out=UKBPPP_CAD_GWAS_signif
$plink2 --pfile $pfile --exclude $hwe_file --score temp_score --out $out $compute --silent

#CKD
score=../sumstats/scores/ckd.gwas_signif.PGS000859.txt.gz
zgrep -v "#" $score > temp_score 
out=UKBPPP_CKD_GWAS_signif
$plink2 --pfile $pfile --exclude $hwe_file --score temp_score --out $out $compute --silent

#diabetic retinopathy
#note there was a duplicated line in score file, just commented out
score=../sumstats/scores/diab_retin.gwas_signif.PGS000862.txt.gz
zgrep -v "#" $score > temp_score 
out=UKBPPP_diab_retin_GWAS_signif
$plink2 --pfile $pfile --exclude $hwe_file --score temp_score --out $out $compute --silent

#T2D
score=../sumstats/scores/t2d.gwas_signif.txt.gz
zgrep -v "rsID" $score | cut -f1,4,6 > temp_score 
out=UKBPPP_T2D_GWAS_signif
$plink2 --pfile $pfile --exclude $hwe_file --score temp_score --out $out $compute --silent


####partioned scores###
scores=$(ls ../sumstats/scores | grep 2018)
for score in $scores; do
out="UKBPPP_$(echo $score | sed 's/.2018.txt//g')"
$plink2 --pfile $pfile --exclude $hwe_file --score ../sumstats/scores/$score --out $out $compute --silent

done 
