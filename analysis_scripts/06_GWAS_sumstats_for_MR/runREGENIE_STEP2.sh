#script for running a GWAS
#get inpputs
pop=$1
chr=$2

###################ukb imputed bgen files ##########################
#bgen=../datasets/ukb/Imputed/ukb_imp_chr${chr}_v3.bgen
#sample=../datasets/ukb/Imputed/ukb_imp_chr${chr}_v3.sample
bgen=../datasets/ukb/Imputed/ukb_c${chr}_b0_v3.bgen
sample=../datasets/ukb/Imputed/ukb_c${chr}_b0_v3.sample
bgi=../datasets/ukb/Imputed/ukb_c${chr}_b0_v3.bgen.bgi
#############################################

#path to regenie
regenie=../software/regenie/regenie

#path to phenotype file

pheno=./data/mr_data/regenie_inputs/UKB_phenotypes_for_regenie.txt

#path to covariate file
covar=./data/mr_data/regenie_inputs/UKB_covariates_for_regenie.txt

#out directory
out=./data/coloc_data/regenie_outputs/

####### run step 2 ########
#prefix
prefix=UKB_${pop}

###BMI###
#path to predictions
pred=./data/mr_data/regenie_outputs/UKB_${pop}_BMIonly_step1_QT_pred.list

$regenie \
--step 2 \
--bgen $bgen \
--sample $sample \
--ref-first \
--minMAC 50 \
--phenoFile $pheno \
--pred $pred \
--phenoCol "BMI" \
--covarFile $covar \
--covarColList age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--catCovarList sex \
--bsize 200 \
--out ${out}/${prefix}_BMIonly_step2_QT_chr${chr}

#compress output
gzip -f ${out}/${prefix}_BMIonly_step2_QT_chr${chr}_BMI.regenie 

#### T2D, BMI adjusted ####
    
#path to predictions
pred=./data/mr_data/regenie_outputs/UKB_${pop}_BMIadj_step1_BT_pred.list
#T2D, BMI adjusted
$regenie \
--step 2 \
--af-cc \
--bgen $bgen \
--ref-first \
--minMAC 50 \
--sample $sample \
--pred $pred \
--phenoFile $pheno \
--phenoCol "T2D" \
--covarFile $covar \
--covarColList age,sex,BMI,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--catCovarList sex \
--bt \
--firth --approx \
--bsize 200 \
--out ${out}/${prefix}_BMIadj_step2_BT_chr${chr}

#compress output
gzip -f ${out}/${prefix}_BMIadj_step2_BT_chr${chr}_T2D.regenie 


#### All other traits #####
pred=./data/mr_data/regenie_outputs/UKB_${pop}_step1_BT_pred.list

$regenie \
--step 2 \
--af-cc \
--bgen $bgen \
--sample $sample \
--ref-first \
--minMAC 50 \
--pred $pred \
--phenoFile $pheno \
--phenoColList "T2D,CAD,CKD,NASH_AST_ALT" \
--covarFile $covar \
--covarColList age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--catCovarList sex \
--bt \
--bsize 200 \
--out ${out}/${prefix}_step2_BT_chr${chr}

#compress output
for trait in T2D CAD CKD NASH_AST_ALT; do 
gzip -f ${out}/${prefix}_step2_BT_chr${chr}_$trait.regenie
done 