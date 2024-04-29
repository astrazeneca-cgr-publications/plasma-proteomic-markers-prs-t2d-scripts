#script for running a GWAS for coloc
#get inpputs

trait=$1
pop=$2
chr=$3
snp=$4

####### snp list for variants +/250 kb of index variant #########
snps=./data/coloc_data/regenie_inputs/${trait}_${pop}_${chr}_${snp}_snplist


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
if [ $pop == "T2D" ]; then
    pheno=./data/mr_data/regenie_inputs/UKB_T2Dsubset_phenotypes_for_regenie.txt
else
    pheno=./data/mr_data/regenie_inputs/UKB_phenotypes_for_regenie.txt
fi

#path to covariate file
covar=./data/mr_data/regenie_inputs/UKB_covariates_for_regenie.txt

#out directory
out=./data/coloc_data/regenie_outputs/

####### run step 2 ########
#prefix
prefix=UKB_${pop}

if [ $trait == "BMI" ]; then
    #path to predictions
    pred=./data/mr_data/regenie_outputs/UKB_${pop}_BMIonly_step1_QT_pred.list

    $regenie \
    --step 2 \
    --bgen $bgen \
    --sample $sample \
    --ref-first \
    --extract $snps \
    --phenoFile $pheno \
    --pred $pred \
    --phenoCol "BMI" \
    --covarFile $covar \
    --covarColList age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --catCovarList sex \
    --bsize 200 \
    --out ${out}/${prefix}_BMIonly_step2_QT_${trait}_${snp}

    #compress output
    gzip -f ${out}/${prefix}_BMIonly_step2_QT_${trait}_${snp}_${trait}.regenie 

elif [ $trait == "T2D" ]; then
    #run both BMI unadjusted and adjusted models
    pred=./data/mr_data/regenie_outputs/UKB_${pop}_step1_BT_pred.list
    
    $regenie \
    --step 2 \
    --af-cc \
    --bgen $bgen \
    --sample $sample \
    --ref-first \
    --extract $snps \
    --pred $pred \
    --phenoFile $pheno \
    --phenoCol "T2D" \
    --covarFile $covar \
    --covarColList age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --catCovarList sex \
    --bt \
    --bsize 200 \
    --out ${out}/${prefix}_step2_BT_${trait}_${snp}

    #compress output 
    gzip -f ${out}/${prefix}_step2_BT_${trait}_${snp}_$trait.regenie
    
  #path to predictions
  pred=./data/mr_data/regenie_outputs/UKB_${pop}_BMIadj_step1_BT_pred.list
    #T2D, BMI adjusted
    $regenie \
    --step 2 \
    --af-cc \
    --bgen $bgen \
    --ref-first \
    --sample $sample \
    --extract $snps \
    --pred $pred \
    --phenoFile $pheno \
    --phenoCol "T2D" \
    --covarFile $covar \
    --covarColList age,sex,BMI,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --catCovarList sex \
    --bt \
    --firth --approx \
    --bsize 200 \
    --out ${out}/${prefix}_BMIadj_step2_BT_${trait}_${snp}

    #compress output
    gzip -f ${out}/${prefix}_BMIadj_step2_BT_${trait}_${snp}_$trait.regenie 
else
  pred=./data/mr_data/regenie_outputs/UKB_${pop}_step1_BT_pred.list

  $regenie \
    --step 2 \
    --af-cc \
    --bgen $bgen \
    --sample $sample \
    --ref-first \
    --extract $snps \
    --pred $pred \
    --phenoFile $pheno \
    --phenoCol $trait \
    --covarFile $covar \
    --covarColList age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --catCovarList sex \
    --bt \
    --bsize 200 \
    --out ${out}/${prefix}_step2_BT_${trait}_${snp}
    #compress output
    gzip -f ${out}/${prefix}_step2_BT_${trait}_${snp}_$trait.regenie
fi