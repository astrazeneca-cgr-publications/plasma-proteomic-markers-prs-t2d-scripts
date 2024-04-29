#script for running step 1 of regenie
#specify inferred ancestry OR T2D status: EAS, AFR, AMR, EUR, SAS, T2D 
pop=$1

#path to regenie
regenie=../software/regenie/regenie

#path to genotype file
geno=./data/mr_data/regenie_inputs/UKB_${pop}

#path to phenotype file
if [ $pop == "T2D" ]; then
    pheno=./data/mr_data/regenie_inputs/UKB_T2Dsubset_phenotypes_for_regenie.txt
else
    pheno=./data/mr_data/regenie_inputs/UKB_phenotypes_for_regenie.txt
fi

#path to covariate file
covar=./data/mr_data/regenie_inputs/UKB_covariates_for_regenie.txt

#out directory
out=./data/mr_data/regenie_outputs/

####### run step 1 ########
#prefix
prefix=UKB_${pop}

#traits different in T2D vs all others
if [ $pop == "T2D" ]; then

  $regenie \
    --step 1 \
    --pgen $geno \
    --phenoFile $pheno \
    --phenoColList "insulin,DR,chronic_renal,acute_renal,CAD,HF,MI,hypertension,stroke,hyperlipidemia,DN,CKD,NASH_AST_ALT,NASH_ICD10,obese" \
    --covarFile $covar \
    --covarColList age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --catCovarList sex \
    --bt \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix ${out}/tmpdir/${prefix}_tmp_preds \
    --out ${out}/${prefix}_step1_BT
else
    $regenie \
    --step 1 \
    --pgen $geno \
    --phenoFile $pheno \
    --phenoColList "T2D,CAD,CKD,NASH_AST_ALT" \
    --covarFile $covar \
    --covarColList age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --catCovarList sex \
    --bt \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix ${out}/tmpdir/${prefix}_tmp_preds \
    --out ${out}/${prefix}_step1_BT

    #T2D, BMI adjusted
    $regenie \
    --step 1 \
    --pgen $geno \
    --phenoFile $pheno \
    --phenoCol "T2D" \
    --covarFile $covar \
    --covarColList age,sex,BMI,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --catCovarList sex \
    --bt \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix ${out}/tmpdir/${prefix}_tmp_preds \
    --out ${out}/${prefix}_BMIadj_step1_BT
fi

#BMI is run in both T2D and pop-specific runs

$regenie \
  --step 1 \
  --pgen $geno \
  --phenoFile $pheno \
  --phenoCol "BMI" \
  --covarFile $covar \
  --covarColList age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --catCovarList sex \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix ${out}/tmpdir/${prefix}_tmp_preds \
  --out ${out}/${prefix}_BMIonly_step1_QT
