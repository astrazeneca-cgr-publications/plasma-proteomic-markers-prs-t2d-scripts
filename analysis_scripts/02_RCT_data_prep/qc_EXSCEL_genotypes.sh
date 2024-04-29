#script for estimating PCs and performing QC on EXSCEL genotype data
trial=EXSCEL

#software
plink2=../software/plink2
plink=../software/plink
king=../software/king
compute="--memory 10000 --threads 4"

#add R path
export PATH="../envs/R/bin:$PATH"

#path to genotype data
bfile=./data/EXSCEL/genotype-quality-control-EXSCEL_merge_update/data/EXSCEL_merge_update


#basic parameters

# Exclusion rules for SNPs:
snp_miss_h=0.02;    # threshold for SNP missing rate
snp_freq_h=0.005;    # threshold for minor allele frequency
snp_hwe_h=1e-15;        # threshold for HWE

# Exclusion rules for individuals:
ind_miss_h=0.05;    # threshold for individual missing rate
ind_het_sd=4;         # standard dev. from mean for outlying heterozygosity
ind_rel=2;      #degree for relationship inference
anc_prob=0.90;    #prob for calling ancestry with KING

#qc file
out=./data/rct_data/EXSCEL
${plink2} --bfile ${bfile} \
      --mind  ${ind_miss_h}\
      --autosome \
      --hwe ${snp_hwe_h} \
      --maf ${snp_freq_h} \
      --geno ${snp_miss_h} \
      --make-bed \
      --out ${out}-IBD --silent \
      --maj-ref \
      $compute;

#missingness and heterozygosity and sex check 
${plink} --bfile ${bfile} \
    --missing \
    --check-sex --het \
    --out ${out} --silent \
    $compute;

#generate ld pruned dataset for ancestry inference
${plink2} --bfile ${out}-IBD \
      --indep-pairwise 1500 150 0.1 \
      --maf 0.01 \
      --out ${out}-LD --silent \
      $compute; 
${plink2} --bfile ${out}-IBD \
      --extract ${out}-LD.prune.in \
      --out ${out}-pruned --silent \
      --make-bed \
      $compute; 

#ancestry estimation from genotype data
ref=./anc_ref/KGref.hg38
${king} -b ${ref}.bed,${out}-pruned.bed --pca --projection --rplot --cpus 8 --prefix ${out}_ancestry

#HWE (ancestry specific)
ancestry=${out}_ancestry_InferredAncestry.txt 
echo -n > ${out}_temp_excl 
for pop in AFR AMR EAS EUR SAS; do
awk -v p=$pop '($5 == p) && ($6 >= 0.9) {print $1, $2}' $ancestry > ${out}_keep 
$plink2 --bfile ${bfile} --keep ${out}_keep --hardy --out ${out}_hwe $compute --silent 
awk '($10 < 1E-15) {print $2}' ${out}_hwe.hardy >> ${out}_temp_excl 
done
rm ${out}_hwe.* 

#prepare final hwe file
sort -u ${out}_temp_excl > ${out}_hwe_excl_genotypes
rm ${out}_temp_excl

#relationship inference
${king} -b ${out}-IBD.bed --unrelated --degree $ind_rel --cpus 8 --prefix ${out}_;

#compile failed variant list
#compile failed individual list - including consent issues!
Rscript compile_qc.R $out $ind_het_sd

failed_snps=${out}_variant_excl
failed_subjects=${out}_ind_excl 

#redo LD pruning 
${plink2} --bfile ${bfile} \
      --mind  ${ind_miss_h}\
      --autosome \
      --exclude $failed_snps \
      --remove $failed_subjects \
      --maf 0.01 \
      --geno ${snp_miss_h} \
      --indep-pairwise 1500 150 0.1 \
      --out ${out}-LD --silent \
      $compute;
#PCA
${plink2} --bfile ${bfile} \
    --remove $failed_subjects \
    --extract ${out}-LD.prune.in \
    --pca \
    --out ${out}-PCA \
    $compute;

#imputation fileset
${plink2} --bfile ${bfile} \
      --mind  ${ind_miss_h}\
      --autosome \
      --exclude $failed_snps \
      --remove $failed_subjects \
      --maf 0.005 \
      --geno ${snp_miss_h} \
      --out ./data/$trial/for_imputation/$trial-CLEAN.RELS-updated-chr${chrname}

#clean up unecessary files
rm ${out}_allsegs.txt  ${out}_ancestry_ancestryplot.ps ${out}-IBD.* ${out}_keep ${out}-pruned.* ${out}.sexcheck ${out}_ancestry_ancestryplot.R ${out}_ancestry_ancestryplot.Rout ${out}-LD.* ${out}.hh ${out}.het ${out}_ancestry_pcplot.R  ${out}.lmiss ${out}.imiss ${out}.log    
#done