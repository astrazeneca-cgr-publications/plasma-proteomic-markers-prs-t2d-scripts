###script for preparping inputs for REGENIE###
#plink parameters
plink2=../software/plink2
compute="--memory 30000 --threads 8"

#genotype file
geno=./ukb/PCA/UKB.genotyped.plink2

#covariate file with ancestry 
covar=./data/phenotypes/UKB.covariates.tsv.gz

#linker file with UKB-PPP information to exclude
linker=../datasets/ukb/olink_sample_map_3k_ukbsamples_forconsort_v2_linked.tsv
cut -f23 $linker | sort -u | awk '{print $1, $1}' > ukbppp_samples

#out directory
out=./data/mr_data/regenie_inputs/

for pop in EAS EUR AFR AMR SAS; do 
echo $pop 
zgrep $pop $covar | awk '($6 >= 0.90) {print $1, $1}' > pop_samples

#LD prune
$plink2 --pfile $geno \
    --remove ukbppp_samples \
    --keep pop_samples \
    --maf 0.01 --mac 100 \
    --geno 0.01 --mind 0.1 \
    --hwe 1E-15 \
    --indep-pairwise 1000 100 0.8 \
    --out ${out}UKB_${pop}_LD \
    $compute --silent

#make file set 
$plink2 --pfile $geno \
    --remove ukbppp_samples \
    --keep pop_samples \
    --extract ${out}UKB_${pop}_LD.prune.in \
    --make-pgen \
    --out ${out}UKB_${pop} \
    $compute --silent
done

#now for T2D cases - European ancestry 
#exclude based on ancestry + UKB-PPP status, include based on T2D status
pop=EUR
#easier to exlude in this case
zgrep $pop $covar | awk '($6 < 0.90) {print $1, $1}' > pop_exclude 
zgrep -v $pop $covar | awk ' {print $1, $1}' >> pop_exclude 
cat ukbppp_samples >> pop_exclude
#T2D status
pheno=./data/phenotypes/UKB.phenotypes.tsv.gz

zcat $pheno | awk ' ($2 == 1) {print $1, $1}' > t2d_include

#LD prune
$plink2 --pfile $geno \
    --remove pop_exclude \
    --keep t2d_include \
    --maf 0.01 --mac 100 \
    --geno 0.01 --mind 0.1 \
    --hwe 1E-15 \
    --indep-pairwise 1000 100 0.8 \
    --out ${out}UKB_T2D_LD \
    $compute --silent

#make file set 
$plink2 --pfile $geno \
    --remove pop_exclude \
    --keep t2d_include \
    --extract ${out}UKB_T2D_LD.prune.in \
    --make-pgen \
    --out ${out}UKB_T2D \
    $compute --silent

#done 
