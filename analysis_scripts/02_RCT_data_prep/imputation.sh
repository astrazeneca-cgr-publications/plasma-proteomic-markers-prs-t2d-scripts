#/bin/bash
### script for imputing genotype data in the clinical trials ###

################################################################################
chrname=$1;
trial=$2;
####software#####
plink2=../genotype_imputation/apps/plink2;
eagle=../genotype_imputation/apps/Eagle_v2.4.1/eagle;
beagle=../genotype_imputation/apps/beagle.27Jan18.7e1.jar;

##target data ###

DATASET_data=./data/$trial/for_imputation/$trial-CLEAN.RELS-updated-chr${chrname}
temp=./data/$trial/imputed/$trail-chr${chrname}_temp;
base=./data/$trial/imputed/$trial-chr${chrname};

##reference data
panel=..//datasets/1000G_ref_panel/1000G_highC_hg38/snpid/chr${chrname}.1000G.high_c.snpid.hg38.vcf.gz

#genetic maps
genetic_map_eagle=../genotype_imputation/apps/eagle_b38/eagle_chr${chrname}_b38.map
genetic_map_beagle=../genotype_imputation/apps/beagle_b38/beagle_chr${chrname}_b38.map;

#cpus
n_cpu=32;
################################################################################
START=$(date)

##### Prep. for imputation ############################################################
${plink2} --bfile ${DATASET_data} --recode vcf-iid bgz --output-chr chrM --out ${temp};
${plink2} --vcf ${temp}.vcf.gz --recode vcf-iid bgz  --vcf-half-call haploid --output-chr chrM --out ${base};
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' ${base}.vcf.gz -Ou | bcftools +fill-tags -Oz -o ${base}_qcd.vcf.gz -- -t AC,AN,AF;


##### phasing for imputation ############################################################
${eagle} \
--vcf ${base}_qcd.vcf.gz \
--chrom ${chrname} \
--geneticMapFile ${genetic_map_eagle} \
--outPrefix ${base}_for_imputation \
--Kpbwt 20000 \
--numThreads $n_cpu


##### imputation ############################################################
java -Xss5m -Xmx32g -jar ${beagle} \
gt=${base}_for_imputation.vcf.gz \
ref=${panel} \
map=${genetic_map_beagle} \
out=${base}_imputed \
niterations=10 \
ne=20000 \
impute=true \
gprobs=true \
window=300000 \
seed=-99999 \
nthreads=$n_cpu

END=$(date)
#echo "Phasing+imputation started at $START"
#echo "Phasing+imputation ended at $END"
echo "START,$START" > chr$chrname.imp_time.txt
echo "END,$END" >> chr$chrname.imp_time.txt

##### Post imputation ############################################################
START=$(date)
bcftools index -t -f ${base}_imputed.vcf.gz --threads $n_cpu;
bcftools +fill-tags ${base}_imputed.vcf.gz -Ou -- -t AF,AC_Hom,AC_Het,HWE | bcftools +impute-info -Oz -o ${base}_imputed_infotags.vcf.gz --threads $n_cpu;

echo -e 'CHR\tSNP\tREF\tALT\tAF\tINFO\tAF_GROUP' > ${base}_varID_AF_INFO_GROUP.txt;
bcftools query -f '%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\t%INFO/INFO\t-\n' ${base}_imputed_infotags.vcf.gz | awk '{if ($5>=0.05 && $5<=0.95) $7=1; else if(($5>=0.005 && $5<0.05) || ($5<=0.995 && $5>0.95)) $7=2; else $7=3} { print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7 }' >> ${base}_varID_AF_INFO_GROUP.txt;

END=$(date)
echo "POST-IMPUTATION" >> chr$chrname.imp_time.txt
echo "START,$START" >> chr$chrname.imp_time.txt
echo "END,$END" >> chr$chrname.imp_time.txt
