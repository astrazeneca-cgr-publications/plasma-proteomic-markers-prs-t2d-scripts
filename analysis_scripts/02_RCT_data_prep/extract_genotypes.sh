##script for extracting genotypes from EXSCEL ##

plink=../software/plink
plink2=../software/plink2
compute="--memory 16000 --threads 4"

#EXSCEL first
echo > n concat.files
out=./data/rct_data
cut -f1 $out/EXSCEL_matched_variants.txt > $out/snp_list

for chr in {1..22}; do
    pfile=../trials/data/EXSCEL_IMP/EXSCEL-chr${chr}_high_impQ_common
    $plink2 --pfile $pfile --extract $out/snp_list --make-pgen --out $out/EXSCEL.chr$chr $compute --silent
    echo $out/EXSCEL.chr$chr >> concat.files
done 

#concatentate
$plink2 --pmerge-list concat.files --make-pgen --out $out/EXSCEL.temp $compute --silent 

#rename variants
$plink2 --pfile $out/EXSCEL.temp --update-name $out/EXSCEL_matched_variants.txt --make-pgen --out EXSCEL.for_scores

#clean up
rm  $out/snp_list concat.files $out/EXSCEL.chr* $out/EXSCEL.temp*


#DECLARE
echo > n concat.files
out=./data/rct_data
cut -f1 $out/DECLARE_matched_variants.txt > $out/snp_list

for chr in {1..22}; do
    pfile=../trials/data/DECLARE_IMP/DECLARE--chr${chr}_high_impQ_common
    $plink2 --pfile $pfile --extract $out/snp_list --make-pgen --out $out/DECLARE.chr$chr $compute --silent
    echo $out/DECLARE.chr$chr >> concat.files
done 

#concatentate
$plink2 --pmerge-list concat.files --make-pgen --out $out/DECLARE.temp $compute --silent 

#rename variants
$plink2 --pfile $out/DECLARE.temp --update-name $out/DECLARE_matched_variants.txt --make-pgen --out DECLARE.for_scores $compute 

#clean up
rm  $out/snp_list concat.files $out/DECLARE.chr* $out/DECLARE.temp*