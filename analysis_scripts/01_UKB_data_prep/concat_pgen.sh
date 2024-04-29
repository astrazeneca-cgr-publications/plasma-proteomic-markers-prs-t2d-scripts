#script for concatenating output of convert_bgen.sh
#plink parameters#
plink2=../software/plink2
compute="--memory 30000 --threads 8"

#generate file lists #
echo -n >  concat.files
echo -n > pqtl.files
for chr in {1..22}; do
echo ./data/genotypes/UKBPPP.chr$chr.pqtls >> pqtl.files
echo ./data/genotypes/UKBPPP.chr$chr.info_filtered >> concat.files
done 



#concatenate pqtl file set#
out=./data/genotypes/UKBPPP
$plink2 --pmerge-list pqtl.files --make-pgen --out $out.pqtls --silent $compute

#concatenate hm3+other variants file set
$plink2 --pmerge-list concat.files --make-pgen --out $out.hm3_plus.info_filtered --silent $compute 


#remove per-chromosome files
rm ./data/genotypes/UKBPPP.chr*.pqtls.*
rm ./data/genotypes/UKBPPP.chr*.info_filtered.*
rm concat.files pqtl.files
rm ./data/genotypes/UKBPPP.hm3_plus.info_filtered-merge.* 
rm ./data/genotypes/UKBPPP.pqtls-merge.*
rm ./data/genotypes/UKBPPP.chr*.*


