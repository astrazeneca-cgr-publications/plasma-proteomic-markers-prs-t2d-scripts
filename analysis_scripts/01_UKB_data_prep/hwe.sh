#script for running remaining QC steps
#first,run hwe on a per-population basis in UKB 
#variants with p-value exact test < 1E-15 in any population will be excluded from score
#plink paramters #
plink2=../software/plink2
compute="--memory 30000 --threads 8"
##
#file specifiying ancestry#
out=./data/genotypes
pfile=$out/UKBPPP.hm3_plus.info_filtered
for pop in AFR EAS EUR SAS AMR; do
zgrep $pop ./data/phenotypes/UKB.covariates.tsv.gz | awk '$6 >= 0.95 {print $1, $1}' > pop.txt

$plink2 --pfile $pfile --keep pop.txt --out $out/$pop $compute --hardy --silent 

#extract variants
awk '($10 < 1E-15) {print $2}' $pop.hardy >> $out/temp_snp_list

#compress
gzip -f $out/$pop.hardy 
done

sort -u $out/temp_snp_list > $out/hwe_fail_snp_list
rm $out/temp_snp_list

