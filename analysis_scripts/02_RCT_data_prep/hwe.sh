#script for running remaining genotype QC steps before estimating genetic scores
#first,run hwe on a per-population basis in EXSCEL
#variants with p-value exact test < 1E-15 in any population will be excluded from score
#plink paramters #
cd ./data/rct_data
plink2=../software/plink2
compute="--memory 15000 --threads 4"
#####

####run in EXCEL#######
out=./data/rct_data
pfile=$out/EXSCEL.for_scores
for pop in AFR EAS EUR SAS AMR; do
    awk -v p=$pop '($5 == p) && ($6 >= 0.95) {print $1, $2}' ../trialstrials/data/ancestry/EXSCEL.hm3.1KGP_ancestry_InferredAncestry.txt > pop.txt

    $plink2 --pfile $pfile --keep pop.txt --out $out/$pop $compute --hardy --silent 

    #extract variants
    awk '($10 < 1E-15) {print $2}' $pop.hardy >> $out/temp_snp_list

    #compress
    gzip -f $out/$pop.hardy 
done

sort -u $out/temp_snp_list > $out/EXSCEL_hwe_fail_snp_list
rm $out/temp_snp_list pop.txt *.hardy.gz *.log

####run in DECLARE ###

pfile=$out/DECLARE.for_scores
for pop in AFR EAS EUR SAS AMR; do
    awk -v p=$pop '($5 == p) && ($6 >= 0.95) {print $1, $2}' ../trialstrials/data/ancestry/DECLARE.hm3.1KGP_ancestry_InferredAncestry.txt > pop.txt

    $plink2 --pfile $pfile --keep pop.txt --out $out/$pop $compute --hardy --silent 
    #extract variants
    awk '($10 < 1E-15) {print $2}' $pop.hardy >> $out/temp_snp_list

    #compress
    gzip -f $out/$pop.hardy 
done

sort -u $out/temp_snp_list > $out/DECLARE_hwe_fail_snp_list
rm $out/temp_snp_list pop.txt *.hardy.gz *.log
