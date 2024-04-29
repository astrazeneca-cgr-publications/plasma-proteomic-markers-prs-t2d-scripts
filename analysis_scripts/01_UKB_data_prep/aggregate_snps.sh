##script for aggregating snps used to estimate PGS##
cd ./data/sumstats
echo -n > temp_snp_list

#partioned scores
cut -f1 ./data/sumstats/scores/beta_cell.2018.txt >> temp_snp_list
cut -f1 ./data/sumstats/scores/liver_lipid.2018.txt >> temp_snp_list
cut -f1 ./data/sumstats/scores/lipodystrophy.2018.txt >> temp_snp_list
cut -f1 ./data/sumstats/scores/obesity.2018.txt >> temp_snp_list
cut -f1 ./data/sumstats/scores/proinsulin.2018.txt >> temp_snp_list

#GWAS-significant scores
#BMI: PGS000841
zgrep -v "#" ./data/sumstats/scores/BMI.gwas_signif.PGS000841.txt.gz | cut -f1 >> temp_snp_list

#CAD: PGS000019
zgrep -v "#" ./data/sumstats/scores/CAD.gwas_signif.PGS000019.txt.gz  | cut -f1 >> temp_snp_list

#CKD: PGS000859
zgrep -v "#" ./data/sumstats/scores/ckd.gwas_signif.PGS000859.txt.gz | cut -f1 >> temp_snp_list

#T2D: PGS00036
#note T2D score is genome-wide, will need to downsample later
zgrep -v "#" ./data/sumstats/scores/t2d.P_T.PGS000036.txt.gz| cut -f1 >> temp_snp_list

#liver scores
cut -f1 ./data/sumstats/scores/anstee_nash_grs.txt  >> temp_snp_list
cut -f1 ./data/sumstats/scores/vujkovic_alt_grs.txt >> temp_snp_list
#diabetic neuropathy PGS000862
zgrep -v "#" ./data/sumstats/scores/diab_retin.gwas_signif.PGS000862.txt.gz | cut -f1 >> temp_snp_list

#the hm3 variants
cut -f2 ./ld_ref/snpinfo_mult_1kg_hm3 >> temp_snp_list 

#the indepdendent pqtls
zcat ./data/pqtl_data/cond_indep_pqtl_list.csv.gz | awk -F, '{print $3}' >> temp_snp_list
#also keep separate pqtl step as this will not be filtered 
zcat ./data/pqtl_data/cond_indep_pqtl_list.csv.gz | awk -F, 'NR>1 {print $3}' | sort -u > pqtl_snp_list

#remove column headers and duplicates
grep -v "rsID" temp_snp_list | grep -v "SNP" | sort -u > hm3_plus_gwas_sgnif_plus_pqtls_snp_list
rm temp_snp_list 