#small script/pipeline for prepare sumstats for prs-csx or prs-csx

export PATH="../envs/R/bin:$PATH"

prefix=specify_prefix
sumstats=path_to_sumstats

#give column names here
a1=effect_allele
a2=other_allele
beta=beta
pval=p_value
N=sample_size #can be column name or just give sample size, e.g. 100000
bp=base_pair
chrom=chromosome
snp=variant
se=standard_error

#this is the path to the listing of variants in the ld reference. Optional, but reduces file sizes
hm3=../ld_ref/snpinfo_mult_1kg_hm3
#hg38 version, if needed:
#hm3=../ld_ref/snpinfo_hm3_hg38

Rscript parse_sumstats.R --prefix $prefix \
--sumstats $sumstats --a1 $a1 --a2 $a2 --beta $beta --pval $pval --snp $snp \
--N $N --se $se --bp $bp --chrom $chrom --format prscs --match $hm3 

#use this option only if making bim from gwas sumstats. Otherwise, need to make a bim file from target genotype data
#bim file format: tab sep file, no header with columns chrom,snp,cm (set to 0), bp, A1, A2
Rscript make_bim.R --prefix $prefix \
--sumstats $sumstats --a1 $a1 --a2 $a2  --snp $snp \
--bp $bp --chrom $chrom --match $hm3 

#put input files into input subdirectory, if desired
mkdir -p inputs
mv $prefix*.bim inputs 
mv $prefix*.parsed.txt inputs 