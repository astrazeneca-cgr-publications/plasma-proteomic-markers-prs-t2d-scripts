###script for parsing GWAS data for reverse MR ##

#log10p threshold for subsetting
thresh <- -log10(5E-08)

#read in ukb GWAS
lambdas <- c()
snps <- c()
chrom <- c()
for(trait in c("BMI", "T2D", "T2D_BMIadj", "CAD", "CKD", "NASH_AST_ALT")){
  pvals <- c()
  for(chr in 1:22){
    
    #read in data
    if(trait == "BMI"){
      f <- paste0("data/coloc_data/regenie_outputs/UKB_EUR_BMIonly_step2_QT_chr", chr, "_BMI.regenie.gz")
    }else if(trait == "T2D_BMIadj"){
      f <- paste0("data/coloc_data/regenie_outputs/UKB_EUR_BMIadj_step2_BT_chr", chr, "_T2D.regenie.gz")
    }else{
      f <- paste0("data/coloc_data/regenie_outputs/UKB_EUR_step2_BT_chr", chr, "_", trait, ".regenie.gz")
    }
    
    if(file.exists(f)){
      foo <- data.table::fread(f)
    }else{
      print(paste("GWAS not done for", trait, "on", chr))
    }
    #save pvals
    pvals <- c(pvals, foo$LOG10P[!is.na(foo$LOG10P)])
    #subset
    foo <- foo[foo$LOG10P > thresh & !is.na(foo$LOG10P),]
    if(chr == 1){
      dat <- foo
    }else{
      dat <- rbind(dat, foo)
    }
  }
  
  pvals <- 10^-pvals
  
  #calculate GC lambda
  z = qnorm(pvals/2)
  lambda = median(z^2) / qchisq(0.5,1)
  print(paste("GC lambda for", trait, "=", lambda))
  lambdas <- c(lambdas, lambda)
  
  #save snp info
  chrom <- c(chrom, dat$CHROM)
  snps <- c(snps, dat$ID)
  
  
  #prepare data so that it has plink assoc format
  dat <- as.data.frame(dat)
  dat$Lxy <- NA
  dat$Uxy <- NA
  dat <- dat[c("CHROM", "ID", "GENPOS", "ALLELE1", "TEST", "N", "BETA", "SE", "Lxy", "Uxy",
               "CHISQ", "LOG10P")]
  colnames(dat) <- c("CHR", "SNP", "BP", "A1", "TEST", "NMISS", "BETA",
                     "SE", "Lxy", "Uxy", "STAT", "P")
  
  dat$P <- 10^-dat$P
  #write out
  data.table::fwrite(dat, paste0("data/reverse_mr_data/sumstats/", trait, "_GWAS_signif_results.tsv.gz"),
                     sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE,
                     na="NA",
                     compress = "gzip")
  remove(dat)
  remove(foo)
  gc()
}

#prepare list of snps for extraction/clumping and write out
snp_list <- data.frame(CHROM=chrom,
                       SNP=snps)
snp_list <- snp_list[!duplicated(snp_list),]

data.table::fwrite(snp_list, "data/reverse_mr_data/sumstats/GWAS_signif_snp_list",
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)
#prepare gc lambda list and write out
gc_lambda <- data.frame(TRAIT=c("BMI", "T2D", "T2D_BMIadj", "CAD", "CKD", "NASH_AST_ALT"),
                        LAMBDA=lambdas)
data.table::fwrite(gc_lambda, "data/reverse_mr_data/sumstats/GWAS_gc_lambda.txt",
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)

