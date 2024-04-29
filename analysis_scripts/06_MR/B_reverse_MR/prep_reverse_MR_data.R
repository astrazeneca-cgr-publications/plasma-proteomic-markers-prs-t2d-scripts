#script for building files for reverse MR
argsv <- commandArgs(trailingOnly = TRUE)
protein <- argsv[1]

#read in combined pqtl sumstats
#file format:ABO_P16442_OID30675_v1_Inflammation_II.tab.gz 
comb <- data.table::fread(paste0("../datasets/ukb//", protein, ".tab.gz"))

#read in disc pqtl sumstats
disc <- data.table::fread(paste0("../datasets/ukb/", protein, ".tab.gz"))

for(trait in c("BMI", "T2D", "T2D_BMIadj", "CAD", "CKD", "NASH_AST_ALT")){
  
  #read in gwas-signiificant data to get chromosomes
  gwas <- data.table::fread(paste0("./data/reverse_mr_data/sumstats/", trait, "_GWAS_signif_results.tsv.gz"))
  
  chroms <- unique(gwas$CHR)
  
  snps <- c()
  for(chr in chroms){
    
    #read in clumping information
    clumped_file <- paste0("./data/reverse_mr_data/clumped/", trait, "_chr", chr, ".clumped")
    if(file.exists(clumped_file)){
      foo <- data.table::fread(clumped_file)
      snps <- c(snps, foo$SNP)
    }else{
      next
    }
  
    
    #read in data
    if(trait == "BMI"){
      f <- paste0("./data/coloc_data/regenie_outputs/UKB_EUR_BMIonly_step2_QT_chr", chr, "_BMI.regenie.gz")
    }else if(trait == "T2D_BMIadj"){
      f <- paste0("./data/coloc_data/regenie_outputs/UKB_EUR_BMIadj_step2_BT_chr", chr, "_T2D.regenie.gz")
    }else{
      f <- paste0("./data/coloc_data/regenie_outputs/UKB_EUR_step2_BT_chr", chr, "_", trait, ".regenie.gz")
    }
    
    if(file.exists(f)){
      foo <- data.table::fread(f)
    }else{
      print(paste("GWAS not done for", trait, "on", chr))
    }
    
    #subset just to clumped snps
 
    foo <- foo[foo$ID %in% snps,]
    if(chr == chroms[1]){
      dat <- foo
    }else{
      dat <- rbind(dat, foo)
    }
  }
  
  #add pqtl sumstats
  dat$IDV1 <- paste(dat$CHR, dat$GENPOS, dat$ALLELE0, dat$ALLELE1, "imp:v1", sep = ":")
  dat$IDV2 <- paste(dat$CHR, dat$GENPOS, dat$ALLELE1, dat$ALLELE0, "imp:v1", sep = ":")
  
  #subset combined pqtl sumstats and merge
  foo <- comb[comb$ID %in% c(dat$IDV1, dat$IDV2),]
  dat$ID_MATCH <- ifelse(dat$IDV1 %in% foo$ID, dat$IDV1, dat$IDV2)
  
  colnames(foo)[c(1:2, 4:ncol(foo))] <- paste0(colnames(foo)[c(1:2, 4:ncol(foo))], "_COMB")
  
  dat <- merge(dat, foo, by.x="ID_MATCH", by.y="ID", all=FALSE, sort=FALSE)
  
  #subset discovery pqtl sumstats and merge
  foo <- disc[disc$ID %in% c(dat$IDV1, dat$IDV2),]
  dat$ID_MATCH <- ifelse(dat$IDV1 %in% foo$ID, dat$IDV1, dat$IDV2)
  
  colnames(foo)[c(1:2, 4:ncol(foo))] <- paste0(colnames(foo)[c(1:2, 4:ncol(foo))], "_DISC")
  
  dat <- merge(dat, foo, by.x="ID_MATCH", by.y="ID", all=FALSE, sort=FALSE)
  
  #align directions of effect
  dat$BETA_COMB <- ifelse(dat$ALLELE1 == dat$ALLELE1_COMB, dat$BETA_COMB, 
                          dat$BETA_COMB*-1)
  dat$BETA_DISC <- ifelse(dat$ALLELE1 == dat$ALLELE1_DISC, dat$BETA_DISC, 
                          dat$BETA_DISC*-1)
  
  #trim down the columns
  dat <- dat[,c(4,1, 2,3, 5:14, 22:24, 26:29, 35:37, 39:42)]
  
  #write out
  data.table::fwrite(dat, paste0("./data/reverse_mr_data/inputs/", 
                                 trait, "_", protein, "_reverse_MR_inputs.tsv"),
                     sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)
  
}


