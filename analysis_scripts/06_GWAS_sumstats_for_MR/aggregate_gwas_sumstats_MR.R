####script for aggregating and prepping summary stats for MR and coloc ######
################## helper functions #############################
#function for meta analysing using fixed effects meta-analysis
meta <- function(b1, b2=NA, b3=NA, b4=NA, b5=NA, 
                 se1, se2=NA, se3=NA, se4=NA, se5=NA){
  #organize betas and standard errors
  bs <- c(b1, b2, b3, b4, b5)
  ses <- c(se1, se2, se3, se4, se5)
  bs <- bs[!is.na(bs) & !is.na(ses)]
  ses <- ses[!is.na(c(b1, b2, b3, b4, b5)) & !is.na(ses)]
  #if only one pop, just return that pop's sum starts
  if(length(bs) ==1){
    
    beta <- bs
    se <- ses
  }else{
    #intermediate statistics
    w <- 1/ses^2
    se <- sqrt(1/sum(w))
    beta <- sum(bs*w)/sum(w)
  }
  #zstatistic and pvals
  z <- beta/se
  pval <- 2*pnorm(abs(-z), lower.tail = FALSE)
  out <- data.frame(B=beta, SE=se, Z=z, P=pval)
  return(out)
}

flip <- function(a){
  a <- ifelse(a == "A", "T",
              ifelse(a == "T", "A",
                     ifelse(a == "C", "G",
                            ifelse(a == "G", "C", "N"))))
  return(a)
}

check_alleles <- function(ref1, alt1, ref2, alt2){
  check <- ifelse(ref1 == ref2 & alt1 == alt2, "SAME",
                  ifelse(ref1 == alt2 & alt1 == ref2, "SWITCH",
                         ifelse(ref1 == flip(ref2) & alt1 == flip(alt2), "FLIP",
                                ifelse(ref1 == flip(alt2) & alt1 == flip(ref2), "SWITCH_FLIP", "FAIL"))))
  return(check)
}



###############read in data #########################
results_dir <- "./data/mr_data/regenie_sumstats/"

########read in BMI-unadjusted results #####

for(trait in c("T2D","CAD","CKD","NASH_AST_ALT")){
  for(pop in c("AFR", "EAS", "EUR", "SAS", "AMR")){
    for(chr in 1:22){
      gwas_file <- paste0(results_dir, "UKB_", pop, "_step2_BT_chr", chr, "_", trait, ".regenie")
      if(file.exists(gwas_file)){
        foo <- data.table::fread(gwas_file, data.table = FALSE)
        foo$TRAIT <- trait
        foo$POP <- pop
        
        if(chr == 1 & pop == "AFR" & trait == "T2D"){
          gwas_dat <- foo
        }else{
          foo <- foo[colnames(gwas_dat)]
          gwas_dat <- rbind(gwas_dat, foo)
        }
      }else{
        print(paste(gwas_file, "does not exist"))
      }
    }
    
  }
}
gwas_dat$BMI_ADJ <- FALSE


####read in BMI-adjusted T2D results #######

trait <- "T2D"
for(pop in c("AFR", "EAS", "EUR", "SAS", "AMR")){
  for(chr in 1:22){
    gwas_file <- paste0(results_dir, "UKB_", pop, "_BMIadj_step2_BT_chr", chr, "_", trait, ".regenie")
    if(file.exists(gwas_file)){
      foo <- data.table::fread(gwas_file, data.table = FALSE)
      foo$TRAIT <- trait
      foo$POP <- pop
      foo$BMI_ADJ <- TRUE
      foo <- foo[colnames(gwas_dat)]
      gwas_dat <- rbind(gwas_dat, foo)
    
    }else{
      print(paste(gwas_file, "does not exist"))
    }
  }
  
}

###read in BMI results
trait <- "BMI"
for(pop in c("AFR", "EAS", "EUR", "SAS", "AMR", "T2D")){
  for(chr in 1:22){
    gwas_file <- paste0(results_dir, "UKB_", pop, "_BMIonly_step2_QT_chr", chr, "_", trait, ".regenie")
    if(file.exists(gwas_file)){
      foo <- data.table::fread(gwas_file, data.table = FALSE)
      foo$TRAIT <- trait
      foo$POP <- pop
      foo$A1FREQ_CASES <- NA
      foo$A1FREQ_CONTROLS <- NA
      foo$BMI_ADJ <- FALSE
      foo <- foo[colnames(gwas_dat)]
      gwas_dat <- rbind(gwas_dat, foo)
    }else{
      print(paste(gwas_file, "does not exist"))
    }
  }
}



#read in T2D-cases results
pop <- "T2D"
for(trait in c("insulin","DR","chronic_renal","acute_renal","CAD","HF","MI","hypertension", 
               "stroke","hyperlipidemia","DN","CKD","NASH_AST_ALT","NASH_ICD10","obese")){
  for(chr in 1:22){
    gwas_file <- paste0(results_dir, "UKB_", pop, "_step2_BT_chr", chr, "_", trait, ".regenie")
    if(file.exists(gwas_file)){
      foo <- data.table::fread(gwas_file, data.table = FALSE)
      foo$TRAIT <- trait
      foo$POP <- pop
      foo$BMI_ADJ <- FALSE
      foo <- foo[colnames(gwas_dat)]
      gwas_dat <- rbind(gwas_dat, foo)
    }else{
      print(paste(gwas_file, "does not exist"))
    }
  }
}


#Done reading in data

########prepare data for coloc ########################
######## Will generate sum stats for coloc for any variant with p-value < 1E-6 #################
thresh <- -log10(1E-6)
coloc_dat <- gwas_dat[gwas_dat$LOG10P > thresh & is.na(gwas_dat$EXTRA),]

#prepare snps list
chroms <- unique(coloc_dat$CHROM)
chroms <- chroms[order(chroms)]

#read in INFO/MAF files to get snp lists
for(chr in chroms){
  #get snp list
  foo <- data.table::fread(paste0("../datasets/ukb/Imputation_genotype/ukb_c", chr, "_b0_v3.mfi.txt"))
  foo$CHR <- chr
  if(chr == chroms[1]){
    snps <- foo
  }else{
    snps <- rbind(snps, foo)
  }
}
colnames(snps) <- c("ID", "rsID", "BP", "REF", "ALT", "MAF", "A1", "INFO", "CHR")

for(i in 1:nrow(coloc_dat)){
  snp <- coloc_dat$ID[i]
  chr <- coloc_dat$CHROM[i]
  bp <- coloc_dat$GENPOS[i]
  
  trait <- coloc_dat$TRAIT[i]
  pop <- coloc_dat$POP[i]
  
  start <- bp - 250000
  end <- bp + 250000
  
  snp_set <- unique(snps$rsID[snps$CHR == chr & snps$BP >= start & snps$BP <= end])
  
  #write out as trait_pop_chr_rsID_snplist
  out <- "./data/coloc_data/regenie_inputs/"
  write.table(snp_set, paste(paste0(out, trait), pop, chr, snp, "snplist", sep="_"),
              sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

}
#quick check
for(i in 1:nrow(coloc_dat)){
  
  snp <- coloc_dat$ID[i]
  chr <- coloc_dat$CHROM[i]
  trait <- coloc_dat$TRAIT[i]
  pop <- coloc_dat$POP[i]
  
  out <- "./data/coloc_data/regenie_inputs/"
  if(file.exists(paste(paste0(out, trait), pop, chr, snp, "snplist", sep="_"))){
    print(paste(paste(paste0(out, trait), pop, chr, snp, "snplist", sep="_"), "exists"))
  }else{
    print(paste(paste(paste0(out, trait), pop, chr, snp, "snplist", sep="_"), "does not exist!"))
  }
}

#subset file to just trait pop chr rsID snplist
coloc_dat <- coloc_dat[c("TRAIT", "POP", "CHROM", "ID")]

write.table(coloc_dat, paste0(out, "index_snplist_for_coloc.txt"),
            sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

#clean up
remove(snps)
gc()

################## prepare data for MR #########################
#read in pqtl list to add additional information 
pqtls <- data.table::fread("./data/pqtl_data/cond_indep_pqtl_list.csv.gz")

gwas_dat <- merge(gwas_dat, pqtls, by.x="ID", by.y = "rsID", all = FALSE, sort=FALSE)

#split out results from T2D cases and write out, no further processing needed
t2d_dat <- gwas_dat[gwas_dat$POP == "T2D",]
data.table::fwrite(t2d_dat, "./data/mr_data/T2D_cases_gwas_sumstats.tsv.gz",
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE, compress = "gzip")




#process data for general population
#BMI-unadjusted results
for(trait in c("CKD", "BMI", "CAD",  "NASH_AST_ALT", "T2D")){
  meta_dat <- as.data.frame(gwas_dat[gwas_dat$POP == "EUR" & gwas_dat$TRAIT == trait & gwas_dat$BMI_ADJ == FALSE,])
  
  for(pop in c("EAS", "AMR", "AFR", "SAS")){
    foo <- as.data.frame(gwas_dat[gwas_dat$POP == pop & gwas_dat$TRAIT == trait  & gwas_dat$BMI_ADJ == FALSE,])
    foo <- foo[c("ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "CHISQ", "LOG10P")]
    colnames(foo)[2:ncol(foo)] <- paste0(colnames(foo)[2:ncol(foo)], "_", pop)
    foo <- foo[!duplicated(foo),]
    
    meta_dat <- merge(meta_dat, foo, by="ID", all.x=TRUE, sort=FALSE)
    
    #match direction of effect
    meta_dat[[paste0("BETA_", pop)]] <- ifelse(meta_dat[[paste0("ALLELE1_", pop)]] == meta_dat$ALLELE1, 
                                               meta_dat[[paste0("BETA_", pop)]], meta_dat[[paste0("BETA_", pop)]]*-1)
  }
  
  #subset to just desired columns
  meta_dat <- meta_dat[c( "TRAIT", "ID", "UKBPPP_ProteinID","Variant_ID", "cis_trans",
                          "CHROM", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", 
                          "N",  "BETA",  "SE", "CHISQ", "LOG10P",
                          "A1FREQ_EAS", "N_EAS", "BETA_EAS", "SE_EAS", "CHISQ_EAS", "LOG10P_EAS",
                          "A1FREQ_AMR", "N_AMR", "BETA_AMR", "SE_AMR", "CHISQ_AMR", "LOG10P_AMR", 
                          "A1FREQ_AFR", "N_AFR", "BETA_AFR", "SE_AFR", "CHISQ_AFR", "LOG10P_AFR",
                          "A1FREQ_SAS", "N_SAS", "BETA_SAS", "SE_SAS", "CHISQ_SAS", "LOG10P_SAS")]
  
  #meta-analyze
  #trans-ancestry meta-analysis without GC correction
  trans.meta <- as.data.frame(t(mapply(meta, b1=meta_dat$BETA, b2=meta_dat$BETA_AFR, b3=meta_dat$BETA_AMR, b4=meta_dat$BETA_EAS, b5=meta_dat$BETA_SAS, 
                                       se1=meta_dat$SE, se2=meta_dat$SE_AFR, se3=meta_dat$SE_AMR, se4=meta_dat$SE_EAS, se5=meta_dat$SE_SAS)))
  #prepare results
  trans.meta$B <- unlist(trans.meta$B)
  trans.meta$SE <- unlist(trans.meta$SE)
  trans.meta$Z <- unlist(trans.meta$Z)
  trans.meta$P <- unlist(trans.meta$P)
  colnames(trans.meta) <- paste0(colnames(trans.meta), "_META")
  #add to data frame
  meta_dat <- cbind(meta_dat, trans.meta)
  
  data.table::fwrite(meta_dat, paste0("./data/mr_data/", trait, "_gwas_sumstats.tsv.gz"),
                     sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE, compress = "gzip")
}



#T2D (BMI-adjusted)
trait <- "T2D"
meta_dat <- as.data.frame(gwas_dat[gwas_dat$POP == "EUR" & gwas_dat$TRAIT == trait & gwas_dat$BMI_ADJ == TRUE,])

for(pop in c("EAS", "AMR", "AFR", "SAS")){
  foo <- as.data.frame(gwas_dat[gwas_dat$POP == pop & gwas_dat$TRAIT == trait  & gwas_dat$BMI_ADJ == TRUE,])
  foo <- foo[c("ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "CHISQ", "LOG10P")]
  colnames(foo)[2:ncol(foo)] <- paste0(colnames(foo)[2:ncol(foo)], "_", pop)
  foo <- foo[!duplicated(foo),]
  
  meta_dat <- merge(meta_dat, foo, by="ID", all.x=TRUE, sort=FALSE)
  
  #match direction of effect
  meta_dat[[paste0("BETA_", pop)]] <- ifelse(meta_dat[[paste0("ALLELE1_", pop)]] == meta_dat$ALLELE1, 
                                             meta_dat[[paste0("BETA_", pop)]], meta_dat[[paste0("BETA_", pop)]]*-1)
}

#subset to just desired columns
meta_dat <- meta_dat[c( "TRAIT", "ID", "UKBPPP_ProteinID","Variant_ID", "cis_trans",
                        "CHROM", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", 
                        "N",  "BETA",  "SE", "CHISQ", "LOG10P",
                        "A1FREQ_EAS", "N_EAS", "BETA_EAS", "SE_EAS", "CHISQ_EAS", "LOG10P_EAS",
                        "A1FREQ_AMR", "N_AMR", "BETA_AMR", "SE_AMR", "CHISQ_AMR", "LOG10P_AMR", 
                        "A1FREQ_AFR", "N_AFR", "BETA_AFR", "SE_AFR", "CHISQ_AFR", "LOG10P_AFR",
                        "A1FREQ_SAS", "N_SAS", "BETA_SAS", "SE_SAS", "CHISQ_SAS", "LOG10P_SAS")]

#meta-analyze
#trans-ancestry meta-analysis without GC correction
trans.meta <- as.data.frame(t(mapply(meta, b1=meta_dat$BETA, b2=meta_dat$BETA_AFR, b3=meta_dat$BETA_AMR, b4=meta_dat$BETA_EAS, b5=meta_dat$BETA_SAS, 
                                     se1=meta_dat$SE, se2=meta_dat$SE_AFR, se3=meta_dat$SE_AMR, se4=meta_dat$SE_EAS, se5=meta_dat$SE_SAS)))
#prepare results
trans.meta$B <- unlist(trans.meta$B)
trans.meta$SE <- unlist(trans.meta$SE)
trans.meta$Z <- unlist(trans.meta$Z)
trans.meta$P <- unlist(trans.meta$P)
colnames(trans.meta) <- paste0(colnames(trans.meta), "_META")
#add to data frame
meta_dat <- cbind(meta_dat, trans.meta)

#change trait to "T2D_BMIadj"
trait <- "T2D_BMIadj"
meta_dat$TRAIT <- trait
data.table::fwrite(meta_dat, paste0("./data/mr_data/", trait, "_gwas_sumstats.tsv.gz"),
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE, compress = "gzip")
