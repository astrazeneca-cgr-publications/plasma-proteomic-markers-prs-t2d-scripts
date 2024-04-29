#script for running 2-sample MR analysis with protein as exposure
#get protein from command line
argsv <- commandArgs(trailingOnly = TRUE)
protein <- argsv[1]
library(MendelianRandomization)

####preparing pqtl data ###########
#read in pqtl list 
pqtls <- data.table::fread("./data/pqtl_data/cond_indep_pqtl_list.csv.gz")
pqtls <- pqtls[pqtls$UKBPPP_ProteinID == protein,]
if(nrow(pqtls) == 0) stop(paste0(protein, "does not have any pQTLs"))
#pqtls$Variant_ID <- gsub("v2", "v1", pqtls$Variant_ID)

#prepare for discovery only analyses
pqtls$A1_DISC <- NA
pqtls$A2_DISC <- NA
pqtls$BETA_DISC <- NA
pqtls$SE_DISC <- NA
disc_files <- list.files("../datasets/ukb/", ".tab.gz")
disc_files <- unique(gsub(".tbi", "", disc_files))

#prepare for combined analyses 
pqtls$A1_COMB <- NA
pqtls$A2_COMB <- NA
pqtls$BETA_COMB <- NA
pqtls$SE_COMB <- NA
comb_files <- list.files("../datasets/ukb//", ".tab.gz")
comb_files <- unique(gsub(".tbi", "", comb_files))

#if data has already been prepared, read in. Otherwise, prepare summary stats
#file name is the same in both combined and discovery analysis
pqtl_file <- disc_files[grep(gsub(":", "_", protein), disc_files)]
all_file <- paste0("./data/mr_data/mr_inputs/", 
                   gsub(".tab.gz", "", basename(pqtl_file)), 
                   ".all_subjects_snp_dat.tsv")

t2d_cases_file <- paste0("./data/mr_data/mr_inputs/", 
                         gsub(".tab.gz", "", basename(pqtl_file)), 
                         ".t2d_cases_snp_dat.tsv")
if(file.exists(all_file) & file.exists(t2d_cases_file)){
  mr_dat <- data.table::fread(all_file, data.table = FALSE)
  t2d_dat <- data.table::fread(t2d_cases_file, data.table = FALSE)
}else{
  #if needed, read in pqtl gwas results and add to pqtl data frame
  #read in discovery results
  
  foo <- data.table::fread(paste0("../datasets/ukb/", pqtl_file))
  foo <- foo[foo$ID %in% pqtls$Variant_ID[pqtls$UKBPPP_ProteinID == protein],]
  gc()
  
  #add to data frame
  for(snp in foo$ID){
    pqtls$A1_DISC[pqtls$Variant_ID == snp & pqtls$UKBPPP_ProteinID == protein] <- foo$ALLELE1[foo$ID == snp]
    pqtls$A2_DISC[pqtls$Variant_ID == snp & pqtls$UKBPPP_ProteinID == protein] <- foo$ALLELE0[foo$ID == snp]
    pqtls$BETA_DISC[pqtls$Variant_ID == snp & pqtls$UKBPPP_ProteinID == protein] <- foo$BETA[foo$ID == snp]
    pqtls$SE_DISC[pqtls$Variant_ID == snp & pqtls$UKBPPP_ProteinID == protein] <- foo$SE[foo$ID == snp]
    
  }
  
  #read in combined results
  pqtl_file <- comb_files[grep(gsub(":", "_", protein), comb_files)]
  foo <- data.table::fread(paste0("../datasets/ukb//", pqtl_file))
  foo <- foo[foo$ID %in% pqtls$Variant_ID[pqtls$UKBPPP_ProteinID == protein],]
  gc()
  
  #add to data frame
  for(snp in foo$ID){
    pqtls$A1_COMB[pqtls$Variant_ID == snp & pqtls$UKBPPP_ProteinID == protein] <- foo$ALLELE1[foo$ID == snp]
    pqtls$A2_COMB[pqtls$Variant_ID == snp & pqtls$UKBPPP_ProteinID == protein] <- foo$ALLELE0[foo$ID == snp]
    pqtls$BETA_COMB[pqtls$Variant_ID == snp & pqtls$UKBPPP_ProteinID == protein] <- foo$BETA[foo$ID == snp]
    pqtls$SE_COMB[pqtls$Variant_ID == snp & pqtls$UKBPPP_ProteinID == protein] <- foo$SE[foo$ID == snp]
    
  }
  
  ####compute f-statistic just using beta^2/se^2 formula#####
  pqtls$FSTAT_DISC <- pqtls$BETA_DISC^2/pqtls$SE_DISC^2
  pqtls$FSTAT_COMB <- pqtls$BETA_COMB^2/pqtls$SE_COMB^2
  
  ########## read in GWAS data ############################
  for(trait in c("CKD", "BMI", "CAD",  "NASH_AST_ALT", "T2D", "T2D_BMIadj")){
    foo <- data.table::fread(paste0("./data/mr_data/", trait, "_gwas_sumstats.tsv.gz"))
    foo <- foo[,-c(3:5)]
    foo <- foo[!duplicated(foo),]
    foo <- foo[foo$ID %in% pqtls$rsID,]
    if(trait == "CKD"){
      mr_dat <- foo
    }else{
      mr_dat <- rbind(mr_dat, foo)
    }
  }
  
  mr_dat <- merge(as.data.frame(mr_dat), as.data.frame(pqtls), 
                  by.x="ID", by.y="rsID", all=FALSE, sort=FALSE)
  
  #align with mr dat
  t2d_dat <- data.table::fread("./data/mr_data/T2D_cases_gwas_sumstats.tsv.gz",
                               data.table = FALSE)
  
  for(varid in colnames(foo)){
    if(varid %in% colnames(t2d_dat)){
      next
    }else{
      t2d_dat[[varid]] <- NA
    }
  }
  for(varid in colnames(t2d_dat)){
    if(varid %in% colnames(foo)){
      next
    }else{
      t2d_dat[[varid]] <- NULL
    }
  }
  t2d_dat <- t2d_dat[t2d_dat$ID %in% pqtls$rsID,]
  t2d_dat <- t2d_dat[!duplicated(t2d_dat),]
  t2d_dat <- merge(t2d_dat, as.data.frame(pqtls), 
                   by.x="ID", by.y="rsID", all=FALSE, sort=FALSE)
  t2d_dat <- t2d_dat[colnames(mr_dat)]
  t2d_dat$TRAIT <- ifelse(t2d_dat$TRAIT %in% c("CKD", "BMI", "CAD",  "NASH_AST_ALT"),
                          paste0(t2d_dat$TRAIT, "_T2Dcases"), t2d_dat$TRAIT)
  mr_dat <- rbind(mr_dat, t2d_dat)
  
  #save mr data (gen pop traits) and t2d data (results from T2D cases gwas)
  
  data.table::fwrite(mr_dat, all_file, sep='\t', quote=FALSE, na="NA", 
                     col.names = TRUE, row.names = FALSE)
  data.table::fwrite(t2d_dat, t2d_cases_file, sep='\t', quote=FALSE, na="NA", 
                     col.names = TRUE, row.names = FALSE)
}


################ run MR ##################################

traits <- unique(mr_dat$TRAIT)
gen_pop_traits <- c("CKD", "BMI", "CAD",  "NASH_AST_ALT", "T2D_BMIadj", "T2D")

#prepare results data frame
mr_results <- expand.grid(PROTEIN=protein,
                          N_SNPS=length(unique(mr_dat$ID)),
                          N_CIS=length(unique(mr_dat$ID[mr_dat$cis_trans == "cis"])),
                          TRAIT=traits, 
                          CIS_ONLY=c(TRUE, FALSE),
                          SUBSET=c("ALL", "EUROPEAN_ONLY"),
                          MAX_GWAS_LOG10P=NA, MEDIAN_FSTAT=NA, 
                          MIN_FSTAT=NA, PROP_FSTAT_10=NA,
                          FSTAT_LOG10P_R=NA, BETAS_R=NA,
                          FSTAT_FILTERED=c(TRUE, FALSE),
                          EST_MEDIAN=NA, SE_MEDIAN=NA, P_MEDIAN=NA,
                          EST_WEIGHTED=NA, SE_WEIGHTED=NA, P_WEIGHTED=NA,
                          EST_IVW=NA, SE_IVW=NA, P_IVW=NA,
                          EST_EGGER=NA, SE_EGGER=NA, P_EGGER=NA,
                          P_INTERCEPT=NA,
                          EST_ENSEMB=NA, SE_ENSEMB=NA, P_ENSEMB=NA,
                          stringsAsFactors = FALSE)

for(trait in traits){
  #running each trait w/ and w/out F-stat filtering
  for(fstat in c(TRUE, FALSE)){
    #cis only or cis+trans 
    for(cis in c(TRUE, FALSE)){
      #subset to either cis only or all
      if(cis){
        foo <- mr_dat[mr_dat$cis_trans == "cis" & mr_dat$TRAIT == trait,]
      }else{
        foo <- mr_dat[mr_dat$TRAIT == trait,]
      } 
      #skip if fewer than 3 variants in subset
   
      if(nrow(foo) < 3) next
      
      #all subjects or european-ancestry only
      for(pop in c("ALL", "EUROPEAN_ONLY")){
        
        if(fstat & pop == "ALL") foo <- foo[foo$FSTAT_COMB > 10,]
        if(fstat & pop == "EUROPEAN_ONLY") foo <- foo[foo$FSTAT_DISC > 10,]
        if(nrow(foo) < 3) next
        
        #add summary level information for pqtl/trait gwas inputs
        if(trait %in% gen_pop_traits){
          mr_results$MAX_GWAS_LOG10P[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- ifelse(pop == "ALL", 
                                                                                                                                  -log10(min(foo$P_META, na.rm = TRUE)), max(foo$LOG10P, na.rm = TRUE))
          mr_results$MEDIAN_FSTAT[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- ifelse(pop == "ALL",
                                                                                                                               median(foo$FSTAT_COMB, na.rm = TRUE), median(foo$FSTAT_DISC, na.rm = TRUE))
          mr_results$MIN_FSTAT[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- ifelse(pop == "ALL",
                                                                                                                            min(foo$FSTAT_COMB, na.rm = TRUE), min(foo$FSTAT_DISC, na.rm = TRUE))
          mr_results$PROP_FSTAT_10[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- ifelse(pop == "ALL",
                                                                                                                                nrow(foo[foo$FSTAT_COMB >= 10,])/nrow(foo),
                                                                                                                                nrow(foo[foo$FSTAT_DISC >= 10,])/nrow(foo))
          
          mr_results$FSTAT_LOG10P_R[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- ifelse(pop == "ALL", cor(foo$FSTAT_COMB, -log10(foo$P_META), use = "complete.obs"),
                                                                                                                                 cor(foo$FSTAT_DISC, foo$LOG10P, use = "complete.obs"))
          mr_results$BETAS_R[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- ifelse(pop == "ALL",
                                                                                                                          cor(foo$BETA_COMB, foo$B_META, use = "complete.obs"), cor(foo$BETA_DISC, foo$BETA, use = "complete.obs"))
        }else{
          mr_results$MAX_GWAS_LOG10P[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait& mr_results$FSTAT_FILTERED == fstat] <- max(foo$LOG10P, na.rm = TRUE)
          mr_results$MEDIAN_FSTAT[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- median(foo$FSTAT_DISC, na.rm = TRUE)
          mr_results$MIN_FSTAT[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <-  min(foo$FSTAT_DISC, na.rm = TRUE)
          mr_results$PROP_FSTAT_10[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- nrow(foo[foo$FSTAT_DISC >= 10,])/nrow(foo)
          
          mr_results$FSTAT_LOG10P_R[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- cor(foo$FSTAT_DISC, foo$LOG10P, use = "complete.obs")
          mr_results$BETAS_R[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- cor(foo$BETA_DISC, foo$BETA, use = "complete.obs")
        }
        
        #prepare gwas and pqtl betas/ses
        if(pop == "ALL"){
          #Combined
          pqtl_betas <- foo$BETA_COMB
          pqtl_se <- foo$SE_COMB
          if(trait %in% gen_pop_traits){
            gwas_betas <- ifelse(foo$ALLELE1 == foo$A1_COMB, foo$B_META, foo$B_META*-1)
            gwas_se <- foo$SE_META
          }else{
            gwas_betas <- ifelse(foo$ALLELE1 == foo$A1_COMB, foo$BETA, foo$BETA*-1)
            gwas_se <- foo$SE
          }
        }else{
          #discovery only
          pqtl_betas <- foo$BETA_DISC
          pqtl_se <- foo$SE_DISC
          gwas_betas <- ifelse(foo$ALLELE1 == foo$A1_DISC, foo$BETA, foo$BETA*-1)
          gwas_se <- foo$SE
        }
        
        #prepare MR input
        dat.input <- mr_input(snps=foo$ID,
                              bx= pqtl_betas,
                              bxse=pqtl_se,
                              by=gwas_betas,
                              byse = gwas_se,
                              exposure=protein,
                              outcome=trait)
        
        #run analysis
        tryCatch({ mr.run <- mr_allmethods(dat.input, iterations=10000, method = "main")},
                 warning=function(w) print(paste("Warning for", protein, "and", trait)))
        
        #parse
        results <- as.data.frame(mr.run@Values)
        #median
        mr_results$EST_MEDIAN[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- results$Estimate[1]
        mr_results$SE_MEDIAN[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- results$`Std Error`[1]
        mr_results$P_MEDIAN[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- results$`P-value`[1]
        #weighted median
        mr_results$EST_WEIGHTED[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- results$Estimate[2]
        mr_results$SE_WEIGHTED[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- results$`Std Error`[2]
        mr_results$P_WEIGHTED[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- results$`P-value`[2]
        #IVW
        mr_results$EST_IVW[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- results$Estimate[3]
        mr_results$SE_IVW[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- results$`Std Error`[3]
        mr_results$P_IVW[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- results$`P-value`[3]
        #Egger
        mr_results$EST_EGGER[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- results$Estimate[4]
        mr_results$SE_EGGER[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- results$`Std Error`[4]
        mr_results$P_EGGER[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- results$`P-value`[4]
        #egger intercept
        mr_results$P_INTERCEPT[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- results$`P-value`[5]
        #ensemble
        mr_results$EST_ENSEMB[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- median(results$Estimate[1:4])
        mr_results$SE_ENSEMB[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- median(results$`Std Error`[1:4])
        mr_results$P_ENSEMB[mr_results$CIS_ONLY == cis & mr_results$SUBSET == pop & mr_results$TRAIT == trait & mr_results$FSTAT_FILTERED == fstat] <- median(results$`P-value`[1:4])
      }
    }

  }
}

#recompute number of pqtls after fstat filtering. Just use T2D as representative trait as it does not differ across traits
mr_results$N_CIS[mr_results$FSTAT_FILTERED == TRUE] <- nrow(mr_dat[mr_dat$cis_trans == "cis" & mr_dat$FSTAT_DISC > 10 & mr_dat$TRAIT == "T2D_BMIadj",])
mr_results$N_SNPS[mr_results$FSTAT_FILTERED == TRUE] <- nrow(mr_dat[mr_dat$FSTAT_DISC > 10 & mr_dat$TRAIT == "T2D_BMIadj",])

#write out results
out.file <- paste0("./results/06_MR/", protein, "_MR_results.tsv")
out.file <- gsub(":", "_", out.file)
data.table::fwrite(mr_results, out.file, 
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)
