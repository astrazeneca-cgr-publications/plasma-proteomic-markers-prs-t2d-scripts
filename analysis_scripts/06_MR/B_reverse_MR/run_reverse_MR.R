# script for running reverse MR with the disease/trait as the exposure and protein as the outcome
#get protein from command line
argsv <- commandArgs(trailingOnly = TRUE)
protein <- argsv[1]
library(MendelianRandomization)

traits <- c("BMI", "T2D", "T2D_BMIadj", "CAD", "CKD", "NASH_AST_ALT")
mr_results <- expand.grid(PROTEIN=protein,
                          TRAIT=traits, 
                          SUBSET=c("ALL", "EUROPEAN_ONLY"),
                          N_SNPS=NA,
                          EST_MEDIAN=NA, SE_MEDIAN=NA, P_MEDIAN=NA,
                          EST_WEIGHTED=NA, SE_WEIGHTED=NA, P_WEIGHTED=NA,
                          EST_IVW=NA, SE_IVW=NA, P_IVW=NA,
                          EST_EGGER=NA, SE_EGGER=NA, P_EGGER=NA,
                          P_INTERCEPT=NA,
                          EST_ENSEMB=NA, SE_ENSEMB=NA, P_ENSEMB=NA,
                          stringsAsFactors = FALSE)

for(trait in traits){
  f <- paste0(trait, "_", protein, "_reverse_MR_inputs.tsv")
  foo <- data.table::fread(paste0("./data/reverse_mr_data/inputs/", f))
  
  for(pop in c("ALL", "EUROPEAN_ONLY")){
    gwas_betas <- foo$BETA
    gwas_se <- foo$SE
    if(pop == "ALL"){
      pqtl_betas <- foo$BETA_COMB
      pqtl_se <- foo$SE_COMB
    }else{
      #discovery only
      pqtl_betas <- foo$BETA_DISC
      pqtl_se <- foo$SE_DISC
    }
    
    #prepare MR input
    dat.input <- mr_input(snps=foo$ID,
                          bx= pqtl_betas,
                          bxse=pqtl_se,
                          by=gwas_betas,
                          byse = gwas_se,
                          exposure=trait,
                          outcome=protein)
    
    tryCatch({ mr.run <- mr_allmethods(dat.input, iterations=10000, method = "main")},
             warning=function(w) print(paste("Warning for", protein, "and", trait)))
    
    #parse
    results <- as.data.frame(mr.run@Values)
    #number of snps included
    mr_results$N_SNPS[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- nrow(foo)
    
    #median
    mr_results$EST_MEDIAN[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- results$Estimate[1]
    mr_results$SE_MEDIAN[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- results$`Std Error`[1]
    mr_results$P_MEDIAN[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- results$`P-value`[1]
    #weighted median
    mr_results$EST_WEIGHTED[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- results$Estimate[2]
    mr_results$SE_WEIGHTED[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- results$`Std Error`[2]
    mr_results$P_WEIGHTED[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- results$`P-value`[2]
    #IVW
    mr_results$EST_IVW[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- results$Estimate[3]
    mr_results$SE_IVW[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- results$`Std Error`[3]
    mr_results$P_IVW[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- results$`P-value`[3]
    #Egger
    mr_results$EST_EGGER[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- results$Estimate[4]
    mr_results$SE_EGGER[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- results$`Std Error`[4]
    mr_results$P_EGGER[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- results$`P-value`[4]
    #egger intercept
    mr_results$P_INTERCEPT[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- results$`P-value`[5]
    #ensemble
    mr_results$EST_ENSEMB[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- median(results$Estimate[1:4])
    mr_results$SE_ENSEMB[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- median(results$`Std Error`[1:4])
    mr_results$P_ENSEMB[mr_results$SUBSET == pop & mr_results$TRAIT == trait] <- median(results$`P-value`[1:4])
  }
}

#save
out.file <- paste0("./results/06_REVERSE_MR/", protein, "_reverse_MR_results.tsv")
out.file <- gsub(":", "_", out.file)
data.table::fwrite(mr_results, out.file, 
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)

