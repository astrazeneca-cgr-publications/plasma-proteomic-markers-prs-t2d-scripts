##### script for running colocalization analysis #########
library(coloc)

#get index
argsv <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(argsv[1])

#get parameters from index file
#index_file <- data.table::fread("./data/coloc_data/regenie_inputs/index_snplist_for_coloc.txt")
index_file <- data.table::fread("./data/coloc_data/regenie_inputs/index_snplist_for_coloc_new_defs.txt")
trait <- index_file$V1[index]
pop <- index_file$V2[index]
chr <- index_file$V3[index]
snp <- index_file$V4[index]

#read in gwas summary stats
gwas_dir <- "./data/coloc_data/regenie_outputs/"
if(trait == "BMI"){
  gwas <- data.table::fread(paste0(gwas_dir, "UKB_", pop, "_BMIonly_step2_QT_", trait, "_", snp, "_", trait, ".regenie.gz")) 
}else if(trait == "T2D"){
  #use bmi adjusted summary stats for T2D
  gwas <- data.table::fread(paste0(gwas_dir, "UKB_", pop, "_BMIadj_step2_BT_", trait, "_", snp, "_", trait, ".regenie.gz"))
  #change trait to reflect BMI adjustment
  trait <- "T2D_BMIadj"
}else{
  gwas <- data.table::fread(paste0(gwas_dir, "UKB_", pop, "_step2_BT_", trait, "_", snp, "_", trait, ".regenie.gz")) 
}

#create id that matches pqtl gwas
gwas$ID2 <- paste(gwas$CHROM, gwas$GENPOS, gwas$ALLELE0, gwas$ALLELE1, "imp:v1", sep=":")


#filter by MAC > 50 and INFO > 0.7
gwas$MAC <- ifelse(gwas$A1FREQ > 0.5, (1-gwas$A1FREQ)*gwas$N*2, gwas$A1FREQ*gwas$N*2)
gwas <- gwas[gwas$MAC > 50 & gwas$INFO > 0.70,]


#prepare data
#prepare for discovery only analyses
gwas$N_DISC <- NA
gwas$A1_DISC <- NA
gwas$A2_DISC <- NA
gwas$AF_DISC <- NA
gwas$BETA_DISC <- NA
gwas$SE_DISC <- NA
disc_files <- list.files("../datasets/ukb/", ".tab.gz")
disc_files <- unique(gsub(".tbi", "", disc_files))

#prepare for combined analyses
gwas$N_COMB <- NA
gwas$A1_COMB <- NA
gwas$A2_COMB <- NA
gwas$AF_COMB <- NA
gwas$BETA_COMB <- NA
gwas$SE_COMB <- NA
comb_files <- list.files("../datasets/ukb//", ".tab.gz")
comb_files <- unique(gsub(".tbi", "", comb_files))


#read in pqtls to find which protein(s) will be tested
pqtls <- data.table::fread("./data/pqtl_data/cond_indep_pqtl_list.csv.gz")
pqtls <- pqtls[pqtls$rsID == snp,]

proteins <- unique(pqtls$UKBPPP_ProteinID)


for(protein in proteins){
  pqtl_file <- disc_files[grep(gsub(":", "_", protein), disc_files)]
  foo <- data.table::fread(paste0("../datasets/ukb/", pqtl_file))
  foo <- foo[foo$`#CHROM` == chr,]
  foo <- foo[foo$ID %in% gwas$ID2,]
  gc()
  
  #add summary stats to file
  for(id in foo$ID){
    gwas$N_DISC[gwas$ID2 == id] <- foo$N[foo$ID == id]
    gwas$A1_DISC[gwas$ID2 == id] <- foo$ALLELE1[foo$ID == id]
    gwas$A2_DISC[gwas$ID2 == id] <- foo$ALLELE0[foo$ID == id]
    gwas$AF_DISC[gwas$ID2 == id] <- foo$A1FREQ[foo$ID == id]
    gwas$BETA_DISC[gwas$ID2 == id] <- foo$BETA[foo$ID == id]
    gwas$SE_DISC[gwas$ID2 == id] <- foo$SE[foo$ID == id]
  }
  
  #quit if index snp is not found due to snp id matching errors
  if(is.na(gwas$BETA_DISC[gwas$ID == snp])) stop("index snp not matched with variant from pqtl gwas")
  
  #read in combined results
  pqtl_file <- comb_files[grep(gsub(":", "_", protein), comb_files)]
  foo <- data.table::fread(paste0("../datasets/ukb//", pqtl_file))
  foo <- foo[foo$`#CHROM` == chr,]
  foo <- foo[foo$ID %in% gwas$ID2,]
  gc()
  
  #add summary stats to file
  for(id in foo$ID){
    gwas$N_COMB[gwas$ID2 == id] <- foo$N[foo$ID == id]
    gwas$A1_COMB[gwas$ID2 == id] <- foo$ALLELE1[foo$ID == id]
    gwas$A2_COMB[gwas$ID2 == id] <- foo$ALLELE0[foo$ID == id]
    gwas$AF_COMB[gwas$ID2 == id] <- foo$A1FREQ[foo$ID == id]
    gwas$BETA_COMB[gwas$ID2 == id] <- foo$BETA[foo$ID == id]
    gwas$SE_COMB[gwas$ID2 == id] <- foo$SE[foo$ID == id]
  }
  
  #quit if index snp is not found due to snp id matching errors
  if(is.na(gwas$BETA_COMB[gwas$ID == snp])) stop("index snp not matched with variant from pqtl gwas")
  
  ####discovery#######
  foo <- gwas[!is.na(gwas$BETA_DISC),]
  #split pqtl and trait gwas
  if(trait == "BMI"){
    d1 <- list(snp=foo$ID, 
               position=foo$GENPOS,
               beta=foo$BETA,
               varbeta=foo$SE^2,
               type="quant",
               N=foo$N,
               MAF=foo$A1FREQ)
  }else{
    d1 <- list(snp=foo$ID, 
               position=foo$GENPOS,
               beta=foo$BETA,
               varbeta=foo$SE^2,
               type="cc")
  }

  
  d2 <- list(snp=foo$ID, 
             position=foo$GENPOS,
             beta=foo$BETA_DISC,
             varbeta=foo$SE_DISC^2,
             type="quant",
             N=foo$N_DISC,
             MAF=foo$AF_DISC)
  test <- coloc.abf(d1, d2)
  
  results <- cbind(data.frame(protein=protein, index_snp=snp, trait=trait, PQTL_SUBSET="DISC", GWAS_POP=pop),
                   as.data.frame(t(test$summary)))
  results$SNP.PP.H4 <- test$results$SNP.PP.H4[test$results$snp == snp]
  ####combined#######
  foo <- gwas[!is.na(gwas$BETA_COMB),]
  #split pqtl and trait gwas
  if(trait == "BMI"){
    d1 <- list(snp=foo$ID, 
               position=foo$GENPOS,
               beta=foo$BETA,
               varbeta=foo$SE^2,
               type="quant",
               N=foo$N,
               MAF=foo$A1FREQ)
  }else{
    d1 <- list(snp=foo$ID, 
               position=foo$GENPOS,
               beta=foo$BETA,
               varbeta=foo$SE^2,
               type="cc")
  }
  
  
  d2 <- list(snp=foo$ID, 
             position=foo$GENPOS,
             beta=foo$BETA_COMB,
             varbeta=foo$SE_COMB^2,
             type="quant",
             N=foo$N_COMB,
             MAF=foo$AF_COMB)
  test <- coloc.abf(d1, d2)
  results_comb <- cbind(data.frame(protein=protein, index_snp=snp, trait=trait, PQTL_SUBSET="COMB", GWAS_POP=pop),
                   as.data.frame(t(test$summary)))
  results_comb$SNP.PP.H4 <- test$results$SNP.PP.H4[test$results$snp == snp]
  
  results <- rbind(results, results_comb)
  cis_trans <- pqtls$cis_trans[pqtls$UKBPPP_ProteinID == protein & pqtls$rsID == snp]
  
  #criteria: PP.H3 + PP.H4 ≥ 0.99 and PP4/PP3 ≥5
  criteria1 <- results$PP.H3.abf + results$PP.H4.abf
  criteria2 <- results$PP.H4.abf/results$PP.H3.abf
  results$SAME_SIGNAL <- ifelse(criteria1 >= 0.99 & criteria2 >=5, TRUE, FALSE)
  
  out.file <- paste0("./results/06_coloc/",
                     gsub(":",  "_", protein), "_", trait, "_", pop, "_", snp, "_",  
                     cis_trans, "_coloc_results.tsv")
  
  data.table::fwrite(results, out.file, 
                     sep = '\t', quote=FALSE, row.names = FALSE, col.names = TRUE)
}


