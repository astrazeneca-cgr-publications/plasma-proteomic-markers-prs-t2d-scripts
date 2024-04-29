####script for performing PRS-npx associations ########
argsv <- commandArgs(trailingOnly = TRUE)
protein <- gsub(":", "_", argsv[1])
protein <- gsub("-", "_", protein)

#####read in supporing files ######
pqtl <- data.table::fread("./data/pqtl_data/cond_indep_pqtl_list.csv.gz")
top_pgs <- data.table::fread("./results/top_PGS_per_trait.tsv")
proteins <- data.table::fread("./results/olink_protein_description.tsv")

####read in files to build analysis file #####
#read in linker file
linker.file <- "../datasets/ukb/olink_sample_map_3k_ukbsamples_forconsort_v2_linked.tsv"

linker <- data.table::fread(linker.file, data.table = FALSE)
linker <- linker[c("userID", "pseudo_ind_id")]
linker <- linker[!duplicated(linker),]

#covariate file
covar.file <- "../datasets/ukb/combined_covars_forconsortium_v1.tsv"
covar <- data.table::fread(covar.file, data.table = FALSE)
covar_rep <- data.table::fread("../datasets/ukb/replication_covars_forconsortium_v1.tsv")
covar$REP <- ifelse(covar$IID %in% covar_rep$IID, 1, 0)

#set factors
covar$Batch <- as.factor(covar$Batch)
covar$array_fct <- as.factor(covar$array_fct)
covar$ukb_centre_fct <- as.factor(covar$ukb_centre_fct)
covar$sample_selection <- as.factor(covar$sample_selection)
covar$sex <- as.factor(covar$sex)

#read in ukb-wide phenotypes and covariate files
p <- data.table::fread("./data/phenotypes/UKB.phenotypes.tsv.gz")
cov_all <- data.table::fread("./data/phenotypes/UKB.covariates.tsv.gz")

#read in relatives file
rels <- data.table::fread("./data/phenotypes/UKB.degree2.drop.tsv")

#read in incident cases file
di <- data.table::fread("./data/phenotypes/UKB.incident_phenotypes.tsv.gz")
incident_cases <- di$userID[di$T2D == 1 | di$CAD == 1 | di$CKD == 1 | di$NASH == 1]

#read in prs files
prs.dir <- "./data/pgs/"
prs.files <- list.files(path = prs.dir, "sscore")

for(prs in prs.files){
  prs.file <- paste0(prs.dir, prs)
  foo <- data.table::fread(prs.file)
  prs_label <- gsub("UKBPPP_", "", gsub(".sscore", "", prs))
  
  if(prs == prs.files[1]){
    prs.dat <- data.frame(userID=foo$IID)
  }
  prs.dat[[prs_label]] <- foo$SCORE1_AVG
}

#specify which PGS models I will use
#for bmi, nash, ckd, cad, and t2d, select top performing model
pgs_models <- top_pgs$PGS[top_pgs$TRAIT %in% c("diabetes_mellitus", "ischaemic_heart_diseases",
                                               "NASH_AST_ALT", "CKD_all_stages", "BMI")]
#also include gwas-significant versions for each (except NASH)
pgs_models <- c(pgs_models, "BMI_GWAS_signif", "CKD_GWAS_signif", "T2D_GWAS_signif", "CAD_GWAS_signif")
#finally, the partioned ps:
pgs_models <- c(pgs_models, "obesity", "liver_lipid", "proinsulin", "beta_cell", "lipodystrophy")

prs.dat <- prs.dat[c("userID", pgs_models)]

#build analysis file
dat <- merge(linker, covar, by.x = "pseudo_ind_id", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, p, by="userID", all=FALSE, sort=FALSE)
dat <- merge(dat, prs.dat, by="userID", all=FALSE, sort=FALSE)

#remove 2nd degree relatives
nrow(dat[dat$userID %in% rels$IID,])

dat <- dat[!dat$userID %in% rels$IID,]

#remove subjects with NASH/CAD/T2D/CKD/heart disease
t2d_terms <- colnames(dat)[grep("diabetes", colnames(dat))]

cmd_indic <- rowSums(dat[c(t2d_terms, "possible_NASH", "NASH_AST_ALT", "CKD_all_stages", "ischaemic_heart_diseases", "heart_failure", "stroke", "acute_myocardial_infarction")])
dat$CMD_FLAG <- ifelse(cmd_indic > 0, 1, 0)
dat <- dat[dat$CMD_FLAG == 0 | dat$userID %in% incident_cases,]


#regress PCs from prs
pcs <- paste0("UKBPC_", 1:10)
for(prs in pgs_models){
  f <- as.formula(paste(prs, "~", paste(pcs, collapse = "+")))
  fit <- lm(f, dat)
  dat[[prs]] <- fit$residuals
}
#read in protein information
comb.file <- "../datasets/ukb/combined_pheno_forconsortium_v1.tsv"
pdat <- data.table::fread(comb.file, data.table = FALSE)

#subset to just desired protein
proteins <- proteins[proteins$ID2 == protein,]
pdat <- pdat[c("pseudo_ind_id", proteins$ID1)]
colnames(pdat)[2] <- proteins$ID2

dat <- merge(dat, pdat, by="pseudo_ind_id", all=FALSE, sort=FALSE)
dat <- dat[!is.na(dat[[proteins$ID2]]),]
####### Done with data prep ###################

####Begin analysis #######
#subsets of cohort being tested
pops <- c("ALL","DISC", "REP", "AFR", "AMR", "EUR", "EAS", "SAS")

results <- expand.grid(PROTEIN=proteins$ID1, GENE_LABEL=proteins$PROTEIN,
                       PGS=pgs_models, SUBSET=pops,
                       N=NA, BETA=NA, SE=NA, P=NA, R2=NA,
                       BETA_BMI=NA, SE_BMI=NA, P_BMI=NA, R2_BMI=NA,
                       SNP=NA, BETA_SNP=NA, SE_SNP=NA, P_SNP=NA,
                       stringsAsFactors = FALSE)

#####tests PGS associations with and without BMI ############

for(pgs in pgs_models){
  for(pop in pops){
    
    #subset accordingly
    if(pop == "ALL"){
      foo <- dat
    }else if(pop == "DISC"){
      foo <- dat[dat$REP == FALSE,]
    }else if(pop == "REP"){
      foo <- dat[dat$REP == TRUE,]
    }else{
      keep <-  cov_all$userID[cov_all$peddy_ancestry_pred == pop & cov_all$peddy_ancestry_prob >= 0.9]
      foo <- dat[dat$userID %in% keep,]
    }
    
    #scale prs
    foo[[pgs]] <- scale(foo[[pgs]])
    
    #specify model
    pcs <- paste0("UKBPC_", 1:20)
    tbms <- paste0(proteins$PANEL, "_tbms")
    
    if(pop %in% c("ALL", "REP", "EUR")){
      b <- as.formula(paste(protein, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                                  pcs, "sample_selection"), collapse = "+")))
      f <- as.formula(paste(protein, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                                  pcs, "sample_selection", pgs), collapse = "+")))
      b_bmi <- as.formula(paste(protein, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                                      pcs,"BMI", "sample_selection"), collapse = "+")))
      f_bmi <- as.formula(paste(protein, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                                      pcs, "BMI", "sample_selection", pgs), collapse = "+")))
    }else if(pop == "DISC"){
      b <- as.formula(paste(protein, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                                  pcs), collapse = "+")))
      f <- as.formula(paste(protein, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                                  pcs, pgs), collapse = "+")))
      b_bmi <- as.formula(paste(protein, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                                      pcs, "BMI"), collapse = "+")))
      f_bmi <- as.formula(paste(protein, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                                      pcs, "BMI", pgs), collapse = "+")))
    }else{
      b <- as.formula(paste(protein, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" , tbms,
                                                  pcs), collapse = "+")))
      f <- as.formula(paste(protein, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" , tbms,
                                                  pcs, pgs), collapse = "+")))
      b_bmi <- as.formula(paste(protein, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct", tbms,
                                                      pcs, "BMI"), collapse = "+")))
      f_bmi <- as.formula(paste(protein, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct", tbms,
                                                      pcs, "BMI", pgs), collapse = "+")))
    }
    
    
    #fit models
    base <- lm(b, foo)
    full <- lm(f, foo)
    base_bmi <- lm(b_bmi, foo)
    full_bmi <- lm(f_bmi, foo)
    
    #r2
    r2 <- summary(full)$r.squared - summary(base)$r.squared
    r2_bmi <- summary(full_bmi)$r.squared - summary(base_bmi)$r.squared
    
    #parse models
    full_sum <- as.data.frame(summary(full)$coefficients)
    bmi_sum <- as.data.frame(summary(full_bmi)$coefficients)
    
    results$N[results$PGS == pgs & results$SUBSET == pop] <- length(full$residuals)
    results$BETA[results$PGS == pgs & results$SUBSET == pop] <- full_sum$Estimate[rownames(full_sum) == pgs]
    results$SE[results$PGS == pgs & results$SUBSET == pop] <- full_sum$`Std. Error`[rownames(full_sum) == pgs]
    results$P[results$PGS == pgs & results$SUBSET == pop] <- full_sum$`Pr(>|t|)`[rownames(full_sum) == pgs]
    results$R2[results$PGS == pgs & results$SUBSET == pop] <- r2
    
    results$BETA_BMI[results$PGS == pgs & results$SUBSET == pop] <- bmi_sum$Estimate[rownames(bmi_sum) == pgs]
    results$SE_BMI[results$PGS == pgs & results$SUBSET == pop] <- bmi_sum$`Std. Error`[rownames(bmi_sum) == pgs]
    results$P_BMI[results$PGS == pgs & results$SUBSET == pop] <- bmi_sum$`Pr(>|t|)`[rownames(bmi_sum) == pgs]
    results$R2_BMI[results$PGS == pgs & results$SUBSET == pop] <- r2_bmi
    
    
  }
}

#######pqtl adjustment procedure ######
plink2 <- "../software/plink2"
out_geno <- "./data/genotypes/"
pfile <- "./data/genotypes/UKBPPP.pqtls"


####create dosages file #####
snps <- unique(c(pqtl$rsID[pqtl$UKBPPP_ProteinID == proteins$PQTL_ID], "rs1260326"))

write.table(snps, paste0(out_geno, protein, "_snps"), sep = '\t', quote=FALSE,
            row.names = FALSE, col.names = FALSE)
system2(plink2, paste("--pfile", pfile, "--export A --memory 4000 --out", paste0(out_geno,protein), "--extract", paste0(out_geno, protein, "_snps"), "--silent"))

pqtl.dat <- data.table::fread(paste0(out_geno, protein, ".raw"))
snps <- colnames(pqtl.dat)[grep("rs", colnames(pqtl.dat))]


#add dosages to data 
dat <- merge(dat, pqtl.dat, by.x = "userID", by.y = "IID", all=FALSE, sort = FALSE)

pops <- c("ALL", "DISC", "REP")

for(pgs in pgs_models){
  snp_results <- expand.grid(SNP=snps, SUBSET=pops, BETA=NA, SE=NA, P=NA,
                             stringsAsFactors = FALSE)
  
  for(snp in snps){
    for(pop in pops){
      
      #subset accordingly
      if(pop == "ALL"){
        foo <- dat
      }else if(pop == "DISC"){
        foo <- dat[dat$REP == FALSE,]
      }else{
        foo <- dat[dat$REP == TRUE,]
      }
      
      #scale prs
      foo[[pgs]] <- scale(foo[[pgs]])
      
      #specify model
      pcs <- paste0("UKBPC_", 1:20)
      tbms <- paste0(proteins$PANEL, "_tbms")
      
      if(pop %in% c("ALL", "REP")){
        f <- as.formula(paste(protein, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                                    pcs, "sample_selection", snp, pgs), collapse = "+")))
        
      }else{
        f <- as.formula(paste(protein, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                                    pcs, snp, pgs), collapse = "+")))
      }
      
      
      #fit models
      full <- lm(f, foo)
      
      
      #parse models
      full_sum <- as.data.frame(summary(full)$coefficients)
      
      snp_results$BETA[snp_results$SNP == snp & snp_results$SUBSET == pop] <- full_sum$Estimate[rownames(full_sum) == pgs]
      snp_results$SE[snp_results$SNP == snp & snp_results$SUBSET == pop] <- full_sum$`Std. Error`[rownames(full_sum) == pgs]
      snp_results$P[snp_results$SNP == snp & snp_results$SUBSET == pop] <- full_sum$`Pr(>|t|)`[rownames(full_sum) == pgs]
      
    }
  }
  
  ##just keep results for snp with largest effect on summary stats##
  snp <- snp_results$SNP[snp_results$SUBSET == "ALL" & snp_results$P == max(snp_results$P[snp_results$SUBSET == "ALL"])][1]
  
  snp_results <- snp_results[snp_results$SNP == snp,]
  
  #add to main data frame
  results$SNP[results$SUBSET %in% c("ALL", "DISC", "REP") & results$PGS == pgs]  <- snp
  results$BETA_SNP[results$SUBSET %in% c("ALL", "DISC", "REP") & results$PGS == pgs] <- snp_results$BETA
  results$SE_SNP[results$SUBSET %in% c("ALL", "DISC", "REP") & results$PGS == pgs] <- snp_results$SE
  results$P_SNP[results$SUBSET %in% c("ALL", "DISC", "REP") & results$PGS == pgs] <- snp_results$P
  
}


##clean up ##
system2("rm", paste0(out_geno,protein, "*"))

#write results
out.file <- paste0("./results/03_pgs_npx/",
                   protein, "_incid_incl_assoc.tsv")

data.table::fwrite(results, out.file, 
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE,
                   na = "NA")
