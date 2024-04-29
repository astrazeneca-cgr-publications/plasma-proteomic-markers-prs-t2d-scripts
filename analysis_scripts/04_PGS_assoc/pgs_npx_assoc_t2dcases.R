####script for performing PRS-npx associations in T2D cases ########
argsv <- commandArgs(trailingOnly = TRUE)
protein <- gsub(":", "_", argsv[1])
protein <- gsub("-", "_", protein)

#####read in supporing files ######
pqtl <- data.table::fread("./data/pqtl_data/cond_indep_pqtl_list.csv.gz")
top_pgs <- data.table::fread("./results/top_PGS_per_trait.tsv")
proteins <- data.table::fread("./results/olink_protein_description.tsv")

####read in files to build analysis file #####
#read in linker file
linker.file <- "../datasets/olink_sample_map_3k_ukbsamples_forconsort_v2_linked.tsv"

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
dat <- dat[!dat$userID %in% rels$IID,]

#regress PCs from prs
pcs <- paste0("UKBPC_", 1:10)
for(prs in pgs_models){
  f <- as.formula(paste(prs, "~", paste(pcs, collapse = "+")))
  fit <- lm(f, dat)
  dat[[prs]] <- fit$residuals
}

#only keep subjects with non-insulin dep diabetes
dat <- dat[dat$non_insulin_dependent_diabetes_mellitus == 1 & dat$insulin_dependent_diabetes_mellitus == 0,]

#read in protein information
comb.file <- "../datasets/ukb/combined_pheno_forconsortium_v1.tsv"
pdat <- data.table::fread(comb.file, data.table = FALSE)


####Begin analysis #######
#subsets of cohort being tested
pops <- c("ALL", "DISC")

#prevalent and all
case_types <- c("PREVALENT", "INCIDENT", "ALL")

#prepare for iteration without changing code too much
d <- dat
protein_list <- proteins

results <- expand.grid(PROTEIN=protein_list$ID1, GENE_LABEL=NA,
                       CASE_STATUS=case_types,
                       PGS=pgs_models, SUBSET=pops,
                       N=NA, BETA=NA, SE=NA, P=NA, R2=NA,
                       BETA_BMI=NA, SE_BMI=NA, P_BMI=NA, R2_BMI=NA,
                       stringsAsFactors = FALSE)


#subset to just desired protein
#for(protein in protein_list$ID2){
  print(protein)
  proteins <- protein_list[protein_list$ID2 == protein,]
  
  #add gene label 
  results$GENE_LABEL[results$PROTEIN == proteins$ID1] <- proteins$PROTEIN
  
  p <- pdat[c("pseudo_ind_id", proteins$ID1)]
  colnames(p)[2] <- proteins$ID2
  
  dat <- merge(d, p, by="pseudo_ind_id", all=FALSE, sort=FALSE)
  dat <- dat[!is.na(dat[[proteins$ID2]]),]
  
  #####tests PGS associations with and without BMI ############
  
  for(pgs in pgs_models){
    for(pop in pops){
      for(case_type in case_types){
        #subset by ancestry
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
        
        #subset by case status
        if(case_type == "PREVALENT"){
          foo <- foo[is.na(foo$T2D_INCIDENT),]
        }else if(case_type == "INCIDENT"){
          foo <- foo[!is.na(foo$T2D_INCIDENT) & foo$T2D_INCIDENT == 1,]
        }else{
          foo <- dat
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
        
        results$N[results$PGS == pgs & results$SUBSET == pop & results$PROTEIN == proteins$ID1 & results$CASE_STATUS == case_type] <- length(full$residuals)
        results$BETA[results$PGS == pgs & results$SUBSET == pop & results$PROTEIN == proteins$ID1 & results$CASE_STATUS == case_type] <- full_sum$Estimate[rownames(full_sum) == pgs]
        results$SE[results$PGS == pgs & results$SUBSET == pop & results$PROTEIN == proteins$ID1 & results$CASE_STATUS == case_type] <- full_sum$`Std. Error`[rownames(full_sum) == pgs]
        results$P[results$PGS == pgs & results$SUBSET == pop & results$PROTEIN == proteins$ID1 & results$CASE_STATUS == case_type] <- full_sum$`Pr(>|t|)`[rownames(full_sum) == pgs]
        results$R2[results$PGS == pgs & results$SUBSET == pop & results$PROTEIN == proteins$ID1 & results$CASE_STATUS == case_type] <- r2
        
        results$BETA_BMI[results$PGS == pgs & results$SUBSET == pop & results$PROTEIN == proteins$ID1 & results$CASE_STATUS == case_type] <- bmi_sum$Estimate[rownames(bmi_sum) == pgs]
        results$SE_BMI[results$PGS == pgs & results$SUBSET == pop & results$PROTEIN == proteins$ID1 & results$CASE_STATUS == case_type] <- bmi_sum$`Std. Error`[rownames(bmi_sum) == pgs]
        results$P_BMI[results$PGS == pgs & results$SUBSET == pop & results$PROTEIN == proteins$ID1 & results$CASE_STATUS == case_type] <- bmi_sum$`Pr(>|t|)`[rownames(bmi_sum) == pgs]
        results$R2_BMI[results$PGS == pgs & results$SUBSET == pop & results$PROTEIN == proteins$ID1 & results$CASE_STATUS == case_type] <- r2_bmi
      }
    }
  }
#}

#remove other proteins 
results <- results[!is.na(results$BETA),]

#write results
out.file <- paste0("./results/03_pgs_npx/",
                   protein, "_T2D_cases_assoc.tsv")

data.table::fwrite(results, out.file, 
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE,
                   na = "NA")
