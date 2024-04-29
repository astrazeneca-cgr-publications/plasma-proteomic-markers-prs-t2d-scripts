########script for running mediation in UKB ####################
library(medflex)

####prepare analyses###
argsv <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(argsv[1])
#####read in supporing files ######
pqtl <- data.table::fread("./data/pqtl_data/cond_indep_pqtl_list.csv.gz") #from UKB-PPP coniditional analysis
top_pgs <- data.table::fread("./results/top_PGS_per_trait.tsv")
proteins <- data.table::fread("./results/olink_protein_description.tsv")

#subset protens file
proteins <- proteins[index,]


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
pi <- data.table::fread("./data/phenotypes/UKB.incident_phenotypes.tsv.gz")
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

#just keep BMI from primary phenotype file and add incident cases
p <- as.data.frame(p)
p <- p[c("userID", "BMI")]
p <- merge(p, pi, by="userID", all=FALSE, sort=FALSE)

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
#read in protein information
comb.file <- "../datasets/ukb/combined_pheno_forconsortium_v1.tsv"
pdat <- data.table::fread(comb.file, data.table = FALSE)

####### Done with data prep ###################

####Begin analysis #######
#subsets of cohort being tested
#to simplify things, just to EUR and ALL
pops <- c("ALL", "EUR")
traits <- c("T2D", "CAD", "CKD", "NASH", "OBESE")
results <- expand.grid(PROTEIN=proteins$ID1, GENE_LABEL=proteins$PROTEIN,
                       PGS=pgs_models, TRAIT=traits, SUBSET=pops,
                       N=NA, N_CASES=NA, BETA_PGS=NA, SE_PGS=NA, P_PGS=NA,
                       BETA_PGS_TRAIT=NA, SE_PGS_TRAIT=NA, P_PGS_TRAIT=NA,
                       BETA_PROTEIN_TRAIT=NA, SE_PROTEIN_TRAIT=NA, P_PROTEIN_TRAIT=NA,
                       EST_DIRECT=NA, SE_DIRECT=NA, P_DIRECT=NA,
                       EST_INDIRECT=NA, SE_INDIRECT=NA, P_INDIRECT=NA,
                       EST_TOTAL=NA, SE_TOTAL=NA, P_TOTAL=NA,
                       stringsAsFactors = FALSE)

#####part1:PGS-trait assocations #######
for(trait in traits){
  for(pgs in pgs_models){
    for(pop in pops){
      if(pop == "ALL"){
        foo <- dat
      }else{
        foo <- dat[dat$userID %in% cov_all$userID[cov_all$peddy_ancestry_pred == "EUR" & cov_all$peddy_ancestry_prob >= 0.9],]
      }
      
      #scale pgs
      foo[[pgs]] <- scale(foo[[pgs]])
      f <- as.formula(paste(trait, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "ukb_centre_fct" ,"array_fct",
                                                pcs, "sample_selection", pgs), collapse = "+")))
      fit <- glm(f, foo, family="binomial")
      #parse results
      fit_sum <- as.data.frame(summary(fit)$coefficients)
      
      #save results
      results$N[results$PGS == pgs & results$TRAIT == trait & results$SUBSET == pop] <- 
      results$BETA_PGS_TRAIT[results$PGS == pgs & results$TRAIT == trait & results$SUBSET == pop] <- fit_sum$Estimate[rownames(fit_sum) == pgs]
      results$SE_PGS_TRAIT[results$PGS == pgs & results$TRAIT == trait & results$SUBSET == pop] <- fit_sum$`Std. Error`[rownames(fit_sum) == pgs]
      results$P_PGS_TRAIT[results$PGS == pgs & results$TRAIT == trait & results$SUBSET == pop] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == pgs]
    }
    
  }
  
}


#test each pgs for association with protein
for(pgs in pgs_models){
    
    for(pop in pops){
      
      #build model
      pcs <- paste0("UKBPC_", 1:20)
      tbms <- paste0(proteins$PANEL, "_tbms")
      #subset data
      foo <- pdat[c("pseudo_ind_id", proteins$ID1)]
      colnames(foo)[2] <- proteins$ID2  
      foo <- merge(foo, dat, by="pseudo_ind_id", all=FALSE, sort = FALSE)
      if(pop == "EUR") foo <- foo[foo$userID %in% cov_all$userID[cov_all$peddy_ancestry_pred == "EUR" & cov_all$peddy_ancestry_prob >= 0.9],]
      
      #scale prs
      foo[[pgs]] <- scale(foo[[pgs]])
      
      f <- as.formula(paste(proteins$ID2, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                                  pcs, "sample_selection", pgs), collapse = "+")))
      
      #fit
      fit <- lm(f, data=foo)
      
      #parse results
      fit_sum <- as.data.frame(summary(fit)$coefficients)
      
      results$N[results$PGS == pgs & results$PROTEIN == proteins$ID1 & results$SUBSET == pop] <- length(fit$residuals)
      results$BETA_PGS[results$PGS == pgs & results$PROTEIN == proteins$ID1  & results$SUBSET == pop] <- fit_sum$Estimate[rownames(fit_sum) == pgs]
      results$SE_PGS[results$PGS == pgs & results$PROTEIN == proteins$ID1  & results$SUBSET == pop] <- fit_sum$`Std. Error`[rownames(fit_sum) == pgs]
      results$P_PGS[results$PGS == pgs & results$PROTEIN == proteins$ID1  & results$SUBSET == pop] <- fit_sum$`Pr(>|t|)`[rownames(fit_sum) == pgs]
    }
    
}


#test protein for association for trait
for(trait in traits){
  for(pop in pops){
    

    #subset data
    foo <- pdat[c("pseudo_ind_id", proteins$ID1)]
    colnames(foo)[2] <- proteins$ID2  
    foo <- merge(foo, dat, by="pseudo_ind_id", all=FALSE, sort = FALSE)
    if(pop == "EUR") foo <- foo[foo$userID %in% cov_all$userID[cov_all$peddy_ancestry_pred == "EUR" & cov_all$peddy_ancestry_prob >= 0.9],]
    
    #build model
    pcs <- paste0("UKBPC_", 1:20)
    tbms <- paste0(proteins$PANEL, "_tbms")
    f <- as.formula(paste(trait, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                                     pcs, "sample_selection", proteins$ID2), collapse = "+")))
    
    #fit
    fit <- glm(f, data=foo, family = "binomial")
    
    #parse results
    fit_sum <- as.data.frame(summary(fit)$coefficients)
    
    results$BETA_PROTEIN_TRAIT[results$TRAIT == trait & results$PROTEIN == proteins$ID1  & results$SUBSET == pop] <- fit_sum$Estimate[rownames(fit_sum) == proteins$ID2]
    results$SE_PROTEIN_TRAIT[results$TRAIT == trait & results$PROTEIN == proteins$ID1  & results$SUBSET == pop] <- fit_sum$`Std. Error`[rownames(fit_sum) == proteins$ID2]
    results$P_PROTEIN_TRAIT[results$TRAIT == trait & results$PROTEIN == proteins$ID1  & results$SUBSET == pop] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == proteins$ID2]
  }
}

#now run mediation#
for(i in 1:nrow(results)){
  if(results$P_PGS[i] < 0.05 & results$P_PGS_TRAIT[i] < 0.05 & results$P_PROTEIN_TRAIT[i] < 0.05){
    print(paste("running mediation for", results$PGS[i], "and", results$PROTEIN[i], "and", results$TRAIT[i], "in", results$SUBSET[i]))
    
    pop <- results$SUBSET[i]
    trait <- results$TRAIT[i]
    pgs <- results$PGS[i]
    protein <- proteins$ID2[i]
 
  
    #subset data
    foo <- pdat[c("pseudo_ind_id", proteins$ID1)]
    colnames(foo)[2] <- proteins$ID2  
    foo <- merge(foo, dat, by="pseudo_ind_id", all=FALSE, sort = FALSE)
    if(pop == "EUR") foo <- foo[foo$userID %in% cov_all$userID[cov_all$peddy_ancestry_pred == "EUR" & cov_all$peddy_ancestry_prob >= 0.9],]
    
    #build model
    pcs <- paste0("UKBPC_", 1:20)
    tbms <- paste0(proteins$PANEL, "_tbms")
    #build model
    pcs <- paste0("UKBPC_", 1:20)
  
    #need to scale age for mediation
    foo$age <- as.numeric(scale(foo$age))
    
    #scale PGS
    foo[[pgs]] <- as.numeric(scale(foo[[pgs]]))
    
    #scale PCs
    pcs <- paste0("UKBPC_", 1:20)
    for(pc in pcs){
      foo[[pc]] <- as.numeric(scale(foo[[pc]]))
    }
    
    #build model
    tbms <- paste0(proteins$PANEL, "_tbms")
    f <- as.formula(paste(trait, "~", paste(c("age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                              pcs, "sample_selection", proteins$ID2), collapse = "+")))
    
    #change pgs column
    colnames(foo)[colnames(foo) == pgs] <- "PGS"
    
    #change protein column
    colnames(foo)[colnames(foo) == proteins$ID2] <- "PROTEIN"
    
    #remove missing subjects
    foo <- foo[complete.cases(foo[[trait]]) & complete.cases(foo$PGS) & complete.cases(foo$PROTEIN),]
    
    
    #run mediation analysis
    f <- as.formula(paste("EVENT~", paste(c("PGS", "PROTEIN", "SEX", "AGE", pcs), collapse = "+")))
    tbms <- paste0(proteins$PANEL, "_tbms")
    f <- as.formula(paste(trait, "~", paste(c("PGS", "PROTEIN", "age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                              pcs, "sample_selection"), collapse = "+")))
    expData <- neImpute(f, family = binomial("logit"), data=foo)
    
    f <- as.formula(paste(trait, "~", paste(c(paste0("PGS", 0), paste0("PGS", 1), "age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                              pcs, "sample_selection"), collapse = "+")))
    nMod <- neModel(f, family = binomial("logit"), expData=expData, se="robust")
    
    lht <- as.data.frame(coef(summary(neLht(nMod, linfct = c("PGS0 = 0",  "PGS1 = 0", "PGS0 + PGS1 = 0")))))
    
    results$EST_DIRECT[i] <- lht$Estimate[1]
    results$SE_DIRECT[i] <- lht$`Std. Error`[1]
    results$P_DIRECT[i] <- lht$`Pr(>|z|)`[1]
    
    results$EST_INDIRECT[i] <- lht$Estimate[2]
    results$SE_INDIRECT[i] <- lht$`Std. Error`[2]
    results$P_INDIRECT[i] <- lht$`Pr(>|z|)`[2]
    
    results$EST_TOTAL[i] <- lht$Estimate[3]
    results$SE_TOTAL[i] <- lht$`Std. Error`[3]
    results$P_TOTAL[i] <- lht$`Pr(>|z|)`[3]
    
    
  }else{
    next
  }
}

#write out
out.file <- paste0("./results/03_mediation/", proteins$ID2, "_", "mediation_results_new_defintions.tsv")

data.table::fwrite(results, out.file, sep='\t', quote=FALSE, row.names = FALSE, 
                   col.names = TRUE, na="NA")



######repeat the above with BMI adjustment ############

results <- expand.grid(PROTEIN=proteins$ID1, GENE_LABEL=proteins$PROTEIN,
                       PGS=pgs_models, TRAIT=traits, SUBSET=pops,
                       N=NA, BETA_PGS=NA, SE_PGS=NA, P_PGS=NA,
                       BETA_PGS_TRAIT=NA, SE_PGS_TRAIT=NA, P_PGS_TRAIT=NA,
                       BETA_PROTEIN_TRAIT=NA, SE_PROTEIN_TRAIT=NA, P_PROTEIN_TRAIT=NA,
                       EST_DIRECT=NA, SE_DIRECT=NA, P_DIRECT=NA,
                       EST_INDIRECT=NA, SE_INDIRECT=NA, P_INDIRECT=NA,
                       EST_TOTAL=NA, SE_TOTAL=NA, P_TOTAL=NA,
                       stringsAsFactors = FALSE)

#test each pgs for association with trait
for(trait in traits){
  for(pgs in pgs_models){
    for(pop in pops){
      if(pop == "ALL"){
        foo <- dat
      }else{
        foo <- dat[dat$userID %in% cov_all$userID[cov_all$peddy_ancestry_pred == "EUR" & cov_all$peddy_ancestry_prob >= 0.9],]
      }
      
      #scale pgs
      foo[[pgs]] <- scale(foo[[pgs]])
      f <- as.formula(paste(trait, "~", paste(c("BMI", "age", "age2", "sex", "age*sex", "age2*sex", "ukb_centre_fct" ,"array_fct",
                                                pcs, "sample_selection", pgs), collapse = "+")))
      fit <- glm(f, foo, family="binomial")
      #parse results
      fit_sum <- as.data.frame(summary(fit)$coefficients)
      
      #save results
      results$BETA_PGS_TRAIT[results$PGS == pgs & results$TRAIT == trait & results$SUBSET == pop] <- fit_sum$Estimate[rownames(fit_sum) == pgs]
      results$SE_PGS_TRAIT[results$PGS == pgs & results$TRAIT == trait & results$SUBSET == pop] <- fit_sum$`Std. Error`[rownames(fit_sum) == pgs]
      results$P_PGS_TRAIT[results$PGS == pgs & results$TRAIT == trait & results$SUBSET == pop] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == pgs]
    }
    
  }
  
}


#test each pgs for association with protein
for(pgs in pgs_models){
  
  for(pop in pops){
    
    #build model
    pcs <- paste0("UKBPC_", 1:20)
    tbms <- paste0(proteins$PANEL, "_tbms")
    #subset data
    foo <- pdat[c("pseudo_ind_id", proteins$ID1)]
    colnames(foo)[2] <- proteins$ID2  
    foo <- merge(foo, dat, by="pseudo_ind_id", all=FALSE, sort = FALSE)
    if(pop == "EUR") foo <- foo[foo$userID %in% cov_all$userID[cov_all$peddy_ancestry_pred == "EUR" & cov_all$peddy_ancestry_prob >= 0.9],]
    
    #scale prs
    foo[[pgs]] <- scale(foo[[pgs]])
    
    f <- as.formula(paste(proteins$ID2, "~", paste(c("BMI", "age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                                     pcs, "sample_selection", pgs), collapse = "+")))
    
    #fit
    fit <- lm(f, data=foo)
    
    #parse results
    fit_sum <- as.data.frame(summary(fit)$coefficients)
    
    results$N[results$PGS == pgs & results$PROTEIN == proteins$ID1 & results$SUBSET == pop] <- length(fit$residuals)
    results$BETA_PGS[results$PGS == pgs & results$PROTEIN == proteins$ID1  & results$SUBSET == pop] <- fit_sum$Estimate[rownames(fit_sum) == pgs]
    results$SE_PGS[results$PGS == pgs & results$PROTEIN == proteins$ID1  & results$SUBSET == pop] <- fit_sum$`Std. Error`[rownames(fit_sum) == pgs]
    results$P_PGS[results$PGS == pgs & results$PROTEIN == proteins$ID1  & results$SUBSET == pop] <- fit_sum$`Pr(>|t|)`[rownames(fit_sum) == pgs]
  }
  
}


#test protein for association for trait
for(trait in traits){
  for(pop in pops){
    
    
    #subset data
    foo <- pdat[c("pseudo_ind_id", proteins$ID1)]
    colnames(foo)[2] <- proteins$ID2  
    foo <- merge(foo, dat, by="pseudo_ind_id", all=FALSE, sort = FALSE)
    if(pop == "EUR") foo <- foo[foo$userID %in% cov_all$userID[cov_all$peddy_ancestry_pred == "EUR" & cov_all$peddy_ancestry_prob >= 0.9],]
    
    #build model
    pcs <- paste0("UKBPC_", 1:20)
    tbms <- paste0(proteins$PANEL, "_tbms")
    f <- as.formula(paste(trait, "~", paste(c("BMI", "age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                              pcs, "sample_selection", proteins$ID2), collapse = "+")))
    
    #fit
    fit <- glm(f, data=foo, family = "binomial")
    
    #parse results
    fit_sum <- as.data.frame(summary(fit)$coefficients)
    
    results$BETA_PROTEIN_TRAIT[results$TRAIT == trait & results$PROTEIN == proteins$ID1  & results$SUBSET == pop] <- fit_sum$Estimate[rownames(fit_sum) == proteins$ID2]
    results$SE_PROTEIN_TRAIT[results$TRAIT == trait & results$PROTEIN == proteins$ID1  & results$SUBSET == pop] <- fit_sum$`Std. Error`[rownames(fit_sum) == proteins$ID2]
    results$P_PROTEIN_TRAIT[results$TRAIT == trait & results$PROTEIN == proteins$ID1  & results$SUBSET == pop] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == proteins$ID2]
  }
}

#now run mediation#
for(i in 1:nrow(results)){
  if(results$P_PGS[i] < 0.05 & results$P_PGS_TRAIT[i] < 0.05 & results$P_PROTEIN_TRAIT[i] < 0.05){
    print(paste("running mediation for", results$PGS[i], "and", results$PROTEIN[i], "and", results$TRAIT[i], "in", results$SUBSET[i]))
    
    pop <- results$SUBSET[i]
    trait <- results$TRAIT[i]
    pgs <- results$PGS[i]
    protein <- proteins$ID2[i]
    
    
    #subset data
    foo <- pdat[c("pseudo_ind_id", proteins$ID1)]
    colnames(foo)[2] <- proteins$ID2  
    foo <- merge(foo, dat, by="pseudo_ind_id", all=FALSE, sort = FALSE)
    if(pop == "EUR") foo <- foo[foo$userID %in% cov_all$userID[cov_all$peddy_ancestry_pred == "EUR" & cov_all$peddy_ancestry_prob >= 0.9],]
    
    #build model
    pcs <- paste0("UKBPC_", 1:20)
    tbms <- paste0(proteins$PANEL, "_tbms")
    #build model
    pcs <- paste0("UKBPC_", 1:20)
    
    #need to scale age for mediation
    foo$age <- as.numeric(scale(foo$age))
    
    #scale PGS
    foo[[pgs]] <- as.numeric(scale(foo[[pgs]]))
    
    #scale PCs
    pcs <- paste0("UKBPC_", 1:20)
    for(pc in pcs){
      foo[[pc]] <- as.numeric(scale(foo[[pc]]))
    }
    
    #change pgs column
    colnames(foo)[colnames(foo) == pgs] <- "PGS"
    
    #change protein column
    colnames(foo)[colnames(foo) == proteins$ID2] <- "PROTEIN"
    
    #remove missing subjects
    foo <- foo[complete.cases(foo[[trait]]) & complete.cases(foo$PGS) & complete.cases(foo$PROTEIN) & complete.cases(foo$BMI),]
    
    
    #run mediation analysis
    tbms <- paste0(proteins$PANEL, "_tbms")
    f <- as.formula(paste(trait, "~", paste(c("PGS", "PROTEIN", "BMI", "age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                              pcs, "sample_selection"), collapse = "+")))
    expData <- neImpute(f, family = binomial("logit"), data=foo)
    
    f <- as.formula(paste(trait, "~", paste(c(paste0("PGS", 0), paste0("PGS", 1), "BMI", "age", "age2", "sex", "age*sex", "age2*sex", "Batch", "ukb_centre_fct" ,"array_fct", tbms,
                                              pcs, "sample_selection"), collapse = "+")))
    nMod <- neModel(f, family = binomial("logit"), expData=expData, se="robust")
    
    lht <- as.data.frame(coef(summary(neLht(nMod, linfct = c("PGS0 = 0",  "PGS1 = 0", "PGS0 + PGS1 = 0")))))
    
    results$EST_DIRECT[i] <- lht$Estimate[1]
    results$SE_DIRECT[i] <- lht$`Std. Error`[1]
    results$P_DIRECT[i] <- lht$`Pr(>|z|)`[1]
    
    results$EST_INDIRECT[i] <- lht$Estimate[2]
    results$SE_INDIRECT[i] <- lht$`Std. Error`[2]
    results$P_INDIRECT[i] <- lht$`Pr(>|z|)`[2]
    
    results$EST_TOTAL[i] <- lht$Estimate[3]
    results$SE_TOTAL[i] <- lht$`Std. Error`[3]
    results$P_TOTAL[i] <- lht$`Pr(>|z|)`[3]
    
    
  }else{
    next
  }
}

#write out
out.file <- paste0("./results/03_mediation/", proteins$ID2, "_", "BMIadj_mediation_results_new_definitions.tsv")

data.table::fwrite(results, out.file, sep='\t', quote=FALSE, row.names = FALSE, 
                   col.names = TRUE, na="NA")


