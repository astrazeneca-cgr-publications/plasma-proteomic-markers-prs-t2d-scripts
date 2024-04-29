#######script for testing PGS associations with traits ##########
#### this evaluates whether score association patterns behave as expected in UKB-PPP ####

#read in linker file to ID UKB-PPP participants
linker.file <- "../datasets/ukb/olink_sample_map_3k_ukbsamples_forconsort_v2_linked.tsv"

linker <- data.table::fread(linker.file, data.table = FALSE)
linker <- linker[c("userID", "pseudo_ind_id")]
linker <- linker[!duplicated(linker),]

#read in UKB-PPP covariate file
covar.file <- "../datasets/ukb/combined_covars_forconsortium_v1.tsv"
covar <- data.table::fread(covar.file, data.table = FALSE)

#set categorial variables as factors
covar$Batch <- as.factor(covar$Batch)
covar$array_fct <- as.factor(covar$array_fct)
covar$ukb_centre_fct <- as.factor(covar$ukb_centre_fct)
covar$sample_selection <- as.factor(covar$sample_selection)
covar$sex <- as.factor(covar$sex)

#read in ukb-wide phenotypes and covariate files
p <- data.table::fread("./data/phenotypes/UKB.phenotypes.tsv.gz")
cov_all <- data.table::fread("./data/phenotypes/UKB.covariates.tsv.gz")

#read in relatives file for seclusions 
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

#build analysis file
dat <- merge(linker, covar, by.x = "pseudo_ind_id", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, p, by="userID", all=FALSE, sort=FALSE)
dat <- merge(dat, prs.dat, by="userID", all=FALSE, sort=FALSE)

#remove 2nd degree relatives
dat <- dat[!dat$userID %in% rels$IID,]

prs_models <- colnames(prs.dat)[2:ncol(prs.dat)]
traits <- colnames(p)[2:ncol(p)]
pops <- c("ALL", "AFR", "AMR", "EUR", "EAS", "SAS")

#first, regress PCs from prs
pcs <- paste0("UKBPC_", 1:10)
for(prs in prs_models){
  f <- as.formula(paste(prs, "~", paste(pcs, collapse = "+")))
  fit <- lm(f, dat)
  dat[[prs]] <- fit$residuals
}

#build results data frame
results <- expand.grid(PGS=prs_models, TRAIT=traits, POP=pops, 
                       N=NA, N_CASES=NA, BETA=NA, SE=NA, P=NA,
                       L_95=NA, U_95=NA,
                       stringsAsFactors = FALSE)


#run associations with overall traits, both all subjects and stratified by ancestry
for(trait in traits){
  
  for(prs in prs_models){
    for(pop in pops){
      if(pop == "ALL"){
        foo <- dat
      }else{
        keep <-  cov_all$userID[cov_all$peddy_ancestry_pred == pop & cov_all$peddy_ancestry_prob >= 0.9]
        foo <- dat[dat$userID %in% keep,]
      }
      
      foo[[prs]] <- scale(foo[[prs]])
      
      #drop array from non-European analyses as they have 4 or less subjects on array 2
      if(pop %in% c("AFR", "EAS", "AMR", "SAS")){
        f <- as.formula(paste(trait, "~", paste(c("age", "sex", "age2", "ukb_centre_fct", pcs, prs), collapse = "+")))
      }else{
        f <- as.formula(paste(trait, "~", paste(c("age", "sex", "age2", "array_fct", "ukb_centre_fct", pcs, prs), collapse = "+")))
      }
      
      if(trait == "BMI"){
        fit <- lm(f, data = foo)
      }else{
        fit <- glm(f, data = foo, family = "binomial")
      }
      results$N[results$PGS == prs & results$POP == pop & results$TRAIT == trait] <- nrow(foo)
      results$N_CASES[results$PGS == prs & results$POP == pop & results$TRAIT == trait] <- nrow(foo[foo[[trait]] == 1,])
      results$BETA[results$PGS == prs & results$POP == pop & results$TRAIT == trait] <- summary(fit)$coefficients[nrow(summary(fit)$coefficients),1]
      results$SE[results$PGS == prs & results$POP == pop & results$TRAIT == trait] <- summary(fit)$coefficients[nrow(summary(fit)$coefficients),2]
      results$P[results$PGS == prs & results$POP == pop & results$TRAIT == trait] <- summary(fit)$coefficients[nrow(summary(fit)$coefficients),4]
      
      #add confidence intervals for linear models
      if(trait == "BMI"){
        results$L_95[results$PGS == prs & results$POP == pop & results$TRAIT == trait] <- confint(fit)[nrow(confint(fit)), 1]
        results$U_95[results$PGS == prs & results$POP == pop & results$TRAIT == trait] <- confint(fit)[nrow(confint(fit)), 2]
      }else{
        beta <- summary(fit)$coefficients[nrow(summary(fit)$coefficients),1]
        se <- summary(fit)$coefficients[nrow(summary(fit)$coefficients),2]
        results$L_95[results$PGS == prs & results$POP == pop & results$TRAIT == trait] <- beta - 1.96*se
        results$U_95[results$PGS == prs & results$POP == pop & results$TRAIT == trait] <- beta + 1.96*se
      }
 
    }
  }
}

#give top performing PGS for trait
top <- data.frame(TRAIT=traits, PGS=NA, PVAL=NA)
for(trait in traits){
  pval <- min(results$P[results$TRAIT == trait & results$POP == "ALL"])
  pgs <- results$PGS[results$TRAIT == trait & results$POP == "ALL" & results$P == pval]
  
  #split ties with beta
  if(length(pgs) >1){
    beta <- max(results$BETA[results$TRAIT == trait & results$POP == "ALL" & results$PGS %in% pgs])
    pgs <- results$PGS[results$TRAIT == trait & results$POP == "ALL" & results$P == pval & results$BETA == beta]
  }

  top$PGS[top$TRAIT == trait] <- pgs 
  top$PVAL[top$TRAIT == trait] <- pval 
}



####write out results####
data.table::fwrite(results, "./results/PGS_trait_assoc.tsv", 
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)
data.table::fwrite(top, "./results/top_PGS_per_trait.tsv", 
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)
