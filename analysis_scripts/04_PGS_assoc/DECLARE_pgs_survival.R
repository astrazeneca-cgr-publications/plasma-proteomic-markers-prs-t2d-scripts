#script for running cox proportional hazards regression in DECLARE using PGS as exposure
library(survival)


#read in covariate file with time to event data
dat <- data.table::fread("./data/rct_data/clinical_outcomes.gz")
dat <- dat[dat$STUDY == "DECLARE",]

#read in risk factors file
covar <- data.table::fread("./data/rct_data/DECLARE_risk_factors.gz")

#read in PCA
pca <- data.table::fread("./data/rct_data/DECLARE-PCA.eigenvec")

#read in ancestry definitions
anc <- data.table::fread("./data/rct_data/DECLARE_ancestry_InferredAncestry.txt")
anc$PC1 <- NULL
anc$PC2 <- NULL 

#build analysis file
dat <- merge(dat, pca, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, anc, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, covar, by=c("EID", "ID"), all=FALSE, sort=FALSE)


#add age2
dat$AGE2 <- dat$AGE^2

#read in PGS
top_pgs <- data.table::fread("./results/top_PGS_per_trait.tsv")
prs.dir <- "./data/rct_data/scores/"
prs.files <- list.files(prs.dir, "sscore")
prs.files <- prs.files[grep("DECLARE", prs.files)]

for(prs in prs.files){
  prs.file <- paste0(prs.dir, prs)
  foo <- data.table::fread(prs.file)
  colnames(foo)[1] <- "IID"
  prs_label <- gsub("DECLARE_", "", gsub(".sscore", "", prs))
  
  if(prs == prs.files[1]){
    prs.dat <- data.frame(ID=foo$IID)
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

prs.dat <- prs.dat[c("ID", pgs_models)]

#add to data frame
dat <- merge(dat, prs.dat, by="ID", all=FALSE, sort=FALSE)

#regress PCs from prs
pcs <- paste0("PC", 1:10)
for(prs in pgs_models){
  f <- as.formula(paste(prs, "~", paste(pcs, collapse = "+")))
  fit <- lm(f, dat)
  dat[[prs]] <- fit$residuals
}

#read in quantiative traits
q <- data.table::fread("./data/rct_data/DECLARE_quant_traits.gz")
q <- q[,-c(2:5)]
dat <- merge(dat, q, by=c("ID", "EID"), all=FALSE, sort=FALSE)

#run analyses
results <- expand.grid(PGS=pgs_models, 
                       OUTCOME=unique(dat$OUTCOME),
                       SUBSET=c("ALL", "EUR"),
                       N=NA, N_EVENT=NA,
                       BETA_LR=NA, SE_LR=NA, P_LR=NA,
                       BETA=NA, HR=NA, SE=NA, P=NA,
                       stringsAsFactors = FALSE)

for(trait in unique(dat$OUTCOME)){
  for(pgs in pgs_models){
    for(pop in c("ALL", "EUR")){
      if(pop == "ALL"){
        foo <- as.data.frame(dat[dat$OUTCOME == trait,])
      }else{
        foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9 & dat$OUTCOME == trait,])
      }
      
      #exclude subjects who had insulin at baseline
      if(trait == "INSULIN") foo <- foo[foo$INSULIN != 1,]
      if(trait == "RENAL") foo <- foo[foo$EGFR_BL > 60,]
      #exclude subjects with 
      #just include necessary columns
      pcs <- paste0("PC", 1:10)
      foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", pcs, pgs)]
      colnames(foo)[ncol(foo)] <- "PGS"
      
      #scale pgs
      foo$PGS <- scale(foo$PGS)
      
      #fit survival object
      sobj <- Surv(foo$TIME, foo$EVENT, type="right")
      foo$SURV <- sobj
      
      #specify formula
      f <- as.formula(paste("EVENT~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "PGS"), collapse = "+")))
      f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "PGS"), collapse = "+")))
      
      #fit logistic model
      fit_glm <- glm(f, foo, family = "binomial")
      #fit cox model
      fit <- coxph(f_cox, foo)
      
      #parse results
      fit_glm_sum <- as.data.frame(summary(fit_glm)$coefficients)
      fit_sum <- as.data.frame(summary(fit)$coefficients)
      
      #summary
      results$N[results$SUBSET == pop & results$PGS == pgs & results$OUTCOME == trait] <- nrow(foo)
      results$N_EVENT[results$SUBSET == pop & results$PGS == pgs ] <- nrow(foo[foo$EVENT == 1,])
      
      #glm results 
      results$BETA_LR[results$SUBSET == pop  & results$PGS == pgs & results$OUTCOME == trait] <- fit_glm_sum$Estimate[rownames(fit_glm_sum) == "PGS"]
      results$SE_LR[results$SUBSET == pop  & results$PGS == pgs & results$OUTCOME == trait] <- fit_glm_sum$`Std. Error`[rownames(fit_glm_sum) == "PGS"]
      results$P_LR[results$SUBSET == pop  & results$PGS == pgs & results$OUTCOME == trait] <- fit_glm_sum$`Pr(>|z|)`[rownames(fit_glm_sum) == "PGS"]
      
      
      #cox results
      results$BETA[results$SUBSET == pop  & results$PGS == pgs & results$OUTCOME == trait] <- fit_sum$coef[rownames(fit_sum) == "PGS"]
      results$HR[results$SUBSET == pop  & results$PGS == pgs & results$OUTCOME == trait] <- fit_sum$`exp(coef)`[rownames(fit_sum) == "PGS"]
      results$SE[results$SUBSET == pop  & results$PGS == pgs & results$OUTCOME == trait] <- fit_sum$`se(coef)`[rownames(fit_sum) == "PGS"]
      results$P[results$SUBSET == pop  & results$PGS == pgs & results$OUTCOME == trait] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PGS"]
      
    }
  }
}

#bonferonni correction
results$P_BONF <- results$P*40
results$P_BONF[results$P_BONF > 1] <- 1
results$P_FDR <- p.adjust(results$P, method="fdr")

#save results
data.table::fwrite(results, "./results/DECLARE_PGS_surival_analysis.tsv",
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)

#for CKD and T2D, test associations with HBA1C, egfr, insulin, DN, DR
traits <- c("EGFR_BL","UACR_BL", "DR", "DN", "INSULIN")
test_models <- c("T2D.V2_EUR_prscs", "CKD_TA_prscs" )

#run analyses
results <- expand.grid(PGS=test_models, 
                       OUTCOME=traits,
                       SUBSET=c("ALL", "EUR"),
                       N=NA, BETA=NA, SE=NA, P=NA,
                       stringsAsFactors = FALSE)

for(trait in traits){
  for(pgs in test_models){
    for(pop in c("ALL", "EUR")){
      if(pop == "ALL"){
        foo <- as.data.frame(dat[dat$OUTCOME == "MACE",])
      }else{
        foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9 & dat$OUTCOME == "MACE",])
      }
      
      #just include necessary columns
      pcs <- paste0("PC", 1:10)
      foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", trait, pcs, pgs)]
      colnames(foo)[ncol(foo)] <- "PGS"
      
      #scale pgs
      foo$PGS <- scale(foo$PGS)
      
      #specify formula
      f <- as.formula(paste(trait, "~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "PGS"), collapse = "+")))

      #fit logistic model
      if(trait == "EGFR_BL" | trait == "UACR_BL"){
        fit <- lm(f, foo)
      }else{
        fit <- glm(f, foo, family = "binomial")
      }
      
      
      #parse results
      fit_sum <- as.data.frame(summary(fit)$coefficients)
      
      #summary
      results$N[results$SUBSET == pop & results$PGS == pgs & results$OUTCOME == trait] <- length(fit$residuals)
      
      #glm results 
      results$BETA[results$SUBSET == pop  & results$PGS == pgs & results$OUTCOME == trait] <- fit_sum$Estimate[rownames(fit_sum) == "PGS"]
      results$SE[results$SUBSET == pop  & results$PGS == pgs & results$OUTCOME == trait] <- fit_sum$`Std. Error`[rownames(fit_sum) == "PGS"]
      if(trait == "EGFR_BL" | trait == "UACR_BL"){
        results$P[results$SUBSET == pop  & results$PGS == pgs & results$OUTCOME == trait] <- fit_sum$`Pr(>|t|)`[rownames(fit_sum) == "PGS"]
      }else{
        results$P[results$SUBSET == pop  & results$PGS == pgs & results$OUTCOME == trait] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PGS"]
      }
      
    }
  }
}

#drop the ckd-t2d ones, not the focus
results <- results[(results$PGS == "CKD_TA_prscs" & results$OUTCOME %in% c("UACR_BL", "EGFR_BL")) | (results$PGS == "T2D.V2_EUR_prscs" & results$OUTCOME %in% c("DN", "DR", "INSULIN")),]
write.table(results, "results/DECLARE_PGS_complications.txt", sep="\t",
            quote=FALSE, row.names = FALSE, col.names = TRUE)

