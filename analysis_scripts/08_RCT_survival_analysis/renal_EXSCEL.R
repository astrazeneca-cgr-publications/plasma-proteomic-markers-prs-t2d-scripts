###co-authors suggested re-running analysis with a different set of clincal covariates for renal outcomes##
#this was the analysis that was included in the manuscript

library(survival)
#read in covariate file with time to event data
dat <- data.table::fread("./data/rct_data/clinical_outcomes.gz")
dat <- dat[dat$STUDY == "EXSCEL",]

#read in risk factors file
covar <- data.table::fread("./data/rct_data/EXSCEL_risk_factors.gz")

#read in PCA
pca <- data.table::fread("./data/rct_data/EXSCEL-PCA.eigenvec")

#read in ancestry definitions
anc <- data.table::fread("./data/rct_data/EXSCEL_ancestry_InferredAncestry.txt")
anc$PC1 <- NULL
anc$PC2 <- NULL 

#build analysis file
dat <- merge(dat, pca, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, anc, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, covar, by=c("EID", "ID"), all=FALSE, sort=FALSE)


#add age2
dat$AGE2 <- dat$AGE^2

#read in proteins measured at baseline and 12 months
p1 <- data.table::fread("./data/rct_data/EXSCEL_SOMASCAN_baseline.tsv",
                        data.table = FALSE)
p2 <- data.table::fread("./data/rct_data/EXSCEL_SOMASCAN_month12.tsv",
                        data.table = FALSE)

proteins <- unique(c(colnames(p1)[2:ncol(p1)], colnames(p2)[2:ncol(p2)]))


#reead in quantiative traits
q <- data.table::fread("./data/rct_data/EXSCEL_quant_traits_no_slope.gz")
q <- q[,-c(2:5)]
dat <- merge(dat, q, by=c("ID", "EID"), all=FALSE, sort=FALSE)


#run analyses
results <- expand.grid(PROTEIN=proteins, 
                       TIMEPOINT=c("BASELINE", "12Mo", "DELTA"),
                       OUTCOME="RENAL",
                       SUBSET=c("ALL", "EUR"),
                       N=NA, N_EVENT=NA,
                       BETA=NA, HR=NA, SE=NA, P=NA,
                       stringsAsFactors = FALSE)


trait <- "RENAL"
for(protein in proteins){
  for(pop in c("ALL", "EUR")){
    if(pop == "ALL"){
      foo <- as.data.frame(dat[dat$OUTCOME == trait,])
    }else{
      foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9 & dat$OUTCOME == trait,])
    }
    
    #add protein data
    foo <- merge(foo, p1, by="ID", all=FALSE, sort=FALSE)
    
    #drop subjecst with insulin
    if(trait == "RENAL") foo <- foo[foo$EGFR_BL >= 30, ]
    
    
    #atherosclerotic cardiovascular disease, heart failure, systolic blood pressure, T2D duration, glycated hemoglobin, eGFR, urine albumin-to-creatinine ratio, and hemoglobin
    #note for heme and UACR, greatly reduces sample size 
    pcs <- paste0("PC", 1:10)
    foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "CVD_ANY", "HF", "DIABDUR", "HBA1C_BL", "EGFR_BL", "SBP_BL", "ALBUMINERIA", pcs, protein)]
    colnames(foo)[ncol(foo)] <- "PROTEIN"
    
    #fit survival object
    sobj <- Surv(foo$TIME, foo$EVENT, type="right")
    foo$SURV <- sobj
    
    #specify formula
    f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "CVD_ANY", "HF", "DIABDUR", "HBA1C_BL", "EGFR_BL", "SBP_BL", "ALBUMINERIA", "PROTEIN"), collapse = "+")))
    
    #fit cox model
    fit <- coxph(f_cox, foo)
    
    #parse results
    fit_sum <- as.data.frame(summary(fit)$coefficients)
    
    #summary
    results$N[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- nrow(foo[!is.na(foo$EVENT),])
    results$N_EVENT[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- nrow(foo[foo$EVENT == 1 & !is.na(foo$EVENT),])
    
    #cox results
    results$BETA[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$coef[rownames(fit_sum) == "PROTEIN"]
    results$HR[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`exp(coef)`[rownames(fit_sum) == "PROTEIN"]
    results$SE[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`se(coef)`[rownames(fit_sum) == "PROTEIN"]
    results$P[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PROTEIN"]
    
  }
}


#####repeat for 12 month timepoint#######
for(protein in proteins){
  for(pop in c("ALL", "EUR")){
    if(pop == "ALL"){
      foo <- as.data.frame(dat[dat$OUTCOME == trait,])
    }else{
      foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9 & dat$OUTCOME == trait,])
    }
    
    #add protein data
    foo <- merge(foo, p2, by="ID", all=FALSE, sort=FALSE)
    
    #drop subjecst with insulin
    if(trait == "RENAL") foo <- foo[foo$EGFR_BL >= 30, ]
    
    
    #atherosclerotic cardiovascular disease, heart failure, systolic blood pressure, T2D duration, glycated hemoglobin, eGFR, urine albumin-to-creatinine ratio, and hemoglobin
    #note for heme and UACR, greatly reduces sample size 
    pcs <- paste0("PC", 1:10)
    foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "CVD_ANY", "HF", "DIABDUR", "HBA1C_BL", "EGFR_BL", "SBP_BL", "ALBUMINERIA", pcs, protein)]
    colnames(foo)[ncol(foo)] <- "PROTEIN"
    
    #fit survival object
    sobj <- Surv(foo$TIME, foo$EVENT, type="right")
    foo$SURV <- sobj
    
    #specify formula
    f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "CVD_ANY", "HF", "DIABDUR", "HBA1C_BL", "EGFR_BL", "SBP_BL", "ALBUMINERIA", "PROTEIN"), collapse = "+")))
    
    #fit cox model
    fit <- coxph(f_cox, foo)
    
    #parse results
    fit_sum <- as.data.frame(summary(fit)$coefficients)
    
    #summary
    results$N[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- nrow(foo)
    results$N_EVENT[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "12Mo"] <- nrow(foo[foo$EVENT == 1 & !is.na(foo$EVENT),])
    
    #cox results
    results$BETA[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$coef[rownames(fit_sum) == "PROTEIN"]
    results$HR[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`exp(coef)`[rownames(fit_sum) == "PROTEIN"]
    results$SE[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`se(coef)`[rownames(fit_sum) == "PROTEIN"]
    results$P[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PROTEIN"]
    
  }
}

#delta
for(protein in proteins){
  for(pop in c("ALL", "EUR")){
    if(pop == "ALL"){
      foo <- as.data.frame(dat[dat$OUTCOME == trait,])
    }else{
      foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9 & dat$OUTCOME == trait,])
    }
    
    #build protein file
    foo1 <- p1[c("ID", protein)]
    foo2 <- p2[c("ID", protein)]
    colnames(foo1)[2] <- "PROTEIN1"
    colnames(foo2)[2] <- "PROTEIN2"
    p <- merge(foo1, foo2, by="ID", all=FALSE, sort=FALSE)
    
    #p$DELTA <- (p$PROTEIN2 - p$PROTEIN1)/p$PROTEIN1
    p$DELTA <- (p$PROTEIN2 - p$PROTEIN1)
    foo <- merge(foo, p, by="ID", all=FALSE, sort=FALSE)
    
    #remove any instances of infinity
    foo <- foo[foo$DELTA != -Inf & foo$DELTA != Inf,]
    
    #drop subjecst with insulin
    if(trait == "RENAL") foo <- foo[foo$EGFR_BL >= 30, ]
    
    
    #atherosclerotic cardiovascular disease, heart failure, systolic blood pressure, T2D duration, glycated hemoglobin, eGFR, urine albumin-to-creatinine ratio, and hemoglobin
    #note for heme and UACR, greatly reduces sample size 
    pcs <- paste0("PC", 1:10)
    foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "CVD_ANY", "HF", "DIABDUR", "HBA1C_BL", "EGFR_BL", "SBP_BL", "ALBUMINERIA", pcs, "DELTA")]
    
    #fit survival object
    sobj <- Surv(foo$TIME, foo$EVENT, type="right")
    foo$SURV <- sobj
    
    #specify formula
    f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "CVD_ANY", "HF", "DIABDUR", "HBA1C_BL", "EGFR_BL", "SBP_BL", "ALBUMINERIA", "DELTA"), collapse = "+")))
    
    #fit cox model
    fit <- coxph(f_cox, foo)
    
    #parse results
    fit_sum <- as.data.frame(summary(fit)$coefficients)
    
    #summary
    results$N[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "DELTA" & results$OUTCOME == trait] <- nrow(foo)
    results$N_EVENT[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "DELTA"] <- nrow(foo[foo$EVENT == 1,])
    
    #cox results
    results$BETA[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "DELTA" & results$OUTCOME == trait] <- fit_sum$coef[rownames(fit_sum) == "DELTA"]
    results$HR[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "DELTA" & results$OUTCOME == trait] <- fit_sum$`exp(coef)`[rownames(fit_sum) == "DELTA"]
    results$SE[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "DELTA" & results$OUTCOME == trait] <- fit_sum$`se(coef)`[rownames(fit_sum) == "DELTA"]
    results$P[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "DELTA" & results$OUTCOME == trait] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "DELTA"]
    
  }
}


#save results
data.table::fwrite(results, "./results/EXSCEL_renal_analysis.tsv.gz",
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE,
                   compress = "gzip")
