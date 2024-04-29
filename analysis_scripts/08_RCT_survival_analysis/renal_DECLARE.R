###co-authors suggested re-running analysis with a different set of clincal covariates for renal outcomes##
#this was the analysis that was included in the manuscript
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

#read in quantiative traits
q <- data.table::fread("./data/rct_data/DECLARE_quant_traits_no_slope.gz")
q <- q[,-c(2:5)]
dat <- merge(dat, q, by=c("ID", "EID"), all=FALSE, sort=FALSE)


#read in proteins measured at baseline
p <- data.table::fread("./data/rct_data/DECLARE_OLINK_visit2.tsv",
                       data.table = FALSE)


proteins <- colnames(p)[2:ncol(p)]

#add proteins to data frame
dat <- merge(dat, p, by="EID", all=FALSE, sort=FALSE)


#run analyses
results <- expand.grid(PROTEIN=proteins, 
                       TIMEPOINT=c("BASELINE", "12Mo", "DELTA"),
                       OUTCOME="RENAL",
                       SUBSET=c("ALL", "EUR"),
                       N=NA, N_EVENT=NA,
                       BETA=NA, HR=NA, SE=NA, P=NA,
                       stringsAsFactors = FALSE)

#atherosclerotic cardiovascular disease, heart failure, systolic blood pressure, T2D duration, glycated hemoglobin, eGFR, urine albumin-to-creatinine ratio, and hemoglobin
trait <- "RENAL"
for(protein in proteins){
  for(pop in c("ALL", "EUR")){
    if(pop == "ALL"){
      foo <- as.data.frame(dat[dat$OUTCOME == trait,])
    }else{
      foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9 & dat$OUTCOME == trait,])
    }
    
    
    #just include necessary columns
    pcs <- paste0("PC", 1:10)
    foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", pcs, "CAD", "HF", "SBP_BL", "DURDIAB", "HBA1C_BL", "EGFR_BL", "HEME_BL", "UACR_BL", protein)]
    colnames(foo)[ncol(foo)] <- "PROTEIN"
    
    #fit survival object
    sobj <- Surv(foo$TIME, foo$EVENT, type="right")
    foo$SURV <- sobj
    
    #specify formula
    #f <- as.formula(paste("EVENT~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", "CAD", "HF", "SBP_BL", "DURDIAB", "HBA1C_BL", "EGFR_BL", "HEME_BL", "UACR_BL", pcs, "PROTEIN"), collapse = "+")))
    f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", "CAD", "HF", "SBP_BL", "DURDIAB", "HBA1C_BL", "EGFR_BL", "HEME_BL", "UACR_BL", pcs, "PROTEIN"), collapse = "+")))
    
    #fit logistic model
    #fit_glm <- glm(f, foo, family = "binomial")
    #fit cox model
    fit <- coxph(f_cox, foo)
    
    #parse results
    fit_glm_sum <- as.data.frame(summary(fit_glm)$coefficients)
    fit_sum <- as.data.frame(summary(fit)$coefficients)
    
    #summary
    results$N[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- nrow(foo)
    results$N_EVENT[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE"  & results$OUTCOME == trait] <- nrow(foo[foo$EVENT == 1,])
    
    #cox results
    results$BETA[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$coef[rownames(fit_sum) == "PROTEIN"]
    results$HR[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`exp(coef)`[rownames(fit_sum) == "PROTEIN"]
    results$SE[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`se(coef)`[rownames(fit_sum) == "PROTEIN"]
    results$P[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PROTEIN"]
    
  }
}

#####repeat for 12 month timepoint#######
#read in covariate file with time to event data
dat <- data.table::fread("./data/rct_data/clinical_outcomes.gz")
dat <- dat[dat$STUDY == "DECLARE",]

#build analysis file
dat <- merge(dat, pca, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, anc, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, covar, by=c("EID", "ID"), all=FALSE, sort=FALSE)

dat <- merge(dat, q, by=c("ID", "EID"), all=FALSE, sort=FALSE)

#add age2
dat$AGE2 <- dat$AGE^2

#read in proteins measured after baseline 
p <- data.table::fread("./data/rct_data/DECLARE_OLINK_visit4.tsv",
                       data.table = FALSE)

proteins <- colnames(p)[2:ncol(p)]

#add proteins to data frame
dat <- merge(dat, p, by="EID", all=FALSE, sort=FALSE)


trait <- "RENAL"
for(protein in proteins){
  for(pop in c("ALL", "EUR")){
    if(pop == "ALL"){
      foo <- as.data.frame(dat[dat$OUTCOME == trait,])
    }else{
      foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9 & dat$OUTCOME == trait,])
    }
    
    
    #just include necessary columns
    pcs <- paste0("PC", 1:10)
    foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", pcs, "CAD", "HF", "SBP_BL", "DURDIAB", "HBA1C_BL", "EGFR_BL", "HEME_BL", "UACR_BL", protein)]
    colnames(foo)[ncol(foo)] <- "PROTEIN"
    
    #fit survival object
    sobj <- Surv(foo$TIME, foo$EVENT, type="right")
    foo$SURV <- sobj
    
    #specify formula
    #f <- as.formula(paste("EVENT~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", "CAD", "HF", "SBP_BL", "DURDIAB", "HBA1C_BL", "EGFR_BL", "HEME_BL", "UACR_BL", pcs, "PROTEIN"), collapse = "+")))
    f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", "CAD", "HF", "SBP_BL", "DURDIAB", "HBA1C_BL", "EGFR_BL", "HEME_BL", "UACR_BL", pcs, "PROTEIN"), collapse = "+")))
    
    #fit logistic model
    #fit_glm <- glm(f, foo, family = "binomial")
    #fit cox model
    fit <- coxph(f_cox, foo)
    
    #parse results
    #fit_glm_sum <- as.data.frame(summary(fit_glm)$coefficients)
    fit_sum <- as.data.frame(summary(fit)$coefficients)
    
    #summary
    results$N[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- nrow(foo)
    results$N_EVENT[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "12Mo"  & results$OUTCOME == trait] <- nrow(foo[foo$EVENT == 1,])
    
    #cox results
    results$BETA[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$coef[rownames(fit_sum) == "PROTEIN"]
    results$HR[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`exp(coef)`[rownames(fit_sum) == "PROTEIN"]
    results$SE[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`se(coef)`[rownames(fit_sum) == "PROTEIN"]
    results$P[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PROTEIN"]
    
  }
}


#delta
#read in covariate file with time to event data
dat <- data.table::fread("./data/rct_data/clinical_outcomes.gz")
dat <- dat[dat$STUDY == "DECLARE",]

#build analysis file
dat <- merge(dat, pca, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, anc, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, covar, by=c("EID", "ID"), all=FALSE, sort=FALSE)
dat <- merge(dat, q, by=c("ID", "EID"), all=FALSE, sort=FALSE)

#add age2
dat$AGE2 <- dat$AGE^2

#read in proteins measured at baseline
p1 <- data.table::fread("./data/rct_data/DECLARE_OLINK_visit2.tsv",
                       data.table = FALSE)

#read in proteins measured after baseline 
p2 <- data.table::fread("./data/rct_data/DECLARE_OLINK_visit4.tsv",
                       data.table = FALSE)


trait <- "RENAL"
for(protein in proteins){
  
  for(pop in c("ALL", "EUR")){
    if(pop == "ALL"){
      foo <- as.data.frame(dat[dat$OUTCOME == trait,])
    }else{
      foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9 & dat$OUTCOME == trait,])
    }
    
    #build protein file
    foo1 <- p1[c("EID", protein)]
    foo2 <- p2[c("EID", protein)]
    colnames(foo1)[2] <- "PROTEIN1"
    colnames(foo2)[2] <- "PROTEIN2"
    p <- merge(foo1, foo2, by="EID", all=FALSE, sort=FALSE)
    #p$DELTA <- (p$PROTEIN2 - p$PROTEIN1)/p$PROTEIN1
    p$DELTA <- (p$PROTEIN2 - p$PROTEIN1)
    
    foo <- merge(foo, p, by="EID", all=FALSE, sort=FALSE)
    #remove any instances of infinity
    foo <- foo[foo$DELTA != -Inf & foo$DELTA != Inf,]
    
    #just include necessary columns
    pcs <- paste0("PC", 1:10)
    foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", pcs, "CAD", "HF", "SBP_BL", "DURDIAB", "HBA1C_BL", "EGFR_BL", "HEME_BL", "UACR_BL", "DELTA")]
    
    #fit survival object
    sobj <- Surv(foo$TIME, foo$EVENT, type="right")
    foo$SURV <- sobj
    
    #specify formula
    #f <- as.formula(paste("EVENT~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", "CAD", "HF", "SBP_BL", "DURDIAB", "HBA1C_BL", "EGFR_BL", "HEME_BL", "UACR_BL", pcs, "DELTA"), collapse = "+")))
    f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", "CAD", "HF", "SBP_BL", "DURDIAB", "HBA1C_BL", "EGFR_BL", "HEME_BL", "UACR_BL", pcs, "DELTA"), collapse = "+")))
    
    #fit logistic model
    #fit_glm <- glm(f, foo, family = "binomial")
    #fit cox model
    fit <- coxph(f_cox, foo)
    
    #parse results
    #fit_glm_sum <- as.data.frame(summary(fit_glm)$coefficients)
    fit_sum <- as.data.frame(summary(fit)$coefficients)
    
    #summary
    results$N[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "DELTA" & results$OUTCOME == trait] <- nrow(foo)
    results$N_EVENT[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "DELTA"  & results$OUTCOME == trait] <- nrow(foo[foo$EVENT == 1,])
    
    #cox results
    results$BETA[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "DELTA" & results$OUTCOME == trait] <- fit_sum$coef[rownames(fit_sum) == "DELTA"]
    results$HR[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "DELTA" & results$OUTCOME == trait] <- fit_sum$`exp(coef)`[rownames(fit_sum) == "DELTA"]
    results$SE[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "DELTA" & results$OUTCOME == trait] <- fit_sum$`se(coef)`[rownames(fit_sum) == "DELTA"]
    results$P[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "DELTA" & results$OUTCOME == trait] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "DELTA"]
    
  }
}

#save results
data.table::fwrite(results, "./results/DECLARE_surival_analysis_renal_v2.tsv.gz",
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE,
                   compress = "gzip")
