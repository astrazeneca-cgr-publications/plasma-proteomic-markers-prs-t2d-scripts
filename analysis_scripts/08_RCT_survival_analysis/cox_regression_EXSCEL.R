#script for running cox proportional hazards regression with proteins as exposure/risk factor
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

#read in proteins measured at baseline
p <- data.table::fread("./data/rct_data/EXSCEL_SOMASCAN_baseline.tsv",
                       data.table = FALSE)

proteins <- colnames(p)[2:ncol(p)]

#add proteins to data frame
dat <- merge(dat, p, by="ID", all=FALSE, sort=FALSE)

#read in quantiative traits
q <- data.table::fread("./data/rct_data/EXSCEL_quant_traits_no_slope.gz")
q <- q[,-c(2:5)]
dat <- merge(dat, q, by=c("ID", "EID"), all=FALSE, sort=FALSE)


#run analyses
results <- expand.grid(PROTEIN=proteins, 
                       TIMEPOINT=c("BASELINE", "12Mo"),
                       OUTCOME=unique(dat$OUTCOME),
                       SUBSET=c("ALL", "EUR"),
                       N=NA, N_EVENT=NA,
                       BETA_LR=NA, SE_LR=NA, P_LR=NA,
                       BETA=NA, HR=NA, SE=NA, P=NA,
                       BETA_M2=NA, HR_M2=NA, SE_M2=NA, P_M2=NA,
                       BETA_M3=NA, HR_M3=NA, SE_M3=NA, P_M3=NA,
                       stringsAsFactors = FALSE)


for(trait in unique(dat$OUTCOME)){
  for(protein in proteins){
    for(pop in c("ALL", "EUR")){
      if(pop == "ALL"){
        foo <- as.data.frame(dat[dat$OUTCOME == trait,])
      }else{
        foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9 & dat$OUTCOME == trait,])
      }
      
      #drop subjecst with insulin
      if(trait == "INSULIN") foo <- foo[foo$INSTFL == "N",]
      if(trait == "RENAL") foo <- foo[foo$EGFR_BL >= 30, ]
      
      #just include necessary columns
      pcs <- paste0("PC", 1:10)
      foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", pcs, protein)]
      colnames(foo)[ncol(foo)] <- "PROTEIN"
      
      #fit survival object
      sobj <- Surv(foo$TIME, foo$EVENT, type="right")
      foo$SURV <- sobj
      
      #specify formula
      f <- as.formula(paste("EVENT~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "PROTEIN"), collapse = "+")))
      f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "PROTEIN"), collapse = "+")))
      
      #fit logistic model
      fit_glm <- glm(f, foo, family = "binomial")
      #fit cox model
      fit <- coxph(f_cox, foo)
      
      #parse results
      fit_glm_sum <- as.data.frame(summary(fit_glm)$coefficients)
      fit_sum <- as.data.frame(summary(fit)$coefficients)
      
      #summary
      results$N[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- nrow(foo)
      results$N_EVENT[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- nrow(foo[foo$EVENT == 1,])
      
      #glm results
      results$BETA_LR[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_glm_sum$Estimate[rownames(fit_glm_sum) == "PROTEIN"]
      results$SE_LR[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_glm_sum$`Std. Error`[rownames(fit_glm_sum) == "PROTEIN"]
      results$P_LR[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_glm_sum$`Pr(>|z|)`[rownames(fit_glm_sum) == "PROTEIN"]
      
      
      #cox results
      results$BETA[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$coef[rownames(fit_sum) == "PROTEIN"]
      results$HR[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`exp(coef)`[rownames(fit_sum) == "PROTEIN"]
      results$SE[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`se(coef)`[rownames(fit_sum) == "PROTEIN"]
      results$P[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PROTEIN"]
      
    }
  }
}

#now, more stringent model

for(trait in unique(dat$OUTCOME)){
  for(protein in proteins){
    for(pop in c("ALL", "EUR")){
      
      #skip if p-value > 0.05
      if(results$P[results$SUBSET == pop & results$OUTCOME == trait & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE"] > 0.05) next
      
      #subset
      if(pop == "ALL"){
        foo <- as.data.frame(dat[dat$OUTCOME == trait,])
      }else{
        foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9 & dat$OUTCOME == trait,])
      }
      
      #drop subjecst with insulin
      if(trait == "INSULIN") foo <- foo[foo$INSTFL == "N",]
      if(trait == "RENAL") foo <- foo[foo$EGFR_BL >= 30, ]
      
      if(trait == "MACE"){
        #hypertension, hypercholesterolemia, obesity (BMI >30 kg/m²), smoking (current, or smoking cessation ≤3 months), positive family history (parent or sibling with CVD before age 65), and atherosclerotic disease (prior MI, PCI/CABG, CVA/TIA, or peripheral arterial disease). 
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "BMI_BL", "CVD_ANY", "CURRENT_SMOKER", "HYPERTENSION", "HYPERLIPIDEMIA", "NPPB_P16860_seq.7655.11", 
                     pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "BMI_BL", "CVD_ANY", "CURRENT_SMOKER", "HYPERTENSION", "HYPERLIPIDEMIA", "PROTEIN"), collapse = "+")))
        f_cox2 <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "BMI_BL", "CVD_ANY", "CURRENT_SMOKER", "HYPERTENSION", "HYPERLIPIDEMIA", "NPPB_P16860_seq.7655.11","PROTEIN"), collapse = "+")))
        
      }else if(trait == "HHF"){
        #coronary artery disease, AFIB, baseline eGFR, baseline UACR, and prior heart failure
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "EGFR_BL", "AFIB",  "CVD_ANY", "HF", "NPPB_P16860_seq.7655.11", 
                     pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "EGFR_BL", "AFIB",  "CVD_ANY", "HF", "PROTEIN"), collapse = "+")))
        f_cox2 <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "EGFR_BL", "AFIB",  "CVD_ANY", "HF", "NPPB_P16860_seq.7655.11","PROTEIN"), collapse = "+")))
      }else if(trait == "RENAL"){
        #eGFR at baseline, UACR at baseline, calcium at baseline, albumin at baseline, and bicarbonate at baseline
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "EGFR_BL", "ALBUMINERIA", pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "EGFR_BL", "ALBUMINERIA", "PROTEIN"), collapse = "+")))
        
      }else{
        #we adjusted for cardiovascular disease, CKD, number of diabetic complications (e.g., diabetic retinopathy, neuropathy, amputations), and duration of diabetes
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "CVD_ANY", "DIABDUR", "ALBUMINERIA" , "NUMBER_COMPLICATIONS", pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "CVD_ANY", "DIABDUR", "ALBUMINERIA" , "NUMBER_COMPLICATIONS", "PROTEIN"), collapse = "+")))
        
      }
      
      #fit survival object
      sobj <- Surv(foo$TIME, foo$EVENT, type="right")
      foo$SURV <- sobj
      
      #fit cox model
      fit <- coxph(f_cox, foo)
      
      #parse results
      fit_sum <- as.data.frame(summary(fit)$coefficients)
      
      #cox results
      results$BETA_M2[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$coef[rownames(fit_sum) == "PROTEIN"]
      results$HR_M2[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`exp(coef)`[rownames(fit_sum) == "PROTEIN"]
      results$SE_M2[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`se(coef)`[rownames(fit_sum) == "PROTEIN"]
      results$P_M2[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PROTEIN"]
      
      if(trait %in% c("MACE", "HHF")){
        #fit cox model with nt-proBNP
        fit <- coxph(f_cox2, foo)
        
        #parse results
        fit_sum <- as.data.frame(summary(fit)$coefficients)
        
        #cox results
        results$BETA_M3[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$coef[rownames(fit_sum) == "PROTEIN"]
        results$HR_M3[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`exp(coef)`[rownames(fit_sum) == "PROTEIN"]
        results$SE_M3[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`se(coef)`[rownames(fit_sum) == "PROTEIN"]
        results$P_M3[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE" & results$OUTCOME == trait] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PROTEIN"]
        
      }
    }
  }
}



#####repeat for 12 month timepoint#######
#read in covariate file with time to event data
dat <- data.table::fread("./data/rct_data/clinical_outcomes.gz")
dat <- dat[dat$STUDY == "EXSCEL",]

#build analysis file
dat <- merge(dat, pca, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, anc, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, covar, by=c("EID", "ID"), all=FALSE, sort=FALSE)

#add age2
dat$AGE2 <- dat$AGE^2

#read in proteins measured at baseline
p <- data.table::fread("./data/rct_data/EXSCEL_SOMASCAN_month12.tsv",
                       data.table = FALSE)

proteins <- colnames(p)[2:ncol(p)]

#add proteins to data frame
dat <- merge(dat, p, by="ID", all=FALSE, sort=FALSE)

#reead in quantiative traits
q <- data.table::fread("./data/rct_data/EXSCEL_quant_traits_no_slope.gz")
q <- q[,-c(2:5)]
dat <- merge(dat, q, by=c("ID", "EID"), all=FALSE, sort=FALSE)


for(trait in unique(dat$OUTCOME)){
  for(protein in proteins){
    for(pop in c("ALL", "EUR")){
      if(pop == "ALL"){
        foo <- as.data.frame(dat[dat$OUTCOME == trait,])
      }else{
        foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9 & dat$OUTCOME == trait,])
      }
      
      #exclude subjects who had insulin at baseline
      if(trait == "INSULIN") foo <- foo[foo$INSTFL !=  "Y",]
      if(trait == "RENAL") foo <- foo[foo$EGFR_BL >= 30, ]
      
      #just include necessary columns
      pcs <- paste0("PC", 1:10)
      foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", pcs, protein)]
      colnames(foo)[ncol(foo)] <- "PROTEIN"
      
      #fit survival object
      sobj <- Surv(foo$TIME, foo$EVENT, type="right")
      foo$SURV <- sobj
      
      #specify formula
      f <- as.formula(paste("EVENT~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "PROTEIN"), collapse = "+")))
      f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "PROTEIN"), collapse = "+")))
      
      #fit logistic model
      fit_glm <- glm(f, foo, family = "binomial")
      #fit cox model
      fit <- coxph(f_cox, foo)
      
      #parse results
      fit_glm_sum <- as.data.frame(summary(fit_glm)$coefficients)
      fit_sum <- as.data.frame(summary(fit)$coefficients)
      
      #summary
      results$N[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- nrow(foo)
      results$N_EVENT[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "12Mo"] <- nrow(foo[foo$EVENT == 1,])
      
      #glm results
      results$BETA_LR[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_glm_sum$Estimate[rownames(fit_glm_sum) == "PROTEIN"]
      results$SE_LR[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_glm_sum$`Std. Error`[rownames(fit_glm_sum) == "PROTEIN"]
      results$P_LR[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_glm_sum$`Pr(>|z|)`[rownames(fit_glm_sum) == "PROTEIN"]
      
      
      #cox results
      results$BETA[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$coef[rownames(fit_sum) == "PROTEIN"]
      results$HR[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`exp(coef)`[rownames(fit_sum) == "PROTEIN"]
      results$SE[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`se(coef)`[rownames(fit_sum) == "PROTEIN"]
      results$P[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PROTEIN"]
      
    }
  }
}

#now, more stringent model

for(trait in unique(dat$OUTCOME)){
  for(protein in proteins){
    for(pop in c("ALL", "EUR")){
      
      #skip if p-value > 0.05
      if(results$P[results$SUBSET == pop & results$OUTCOME == trait & results$PROTEIN == protein & results$TIMEPOINT == "12Mo"] > 0.05) next
      
      #subset
      if(pop == "ALL"){
        foo <- as.data.frame(dat[dat$OUTCOME == trait,])
      }else{
        foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9 & dat$OUTCOME == trait,])
      }
      
      #exclude subjects who had insulin at baseline
      if(trait == "INSULIN") foo <- foo[foo$INSTFL !=  "Y",]
      if(trait == "RENAL") foo <- foo[foo$EGFR_BL >= 30, ]
      
      if(trait == "MACE"){
        #hypertension, hypercholesterolemia, obesity (BMI >30 kg/m²), smoking (current, or smoking cessation ≤3 months), positive family history (parent or sibling with CVD before age 65), and atherosclerotic disease (prior MI, PCI/CABG, CVA/TIA, or peripheral arterial disease). 
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "BMI_BL", "CVD_ANY", "CURRENT_SMOKER", "HYPERTENSION", "HYPERLIPIDEMIA", "NPPB_P16860_seq.7655.11", 
                     pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "BMI_BL", "CVD_ANY", "CURRENT_SMOKER", "HYPERTENSION", "HYPERLIPIDEMIA", "PROTEIN"), collapse = "+")))
        f_cox2 <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "BMI_BL", "CVD_ANY", "CURRENT_SMOKER", "HYPERTENSION", "HYPERLIPIDEMIA", "NPPB_P16860_seq.7655.11","PROTEIN"), collapse = "+")))
        
      }else if(trait == "HHF"){
        #coronary artery disease, AFIB, baseline eGFR, baseline UACR, and prior heart failure
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "EGFR_BL", "AFIB",  "CVD_ANY", "HF", "NPPB_P16860_seq.7655.11", 
                     pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "EGFR_BL", "AFIB",  "CVD_ANY", "HF", "PROTEIN"), collapse = "+")))
        f_cox2 <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "EGFR_BL", "AFIB",  "CVD_ANY", "HF", "NPPB_P16860_seq.7655.11","PROTEIN"), collapse = "+")))
      }else if(trait == "RENAL"){
        #eGFR at baseline, UACR at baseline, calcium at baseline, albumin at baseline, and bicarbonate at baseline
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "EGFR_BL", "ALBUMINERIA", pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "EGFR_BL", "ALBUMINERIA", "PROTEIN"), collapse = "+")))
        
      }else{
        #we adjusted for cardiovascular disease, CKD, number of diabetic complications (e.g., diabetic retinopathy, neuropathy, amputations), and duration of diabetes
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "CVD_ANY", "DIABDUR", "ALBUMINERIA" , "NUMBER_COMPLICATIONS", pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "CVD_ANY", "DIABDUR", "ALBUMINERIA" , "NUMBER_COMPLICATIONS", "PROTEIN"), collapse = "+")))
        
      }
      
      #fit survival object
      sobj <- Surv(foo$TIME, foo$EVENT, type="right")
      foo$SURV <- sobj
      
      #fit cox model
      fit <- coxph(f_cox, foo)
      
      #parse results
      fit_sum <- as.data.frame(summary(fit)$coefficients)
      
      #cox results
      results$BETA_M2[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$coef[rownames(fit_sum) == "PROTEIN"]
      results$HR_M2[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`exp(coef)`[rownames(fit_sum) == "PROTEIN"]
      results$SE_M2[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`se(coef)`[rownames(fit_sum) == "PROTEIN"]
      results$P_M2[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PROTEIN"]
      
      if(trait %in% c("MACE", "HHF")){
        #fit cox model with nt-proBNP
        fit <- coxph(f_cox2, foo)
        
        #parse results
        fit_sum <- as.data.frame(summary(fit)$coefficients)
        
        #cox results
        results$BETA_M3[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$coef[rownames(fit_sum) == "PROTEIN"]
        results$HR_M3[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`exp(coef)`[rownames(fit_sum) == "PROTEIN"]
        results$SE_M3[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`se(coef)`[rownames(fit_sum) == "PROTEIN"]
        results$P_M3[results$SUBSET == pop  & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PROTEIN"]
        
      }
    }
  }
}


#save results
data.table::fwrite(results, "./results/EXSCEL_surival_analysis.tsv.gz",
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE,
                   compress = "gzip")
