#script for running cox proportional hazards regression in DECLARE
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

#read in proteins measured at baseline
p <- data.table::fread("./data/rct_data/DECLARE_OLINK_visit2.tsv",
                       data.table = FALSE)

proteins <- colnames(p)[2:ncol(p)]

#add proteins to data frame
dat <- merge(dat, p, by="EID", all=FALSE, sort=FALSE)

#read in quantiative traits
q <- data.table::fread("./data/rct_data/DECLARE_quant_traits.gz")
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
      
      #exclude subjects with insulin at baseline for insulin analyses
      if(trait == "INSULIN") foo <- foo[foo$INSULIN == 0,]
      
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
      results$N_EVENT[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "BASELINE"  & results$OUTCOME == trait] <- nrow(foo[foo$EVENT == 1,])
      
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
      
      #exclude subjects with insulin at baseline for insulin analyses
      if(trait == "INSULIN") foo <- foo[foo$INSULIN == 0,]
      
      if(trait == "MACE"){
        #for MACE, multiple risk factors flag covers everything except BMI
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "ECVDFL", "BMI_BL", "MRFFL", "NTproBNP_NA_OID00131", 
                     pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "ECVDFL", "BMI_BL", "MRFFL", "PROTEIN"), collapse = "+")))
        f_cox2 <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "ECVDFL", "BMI_BL", "MRFFL", "NTproBNP_NA_OID00131","PROTEIN"), collapse = "+")))
        
      }else if(trait == "HHF"){
        #coronary artery disease, AFIB, baseline eGFR, baseline UACR, and prior heart failure
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "UACR_BL", "EGFR_BL", "AFIB",  "ECVDFL", "HF", "NTproBNP_NA_OID00131", 
                     pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs,"UACR_BL", "EGFR_BL", "AFIB",  "ECVDFL", "HF", "PROTEIN"), collapse = "+")))
        f_cox2 <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "UACR_BL", "EGFR_BL", "AFIB",  "ECVDFL", "HF", "NTproBNP_NA_OID00131","PROTEIN"), collapse = "+")))
      }else if(trait == "RENAL"){
        #eGFR at baseline, UACR at baseline, calcium at baseline, albumin at baseline, and bicarbonate at baseline
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "UACR_BL", "EGFR_BL", "CALCIUM_BL", "ALBUMIN_BL", "BICARB_BL", pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "UACR_BL", "EGFR_BL", "CALCIUM_BL", "ALBUMIN_BL", "BICARB_BL", "PROTEIN"), collapse = "+")))
        
      }else{
        #we adjusted for cardiovascular disease, CKD, number of diabetic complications (e.g., diabetic retinopathy, neuropathy, amputations), and duration of diabetes
        foo$N_COMPS <- foo$AMPUTATION + foo$DN + foo$DR
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "ECVDFL", "DURDIAB", "CKD" , "N_COMPS", pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "ECVDFL", "DURDIAB", "CKD" , "N_COMPS", "PROTEIN"), collapse = "+")))
        
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
dat <- dat[dat$STUDY == "DECLARE",]

#build analysis file
dat <- merge(dat, pca, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, anc, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, covar, by=c("EID", "ID"), all=FALSE, sort=FALSE)

#add age2
dat$AGE2 <- dat$AGE^2

#read in proteins measured at baseline
p <- data.table::fread("./data/rct_data/DECLARE_OLINK_visit4.tsv",
                       data.table = FALSE)

proteins <- colnames(p)[2:ncol(p)]

#add proteins to data frame
dat <- merge(dat, p, by="EID", all=FALSE, sort=FALSE)

#reead in quantiative traits
q <- data.table::fread("./data/rct_data/DECLARE_quant_traits.gz")
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
      
      #exclude subjects with insulin at baseline for insulin analyses
      if(trait == "INSULIN") foo <- foo[foo$INSULIN == 0,]
      
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
      results$N_EVENT[results$SUBSET == pop & results$PROTEIN == protein & results$TIMEPOINT == "12Mo" & results$OUTCOME == trait] <- nrow(foo[foo$EVENT == 1,])
      
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
      
      #exclude subjects with insulin at baseline for insulin analyses
      if(trait == "INSULIN") foo <- foo[foo$INSULIN == 0,]
      
      if(trait == "MACE"){
        #for MACE, multiple risk factors flag covers everything except BMI
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "ECVDFL", "BMI_BL", "MRFFL", "NTproBNP_NA_OID00131", 
                     pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "ECVDFL", "BMI_BL", "MRFFL", "PROTEIN"), collapse = "+")))
        f_cox2 <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "ECVDFL", "BMI_BL", "MRFFL", "NTproBNP_NA_OID00131","PROTEIN"), collapse = "+")))
        
      }else if(trait == "HHF"){
        #coronary artery disease, AFIB, baseline eGFR, baseline UACR, and prior heart failure
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "UACR_BL", "EGFR_BL", "AFIB",  "ECVDFL", "HF", "NTproBNP_NA_OID00131", 
                     pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs,"UACR_BL", "EGFR_BL", "AFIB",  "ECVDFL", "HF", "PROTEIN"), collapse = "+")))
        f_cox2 <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "UACR_BL", "EGFR_BL", "AFIB",  "ECVDFL", "HF", "NTproBNP_NA_OID00131","PROTEIN"), collapse = "+")))
      }else if(trait == "RENAL"){
        #eGFR at baseline, UACR at baseline, calcium at baseline, albumin at baseline, and bicarbonate at baseline
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "UACR_BL", "EGFR_BL", "CALCIUM_BL", "ALBUMIN_BL", "BICARB_BL", pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "UACR_BL", "EGFR_BL", "CALCIUM_BL", "ALBUMIN_BL", "BICARB_BL", "PROTEIN"), collapse = "+")))
        
      }else{
        #we adjusted for cardiovascular disease, CKD, number of diabetic complications (e.g., diabetic retinopathy, neuropathy, amputations), and duration of diabetes
        foo$N_COMPS <- foo$AMPUTATION + foo$DN + foo$DR
        pcs <- paste0("PC", 1:10)
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", "ECVDFL", "DURDIAB", "CKD" , "N_COMPS", pcs, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        #specify formula
        f_cox <- as.formula(paste("SURV~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "ECVDFL", "DURDIAB", "CKD" , "N_COMPS", "PROTEIN"), collapse = "+")))
        
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
data.table::fwrite(results, "./results/DECLARE_surival_analysis.tsv.gz",
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE,
                   compress = "gzip")
