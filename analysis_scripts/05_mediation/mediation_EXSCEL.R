#script for running mediation analyses in EXSCEL using medflex

library(medflex)

#read in file with time to event data
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

#read in PGS
top_pgs <- data.table::fread("./results/top_PGS_per_trait.tsv")
prs.dir <- "./data/rct_data/scores/"
prs.files <- list.files(prs.dir, "sscore")
prs.files <- prs.files[grep("EXSCEL", prs.files)]

for(prs in prs.files){
  prs.file <- paste0(prs.dir, prs)
  foo <- data.table::fread(prs.file)
  colnames(foo)[1] <- "IID"
  prs_label <- gsub("EXSCEL_", "", gsub(".sscore", "", prs))
  
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

#build analysis file
dat <- merge(dat, pca, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, anc, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, covar, by=c("EID", "ID"), all=FALSE, sort=FALSE)
dat <- merge(dat, prs.dat, by="ID", all=FALSE, sort=FALSE)

#add age2
dat$AGE2 <- dat$AGE^2

#regress PCs from prs
pcs <- paste0("PC", 1:10)
for(prs in pgs_models){
  f <- as.formula(paste(prs, "~", paste(pcs, collapse = "+")))
  fit <- lm(f, dat)
  dat[[prs]] <- fit$residuals
}

#read in proteins measured at baseline
p <- data.table::fread("./data/rct_data/EXSCEL_SOMASCAN_baseline.tsv",
                       data.table = FALSE)

proteins <- colnames(p)[2:ncol(p)]

#add proteins to data frame
dat <- merge(dat, p, by="ID", all=FALSE, sort=FALSE)

#reead in quantiative traits
q <- data.table::fread("./data/rct_data/EXSCEL_quant_traits.gz")
q <- q[,-c(2:5)]
dat <- merge(dat, q, by=c("ID", "EID"), all=FALSE, sort=FALSE)

type <- "event"

if(type == "outcome"){
  #prepare analysis
  traits <- unique(dat$OUTCOME)
  
  results <- expand.grid(TRAIT=traits, PGS=pgs_models, PROTEIN=proteins,
                         SUBSET=c("ALL", "EUR"),
                         BETA_PGS=NA, SE_PGS=NA, P_PGS=NA,
                         OR_PGS_TRAIT=NA, SE_PGS_TRAIT=NA, P_PGS_TRAIT=NA,
                         OR=NA, SE=NA, P=NA,
                         OR_DIRECT=NA, SE_DIRECT=NA, P_DIRECT=NA,
                         OR_INDIRECT=NA, SE_INDIRECT=NA, P_INDIRECT=NA,
                         OR_TOTAL=NA, SE_TOTAL=NA, P_TOTAL=NA)
  
  #start iteration on trait
  for(trait in traits){
    #protein
    for(protein in proteins){
      
      #subset
      for(pop in c("ALL", "EUR")){
        if(pop == "ALL"){
          foo <- as.data.frame(dat[dat$OUTCOME == trait,])
        }else{
          foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9 & dat$OUTCOME == trait,])
        }
        
        #exclude subjects who had insulin at baseline
        if(trait == "INSULIN") foo <- foo[foo$INSTFL !=  "Y",]
        if(trait == "RENAL") foo <- foo[foo$EGFR_BL >= 30, ]
        
        #test protein for association with a trait
        pcs <- paste0("PC", 1:10)
        f <- as.formula(paste("EVENT~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "PROTEIN"), collapse = "+")))
        
        #just include necessary columns
        foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", "TIME", pcs, pgs_models, protein)]
        colnames(foo)[ncol(foo)] <- "PROTEIN"
        
        #fit logistic model
        fit <- glm(f, foo, family = "binomial")
        fit_sum <- as.data.frame(summary(fit)$coefficients)
        results$OR[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein] <- exp(fit_sum$Estimate[rownames(fit_sum) == "PROTEIN"])
        results$SE[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein] <- fit_sum$`Std. Error`[rownames(fit_sum) == "PROTEIN"]
        results$P[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PROTEIN"]
        
        #skip if not nominally signiificant
        if(fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PROTEIN"] > 0.05) next
        #pgs in pgs models
        for(pgs in pgs_models){
          #scale PGS
          foo[[pgs]] <- scale(foo[[pgs]])
          
          pcs <- paste0("PC", 1:10)
          f <- as.formula(paste("PROTEIN~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX", pcs, pgs), collapse = "+")))
          #fit pgs-npx model
          fit <- lm(f, foo)
          fit_sum <- as.data.frame(summary(fit)$coefficients)
          results$BETA_PGS[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- fit_sum$Estimate[rownames(fit_sum) == pgs]
          results$SE_PGS[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- fit_sum$`Std. Error`[rownames(fit_sum) == pgs]
          results$P_PGS[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- fit_sum$`Pr(>|t|)`[rownames(fit_sum) == pgs]
          if(fit_sum$`Pr(>|t|)`[rownames(fit_sum) == pgs] > 0.05) next
          
          
          #test if pgs is associated with outcome
          f <- as.formula(paste("EVENT~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX", pcs, pgs), collapse = "+")))
          fit <- glm(f, family = "binomial", data=foo)
          
          fit_sum <- as.data.frame(summary(fit)$coefficients)
          results$OR_PGS_TRAIT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- exp(fit_sum$Estimate[rownames(fit_sum) == pgs])
          results$SE_PGS_TRAIT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- fit_sum$`Std. Error`[rownames(fit_sum) == pgs]
          results$P_PGS_TRAIT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == pgs]
          
          if(fit_sum$`Pr(>|z|)`[rownames(fit_sum) == pgs] > 0.05) next
          
          print(paste("running mediation for", pgs, "and", trait, "and", protein))
          #need to scale age for mediation
          foo$AGE <- scale(foo$AGE)
          foo$AGE <- as.numeric(foo$AGE)
          
          #scale PCs
          for(pc in pcs){
            foo[[pc]] <- as.numeric(scale(foo[[pc]]))
          }
          
          #change pgs column
          colnames(foo)[colnames(foo) == pgs] <- "PGS"
          foo$PGS <- as.numeric(foo$PGS)
          
          #remove missing subjects
          foo <- foo[complete.cases(foo),]
          
          #run mediation analysis
          f <- as.formula(paste("EVENT~", paste(c("PGS", "PROTEIN", "SEX", "AGE", pcs), collapse = "+")))
          expData <- neImpute(f, family = binomial("logit"), data=foo)
          f <- as.formula(paste("EVENT~", paste(c(paste0("PGS", "0"), paste0("PGS", "1"), "SEX", "AGE", pcs), collapse = "+")))
          nMod <- neModel(f, family = binomial("logit"), expData=expData, se="robust")
          #this does not like my variable names, I think
          #cf <- coef(summary(neEffdecomp(nMod)))
          lht <- as.data.frame(coef(summary(neLht(nMod, linfct = c("PGS0 = 0",  "PGS1 = 0", "PGS0 + PGS1 = 0")))))
          
          results$OR_DIRECT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- exp(lht$Estimate)[1]
          results$SE_DIRECT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- lht$`Std. Error`[1]
          results$P_DIRECT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- lht$`Pr(>|z|)`[1]
          
          results$OR_INDIRECT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- exp(lht$Estimate)[2]
          results$SE_INDIRECT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- lht$`Std. Error`[2]
          results$P_INDIRECT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- lht$`Pr(>|z|)`[2]
          
          results$OR_TOTAL[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- exp(lht$Estimate)[3]
          results$SE_TOTAL[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- lht$`Std. Error`[3]
          results$P_TOTAL[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- lht$`Pr(>|z|)`[3]
        }
        
      }
    }
  } 
  
  data.table::fwrite(results, "./results/EXSCEL_clinical_outcomes_mediation.tsv",
                     sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE, na="NA")
  
  remove(results)
  gc()
}else{
  print("Skipping clinical outcomes")
}



#do the same analysis for events
#read in file with time to event data
dat <- data.table::fread("./data/rct_data/EXSCEL_med_events.gz")
traits <- unique(dat$PARAM)

#build analysis file
dat <- merge(dat, pca, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, anc, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, covar, by=c("EID", "ID"), all=FALSE, sort=FALSE)
dat <- merge(dat, prs.dat, by="ID", all=FALSE, sort=FALSE)

#add age2
dat$AGE2 <- dat$AGE^2

#regress PCs from prs
pcs <- paste0("PC", 1:10)
for(prs in pgs_models){
  f <- as.formula(paste(prs, "~", paste(pcs, collapse = "+")))
  fit <- lm(f, dat)
  dat[[prs]] <- fit$residuals
}

#add proteins to data frame
dat <- merge(dat, p, by="ID", all=FALSE, sort=FALSE)

#prepare to run analyses
results <- expand.grid(TRAIT=traits, PGS=pgs_models, PROTEIN=proteins,
                       SUBSET=c("ALL", "EUR"),
                       BETA_PGS=NA, SE_PGS=NA, P_PGS=NA,
                       OR_PGS_TRAIT=NA, SE_PGS_TRAIT=NA, P_PGS_TRAIT=NA,
                       OR=NA, SE=NA, P=NA,
                       OR_DIRECT=NA, SE_DIRECT=NA, P_DIRECT=NA,
                       OR_INDIRECT=NA, SE_INDIRECT=NA, P_INDIRECT=NA,
                       OR_TOTAL=NA, SE_TOTAL=NA, P_TOTAL=NA)

#start iteration on trait
for(trait in traits){
  #protein
  for(protein in proteins){
    
    #subset
    for(pop in c("ALL", "EUR")){
      if(pop == "ALL"){
        foo <- as.data.frame(dat[dat$PARAM == trait,])
      }else{
        foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9 & dat$PARAM == trait,])
      }
      
      #exclude insulin
      if(trait == "INSULIN") foo <- foo[foo$INSTFL == "N",]
      if(trait == "RENAL") foo <- foo[foo$EGFR_BL >= 30, ]
      
      #test protein for association with a trait
      pcs <- paste0("PC", 1:10)
      f <- as.formula(paste("EVENT~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX","TRT", pcs, "PROTEIN"), collapse = "+")))
      
      #just include necessary columns
      foo <- foo[c("AGE","AGE2", "SEX","TRT", "EVENT", pcs, pgs_models, protein)]
      colnames(foo)[ncol(foo)] <- "PROTEIN"
      
      #fit logistic model
      fit <- glm(f, foo, family = "binomial")
      fit_sum <- as.data.frame(summary(fit)$coefficients)
      results$OR[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein] <- exp(fit_sum$Estimate[rownames(fit_sum) == "PROTEIN"])
      results$SE[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein] <- fit_sum$`Std. Error`[rownames(fit_sum) == "PROTEIN"]
      results$P[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PROTEIN"]
      
      #skip if not nominally signiificant
      if(fit_sum$`Pr(>|z|)`[rownames(fit_sum) == "PROTEIN"] > 0.05) next
      #pgs in pgs models
      for(pgs in pgs_models){
        #scale PGS
        foo[[pgs]] <- scale(foo[[pgs]])
        
        pcs <- paste0("PC", 1:10)
        f <- as.formula(paste("PROTEIN~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX", pcs, pgs), collapse = "+")))
        #fit pgs-npx model
        fit <- lm(f, foo)
        fit_sum <- as.data.frame(summary(fit)$coefficients)
        results$BETA_PGS[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- fit_sum$Estimate[rownames(fit_sum) == pgs]
        results$SE_PGS[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- fit_sum$`Std. Error`[rownames(fit_sum) == pgs]
        results$P_PGS[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- fit_sum$`Pr(>|t|)`[rownames(fit_sum) == pgs]
        if(fit_sum$`Pr(>|t|)`[rownames(fit_sum) == pgs] > 0.05) next
        
        
        #test if pgs is associated with outcome
        f <- as.formula(paste("EVENT~", paste(c("AGE", "AGE2", "AGE*SEX", "AGE2*SEX", pcs, pgs), collapse = "+")))
        fit <- glm(f, family = "binomial", data=foo)
        
        fit_sum <- as.data.frame(summary(fit)$coefficients)
        results$OR_PGS_TRAIT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- exp(fit_sum$Estimate[rownames(fit_sum) == pgs])
        results$SE_PGS_TRAIT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- fit_sum$`Std. Error`[rownames(fit_sum) == pgs]
        results$P_PGS_TRAIT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- fit_sum$`Pr(>|z|)`[rownames(fit_sum) == pgs]
        
        if(fit_sum$`Pr(>|z|)`[rownames(fit_sum) == pgs] > 0.05) next
        
        print(paste("running mediation for", pgs, "and", trait, "and", protein))
        #need to scale age for mediation
        foo$AGE <- scale(foo$AGE)
        foo$AGE <- as.numeric(foo$AGE)
        
        #scale PCs
        for(pc in pcs){
          foo[[pc]] <- as.numeric(scale(foo[[pc]]))
        }
        
        #change pgs column
        colnames(foo)[colnames(foo) == pgs] <- "PGS"
        foo$PGS <- as.numeric(foo$PGS)
        
        #remove missing subjects
        foo <- foo[complete.cases(foo),]
        
        #run mediation analysis
        f <- as.formula(paste("EVENT~", paste(c("PGS", "PROTEIN", "SEX", "AGE", pcs), collapse = "+")))
        expData <- neImpute(f, family = binomial("logit"), data=foo)
        f <- as.formula(paste("EVENT~", paste(c(paste0("PGS", "0"), paste0("PGS", "1"), "SEX", "AGE", pcs), collapse = "+")))
        nMod <- neModel(f, family = binomial("logit"), expData=expData, se="robust")
        #this does not like my variable names, I think
        #cf <- coef(summary(neEffdecomp(nMod)))
        lht <- as.data.frame(coef(summary(neLht(nMod, linfct = c("PGS0 = 0",  "PGS1 = 0", "PGS0 + PGS1 = 0")))))
        
        results$OR_DIRECT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- exp(lht$Estimate)[1]
        results$SE_DIRECT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- lht$`Std. Error`[1]
        results$P_DIRECT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- lht$`Pr(>|z|)`[1]
        
        results$OR_INDIRECT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- exp(lht$Estimate)[2]
        results$SE_INDIRECT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- lht$`Std. Error`[2]
        results$P_INDIRECT[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- lht$`Pr(>|z|)`[2]
        
        results$OR_TOTAL[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- exp(lht$Estimate)[3]
        results$SE_TOTAL[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- lht$`Std. Error`[3]
        results$P_TOTAL[results$TRAIT == trait & results$SUBSET == pop & results$PROTEIN == protein & results$PGS == pgs] <- lht$`Pr(>|z|)`[3]
      }
      
    }
  }
} 


data.table::fwrite(results, "./results/EXSCEL_clinical_events_mediation.tsv",
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE, na="NA")
