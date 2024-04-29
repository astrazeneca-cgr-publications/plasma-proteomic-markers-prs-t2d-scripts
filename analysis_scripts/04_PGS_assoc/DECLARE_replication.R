#script for running PGS-npx replication in DECLARE
#read in covariate file
covar <- data.table::fread("./clinical_data/COMBINED.MACE.tsv")
covar <- covar[covar$STUDY == "DECLARE",]

#read in PCA
pca <- data.table::fread("./data/rct_data/DECLARE-PCA.eigenvec")

#read in ancestry definitions
anc <- data.table::fread("./data/rct_data/DECLARE_ancestry_InferredAncestry.txt")
anc$PC1 <- NULL
anc$PC2 <- NULL 

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


#build analysis file
dat <- merge(covar, pca, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
dat <- merge(dat, anc, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)
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
p <- data.table::fread("./data/rct_data/DECLARE_OLINK_visit2.tsv",
                       data.table = FALSE)
proteins <- colnames(p)[2:ncol(p)]

#add proteins to data frame
dat <- merge(dat, p, by="EID", all=FALSE, sort=FALSE)

#read in quantiative traits
q <- data.table::fread("./data/rct_data/DECLARE_quant_traits.gz")
q <- q[,-c(2:5)]
dat <- merge(dat, q, by=c("ID", "EID"), all.x=TRUE, sort=FALSE)

results <- expand.grid(PROTEIN=proteins,
                       PGS=pgs_models,
                       SUBSET=c("ALL", "EUR"),
                       N=NA, BETA=NA, SE=NA, P=NA,
                       BETA_BMI=NA, SE_BMI=NA, P_BMI=NA)

###much smaller than UKB-PPP, so will run linearly instead of in parallel
for(protein in proteins){
  for(pgs in pgs_models){
    for(pop in c("ALL", "EUR")){
      if(pop == "ALL"){
        foo <- as.data.frame(dat)
      }else{
        foo <- as.data.frame(dat[dat$Anc_1st == "EUR" & dat$Pr_1st >= 0.9,])
      }
      
      #just include necessary columns
      foo <- foo[c("AGE","AGE2", "SEX", "BMI_BL", pcs, protein, pgs)]
      
      #scale prs
      foo[[pgs]] <- scale(foo[[pgs]])
      
      f <- as.formula(paste(protein, "~", paste(c("AGE", "SEX", "AGE2", "AGE*SEX", "AGE2*SEX", pcs, pgs), collapse = "+")))
      fit <- lm(f, foo)
      
      #parse results
      fit_sum <- as.data.frame(summary(fit)$coefficients)
      
      results$N[results$PGS == pgs & results$SUBSET == pop & results$PROTEIN == protein] <- length(fit$residuals)
      results$BETA[results$PGS == pgs & results$SUBSET == pop  & results$PROTEIN == protein] <- fit_sum$Estimate[rownames(fit_sum) == pgs]
      results$SE[results$PGS == pgs & results$SUBSET == pop  & results$PROTEIN == protein] <- fit_sum$`Std. Error`[rownames(fit_sum) == pgs]
      results$P[results$PGS == pgs & results$SUBSET == pop  & results$PROTEIN == protein] <- fit_sum$`Pr(>|t|)`[rownames(fit_sum) == pgs]
      
      #BMI adjusted model
      f <- as.formula(paste(protein, "~", paste(c("AGE", "SEX", "AGE2", "AGE*SEX", "AGE2*SEX", "BMI_BL", pcs, pgs), collapse = "+")))
      fit <- lm(f, foo)
      
      #parse results
      fit_sum <- as.data.frame(summary(fit)$coefficients)
      
      results$BETA_BMI[results$PGS == pgs & results$SUBSET == pop  & results$PROTEIN == protein] <- fit_sum$Estimate[rownames(fit_sum) == pgs]
      results$SE_BMI[results$PGS == pgs & results$SUBSET == pop  & results$PROTEIN == protein] <- fit_sum$`Std. Error`[rownames(fit_sum) == pgs]
      results$P_BMI[results$PGS == pgs & results$SUBSET == pop  & results$PROTEIN == protein] <- fit_sum$`Pr(>|t|)`[rownames(fit_sum) == pgs]
    }

  }
}

#write out results
data.table::fwrite(results, "./results/DECLARE_pgs_npx_replication.tsv",
                   sep='\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
