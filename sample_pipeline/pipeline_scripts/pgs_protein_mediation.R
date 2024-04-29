#script for running association analysis between PGS and proteins#
#if (!require("pacman")) install.packages("pacman")
#pacman::p_load(package1, package2, package_n) #could also use pacman 

if (!require("argparser")) install.packages("argparser")
require(argparser)

if (!require("data.table")) install.packages("data.table")
if (!require("medflex")) install.packages("medflex")

########set argument parser##############
argp <- arg_parser("run mediation analysis with PGS as direct exposure and protein as indirect/mediator")

##inputs
#files
argp <- add_argument(argp, "--pheno_file", help="file with phenotype and covariate data", 
                     type = "character", default="./test_data/simulated_pheno_data.txt")
argp <- add_argument(argp, "--protein_file", help="file with protein data", 
                     type = "character", default="./test_data/simulated_protein_data.txt")
#column names
argp <- add_argument(argp, "--id", help="name of ID column, which needs to be the same across files", 
                     type="character", default = "ID")
argp <- add_argument(argp, "--prev", help="name of phenotype column to exclude prevalent cases", 
                     type="character", default = "PHENO")
argp <- add_argument(argp, "--incid", help="name of incident cases to add back to analysis", 
                     type="character", default = "INCIDENT")
argp <- add_argument(argp, "--qcovars", help="column names for quantitative covars, seperated by commas", 
                     type="character", default = "AGE,PC1,PC2")
argp <- add_argument(argp, "--cat", help="column names for categorical covars, seperated by commas", 
                     type="character", default = "SEX")
argp <- add_argument(argp, "--sensitivity", help="column name of quantitative covariate to run with and without, e.g. T2D + BMI", 
                     type="character", default = "BMI")
argp <- add_argument(argp, "--pgs", help="column name for polygenic score", 
                     type="character", default = "PRS")
argp <- add_argument(argp, "--protein", help="column name for protein", 
                     type="character", default = "PROTEIN")

#output file
argp <- add_argument(argp, "--output", help="file name for output", 
                     type="character")
argp <- add_argument(argp, "--wide", help="wide or long format", 
                     type="logical", default=TRUE)
#parse arguments
argsv <- parse_args(argp)

######### read in data ###########
pheno <- data.table::fread(argsv$pheno_file, data.table = FALSE)
proteins <- data.table::fread(argsv$protein_file, data.table = FALSE)

#prep data#
argsv$qcovars <- unlist(unname(strsplit(argsv$qcovars, split=",")))
argsv$cat <- unlist(unname(strsplit(argsv$cat, split=",")))

proteins <- proteins[c(argsv$id, argsv$protein)]

dat <- merge(pheno, proteins, by=argsv$id, all=FALSE, sort=FALSE)
dat <- merge(dat, dosage, by=argsv$id, all=FALSE, sort=FALSE)

#scale quantitative covarites
for(q in c(argsv$qcovars, argsv$sensitivity, argsv$pgs, argsv$protein)){
  dat[[q]] <- scale(dat[[q]])
}
#set categorical covariates to factors
for(covar in argsv$cat){
  dat[[covar]] <- as.factor(dat[[covar]])
}


#check if genetic PCs are included. If they are, adjust PGS for the PCs to account for pop structure
pcs <- colnames(dat)[grep("PC", colnames(dat))]

if(length(pcs) > 0){
  f <- as.formula(paste(argsv$pgs, "~", paste(pcs, collapse = "+")))
  fit <- lm(f, data=dat)
  dat[[argsv$pgs]][!is.na(dat[[argsv$pgs]])] <- fit$residuals
}

####exclude prevalent cases, but keep incident cases ######
dat <- dat[!(dat[[argsv$prev]] == 1 & dat[[argsv$incid]] == 0),]

####### prepare results data.frame #########
results <- expand.grid(PROTEIN=argsv$protein,
                      PGS=argsv$pgs,
                      TRAIT=argsv$prev, N_INCIDENT=NA,
                      SENSITIVITY=c(FALSE, TRUE),
                      SENS_COVAR=argsv$sensitivity,
                      BETA_PGS=NA, SE_PGS=NA, P_PGS=NA,
                      BETA_PROTEIN=NA, SE_PROTEIN=NA, P_PROTEIN=NA,
                      EST_DIRECT=NA, SE_DIRECT=NA, P_DIRECT=NA,
                      EST_INDIRECT=NA, SE_INDIRECT=NA, P_INDIRECT=NA,
                      EST_TOTAL=NA, SE_TOTAL=NA, P_TOTAL=NA,
                      stringsAsFactors = FALSE)


#first run with out BMI adjustment (or any other covariate run in a +/- manner)

f1 <- as.formula(paste(argsv$incid, "~", 
                       paste(argsv$qcovars, collapse = "+"),
                       "+",
                       paste(argsv$cat, collapse = "+"),
                       "+", argsv$pgs))

f2 <- as.formula(paste(argsv$incid, "~", 
                       paste(argsv$qcovars, collapse = "+"),
                       "+",
                       paste(argsv$cat, collapse = "+"),
                       "+", argsv$protein))


#fit the two models
fit <- glm(f1, data=dat, family="binomial")
fit2 <- glm(f2, data=dat, family = "binomial")

#sample size
results$N_INCIDENT <- nrow(dat[dat[[argsv$incid]] == 1,])

#parse the fitted models
#first, parse the model testing PGS association with incident trait
p <- as.data.frame(summary(fit)$coefficients)
results$BETA_PGS[1] <- p$Estimate[nrow(p)]
results$SE_PGS[1] <- p$`Std. Error`[nrow(p)]
results$P_PGS[1] <- p$`Pr(>|z|)`[nrow(p)]

#then, parse the model testing PGS association with incident trait
p <- as.data.frame(summary(fit2)$coefficients)
results$BETA_PROTEIN[1] <- p$Estimate[nrow(p)]
results$SE_PROTEIN[1] <- p$`Std. Error`[nrow(p)]
results$P_PROTEIN[1] <- p$`Pr(>|z|)`[nrow(p)]


#only run mediation if p-values are < 0.05
if(results$P_PGS[1] < 0.05 & results$P_PROTEIN[1] < 0.05){
  
  f <- as.formula(paste(argsv$incid, "~",
                         argsv$pgs, "+",
                         argsv$protein, "+",
                         paste(argsv$qcovars, collapse = "+"),
                         "+",
                         paste(argsv$cat, collapse = "+")))
  dat[[argsv$pgs]] <- as.numeric(dat[[argsv$pgs]])
  expData <- neImpute(f, family = binomial("logit"), data=dat)
  
  f <- as.formula(paste(argsv$incid, "~",
                        paste0(argsv$pgs, 0), "+",
                        paste0(argsv$pgs, 1), "+",
                        paste(argsv$qcovars, collapse = "+"),
                        "+",
                        paste(argsv$cat, collapse = "+")))
  
  nMod <- neModel(f, family = binomial("logit"), expData=expData, se="robust")
  
  lht <- as.data.frame(coef(summary(neLht(nMod, linfct = c(paste0(argsv$pgs, "0 = 0"),  
                                                           paste0(argsv$pgs, "1 = 0"), 
                                                           paste0(argsv$pgs, "0 + ", argsv$pgs, "1 = 0")
                                                           )))))
  
  results$EST_DIRECT[1] <- lht$Estimate[1]
  results$SE_DIRECT[1] <- lht$`Std. Error`[1]
  results$P_DIRECT[1] <- lht$`Pr(>|z|)`[1]
  
  results$EST_INDIRECT[1] <- lht$Estimate[2]
  results$SE_INDIRECT[1] <- lht$`Std. Error`[2]
  results$P_INDIRECT[1] <- lht$`Pr(>|z|)`[2]
  
  results$EST_TOTAL[1] <- lht$Estimate[3]
  results$SE_TOTAL[1] <- lht$`Std. Error`[3]
  results$P_TOTAL[1] <- lht$`Pr(>|z|)`[3]
  
}

#repeat the above adjusting for the additional covariate

f1 <- as.formula(paste(argsv$incid, "~", 
                       paste(argsv$qcovars, collapse = "+"),
                       "+",
                       paste(argsv$sensitivity), "+",
                       paste(argsv$cat, collapse = "+"),
                       "+", argsv$pgs))

f2 <- as.formula(paste(argsv$incid, "~", 
                       paste(argsv$qcovars, collapse = "+"),
                       "+",
                       paste(argsv$sensitivity), "+",
                       paste(argsv$cat, collapse = "+"),
                       "+", argsv$protein))


#fit the two models
fit <- glm(f1, data=dat, family="binomial")
fit2 <- glm(f2, data=dat, family = "binomial")

#parse the fitted models
#first, parse the model testing PGS association with incident trait
p <- as.data.frame(summary(fit)$coefficients)
results$BETA_PGS[2] <- p$Estimate[nrow(p)]
results$SE_PGS[2] <- p$`Std. Error`[nrow(p)]
results$P_PGS[2] <- p$`Pr(>|z|)`[nrow(p)]

#then, parse the model testing PGS association with incident trait
p <- as.data.frame(summary(fit2)$coefficients)
results$BETA_PROTEIN[2] <- p$Estimate[nrow(p)]
results$SE_PROTEIN[2] <- p$`Std. Error`[nrow(p)]
results$P_PROTEIN[2] <- p$`Pr(>|z|)`[nrow(p)]


#only run mediation if p-values are < 0.05
if(results$P_PGS[1] < 0.05 & results$P_PROTEIN[1] < 0.05){
  
  f <- as.formula(paste(argsv$incid, "~",
                        argsv$pgs, "+",
                        argsv$protein, "+",
                        paste(argsv$qcovars, collapse = "+"),
                        "+",
                        paste(argsv$sensitivity), "+",
                        paste(argsv$cat, collapse = "+")))
  dat[[argsv$pgs]] <- as.numeric(dat[[argsv$pgs]])
  expData <- neImpute(f, family = binomial("logit"), data=dat)
  
  f <- as.formula(paste(argsv$incid, "~",
                        paste0(argsv$pgs, 0), "+",
                        paste0(argsv$pgs, 1), "+",
                        paste(argsv$qcovars, collapse = "+"),
                        "+",
                        paste(argsv$sensitivity), "+",
                        paste(argsv$cat, collapse = "+")))
  
  nMod <- neModel(f, family = binomial("logit"), expData=expData, se="robust")
  
  lht <- as.data.frame(coef(summary(neLht(nMod, linfct = c(paste0(argsv$pgs, "0 = 0"),  
                                                           paste0(argsv$pgs, "1 = 0"), 
                                                           paste0(argsv$pgs, "0 + ", argsv$pgs, "1 = 0")
  )))))
  
  results$EST_DIRECT[2] <- lht$Estimate[1]
  results$SE_DIRECT[2] <- lht$`Std. Error`[1]
  results$P_DIRECT[2] <- lht$`Pr(>|z|)`[1]
  
  results$EST_INDIRECT[2] <- lht$Estimate[2]
  results$SE_INDIRECT[2] <- lht$`Std. Error`[2]
  results$P_INDIRECT[2] <- lht$`Pr(>|z|)`[2]
  
  results$EST_TOTAL[2] <- lht$Estimate[3]
  results$SE_TOTAL[2] <- lht$`Std. Error`[3]
  results$P_TOTAL[2] <- lht$`Pr(>|z|)`[3]
  
}


#done, write out
#create outfile name if none given from pgs and protein names
outfile <- ifelse(is.na(argsv$output), 
                  paste0(argsv$protein, "_", argsv$pgs, "_mediation_results.tsv"),
                  argsv$output)

if(argsv$wide){
  data.table::fwrite(results, outfile, sep='\t')
}else{
  results <- data.frame(TERM=rownames(t(results)),
                        VALUE=unname(t(results)[,1]),
                        VALUE_SENS=unname(t(results)[,2]))
  data.table::fwrite(results, outfile, sep='\t')
}
