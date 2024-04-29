#script for running association analysis between PGS and proteins#
#if (!require("pacman")) install.packages("pacman")
#pacman::p_load(package1, package2, package_n) #could also use pacman 

if (!require("argparser")) install.packages("argparser")
require(argparser)

if (!require("data.table")) install.packages("data.table")

########set argument parser##############
argp <- arg_parser("test PGS for association with protein levels")

##inputs
#files
argp <- add_argument(argp, "--pheno_file", help="file with phenotype and covariate data", 
                     type = "character", default="./test_data/simulated_pheno_data.txt")
argp <- add_argument(argp, "--protein_file", help="file with protein data", 
                     type = "character", default="./test_data/simulated_protein_data.txt")
argp <- add_argument(argp, "--dosage_file", help="file with genotype dosages of pQTLs", 
                     type = "character", default="./test_data/simulated_dosage.txt")
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
dosage <- data.table::fread(argsv$dosage_file, data.table = FALSE)

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
results <- data.frame(PROTEIN=argsv$protein,
                      PGS=argsv$pgs,
                      EXCL=argsv$prev, N=NA,
                      BETA=NA, SE=NA, P=NA, R2=NA,
                      BETA_SENS=NA, SE_SENS=NA, P_SENS=NA,
                      PQTL=NA, BETA_PQTL=NA, SE_PQTL=NA, P_PQTL=NA)

b <- as.formula(paste(argsv$protein, "~", 
                      paste(argsv$qcovars, collapse = "+"),
                      "+",
                      paste(argsv$cat, collapse = "+")))
f1 <- as.formula(paste(argsv$protein, "~", 
                      paste(argsv$qcovars, collapse = "+"),
                      "+",
                      paste(argsv$cat, collapse = "+"),
                      "+", argsv$pgs))

f2 <- as.formula(paste(argsv$protein, "~", 
                       paste(argsv$qcovars, collapse = "+"),
                       "+",
                       paste(argsv$cat, collapse = "+"),
                       "+", argsv$sensitivity,
                       "+", argsv$pgs))


#fit the three models
base <- lm(b, data=dat)
fit <- lm(f1, data=dat)
fit2 <- lm(f2, data=dat)

#sample size
results$N <- length(fit$residuals)

#get r2 of pgs
results$R2 <- summary(fit)$r.squared - summary(base)$r.squared

#parse the fitted models
p <- as.data.frame(summary(fit)$coefficients)
results$BETA <- p$Estimate[nrow(p)]
results$SE <- p$`Std. Error`[nrow(p)]
results$P <- p$`Pr(>|t|)`[nrow(p)]


p <- as.data.frame(summary(fit2)$coefficients)
results$BETA_SENS <- p$Estimate[nrow(p)]
results$SE_SENS <- p$`Std. Error`[nrow(p)]
results$P_SENS <- p$`Pr(>|t|)`[nrow(p)]


#now, test each pqtl, just save pvalues initially
pqtls <- colnames(dosage)[colnames(dosage) != argsv$id]
pvals <- c()
for(pqtl in pqtls){
  f1 <- as.formula(paste(argsv$protein, "~", 
                         paste(argsv$qcovars, collapse = "+"),
                         "+",
                         paste(argsv$cat, collapse = "+"),
                         "+", pqtl,
                         "+", argsv$pgs))
  fit <- lm(f1, dat)
  p <- as.data.frame(summary(fit)$coefficients)
  pvals <- c(pvals, p$`Pr(>|t|)`[nrow(p)])
}

#just keep pqtl with largest impact, e.g. worst pvalue for pgs
pqtl <- pqtls[pvals == max(pvals)][1] # just keep first one in case of ties

results$PQTL <- pqtl

#refit
f1 <- as.formula(paste(argsv$protein, "~", 
                       paste(argsv$qcovars, collapse = "+"),
                       "+",
                       paste(argsv$cat, collapse = "+"),
                       "+", pqtl,
                       "+", argsv$pgs))
fit <- lm(f1, dat)
p <- as.data.frame(summary(fit)$coefficients)
results$BETA_PQTL <- p$Estimate[nrow(p)]
results$SE_PQTL <- p$`Std. Error`[nrow(p)]
results$P_PQTL <- p$`Pr(>|t|)`[nrow(p)]

#done, write out
#create outfile name if none given from pgs and protein names
outfile <- ifelse(is.na(argsv$output), 
              paste0(argsv$protein, "_", argsv$pgs, "_assoc_results.tsv"),
              argsv$output)

if(argsv$wide){
  data.table::fwrite(results, outfile, sep='\t')
}else{
  results <- data.frame(TERM=rownames(t(results)),
                        VALUE=unname(t(results)[,1]))
  data.table::fwrite(results, outfile, sep='\t')
}

