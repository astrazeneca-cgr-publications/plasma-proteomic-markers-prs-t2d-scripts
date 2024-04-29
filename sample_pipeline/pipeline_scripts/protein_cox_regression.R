#script for running time to event analysis with proteins as exposure#
#if (!require("pacman")) install.packages("pacman")
#pacman::p_load(package1, package2, package_n) #could also use pacman 

if (!require("argparser")) install.packages("argparser")
require(argparser)

if (!require("data.table")) install.packages("data.table")

if (!require("survival")) install.packages("survival")

########set argument parser##############
argp <- arg_parser("test protein for association with time to event")

##inputs
#files
argp <- add_argument(argp, "--pheno_file", help="file with phenotype and covariate data", 
                     type = "character", default="./test_data/simulated_pheno_data.txt")
argp <- add_argument(argp, "--protein_file", help="file with protein data", 
                     type = "character", default="./test_data/simulated_protein_data.txt")
argp <- add_argument(argp, "--event_file", help="file with time-to-event data", 
                     type = "character", default="./test_data/simulated_time_to_event.txt")
argp <- add_argument(argp, "--exclude_file", 
                     help="file specifing sample exclusions, e.g. already experienced event, medication, want to exclude diseaes, etc.", 
                     type = "character", default="./test_data/simulated_time_to_event_exclusions.txt")
#column names
argp <- add_argument(argp, "--id", help="name of ID column, which needs to be the same across files", 
                     type="character", default = "ID")
argp <- add_argument(argp, "--qcovars", help="column names for quantitative covars, seperated by commas", 
                     type="character", default = "AGE,PC1,PC2")
argp <- add_argument(argp, "--cat", help="column names for categorical covars, seperated by commas", 
                     type="character", default = "SEX")
argp <- add_argument(argp, "--clinical_covars", help="column names for clinical risk factors", 
                     type="character", default = "BMI")
argp <- add_argument(argp, "--biomarkers", help="column names for known biomakers", 
                     type="character", default = NA)
argp <- add_argument(argp, "--protein", help="column name for protein", 
                     type="character", default = "PROTEIN")
argp <- add_argument(argp, "--event", help="column name for event", 
                     type="character", default = "EVENT")
argp <- add_argument(argp, "--time", help="column name for time to event", 
                     type="character", default = "TIME")

#output file
argp <- add_argument(argp, "--output", help="file name for output", 
                     type="character")
argp <- add_argument(argp, "--wide", help="wide or long format", 
                     type="logical", default=TRUE)
#parse arguments
argsv <- parse_args(argp)

######### read in data ###########
pheno <- data.table::fread(argsv$pheno_file, data.table = FALSE)
event <- data.table::fread(argsv$event_file, data.table = FALSE)
proteins <- data.table::fread(argsv$protein_file, data.table = FALSE)
exclude <- data.table::fread(argsv$exclude_file, data.table = FALSE,
                             header = FALSE)

#prep data#
argsv$qcovars <- unlist(unname(strsplit(argsv$qcovars, split=",")))
argsv$cat <- unlist(unname(strsplit(argsv$cat, split=",")))
argsv$clinical_covars <- unlist(unname(strsplit(argsv$clinical_covars, split=",")))
argsv$biomarkers <- unlist(unname(strsplit(argsv$biomarkers, split=",")))

proteins <- proteins[c(argsv$id, argsv$protein)]

dat <- merge(pheno, proteins, by=argsv$id, all=FALSE, sort=FALSE)
dat <- merge(dat, event, by=argsv$id, all=FALSE, sort=FALSE)

#remove exclusions
dat <- dat[!dat$ID %in% exclude$V1,]

#scale quantitative covarites
for(q in c(argsv$qcovars, argsv$biomarkers, argsv$protein)){
  dat[[q]] <- scale(dat[[q]])
}
#set categorical covariates to factors
for(covar in argsv$cat){
  dat[[covar]] <- as.factor(dat[[covar]])
}

#clinical risk factors can be categorical or quantitative
for(covar in argsv$clinical_covars){
  if(class(dat[[covar]]) %in% c("numeric", "integer")){
    
    #sometimes, a categorical factor can appear nummeric, e.g. categores are 1,2,3...
    if(length(unique(dat[[covar]])) < 10){
      dat[[covar]] <- as.factor(dat[[covar]])
    }else{
      dat[[covar]] <- scale(dat[[covar]])
    }
    
  }else{
    dat[[covar]] <- as.factor(dat[[covar]])
  }
}




####### prepare results data.frame #########
results <- data.frame(PROTEIN=argsv$protein,
                      EVENT=argsv$event,
                      N_EXCL=nrow(exclude), 
                      N_EVENT=nrow(dat[dat[[argsv$event]] == 1,]),
                      BETA_M1=NA, HR_M1=NA, SE_M1=NA, P_M1=NA,
                      BETA_M2=NA, HR_M2=NA, SE_M2=NA, P_M2=NA,
                      BETA_M3=NA, HR_M3=NA, SE_M3=NA, P_M3=NA)


sobj <- Surv(dat[[argsv$time]], dat[[argsv$event]], type="right")

foo <- dat
foo$SURV <- sobj

#model 1, only adjust for "standard" covariates
f1 <- as.formula(paste("SURV ~", 
                       paste(argsv$qcovars, collapse = "+"),
                       "+",
                       paste(argsv$cat, collapse = "+"),
                       "+", argsv$protein))
#model 2, also include clinical risk factors
f2 <- as.formula(paste("SURV ~", 
                       paste(argsv$qcovars, collapse = "+"),
                       "+",
                       paste(argsv$cat, collapse = "+"),
                       "+", paste(argsv$clinical_covars, collapse = "+"),
                       "+", argsv$protein))

#model 3, if provided, also include known biomarkers
f3 <- as.formula(paste( "SURV ~", 
                       paste(argsv$qcovars, collapse = "+"),
                       "+",
                       paste(argsv$cat, collapse = "+"),
                       "+", paste(argsv$clinical_covars, collapse = "+"),
                       "+",
                       paste(argsv$biomarkers, collapse = "+"),
                       "+", argsv$protein))

#fit the three models
fit1 <- coxph(f1, data=foo)
fit2 <- coxph(f2, data=foo)
if(!is.na(argsv$biomarkers)){
  fit3 <- coxph(f3, data=foo)
}


#parse the fitted models
p <- as.data.frame(summary(fit1)$coefficients)
results$BETA_M1 <- p$coef[nrow(p)]
results$HR_M1 <- p$`exp(coef)`[nrow(p)]
results$SE_M1 <- p$`se(coef)`[nrow(p)]
results$P_M1 <- p$`Pr(>|z|)`[nrow(p)]


p <- as.data.frame(summary(fit2)$coefficients)
results$BETA_M2 <- p$coef[nrow(p)]
results$HR_M2 <- p$`exp(coef)`[nrow(p)]
results$SE_M2 <- p$`se(coef)`[nrow(p)]
results$P_M2 <- p$`Pr(>|z|)`[nrow(p)]

if(!is.na(argsv$biomarkers)){
  p <- as.data.frame(summary(fit3)$coefficients)
  results$BETA_M3 <- p$coef[nrow(p)]
  results$HR_M3 <- p$`exp(coef)`[nrow(p)]
  results$SE_M3 <- p$`se(coef)`[nrow(p)]
  results$P_M3 <- p$`Pr(>|z|)`[nrow(p)]
}


#done, write out
#create outfile name if none given from pgs and protein names
outfile <- ifelse(is.na(argsv$output), 
                  paste0(argsv$protein, "_", argsv$event, "_cox_results.tsv"),
                  argsv$output)

if(argsv$wide){
  data.table::fwrite(results, outfile, sep='\t')
}else{
  results <- data.frame(TERM=rownames(t(results)),
                        VALUE=unname(t(results)[,1]))
  data.table::fwrite(results, outfile, sep='\t')
}
