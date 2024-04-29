#script for running MR analysis#
#if (!require("pacman")) install.packages("pacman")
#pacman::p_load(package1, package2, package_n) #could also use pacman 

if (!require("argparser")) install.packages("argparser")
require(argparser)

if (!require("data.table")) install.packages("data.table")

if (!require("MendelianRandomization")) install.packages("MendelianRandomization")

########set argument parser##############
argp <- arg_parser("runs MR analysis testing if protein isa causal for trait")

##inputs
#files
argp <- add_argument(argp, "--trait_gwas", help="file with trait gwas data", 
                     type = "character", default="./test_data/TRAIT_GWAS.txt")
argp <- add_argument(argp, "--exposure_gwas", help="file with exposure GWAS, e.g. protein data", 
                     type = "character", default="./test_data/PROTEIN_GWAS.txt")
argp <- add_argument(argp, "--instruments", help="file specifying which instruments to use", 
                     type = "character", default="./test_data/pqtl_list")

#optional filtering step
argp <- add_argument(argp, "--fstat", help = "filter using F-statistic",
                     type="logical", default=FALSE)

#information for output
argp <- add_argument(argp, "--exposure", help="specify exposure", 
                     type="character", default = "PROTEIN")
argp <- add_argument(argp, "--trait", help="specify trait", 
                     type="character", default = "TRAIT")


#output file
argp <- add_argument(argp, "--output", help="file name for output", 
                     type="character")
#parse arguments
argsv <- parse_args(argp)


##### read in data ##########
exp_gwas <- data.table::fread(argsv$exposure_gwas)
trait_gwas <- data.table::fread(argsv$trait_gwas)
pqtls <- data.table::fread(argsv$instruments, header=FALSE)

####### helper functions ######

#flip alleles
flip <- function(x){
  y <- ifelse(x == "A","T",
              ifelse(x == "T", "A",
                     ifelse(x == "C", "G",
                            ifelse(x == "G", "C", "NA"))))
  if(y == "NA") warning("variant not a SNP")
  return(y)
}

#check alleles
checkAlleles <- function(ref1, alt1, ref2, alt2){
  if(ref1 == ref2 & alt1 == alt2){
    y <- "SAME"
  }else if(ref1 == alt2 & alt1 == ref2){
    y <- "SWITCH"
  }else if(ref1 == flip(ref2) & alt1 == flip(alt2)){
    y <- "FLIP"
  }else if(ref1 == flip(alt2) & alt1 == flip(ref2)){
    y <- "SWITCH_FLIP"
  }else{
    y <- "MISMATCH"
  }
  return(y)
}
###########
#prepare data for merging
colnames(exp_gwas)[2:ncol(exp_gwas)] <- paste0(colnames(exp_gwas)[2:ncol(exp_gwas)], "_exp")
colnames(trait_gwas)[2:ncol(trait_gwas)] <- paste0(colnames(trait_gwas)[2:ncol(trait_gwas)], "_trait")

dat <- merge(exp_gwas, trait_gwas, by="id", all=FALSE, sort=FALSE)

#check alleles
dat$CHECK <- mapply(checkAlleles, dat$effect_allele_exp, dat$other_allele_exp,
                    dat$effect_allele_trait, dat$other_allele_trait, 
                    USE.NAMES = FALSE)

#harmonize
dat <- dat[dat$CHECK != "MISMATCH",]
dat$beta_trait <- ifelse(dat$CHECK %in% c("SWITCH", "SWITCH_FLIP"),
                         dat$beta_trait*-1, dat$beta_trait)

#filter by fstatistic if needed
if(argsv$fstat){
  dat$fstat <- dat$beta_exp^2/dat$se_exp^2
  dat <- dat[dat$fstat >= 10,]
}


#subset to just user-provided instruments
dat <- dat[dat$id %in% pqtls$V1,]

#exit if there are not 3+ pQTLs
if(nrow(dat) < 3){
  stop("Less than the required 3 pQTLs")
}

dat.input <- mr_input(snps=dat$id,
                      bx= dat$beta_exp,
                      bxse=dat$se_exp,
                      by=dat$beta_trait,
                      byse = dat$se_trait,
                      exposure=argsv$exposure,
                      outcome=argsv$trait)
#run analysis
tryCatch({ mr.run <- mr_allmethods(dat.input, iterations=10000, method = "main")},
         warning=function(w) print(paste("Warning for", argsv$exposure, "and", argsv$trait)))

#parse
results <- as.data.frame(mr.run@Values)
colnames(results) <- c("Method", "Estimate", "SE", "Lower_95", "Upper_95", "Pvalue")


#add ensemble estimates by taking the median of all methods (but not intercept)
results <- rbind(results,
                 data.frame(Method="Ensemble",
                            Estimate = median(results$Estimate[1:4]),
                            SE=median(results$SE[1:4]),
                            Lower_95 = median(results$Lower_95[1:4]),
                            Upper_95 = median(results$Upper_95[1:4]),
                            Pvalue = median(results$Pvalue[1:4]))
                 )
results$Exposure <- argsv$exposure
results$Trait <- argsv$trait

outfile <- ifelse(is.na(argsv$output), 
                  paste0(argsv$exposure, "_", argsv$trait, "_MR_results.tsv"),
                  argsv$output)
data.table::fwrite(results, outfile, sep='\t')

  