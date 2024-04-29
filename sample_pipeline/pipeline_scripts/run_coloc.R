#script for running coloc analysis#
#if (!require("pacman")) install.packages("pacman")
#pacman::p_load(package1, package2, package_n) #could also use pacman 

if (!require("argparser")) install.packages("argparser")
require(argparser)

if (!require("data.table")) install.packages("data.table")

if (!require("coloc")) install.packages("coloc")

########set argument parser##############
argp <- arg_parser("runs colocalization analysis for two traits")

##inputs
#files
argp <- add_argument(argp, "--trait_gwas", help="file with trait gwas data", 
                     type = "character", default="./test_data/TRAIT_GWAS.txt")
argp <- add_argument(argp, "--exposure_gwas", help="file with exposure GWAS, e.g. protein data", 
                     type = "character", default="./test_data/PROTEIN_GWAS.txt")
argp <- add_argument(argp, "--index_snp", help="specify index snp", 
                     type = "character", default="rs505922")


#information for output
argp <- add_argument(argp, "--exposure", help="specify exposure", 
                     type="character", default = "PROTEIN")
argp <- add_argument(argp, "--trait", help="specify trait", 
                     type="character", default = "TRAIT")
argp <- add_argument(argp, "--quantitative_trait", help="state if trait is quantitative", 
                     type = "logical", default=FALSE)

#output file
argp <- add_argument(argp, "--output", help="file name for output", 
                     type="character")
#parse arguments
argsv <- parse_args(argp)

##### read in data ##########
exp_gwas <- data.table::fread(argsv$exposure_gwas)
trait_gwas <- data.table::fread(argsv$trait_gwas)

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


#center on index snp, +/- 200 kb
start <- dat$bp_trait[dat$id == argsv$index_snp] - 2E5
end <- dat$bp_trait[dat$id == argsv$index_snp] + 2E5

dat <- dat[dat$bp_trait >= start & dat$bp_trait <= end,]
dat <- dat[complete.cases(dat),]
if(argsv$quantitative_trait){
  d1 <- list(snp=dat$id, 
             position=dat$bp_trait,
             beta=dat$beta_trait,
             varbeta=dat$se_trait^2,
             type="quant",
             N=dat$N_trait,
             MAF=dat$af_trait)
} else{
  d1 <- list(snp=dat$id, 
             position=dat$bp_trait,
             beta=dat$beta_trait,
             varbeta=dat$se_trait^2,
             type="cc")
}


d2 <- list(snp=dat$id, 
           position=dat$bp_trait, #use trait base position as matched by id
           beta=dat$beta_exp,
           varbeta=dat$se_exp^2,
           type="quant",
           N=dat$N_exp,
           MAF=dat$af_exp)
test <- coloc.abf(d1, d2)

results <- cbind(data.frame(protein=argsv$exposure, index_snp=argsv$index_snp, trait=argsv$trait),
                 as.data.frame(t(test$summary)))
results$SNP.PP.H4 <- test$results$SNP.PP.H4[test$results$snp == argsv$index_snp]

#criteria: PP.H3 + PP.H4 ≥ 0.99 and PP4/PP3 ≥5
criteria1 <- results$PP.H3.abf + results$PP.H4.abf
criteria2 <- results$PP.H4.abf/results$PP.H3.abf
results$SAME_SIGNAL <- ifelse(criteria1 >= 0.99 & criteria2 >=5, TRUE, FALSE)

outfile <- ifelse(is.na(argsv$output), 
                  paste0(argsv$exposure, "_", argsv$trait, "_coloc_results.tsv"),
                  argsv$output)
data.table::fwrite(results, outfile, sep='\t')