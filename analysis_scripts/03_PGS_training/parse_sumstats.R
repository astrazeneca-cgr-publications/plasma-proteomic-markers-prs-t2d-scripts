# script for preparing GWAS sumstats for PRS estimation #
library(argparser)

########set argument parser##############
argp <- arg_parser("parse GWAS sumstats for PRS-CS, LDSC, or PRSice2")
##inputs
argp <- add_argument(argp, "--prefix", help="prefix for output", 
                     type = "character")
argp <- add_argument(argp, "--sumstats", help="GWAS summary stats file", 
                     type="character")
argp <- add_argument(argp, "--snp", help="snp column", 
                     type="character")
argp <- add_argument(argp, "--bp", help="snp basepair column", 
                     type="character")
argp <- add_argument(argp, "--chrom", help="chromosome column", 
                     type="character")
argp <- add_argument(argp, "--a1", help="A1 column", 
                     type="character")
argp <- add_argument(argp, "--a2", help="A2 column", type="character")
argp <- add_argument(argp, "--pval", help="Pval column", 
                     type="character")
argp <- add_argument(argp, "--beta", help="effect size column", 
                     type="character")
argp <- add_argument(argp, "--OR", help="Are effect sizes odds ratios? Default FALSE", 
                     type="logical", default=FALSE)
argp <- add_argument(argp, "--N", help="Sample size column or overall sample size", 
                     type="character", default=NA)
argp <- add_argument(argp, "--se", help="standard error column", 
                     type="character", default=NA)
argp <- add_argument(argp, "--format", help="Desired format: prscs, ldsc, or prsice2", 
                     type="character", default="prscs")
argp <- add_argument(argp, "--match", help="file path to hm3 variants provided by PRS-CS", 
                     type="character", default=NA)
#parse arguments
argsv <- parse_args(argp)

#########
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

############### begin parsing ##################
#read in sumstats
if(!file.exists(argsv$sumstats)) stop("Specifify path to sumstats file")

print("reading in sumstats")
sumstats <- data.table::fread(argsv$sumstats)

#sort by chromosome
sumstats <- sumstats[order(sumstats[[argsv$chrom]], sumstats[[argsv$bp]]),]

#subset just to set of snps
if(!is.na(argsv$match)){
  print("subsetting to hm3 variants")
  if(!file.exists(argsv$match)) stop("Specifify path to prs-cs hm3 file")
  
  hm3 <- data.table::fread(argsv$match)
  
  #try to match first by snp
  if(nrow(sumstats[sumstats[[argsv$snp]] %in% hm3$SNP,]) > 0.70*nrow(hm3)){
    sumstats <- sumstats[sumstats[[argsv$snp]] %in% hm3$SNP,]
  }else{
    
     #look up data.table method for merging and replace this
    sumstats <- as.data.frame(sumstats)
    hm3 <- as.data.frame(hm3)
    colnames(hm3) <- paste0(colnames(hm3), ".hm3")

    sumstats[[argsv$chrom]] <- as.numeric(sumstats[[argsv$chrom]])
    sumstats <- sumstats[!is.na(sumstats[[argsv$chrom]]),]

    sumstats <- merge(sumstats, hm3, by.x=c(argsv$chrom, argsv$bp),
                      by.y=c("CHR.hm3", "BP.hm3"), all=FALSE, sort=FALSE)
    
    
    sumstats$CHECK <- mapply(checkAlleles, sumstats[[argsv$a2]], sumstats[[argsv$a1]], 
                             sumstats$A2.hm3, sumstats$A1.hm3)
    sumstats <- sumstats[sumstats$CHECK != "MISMATCH",]
    sumstats[[argsv$snp]] <- sumstats$SNP.hm3 
    
  }

  
  if(nrow(sumstats) == 0) stop("No snps remainig; please provide hm3 variant list provided by PRS-CS")
}

print("done reading in data")
#convert odds ratios to betas, if needed
if(argsv$OR) sumstats[[argsv$beta]] <- log(sumstats[[argsv$beta]])

#determine if N is single value or clolumn
if(argsv$format == "prscs" | argsv$format == "ldsc"){
  if(is.na(argsv$N)) stop("provide column for N or overall sample size")
  
  if(argsv$N %in% colnames(sumstats)){
    N <- as.numeric(sumstats[[argsv$N]])
  }else{
    N <- as.numeric(argsv$N)
  }
}

#parse
print("parsing")
if(argsv$format == "prscs"){
  dat <- data.frame(SNP=sumstats[[argsv$snp]],
                    A1=toupper(sumstats[[argsv$a1]]),
                    A2=toupper(sumstats[[argsv$a2]]),
                    BETA=sumstats[[argsv$beta]],
                    P=sumstats[[argsv$pval]],
                    stringsAsFactors = FALSE)
  print(paste("sample size:", mean(N)))
  write.table(mean(N), paste0(argsv$prefix, ".N.txt"), 
              sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)
}else if(argsv$format == "ldsc"){
  if(is.na(argsv$se)) stop("Provide standard error col name")
  Z <- sumstats[[argsv$beta]]/sumstats[[argsv$se]]
  dat <- data.frame(SNP=sumstats[[argsv$snp]],
                    Z=Z,
                    N=N,
                    A1=toupper(sumstats[[argsv$a1]]),
                    A2=toupper(sumstats[[argsv$a2]]),
                    stringsAsFactors = FALSE)
}else if(argsv$format == "prsice2"){
  if(is.na(argsv$chrom)) stop("Provide chrom col name")
  if(is.na(argsv$bp)) stop("Provide basepair col name")
  dat <- data.frame(SNP=sumstats[[argsv$snp]],
                    CHROM=sumstats[[argsv$chrom]],
                    BP=sumstats[[argsv$bp]],
                    A1=toupper(sumstats[[argsv$a1]]),
                    A2=toupper(sumstats[[argsv$a2]]),
                    BETA=sumstats[[argsv$beta]],
                    P=sumstats[[argsv$pval]],
                    stringsAsFactors = FALSE)
  
}

#write-out
print("writing out")

data.table::fwrite(dat, paste0(argsv$prefix, ".", argsv$format, ".parsed.txt"),
                   sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)


print("done")
