# script for making a bim file that is needed for PRS-CS #
library(argparser)

#############set argument parser #####################
argp <- arg_parser("make bim file from GWAS sumstats")
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
argp <- add_argument(argp, "--a2", help="A2 column", 
                     type="character")
argp <- add_argument(argp, "--match", 
                     help="file path to hm3 variants provided by PRS-CS", 
                     type="character", default=NA)

#parse arguments
argsv <- parse_args(argp)

#######

########### helper functions ###################
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

##########
####### begin parsing ##########
#read in sumstats
if(!file.exists(argsv$sumstats)) stop("Specifify path to sumstats file")

print("reading in sumstats")
sumstats <- data.table::fread(argsv$sumstats)

#sort by chromosome
sumstats <- sumstats[order(sumstats[[argsv$chrom]], sumstats[[argsv$bp]]),]

#subset just to set of snps used by prs-cs (hm3)
if(!is.na(argsv$match)){
  print("subsetting to hm3 variants")
  if(!file.exists(argsv$match)) stop("Specifify path to prs-cs hm3 file")
  
  hm3 <- data.table::fread(argsv$match)
  
  #try to match first by snp
  if(nrow(sumstats[sumstats[[argsv$snp]] %in% hm3$SNP,]) > 0.70*nrow(hm3)){
    sumstats <- sumstats[sumstats[[argsv$snp]] %in% hm3$SNP,]
  }else{
    
    #merge by chrom and bp instead
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


#parse
print("parsing")

dat <- data.frame(CHROM=sumstats[[argsv$chrom]],
                  SNP=sumstats[[argsv$snp]],
                  CM=0,
                  BP=sumstats[[argsv$bp]],
                  A1=toupper(sumstats[[argsv$a1]]),
                  A2=toupper(sumstats[[argsv$a2]]),
                  stringsAsFactors = FALSE)


#write-out
print("writing out")

data.table::fwrite(dat, paste0(argsv$prefix, ".bim"),
                   sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)


print("done")