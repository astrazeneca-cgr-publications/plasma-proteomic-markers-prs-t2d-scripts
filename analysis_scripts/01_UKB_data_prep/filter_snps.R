#script for filtering variants in ukb plink file for PGS estimating
argsv <- commandArgs(trailingOnly = TRUE)

pvar.file <- argsv[1]
pvar <- data.table::fread(pvar.file)
pvar.cols <- colnames(pvar)
n <- nrow(pvar)
snps <- pvar$ID
#get chrom from pvar data
chr <- unique(pvar$`#CHROM`)

#flag duplicates
#append "_DROP" on duplicated variants and overwrite pvar file
#read in info data
info <- data.table::fread(paste0("../datasets/ukb/Imputed/ukb_c", chr,"_b0_v3.mfi.txt"))
colnames(info) <- c("ID2", "ID", "POS", "REF", "ALT", "FREQ", "A1", "INFO")

pvar <- merge(pvar, info, by=c("ID", "POS", "REF", "ALT"), all=FALSE, sort=FALSE)

#run some checks
#if(nrow(pvar) != n) stop("Error when adding info score to pvar file")
#if(nrow(pvar) != nrow(pvar[pvar$ID == snps,])) stop("Error when adding info score to pvar file")

dupes <- pvar$ID[duplicated(pvar$ID)]
for(dupe in dupes){
  pvar$ID <- ifelse(pvar$ID != dupe, pvar$ID, 
                    ifelse(pvar$ID == dupe & pvar$INFO == max(pvar$INFO[pvar$ID == dupe]), 
                           pvar$ID, paste0(pvar$ID, "_DROP")))
}
#overwrite
pvar.out <- as.data.frame(pvar)
pvar.out <- pvar.out[pvar.cols]

data.table::fwrite(pvar.out, pvar.file, sep='\t', quote=FALSE, 
                   row.names = FALSE, col.names = TRUE)

#flag cgat variants
get_cgat <- function(a1, a2){
  a1 <- toupper(a1)
  a2 <- toupper(a2)
  
  is.cgat <- ifelse(a1 == "C" & a2 == "G", TRUE,
                    ifelse(a1 == "G" & a2 == "C", TRUE,
                           ifelse(a1 == "A" & a2 == "T", TRUE,
                                  ifelse(a1 == "T" & a2 == "A", TRUE, FALSE))))
  return(is.cgat)
}

pvar$CGAT <- mapply(get_cgat, pvar$REF, pvar$ALT)


#write out snp list

snp_list <- pvar$ID[ pvar$INFO >= 0.5 & pvar$CGAT == FALSE]

#write out snp lists
out.prefix <- gsub(".pvar", "", pvar.file)
write.table(unique(snp_list), paste0(out.prefix, ".info_filter_snp_list"),
            sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)