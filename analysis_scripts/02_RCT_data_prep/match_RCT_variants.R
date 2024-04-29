####script for QC/matching variants used for PGS scoring in UKB with the clinical trials  ######
#helper function
flip <- function(a){
  a <- ifelse(a == "A", "T",
              ifelse(a == "T", "A",
                     ifelse(a == "C", "G",
                            ifelse(a == "G", "C", "N"))))
  return(a)
}

check_alleles <- function(ref1, alt1, ref2, alt2){
  check <- ifelse(ref1 == ref2 & alt1 == alt2, "SAME",
                  ifelse(ref1 == alt2 & alt1 == ref2, "SWITCH",
                         ifelse(ref1 == flip(ref2) & alt1 == flip(alt2), "FLIP",
                                ifelse(ref1 == flip(alt2) & alt1 == flip(ref2), "SWITCH_FLIP", "FAIL"))))
  return(check)
}

#read in variant lists
ukb_variants <- data.table::fread("./data/genotypes/UKBPPP.hm3_plus.info_filtered.pvar")


exscel_variants <- data.table::fread("../trials/data/EXSCEL_IMP/high_impQ_common_variants.pos.txt.gz")
declare_variants <- data.table::fread("../trials/data/DECLARE_IMP/high_impQ_common_variants.pos.txt.gz")

#read in snp map from UKB-PPP consortia
rsid <- data.table::fread("../datasets/ukb/olink_rsid_map_ALLAUTOSOMAL.tsv.gz",
                          data.table=FALSE)
rsid_missing <- ukb_variants$ID[!ukb_variants$ID %in% rsid$rsid]
rsid <- rsid[rsid$rsid %in% ukb_variants$ID,]

gc()

#look at duplicated snps
dupe_variants <- unique(rsid$rsid[duplicated(rsid$rsid)])

d <- rsid[rsid$rsid %in% dupe_variants,]

for(i in 1:nrow(d)){
  parsed <- unlist(unname(strsplit(d$ID[i], split = ":")))
  d$CHR[i] <- parsed[1]
  d$REF[i] <- parsed[3]
  d$ALT[i] <- parsed[4]
}

d <- merge(d, ukb_variants, by.x = "rsid", by.y= "ID", all = FALSE, sort=FALSE)
d$MATCH <- mapply(check_alleles, d$REF.x, d$ALT.x, d$REF.y, d$ALT.y)
drop <- d$ID[d$MATCH == "FAIL"]

rsid <- rsid[!rsid$ID %in% drop,]

hg38 <- data.frame(ID=rsid$rsid, POS38=rsid$POS38)

ukb_variants <- merge(ukb_variants, hg38, by="ID", all = FALSE, sort = FALSE)

#merge with excel and match
ukb_variants$ID_MATCH <- paste0(ukb_variants$`#CHROM`, ":", ukb_variants$POS38)

exscel_variants$ID_MATCH <- paste0(exscel_variants$V2, ":", exscel_variants$V3)
d <- merge(ukb_variants, exscel_variants, by="ID_MATCH", all = FALSE, sort = FALSE)
d$MATCH <- mapply(check_alleles, d$REF, d$ALT, d$V4, d$V5)
d <- d[d$MATCH != "FAIL",]

out <- data.frame(ID1=d$V1, ID2=d$ID)

data.table::fwrite(out, "./results/EXSCEL_matched_variants.txt",
                   sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

#merge with declare and match

declare_variants$ID_MATCH <- paste0(declare_variants$V2, ":", declare_variants$V3)
d <- merge(ukb_variants, declare_variants, by="ID_MATCH", all = FALSE, sort = FALSE)
d$MATCH <- mapply(check_alleles, d$REF, d$ALT, d$V4, d$V5)
d <- d[d$MATCH != "FAIL",]

out <- data.frame(ID1=d$V1, ID2=d$ID)

data.table::fwrite(out, "./results/DECLARE_matched_variants.txt",
                   sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

