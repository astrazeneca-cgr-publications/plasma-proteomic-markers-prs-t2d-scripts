#script for resolving relatives given king-generated kinship coefficients
rel.dat <- data.table::fread("../datasets/ukb/Relatedness/ukb_rel.dat")

#specify degree
degree <- 2
if(degree == 2){
  d <- 0.0884
}

if(degree == 3){
  d <- 0.0442 
}


#subset matrix down to those with kinship > than cutoff
rels <- rel.dat[rel.dat$Kinship > d,]
counts <- data.frame(ID=unique(c(rels$ID1, rels$ID2)),
                     stringsAsFactors = FALSE)
ids <- c(rels$ID1, rels$ID2)


#count number of occurences each subject appears as a proxy for number of relatives
getCount <- function(x){
  y <- length(ids[ids == x])
  return(y)
}

counts$N <- sapply(X = counts$ID, FUN = getCount)

#subset to subjects who appear > 1 time
counts <- counts[counts$N > 1,]
counts <- counts[order(counts$N, decreasing = TRUE),]
drop <- c()
for(i in 1:nrow(counts)){
  id <- counts$ID[i]
  #if subject has not already been resolved, remove
  if(id %in% rels$ID1 | id %in% rels$ID2){
    rels <- rels[rels$ID1 != id & rels$ID2 != id,]
    drop <- c(drop, id)
  }else{
    next
  }
  
}

#remaining relative pairs need to be resolved
#read in phenotype data
p <- data.table::fread("./data/phenotypes/UKB.phenotypes.tsv.gz")

cm_cases <- p$userID[p$diabetes_mellitus == 1 | p$CKD_all_stages == 1 | p$ischaemic_heart_diseases == 1]
rels$CM1 <- ifelse(rels$ID1 %in% cm_cases, 1, 0)
rels$CM2 <- ifelse(rels$ID2 %in% cm_cases, 1, 0)

drop_step2 <- c()
for(i in 1:nrow(rels)){
  #prioritise CVM subjects
  id <- ifelse(rels$CM1[i] == 1 & rels$CM2[i] == 0, rels$ID2[i], 
               ifelse(rels$CM1[i] == 0 & rels$CM2[i] == 1, rels$ID1[i], 
                      sample(c(rels$ID1[i], rels$ID2[i]), 1)) )
  
  drop_step2 <- c(drop_step2, id)
}

#get final list
drop <- c(drop, drop_step2)

out.dat <- data.frame(FID=drop, IID=drop)
out <- paste0("./data/phenotypes/UKB.degree",
              degree, ".drop.tsv")
write.table(out.dat, out, sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE)
