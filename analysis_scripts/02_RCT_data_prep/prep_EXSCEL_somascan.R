### script for preparing EXSCEL SomaScan protein data ##
dat <- data.table::fread("./SOMASCAN/SOMASCAN_proteins.gct", 
                         data.table = FALSE)

#meta data
meta <- data.table::fread("./SOMASCAN/SOMASCAN_samples_idap.tsv__metadata.csv")

#remove any sample with an assay note
meta <- meta[is.na(meta$AssayNotes),]

#split baseline vs 12 months
mbl <- meta$Name[meta$`Sample Group` == "SomaSCAN_BL"]
m12 <- meta$Name[meta$`Sample Group` == "SomaSCAN_12mo"]
dat12 <- dat[c("Name","Description", "UniProt", "EntrezGeneSymbol", "Target", "TargetFullName", m12)]
dat <- dat[c("Name","Description", "UniProt", "EntrezGeneSymbol", "Target", "TargetFullName", mbl)]


#baseline first
#identify outliers
qc <- data.frame(ID=colnames(dat)[7:ncol(dat)],
                 stringsAsFactors = FALSE)
for(id in qc$ID){
  qc$MEDIAN[qc$ID == id] <- median(dat[[id]], na.rm = TRUE)
  qc$IQR[qc$ID == id] <- IQR(dat[[id]], na.rm = TRUE)
}

#flag samples with PC1, PC2, IQR or median > 5 sd from the mean 
qc$M_FLAG <- ifelse(qc$MEDIAN > (5*sd(qc$MEDIAN)+mean(qc$MEDIAN)), "FAIL", "PASS")
qc$IQR_FLAG <- ifelse(qc$IQR > (5*sd(qc$IQR)+mean(qc$IQR)), "FAIL", "PASS")

#create matrix for pca
mat <- dat[,7:ncol(dat)]
mat <- t(mat)
colnames(mat) <- dat$Description

#pca
pca <- prcomp(mat, center = TRUE, scale. = TRUE)
pcs <- pca$x

qc$PC1 <- pcs[,1]
qc$PC2 <- pcs[,2]

qc$PC1_FLAG <- ifelse(qc$PC1 > (5*sd(qc$PC1)+mean(qc$PC1)), "FAIL", 
                      ifelse(qc$PC1 < (mean(qc$PC1)-5*sd(qc$PC1)), 
                             "FAIL", "PASS"))
qc$PC2_FLAG <- ifelse(qc$PC2 > (5*sd(qc$PC2)+mean(qc$PC2)), "FAIL", 
                      ifelse(qc$PC2 < (mean(qc$PC2)-5*sd(qc$PC2)), 
                             "FAIL", "PASS"))

drop <- unique(c(qc$ID[qc$M_FLAG == "FAIL"], qc$ID[qc$IQR_FLAG == "FAIL"],
                 qc$ID[qc$PC1_FLAG == "FAIL"], qc$ID[qc$PC2_FLAG == "FAIL"]))

#prepare protein data matrix
mat <- as.data.frame(mat)
labels <- colnames(mat)
mat$ID <- rownames(mat)
ids <- data.frame(ID=meta$Name,
                  EID=meta$`Subject ID`)
mat <- merge(mat, ids, by="ID", all=FALSE, sort=FALSE)
mat <- mat[c("ID","EID", labels)]
mat <- mat[!mat$ID %in% drop,]

#following steps in Sun et al. 2018
#natural log, regress age, sex, PCs1-3, and time from data to plate running, inverse-normalize

#read in covariates file
covar <- data.table::fread("./clinical_data/COMBINED.MACE.tsv")
covar <- covar[covar$STUDY == "EXSCEL",]

#read in pca
pca <- data.table::fread("./data/rct_data/EXSCEL-PCA.eigenvec")

#prepare analysis data set
covar <- merge(covar, pca, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)

#create tbms variable
d1 <- as.Date(meta$PlateRunDate, format="%m/%d/%Y")
d2 <- as.Date(meta$`Sample Collection Day`, format="%d-%b-%Y")
meta$tbms <- as.numeric(d1 - d2)
tbms <- data.frame(ID=meta$`Subject ID`[meta$Name %in% mbl], tmbs=meta$tbms[meta$Name %in% mbl])
covar <- merge(covar, tbms, by.x="EID", by.y="ID", all=FALSE, sort=FALSE)

pdat <- data.frame(ID=covar$ID)

for(apt in colnames(mat)[3:ncol(mat)]){
  print(apt)
  foo <- mat[c("ID", "EID", apt)]
  colnames(foo) <- c("PID", "EID", "APT")
  
  #natrual log transform
  foo$APT <- log(foo$APT)
  
  #combine data
  foo <- merge(foo, covar, by="EID", all=FALSE, sort=FALSE)
  foo <- foo[complete.cases(foo),]
  
  #regress
  f <- as.formula("APT~AGE+SEX+tmbs+PC1+PC2+PC3")
  fit <- lm(f, foo)
  
  #inverse normalize residuals
  out <- data.frame(ID=foo$ID,
                    NPX=qnorm((rank(fit$residuals,na.last="keep")-0.5)/sum(!is.na(fit$residuals))))
  colnames(out)[2] <- apt
  
  #add to data frame
  pdat <- merge(pdat, out, by="ID", all.x=TRUE, sort=FALSE)
  
}

#update labels: gene name + uniprot + somascan aptamber "seq" number
labels <- data.frame(OLD=colnames(pdat)[2:ncol(pdat)],
                     NEW=NA)

for(label in labels$OLD){
  
  gene <- dat$EntrezGeneSymbol[dat$Description == label]
  if(gene == "" | is.na(gene)) gene <- dat$Target[dat$Description == label]
  labels$NEW[labels$OLD == label] <- paste0(gene, "_", 
                                            dat12$UniProt[dat$Description == label], "_", 
                                            label)
}
colnames(pdat)[2:ncol(pdat)] <- labels$NEW

#remove spaces and commas
colnames(pdat) <- gsub(" ", "_", colnames(pdat))
colnames(pdat) <- gsub(",", "_", colnames(pdat))
colnames(pdat) <- gsub("__", "_", colnames(pdat))
colnames(pdat) <- gsub("__", "_", colnames(pdat))
colnames(pdat) <- gsub("-", "", colnames(pdat))

#write out
data.table::fwrite(pdat, "./data/rct_data/EXSCEL_SOMASCAN_baseline.tsv",
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE, na="NA")

#12 month
#identify outliers
qc <- data.frame(ID=colnames(dat12)[7:ncol(dat12)],
                 stringsAsFactors = FALSE)
for(id in qc$ID){
  qc$MEDIAN[qc$ID == id] <- median(dat12[[id]], na.rm = TRUE)
  qc$IQR[qc$ID == id] <- IQR(dat12[[id]], na.rm = TRUE)
}

#flag samples with PC1, PC2, IQR or median > 5 sd from the mean 
qc$M_FLAG <- ifelse(qc$MEDIAN > (5*sd(qc$MEDIAN)+mean(qc$MEDIAN)), "FAIL", "PASS")
qc$IQR_FLAG <- ifelse(qc$IQR > (5*sd(qc$IQR)+mean(qc$IQR)), "FAIL", "PASS")

#create matrix for pca
mat <- dat12[,7:ncol(dat12)]
mat <- t(mat)
colnames(mat) <- dat12$Description

#pca
pca <- prcomp(mat, center = TRUE, scale. = TRUE)
pcs <- pca$x

qc$PC1 <- pcs[,1]
qc$PC2 <- pcs[,2]

qc$PC1_FLAG <- ifelse(qc$PC1 > (5*sd(qc$PC1)+mean(qc$PC1)), "FAIL", 
                      ifelse(qc$PC1 < (mean(qc$PC1)-5*sd(qc$PC1)), 
                             "FAIL", "PASS"))
qc$PC2_FLAG <- ifelse(qc$PC2 > (5*sd(qc$PC2)+mean(qc$PC2)), "FAIL", 
                      ifelse(qc$PC2 < (mean(qc$PC2)-5*sd(qc$PC2)), 
                             "FAIL", "PASS"))

drop <- unique(c(qc$ID[qc$M_FLAG == "FAIL"], qc$ID[qc$IQR_FLAG == "FAIL"],
                 qc$ID[qc$PC1_FLAG == "FAIL"], qc$ID[qc$PC2_FLAG == "FAIL"]))

#prepare protein data matrix
mat <- as.data.frame(mat)
labels <- colnames(mat)
mat$ID <- rownames(mat)
ids <- data.frame(ID=meta$Name,
                  EID=meta$`Subject ID`)
mat <- merge(mat, ids, by="ID", all=FALSE, sort=FALSE)
mat <- mat[c("ID","EID", labels)]
mat <- mat[!mat$ID %in% drop,]

#following steps in Sun et al. 2018
#natural log, regress age, sex, PCs1-3, and time from data to plate running, inverse-normalize

#read in covariates file
covar <- data.table::fread("./clinical_data/COMBINED.MACE.tsv")
covar <- covar[covar$STUDY == "EXSCEL",]

#read in pca
pca <- data.table::fread("./data/rct_data/EXSCEL-PCA.eigenvec")

#prepare analysis data set
covar <- merge(covar, pca, by.x="ID", by.y="IID", all=FALSE, sort=FALSE)

#create tbms variable
d1 <- as.Date(meta$PlateRunDate, format="%m/%d/%Y")
d2 <- as.Date(meta$`Sample Collection Day`, format="%d-%b-%Y")
meta$tbms <- as.numeric(d1 - d2)
tbms <- data.frame(ID=meta$`Subject ID`[meta$Name %in% m12], tmbs=meta$tbms[meta$Name %in% m12])
covar <- merge(covar, tbms, by.x="EID", by.y="ID", all=FALSE, sort=FALSE)

pdat <- data.frame(ID=covar$ID)

for(apt in colnames(mat)[3:ncol(mat)]){
  print(apt)
  foo <- mat[c("ID", "EID", apt)]
  colnames(foo) <- c("PID", "EID", "APT")
  
  #natrual log transform
  foo$APT <- log(foo$APT)
  
  #combine data
  foo <- merge(foo, covar, by="EID", all=FALSE, sort=FALSE)
  foo <- foo[complete.cases(foo),]
  
  #regress
  f <- as.formula("APT~AGE+SEX+tmbs+PC1+PC2+PC3")
  fit <- lm(f, foo)
  
  #inverse normalize residuals
  out <- data.frame(ID=foo$ID,
                    NPX=qnorm((rank(fit$residuals,na.last="keep")-0.5)/sum(!is.na(fit$residuals))))
  colnames(out)[2] <- apt
  
  #add to data frame
  pdat <- merge(pdat, out, by="ID", all.x=TRUE, sort=FALSE)
  
}

#update labels: gene name + uniprot + somascan aptamber "seq" number
labels <- data.frame(OLD=colnames(pdat)[2:ncol(pdat)],
                     NEW=NA)

for(label in labels$OLD){
  
  gene <- dat$EntrezGeneSymbol[dat$Description == label]
  if(gene == "" | is.na(gene)) gene <- dat$Target[dat$Description == label]
  labels$NEW[labels$OLD == label] <- paste0(gene, "_", 
                                            dat12$UniProt[dat$Description == label], "_", 
                                            label)
}
colnames(pdat)[2:ncol(pdat)] <- labels$NEW

#remove spaces and commas
colnames(pdat) <- gsub(" ", "_", colnames(pdat))
colnames(pdat) <- gsub(",", "_", colnames(pdat))
colnames(pdat) <- gsub("__", "_", colnames(pdat))
colnames(pdat) <- gsub("__", "_", colnames(pdat))
colnames(pdat) <- gsub("-", "", colnames(pdat))
#write out
data.table::fwrite(pdat, "./data/rct_data/EXSCEL_SOMASCAN_month12.tsv",
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE, na="NA")
