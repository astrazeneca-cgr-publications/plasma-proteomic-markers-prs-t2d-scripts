#script for preparing DECLARE proteomics data for analyses
#treat it similar as UKB-PPP manuscript

#read in data
dat <- data.table::fread("./OLINK/deid_olinkresult_timi.csv")
dat$SampleId <- as.character(dat$SampleId)

#make linker
linker <- data.frame(ID=dat$SampleId, EID=dat$subject)
linker <- linker[!duplicated(linker),]

#remove QC warning
dat <- dat[dat$QCWarning == "Pass",]

#split 2 and 4 visits and treat separately
dat4 <- dat[dat$visit == 4,]
dat <- dat[dat$visit == 2,]


######run for visit 2 (baseline) ###########
qc <- data.frame(ID=unique(dat$SampleId),
                 stringsAsFactors = FALSE)

mat <- qc 

for(i in qc$ID){
  #calculate median and IQR
  qc$MEDIAN[qc$ID == i] <- median(dat$LBORRES[dat$SampleId == i], na.rm = TRUE)
  qc$IQR[qc$ID == i] <- IQR(dat$LBORRES[dat$SampleId == i], na.rm = TRUE)
  
  #construct matrix of NPX values
  for(u in unique(dat$UNPROTID)){
    if(length(dat$LBORRES[dat$UNPROTID == u & dat$SampleId == i]) == 0){
      mat[[u]][mat$ID == i] <- NA
    }else{
      mat[[u]][mat$ID == i] <- dat$LBORRES[dat$UNPROTID == u & dat$SampleId == i]
    }
    
  }
}

#flag samples with PC1, PC2, IQR or median > 5 sd from the mean 
qc$M_FLAG <- ifelse(qc$MEDIAN > (5*sd(qc$MEDIAN)+mean(qc$MEDIAN)), "FAIL", "PASS")
qc$IQR_FLAG <- ifelse(qc$IQR > (5*sd(qc$IQR)+mean(qc$IQR)), "FAIL", "PASS")

#mean impute missing values for pca
mat.imp <- mat
for(j in 2:ncol(mat.imp)){
  mat.imp[,j][is.na(mat.imp[,j])] <- mean(mat.imp[,j], na.rm = TRUE)
}

pca <- prcomp(mat.imp[2:ncol(mat.imp)], center = TRUE, scale. = TRUE)
qc <- cbind(qc, as.data.frame(pca$x[,1:2]))


qc$PC1_FLAG <- ifelse(qc$PC1 > (5*sd(qc$PC1)+mean(qc$PC1)), "FAIL", 
                      ifelse(qc$PC1 < (mean(qc$PC1)-5*sd(qc$PC1)), 
                             "FAIL", "PASS"))
qc$PC2_FLAG <- ifelse(qc$PC2 > (5*sd(qc$PC2)+mean(qc$PC2)), "FAIL", 
                      ifelse(qc$PC2 < (mean(qc$PC2)-5*sd(qc$PC2)), 
                             "FAIL", "PASS"))

drop <- unique(c(qc$ID[qc$M_FLAG == "FAIL"], qc$ID[qc$IQR_FLAG == "FAIL"],
                 qc$ID[qc$PC1_FLAG == "FAIL"], qc$ID[qc$PC2_FLAG == "FAIL"]))
mat.qc <- mat[! mat$ID %in% drop, ]

#inverse normalize each column
for(u in colnames(mat.qc)[2:ncol(mat.qc)]){
  mat.qc[[u]] <- qnorm((rank(mat.qc[[u]],na.last="keep")-0.5)/sum(!is.na(mat.qc[[u]])))
}

#prepare column names
#protein_uniprot_olinki
labels <- data.frame(UNIPROT=colnames(mat.qc)[2:ncol(mat.qc)])
for(label in labels){
  labels$NEW[labels$UNIPROT == label] <- paste0(unique(dat$BMCODE[dat$UNPROTID == label]), "_", label, "_", unique(dat$olinkid[dat$UNPROTID == label]))
}

#remove problematic characeters
labels$NEW <- gsub(",", "_", labels$NEW)
labels$NEW <- gsub("-", "", labels$NEW)
labels$NEW <- gsub(" ", "_", labels$NEW)
colnames(mat.qc)[2:ncol(mat.qc)] <- labels$NEW

#add subject id
proteins <- colnames(mat.qc)[2:ncol(mat.qc)]
mat.qc <- merge(mat.qc, linker, by="ID", all=FALSE, sort=FALSE)
mat.qc <- mat.qc[c("EID", proteins)]

#write out
data.table::fwrite(mat.qc, "./data/rct_data/DECLARE_OLINK_visit2.tsv",
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE, na="NA")


######also run for visit 4 ###########
qc <- data.frame(ID=unique(dat4$SampleId),
                 stringsAsFactors = FALSE)

mat <- qc 

for(i in qc$ID){
  #calculate median and IQR
  qc$MEDIAN[qc$ID == i] <- median(dat4$LBORRES[dat4$SampleId == i], na.rm = TRUE)
  qc$IQR[qc$ID == i] <- IQR(dat4$LBORRES[dat4$SampleId == i], na.rm = TRUE)
  
  #construct matrix of NPX values
  for(u in unique(dat4$UNPROTID)){
    if(length(dat4$LBORRES[dat4$UNPROTID == u & dat4$SampleId == i]) == 0){
      mat[[u]][mat$ID == i] <- NA
    }else{
      mat[[u]][mat$ID == i] <- dat4$LBORRES[dat4$UNPROTID == u & dat4$SampleId == i]
    }
    
  }
}

#flag samples with PC1, PC2, IQR or median > 5 sd from the mean 
qc$M_FLAG <- ifelse(qc$MEDIAN > (5*sd(qc$MEDIAN)+mean(qc$MEDIAN)), "FAIL", "PASS")
qc$IQR_FLAG <- ifelse(qc$IQR > (5*sd(qc$IQR)+mean(qc$IQR)), "FAIL", "PASS")

#mean impute missing values for pca
mat.imp <- mat
for(j in 2:ncol(mat.imp)){
  mat.imp[,j][is.na(mat.imp[,j])] <- mean(mat.imp[,j], na.rm = TRUE)
}

pca <- prcomp(mat.imp[2:ncol(mat.imp)], center = TRUE, scale. = TRUE)
qc <- cbind(qc, as.data.frame(pca$x[,1:2]))


qc$PC1_FLAG <- ifelse(qc$PC1 > (5*sd(qc$PC1)+mean(qc$PC1)), "FAIL", 
                      ifelse(qc$PC1 < (mean(qc$PC1)-5*sd(qc$PC1)), 
                             "FAIL", "PASS"))
qc$PC2_FLAG <- ifelse(qc$PC2 > (5*sd(qc$PC2)+mean(qc$PC2)), "FAIL", 
                      ifelse(qc$PC2 < (mean(qc$PC2)-5*sd(qc$PC2)), 
                             "FAIL", "PASS"))

drop <- unique(c(qc$ID[qc$M_FLAG == "FAIL"], qc$ID[qc$IQR_FLAG == "FAIL"],
                 qc$ID[qc$PC1_FLAG == "FAIL"], qc$ID[qc$PC2_FLAG == "FAIL"]))
mat.qc <- mat[! mat$ID %in% drop, ]

#inverse normalize each column
for(u in colnames(mat.qc)[2:ncol(mat.qc)]){
  mat.qc[[u]] <- qnorm((rank(mat.qc[[u]],na.last="keep")-0.5)/sum(!is.na(mat.qc[[u]])))
}

#prepare column names
#protein_uniprot_olinki
labels <- data.frame(UNIPROT=colnames(mat.qc)[2:ncol(mat.qc)])
for(label in labels){
  labels$NEW[labels$UNIPROT == label] <- paste0(unique(dat4$BMCODE[dat4$UNPROTID == label]), "_", label, "_", unique(dat4$olinkid[dat4$UNPROTID == label]))
}

#remove problematic characeters
labels$NEW <- gsub(",", "_", labels$NEW)
labels$NEW <- gsub("-", "", labels$NEW)
labels$NEW <- gsub(" ", "_", labels$NEW)
colnames(mat.qc)[2:ncol(mat.qc)] <- labels$NEW

#add subject id
proteins <- colnames(mat.qc)[2:ncol(mat.qc)]
mat.qc <- merge(mat.qc, linker, by="ID", all=FALSE, sort=FALSE)
mat.qc <- mat.qc[c("EID", proteins)]

#write out
data.table::fwrite(mat.qc, "./data/rct_data/DECLARE_OLINK_visit4.tsv",
                   sep='\t', quote=FALSE, row.names = FALSE, col.names = TRUE, na="NA")
