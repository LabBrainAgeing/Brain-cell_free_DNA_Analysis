rm(list=ls())
setwd("/path/to/working/directory")

#1. Chromatin Analysis ####
load('/path/to/cfDNA_chromatinstructure_roadmap_EPICv450k.RData')
ROAD_ss <- read.table("/path/to/ROADMAP_Samplesheet.txt", header=T, sep="\t")

##Genera Dataset Brain
brain_match <- ROAD_ss[ROAD_ss$ANATOMY=="BRAIN",]
brain_match <- brain_match[-grep("Fetal", brain_match$TISSUE),]
brain_match <- brain_match[brain_match$TYPE=="PrimaryTissue",]
brain_match <- brain_match[brain_match$AGE=="Adult",]
brain_data <- full_data[,unique(grep(paste(brain_match$ID,collapse="|"),colnames(full_data),value=T))]
colnames(brain_data) <- paste(brain_match$ID,brain_match$TISSUE,sep="_")
brain_data <- cbind(full_data[,1:6],brain_data)
head(brain_data)
write.table(brain_data, file="cfDNA_braindata_AGE.txt", sep="\t", col.names=T)

#Genera Dataset altri tessuti
other_match <- ROAD_ss[!(ROAD_ss$ANATOMY=="IPSC"|ROAD_ss$ANATOMY=="ESC"|ROAD_ss$ANATOMY=="ESC_DERIVED"|ROAD_ss$ANATOMY=="BRAIN"|ROAD_ss$ANATOMY=="CERVIX"),]
other_match <- other_match[-grep("Fetal", other_match$TISSUE),]
other_match <- other_match[other_match$TYPE=="PrimaryTissue"|other_match$TYPE=="PrimaryCell",]
other_match <- other_match[other_match$AGE=="Adult",]
other_match <- na.omit(other_match)
other_data <- full_data[,unique(grep(paste(other_match$ID,collapse="|"),colnames(full_data),value=T))]
colnames(other_data) <- paste(other_match$ID,other_match$TISSUE,sep="_")
other_data <- cbind(full_data[,1:6],other_data)
head(other_data)
write.table(other_data, file="cfDNA_otherdata_AGE.txt", sep="\t", col.names=T)

#Interseca
v0 <- matrix(nrow=nrow(full_data),ncol=3)
v0[,1] <- brain_data$posnames
v0[,2] <- brain_data$probeID
for (i in 1:nrow(full_data)){
  v1 <-c(as.character(brain_data[i,7:ncol(brain_data)])) 
  v2 <-c(as.character(other_data[i,7:ncol(other_data)])) 
  vi <- intersect(v1,v2)
  if (length(vi)==0) {
    v0[i,3] <- "1"
  } else {v0[i,3] <- "0"}
}

unique <- v0[v0[,3]=="1",]
dim(unique)
write.table(unique, file="cfDNA_BrainVSOther_AGE_Full_Data.txt", sep="\t", quote=F)


#2. DNAm Analysis ####
target <- read.table("/path/to/cfDNA_BrainVSOther_AGE_Full_Data.txt",sep="\t",header=T,quote="")
dim(target)
target <- target[,-3]
colnames(target) <- c("Localization","ID_REF")
ss <- read.table("/path/to/GPL13534_target_cfDNA_BrainVSOther_AGE_Full_Data_SampleSheet.txt",sep="\t",header=T,quote="")
dim(ss)

table(ss$OXBS)
ss <- ss[ss$OXBS != "X",]
dim(ss)
ss <- droplevels(ss)

list <- list.files(pattern="_target_cfDNA_BrainVSOther_AGE_Full_Data.txt")
length(list)
for (i in 1:length(list)) {
  print(i)
  temp <- read.table(list[i],sep="\t",header=T,quote="")
  print(dim(temp))
  GSE_temp <- strsplit(list[i],split="_")[[1]][1]
  colnames(temp)[2:ncol(temp)] <- paste(GSE_temp,colnames(temp),sep="_")
  target <- merge(target,temp,by="ID_REF",all.x=T)
  rm(temp)
}
dim(target)
colnames(target)
table(colnames(target) %in% ss$SampleID)
table(ss$SampleID %in% colnames(target))
!ss$SampleID %in% colnames(target)
ss$SampleID[!ss$SampleID %in% colnames(target)]

rownames(target) <- target$ID_REF
target <- target[,colnames(target) %in% ss$SampleID]
dim(target)
target <- target[rowSums(is.na(target))<(nrow(target)*0.1),]
dim(target)
# 1378

ss <- ss[ss$SampleID %in% colnames(target),]
dim(ss)

ss <- ss[match(colnames(target),ss$SampleID),]
table(colnames(target)==ss$SampleID)

table(ss$Tissue_macro)
length(unique(ss$series_id))

ss_per_SupplementaryFile1 <- data.frame(ss,paste0(ss$series_id,"_",ss$Tissue_micro))
colnames(ss_per_SupplementaryFile1)[ncol(ss_per_SupplementaryFile1)] <- "SeriesID_Tissue"
ss_per_SupplementaryFile1$SeriesID_Tissue <- as.factor(ss_per_SupplementaryFile1$SeriesID_Tissue)
matrix_ss_per_SupplementaryFile1 <- matrix(nrow=nlevels(ss_per_SupplementaryFile1$SeriesID_Tissue),ncol=5)

for (i in 1:nlevels(ss_per_SupplementaryFile1$SeriesID_Tissue)) {
  matrix_ss_per_SupplementaryFile1[i,1] <- levels(ss_per_SupplementaryFile1$SeriesID_Tissue)[i]
  matrix_ss_per_SupplementaryFile1[i,2] <- ss_per_SupplementaryFile1[ss_per_SupplementaryFile1$SeriesID_Tissue==levels(ss_per_SupplementaryFile1$SeriesID_Tissue)[i],"series_id"][1]
  matrix_ss_per_SupplementaryFile1[i,3] <- ss_per_SupplementaryFile1[ss_per_SupplementaryFile1$SeriesID_Tissue==levels(ss_per_SupplementaryFile1$SeriesID_Tissue)[i],"Tissue_micro"][1]
  matrix_ss_per_SupplementaryFile1[i,4] <- ss_per_SupplementaryFile1[ss_per_SupplementaryFile1$SeriesID_Tissue==levels(ss_per_SupplementaryFile1$SeriesID_Tissue)[i],"Tissue_macro"][1]
  matrix_ss_per_SupplementaryFile1[i,5] <- nrow(ss_per_SupplementaryFile1[ss_per_SupplementaryFile1$SeriesID_Tissue==levels(ss_per_SupplementaryFile1$SeriesID_Tissue)[i],])
}
write.table(matrix_ss_per_SupplementaryFile1,file="matrix_ss_per_SupplementaryFile1.txt",sep="\t",quote=F,row.names=F)

library(limma)
design <- model.matrix(~0+ss$Tissue_macro)
colnames(design) <- c("Brain","Other")
head(design)
fit <- lmFit(target,design=design)
contrast.matrix <- makeContrasts(Brain-Other,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topl1 <- topTable(fit2, num=Inf)
head(topl1)
topl1 <- data.frame(rownames(topl1),topl1)
colnames(topl1)[1] <- "ID_REF"
head(topl1)
tail(topl1)
write.table(topl1,file="~/Dropbox/cfDNA/cfDNA_chromatinAnalysis/cfDNA_BrainVSOther_AGE_Full_Data_limma_results.txt",sep="\t",row.names=F)
topl1_sign <- topl1[topl1$adj.P.Val<0.01,]
dim(topl1_sign)
summary(topl1_sign$logFC)
topl1_sign <- topl1_sign[order(topl1_sign$logFC),]
head(topl1_sign)
topl1_sign_top20_ipo <- topl1_sign[1:20,]
topl1_sign <- topl1_sign[order(-topl1_sign$logFC),]
head(topl1_sign)
topl1_sign_top20_iper <- topl1_sign[1:20,]
topl1_sign_top20 <- rbind(topl1_sign_top20_ipo,topl1_sign_top20_iper)

target_topl1_sign_top20 <- target[rownames(target) %in% topl1_sign_top20$ID_REF,]
dim(target_topl1_sign_top20)
target_topl1_sign_top20 <- target_topl1_sign_top20[match(as.character(topl1_sign_top20$ID_REF),rownames(target_topl1_sign_top20)),]

ss <- ss[with(ss,order(ss$Tissue_macro,ss$Tissue_micro)),]
target_topl1_sign_top20 <- target_topl1_sign_top20[as.character(ss$SampleID)]
table(colnames(target_topl1_sign_top20)==ss$SampleID)

datasets <- paste(ss$series_id,ss$Tissue_micro,ss$Tissue_macro,sep="_")
length(datasets)
datasets_levels <- datasets[!duplicated(datasets)]  
length(datasets_levels)
library(stringr)
datasets_levels_tissue <- str_sub(datasets_levels, -5)
datasets_levels_tissue <- gsub("Brain","yellow",datasets_levels_tissue)
datasets_levels_tissue <- gsub("Other","blue",datasets_levels_tissue)

ss <- data.frame(ss,datasets)
ss$datasets <- factor(ss$datasets,levels=datasets_levels)
levels(ss$datasets)

load("/path/to/annotations450k.RData")
annotations450k <- annotations450k[annotations450k$ID_REF %in% rownames(target_topl1_sign_top20),]
dim(annotations450k)
annotations450k <- annotations450k[match(rownames(target_topl1_sign_top20),annotations450k$ID_REF),]
table(rownames(target_topl1_sign_top20)==annotations450k$ID_REF)

sign20_ipoiper <- annotations450k[,c(2,3)]
write.table(sign20_ipoiper, file="sign20_ipoiper.bed", row.names=F, col.names=F, sep="\t", quote=F)

pdf("cfDNA_BrainVSOther_AGE_Full_Data_limma_sign20ipo_iper.pdf",height=5,width=20)
for (i in 1:nrow(target_topl1_sign_top20)) {
boxplot(as.numeric(as.character(target_topl1_sign_top20[i,]))~ss$datasets,ylim=c(0,1),las=2,cex.axis=0.4,col=datasets_levels_tissue,xlab="",main=paste0(rownames(target_topl1_sign_top20)[i]," ",annotations450k[i,]$CHR," ",annotations450k[i,]$MAPINFO," ",annotations450k[i,]$UCSC_REFGENE_NAME," ",annotations450k[i,]$RELATION_TO_UCSC_CPG_ISLAND))
}
dev.off()

