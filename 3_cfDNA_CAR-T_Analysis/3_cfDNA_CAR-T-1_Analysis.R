rm(list=ls())
library(stringr)
setwd("/path/to/output/folder/")

experiment_name <- "your_experiment" #Input Experiment name
target_amplicons <- c("APC2","FAM123","PACRG","cg09787504","cg02619656") #Input list of target amplicons

ss <- read.table("/path/to/cfDNACAR-T-1_SampleSheet_Analysis.txt",sep="\t",header=T,stringsAsFactors = T)
ss <- ss[ss$Group%in%c("CTR","PosTox"),]
ss$Group<-factor(ss$Group,levels=c("CTR","PosTox"))
str(ss)
rownames(ss) <- ss$SampleID

for (i in 1:length(target_amplicons)) {
  pdf(paste0(target_amplicons[i],"_Boxplot.pdf"),height=8,width=10)
  print(target_amplicons[i])
  data <- read.csv(paste0("/path/to/AMPLIMETH/folder/",target_amplicons[i],"/",target_amplicons[i],"_CovMeth_Amplimeth.csv"),header=T)
  data <- data.frame(rownames(data),data)
  colnames(data)[1] <- "SampleID"
  data$SampleID <- gsub(paste0(target_amplicons[i],"_"),"",data$SampleID)
  data_annotated <- merge(ss,data,by="SampleID",all.x=T)
  for (j in 5:ncol(data_annotated)) {
    boxplot(data_annotated[,j]~data_annotated$Group,colnames=colnames(data_annotated[j]), ylim=c(0,1))
    stripchart(data_annotated[,j]~data_annotated$Group, col="blue",vertical=TRUE,  method="jitter", add=TRUE, pch=21)
    
  }
  dev.off()
}

##APC2 ####
setwd("/path/to/AMPLIMETH/APC2/folder/")
amplicon_name <- "APC2"
ls_stat <- list.files("/path/to/AMPLIMETH/APC2/folder",pattern=".csv$", recursive=TRUE, full.names=TRUE)
ls_stat <- ls_stat[-grep("Stat|CovMeth",ls_stat)]
cat("There are ",length(ls_stat)," samples analyzed with target amplicon ",amplicon_name,"\n","\n")
refseq_fasta <- read.table(list.files("/path/to/AMPLIMETH/folder", pattern=paste0(amplicon_name,"_ReferenceSeq.fasta"), recursive=TRUE, full.names=TRUE),header=T,stringsAsFactors=F)
n <- str_count(refseq_fasta, pattern="CG")
n <- n-1

matrix <- data.frame(matrix(nrow=length(ls_stat),ncol=n+1))

i=1
epialleles <- read.csv(ls_stat[i],header=T,sep="\t",colClasses = c("character","character","numeric","numeric"))
epialleles$profile <- gsub('^.', '', epialleles$profile)
table(nchar(epialleles$profile))
for (w in 1:nrow(epialleles)) {
  epialleles[w,]$n_meths <- str_count(epialleles[w,]$profile, "1")
}

for (z in c(0:n)) {
  matrix[i,z+1] <- sum(epialleles[epialleles$n_meths==z,]$count)/sum(epialleles$count)
}
sum(matrix[i,])
rownames(matrix)[i] <- epialleles[1,1]
colnames(epialleles)[3:4] <- paste0(epialleles[1,1],"_",colnames(epialleles)[3:4])
epialleles <- epialleles[,-1]

for (i in 2:(length(ls_stat))){
  input <- read.csv(ls_stat[i],header=T,sep="\t",colClasses = c("character","character","numeric","numeric"))
  input$profile <- gsub('^.', '', input$profile)
  print(table(nchar(input$profile)))
  for (w in 1:nrow(input)) {
    input[w,]$n_meths <- str_count(input[w,]$profile, "1")
  }
  
  for (z in c(0:n)) {
    matrix[i,z+1] <- sum(input[input$n_meths==z,]$count)/sum(input$count)
  }
  print(sum(matrix[i,]))
  rownames(matrix)[i] <- input[1,1]
  #  colnames(input)[3:4] <- paste0(input[1,1],"_",colnames(input)[3:4])
  #  input <- input[,-1]
  #  epialleles <- merge(epialleles,input,by="profile",all.x=T,all.y=T)
}
head(matrix)
colnames(matrix) <- c(0:n)
matrix_annotated <- merge(ss,matrix,by="row.names")

setwd("/path/to/result/folder/")
pdf(paste0(amplicon_name,"_boxplot_number_meth.pdf"))
for (i in 5:ncol(matrix_annotated)) {
  boxplot(matrix_annotated[,i]~matrix_annotated$Group,main=colnames(matrix_annotated)[i])
  stripchart(matrix_annotated[,i]~matrix_annotated$Group, col="blue",vertical=TRUE,  method="jitter", add=TRUE, pch=21)
  
}
dev.off()
write.table(matrix_annotated,file=paste0(amplicon_name,"_number_meth.txt"),sep="\t",row.names=F,quote=F)

## cg02619656 ####
setwd("/path/to/AMPLIMETH/cg02619656/folder/")
amplicon_name <- "cg02619656"
ls_stat <- list.files("/path/to/AMPLIMETH/cg02619656/folder", pattern=".csv$", recursive=TRUE, full.names=TRUE)
ls_stat <- ls_stat[-grep("Stat|CovMeth",ls_stat)]
cat("There are ",length(ls_stat)," samples analyzed with target amplicon ",amplicon_name,"\n","\n")
refseq_fasta <- read.table(list.files("/path/to/AMPLIMETH/folder", pattern=paste0(amplicon_name,"_ReferenceSeq.fasta"), recursive=TRUE, full.names=TRUE),header=T,stringsAsFactors=F)
n <- str_count(refseq_fasta, pattern="CG")
n <- n-2

matrix <- data.frame(matrix(nrow=length(ls_stat),ncol=n+1))

i=1
epialleles <- read.csv(ls_stat[i],header=T,sep="\t",colClasses = c("character","character","numeric","numeric"))
epialleles$profile <- gsub('^.', '', epialleles$profile)
epialleles$profile <- gsub('.$', '', epialleles$profile)

table(nchar(epialleles$profile))
for (w in 1:nrow(epialleles)) {
  epialleles[w,]$n_meths <- str_count(epialleles[w,]$profile, "1")
}

for (z in c(0:n)) {
  matrix[i,z+1] <- sum(epialleles[epialleles$n_meths==z,]$count)/sum(epialleles$count)
}
sum(matrix[i,])
rownames(matrix)[i] <- epialleles[1,1]
colnames(epialleles)[3:4] <- paste0(epialleles[1,1],"_",colnames(epialleles)[3:4])
epialleles <- epialleles[,-1]

for (i in 2:(length(ls_stat))){
  input <- read.csv(ls_stat[i],header=T,sep="\t",colClasses = c("character","character","numeric","numeric"))
  input$profile <- gsub('^.', '', input$profile)
  input$profile <- gsub('.$', '', input$profile)
  
  print(table(nchar(input$profile)))
  for (w in 1:nrow(input)) {
    input[w,]$n_meths <- str_count(input[w,]$profile, "1")
  }
  for (z in c(0:n)) {
    matrix[i,z+1] <- sum(input[input$n_meths==z,]$count)/sum(input$count)
  }
  print(sum(matrix[i,]))
  rownames(matrix)[i] <- input[1,1]
  #  colnames(input)[3:4] <- paste0(input[1,1],"_",colnames(input)[3:4])
  #  input <- input[,-1]
  #  epialleles <- merge(epialleles,input,by="profile",all.x=T,all.y=T)
}
head(matrix)
colnames(matrix) <- c(0:n)
matrix_annotated <- merge(ss,matrix,by="row.names")
setwd("/path/to/output/folder/")
pdf(paste0(amplicon_name,"_boxplot_number_meth.pdf"))
for (i in 5:ncol(matrix_annotated)) {
  boxplot(matrix_annotated[,i]~matrix_annotated$Group,main=colnames(matrix_annotated)[i])
  stripchart(matrix_annotated[,i]~matrix_annotated$Group, col="blue",vertical=TRUE,  method="jitter", add=TRUE, pch=21)
  
}
dev.off()
write.table(matrix_annotated,file=paste0(amplicon_name,"_number_meth.txt"),sep="\t",row.names=F,quote=F)


## cg09787504 ####
setwd("/path/to/AMPLIMETH/cg09787504/folder/")
amplicon_name <- "cg09787504"
ls_stat <- list.files("/path/to/AMPLIMETH/cg09787504/folder", pattern=".csv$", recursive=TRUE, full.names=TRUE)
ls_stat <- ls_stat[-grep("Stat|CovMeth",ls_stat)]
cat("There are ",length(ls_stat)," samples analyzed with target amplicon ",amplicon_name,"\n","\n")
refseq_fasta <- read.table(list.files("/path/to/AMPLIMETH/folder", pattern=paste0(amplicon_name,"_ReferenceSeq.fasta"), recursive=TRUE, full.names=TRUE),header=T,stringsAsFactors=F)
n <- str_count(refseq_fasta, pattern="CG")
n <- n-1

matrix <- data.frame(matrix(nrow=length(ls_stat),ncol=n+1))

i=1
epialleles <- read.csv(ls_stat[i],header=T,sep="\t",colClasses = c("character","character","numeric","numeric"))
epialleles$profile <- gsub('^.', '', epialleles$profile)
table(nchar(epialleles$profile))
for (w in 1:nrow(epialleles)) {
  epialleles[w,]$n_meths <- str_count(epialleles[w,]$profile, "1")
}

for (z in c(0:n)) {
  matrix[i,z+1] <- sum(epialleles[epialleles$n_meths==z,]$count)/sum(epialleles$count)
}
sum(matrix[i,])
rownames(matrix)[i] <- epialleles[1,1]
colnames(epialleles)[3:4] <- paste0(epialleles[1,1],"_",colnames(epialleles)[3:4])
epialleles <- epialleles[,-1]

for (i in 2:(length(ls_stat))){
  input <- read.csv(ls_stat[i],header=T,sep="\t",colClasses = c("character","character","numeric","numeric"))
  input$profile <- gsub('^.', '', input$profile)
  print(table(nchar(input$profile)))
  for (w in 1:nrow(input)) {
    input[w,]$n_meths <- str_count(input[w,]$profile, "1")
  }
  
  for (z in c(0:n)) {
    matrix[i,z+1] <- sum(input[input$n_meths==z,]$count)/sum(input$count)
  }
  print(sum(matrix[i,]))
  rownames(matrix)[i] <- input[1,1]
  #  colnames(input)[3:4] <- paste0(input[1,1],"_",colnames(input)[3:4])
  #  input <- input[,-1]
  #  epialleles <- merge(epialleles,input,by="profile",all.x=T,all.y=T)
}
head(matrix)
colnames(matrix) <- c(0:n)
matrix_annotated <- merge(ss,matrix,by="row.names")
setwd("/path/to/output/folder/")
pdf(paste0(amplicon_name,"_boxplot_number_meth.pdf"))
for (i in 5:ncol(matrix_annotated)) {
  boxplot(matrix_annotated[,i]~matrix_annotated$Group,main=colnames(matrix_annotated)[i])
  stripchart(matrix_annotated[,i]~matrix_annotated$Group, col="blue",vertical=TRUE,  method="jitter", add=TRUE, pch=21)
  
}
dev.off()
write.table(matrix_annotated,file=paste0(amplicon_name,"_number_meth.txt"),sep="\t",row.names=F,quote=F)

## FAM123 ####
setwd("/path/to/AMPLIMETH/FAM123/folder/")
amplicon_name <- "FAM123"
ls_stat <- list.files("/path/to/AMPLIMETH/FAM123/folder", pattern=".csv$", recursive=TRUE, full.names=TRUE)
ls_stat <- ls_stat[-grep("Stat|CovMeth",ls_stat)]
cat("There are ",length(ls_stat)," samples analyzed with target amplicon ",amplicon_name,"\n","\n")
refseq_fasta <- read.table(list.files("/path/to/AMPLIMETH/folder", pattern=paste0(amplicon_name,"_ReferenceSeq.fasta"), recursive=TRUE, full.names=TRUE),header=T,stringsAsFactors=F)
n <- str_count(refseq_fasta, pattern="CG")

matrix <- data.frame(matrix(nrow=length(ls_stat),ncol=n+1))

i=1
epialleles <- read.csv(ls_stat[i],header=T,sep="\t",colClasses = c("character","character","numeric","numeric"))
for (w in 1:nrow(epialleles)) {
  epialleles[w,]$n_meths <- str_count(epialleles[w,]$profile, "1")
}

table(nchar(epialleles$profile))
for (z in c(0:n)) {
  matrix[i,z+1] <- sum(epialleles[epialleles$n_meths==z,]$count)/sum(epialleles$count)
}
sum(matrix[i,])
rownames(matrix)[i] <- epialleles[1,1]
colnames(epialleles)[3:4] <- paste0(epialleles[1,1],"_",colnames(epialleles)[3:4])
epialleles <- epialleles[,-1]

for (i in 2:(length(ls_stat))){
  input <- read.csv(ls_stat[i],header=T,sep="\t",colClasses = c("character","character","numeric","numeric"))
  print(table(nchar(input$profile)))
  for (w in 1:nrow(input)) {
    input[w,]$n_meths <- str_count(input[w,]$profile, "1")
  }
  
  for (z in c(0:n)) {
    matrix[i,z+1] <- sum(input[input$n_meths==z,]$count)/sum(input$count)
  }
  print(sum(matrix[i,]))
  rownames(matrix)[i] <- input[1,1]
  #  colnames(input)[3:4] <- paste0(input[1,1],"_",colnames(input)[3:4])
  #  input <- input[,-1]
  #  epialleles <- merge(epialleles,input,by="profile",all.x=T,all.y=T)
}
head(matrix)
colnames(matrix) <- c(0:n)
matrix_annotated <- merge(ss,matrix,by="row.names")
setwd("/path/to/output/folder/")
pdf(paste0(amplicon_name,"_boxplot_number_meth.pdf"))
for (i in 5:ncol(matrix_annotated)) {
  boxplot(matrix_annotated[,i]~matrix_annotated$Group,main=colnames(matrix_annotated)[i])
  stripchart(matrix_annotated[,i]~matrix_annotated$Group, col="blue",vertical=TRUE,  method="jitter", add=TRUE, pch=21)
  
}
dev.off()
write.table(matrix_annotated,file=paste0(amplicon_name,"_number_meth.txt"),sep="\t",row.names=F,quote=F)


## PACRG ####
setwd("/path/to/AMPLIMETH/PACRG/folder/")
amplicon_name <- "PACRG"
ls_stat <- list.files("/path/to/AMPLIMETH/PACRG/folder", pattern=".csv$", recursive=TRUE, full.names=TRUE)
ls_stat <- ls_stat[-grep("Stat|CovMeth",ls_stat)]
cat("There are ",length(ls_stat)," samples analyzed with target amplicon ",amplicon_name,"\n","\n")
refseq_fasta <- read.table(list.files("/path/to/AMPLIMETH/folder", pattern=paste0(amplicon_name,"_ReferenceSeq.fasta"), recursive=TRUE, full.names=TRUE),header=T,stringsAsFactors=F)
n <- str_count(refseq_fasta, pattern="CG")
n <- n-1

matrix <- data.frame(matrix(nrow=length(ls_stat),ncol=n+1))

i=1
epialleles <- read.csv(ls_stat[i],header=T,sep="\t",colClasses = c("character","character","numeric","numeric"))
epialleles$profile <- gsub('.$', '', epialleles$profile)
table(nchar(epialleles$profile))
for (w in 1:nrow(epialleles)) {
  epialleles[w,]$n_meths <- str_count(epialleles[w,]$profile, "1")
}
for (z in c(0:n)) {
  matrix[i,z+1] <- sum(epialleles[epialleles$n_meths==z,]$count)/sum(epialleles$count)
}
sum(matrix[i,])
rownames(matrix)[i] <- epialleles[1,1]
colnames(epialleles)[3:4] <- paste0(epialleles[1,1],"_",colnames(epialleles)[3:4])
epialleles <- epialleles[,-1]

for (i in 2:(length(ls_stat))){
  input <- read.csv(ls_stat[i],header=T,sep="\t",colClasses = c("character","character","numeric","numeric"))
  input$profile <- gsub('.$', '', input$profile)
  print(table(nchar(input$profile)))
  for (w in 1:nrow(input)) {
    input[w,]$n_meths <- str_count(input[w,]$profile, "1")
  }
  
  for (z in c(0:n)) {
    matrix[i,z+1] <- sum(input[input$n_meths==z,]$count)/sum(input$count)
  }
  print(sum(matrix[i,]))
  rownames(matrix)[i] <- input[1,1]
  #  colnames(input)[3:4] <- paste0(input[1,1],"_",colnames(input)[3:4])
  #  input <- input[,-1]
  #  epialleles <- merge(epialleles,input,by="profile",all.x=T,all.y=T)
}
head(matrix)
colnames(matrix) <- c(0:n)
matrix_annotated <- merge(ss,matrix,by="row.names")
setwd("/path/to/output/folder")
pdf(paste0(amplicon_name,"_boxplot_number_meth.pdf"))
for (i in 5:ncol(matrix_annotated)) {
  boxplot(matrix_annotated[,i]~matrix_annotated$Group,main=colnames(matrix_annotated)[i])
  stripchart(matrix_annotated[,i]~matrix_annotated$Group, col="blue",vertical=TRUE,  method="jitter", add=TRUE, pch=21)
  
}
dev.off()
write.table(matrix_annotated,file=paste0(amplicon_name,"_number_meth.txt"),sep="\t",row.names=F,quote=F)

# bcfDNA contribution calculation ####
rm(list=setdiff(ls(), "experiment_name"))

APC2short2 <- read.table("APC2_number_meth.txt",sep="\t",header=T)
APC2short2 <- APC2short2[,c(2,ncol(APC2short2)),drop=F]
colnames(APC2short2)[2] <- "APC2short2"

FAM123 <- read.table("FAM123_number_meth.txt",sep="\t",header=T)
rownames(FAM123) <- FAM123$SampleID
FAM123 <- FAM123[,c(2,5),drop=F]
colnames(FAM123)[2] <- "FAM123"

PACRG <- read.table("PACRG_number_meth.txt",sep="\t",header=T)
rownames(PACRG) <- PACRG$SampleID
PACRG <- PACRG[,c(2,5),drop=F]
colnames(PACRG)[2] <- "PACRG"

cg09787504 <- read.table("cg09787504_number_meth.txt",sep="\t",header=T)
rownames(cg09787504) <- cg09787504$SampleID
cg09787504 <- cg09787504[,c(2,5),drop=F]
colnames(cg09787504)[2] <- "cg09787504"

cg02619656 <- read.table("cg02619656_number_meth.txt",sep="\t",header=T)
rownames(cg02619656) <- cg02619656$SampleID
cg02619656 <- cg02619656[,c(2,5),drop=F]
colnames(cg02619656)[2] <- "cg02619656"

tot <- merge(APC2short2,FAM123,by="SampleID")
tot <- merge(tot,PACRG,by="SampleID")
tot <- merge(tot,cg09787504,by="SampleID")
tot <- merge(tot,cg02619656,by="SampleID")

ss <- read.table("~/Desktop/20230912-MMISEQ-cfDNACAR-T-1_AMPLIMETH/20230912-MMISEQ-cfDNACAR-T-1_SampleSheet_Analysis.txt",sep="\t",header=T,stringsAsFactors = T)
ss$SampleID <- as.character(ss$SampleID)
ss <- ss[ss$Group%in%c("CTR","PosTox"),]
ss$Group<-factor(ss$Group,levels=c("CTR","PosTox"))
str(ss)

tot <- merge(ss,tot,by="SampleID")
rownames(tot) <- tot$SampleID

tot$Group <- gsub("CTR","NO_ICANS",tot$Group)
tot$Group <- gsub("PosTox","ICANS",tot$Group)
tot$Group <- factor(tot$Group,levels=c("NO_ICANS","ICANS"))

tot <- tot[order(tot$Group),]
tot$SampleID <- c("NO_ICANS 1","NO_ICANS 2","ICANS 1","ICANS 2","ICANS 3")
colnames(tot)[4:ncol(tot)] <- c("APC2","FAM123","PACRG","cg09787504","cg02619656")

bcfDNA_sum <- apply(tot[,4:ncol(tot)],1,sum)

pdf(paste0(experiment_name,"_bcfDNA_content.pdf"))
par(mar = c(4, 6, 2, 2))
boxplot(bcfDNA_sum~tot$Group,ylim=c(0,0.1),xlab="",ylab="Sum of frequencies of \nbrain-specific epihaplotypes",cex=2,col=c("#addd8e","#fec44f"))
stripchart(bcfDNA_sum~tot$Group, vertical=TRUE,  method="jitter", add=TRUE, pch=21,las=2,col="black",las=2)
dev.off()

# Stastistical Analysis ####

NOICANS<-bcfDNA_sum[which(tot$Group=="NO_ICANS")]
ICANS<-bcfDNA_sum[which(tot$Group=="ICANS")]
t_contrast<-t.test(ICANS,NOICANS)

results<-data.frame(matrix(ncol=2,nrow=1))
colnames(results)<-c("Delta","pValue")
rownames(results)<-c("t.test")
results[1,1]<-mean(ICANS)-mean(NOICANS)
results[1,2]<-t_contrast[3]
write.table(results, file=paste0(experiment_name,"_Statistics.txt"), sep="\t")

#bcfDNA Genome Equivalent calculation ####
qtot<-cbind(tot,bcfDNA_sum)
cfDNAq<-read.table("/path/to/cfDNA/quantitation/measures.txt", sep="\t", header=T, stringsAsFactors=T)
table(tot$SampleID==cfDNAq$SampleID)
qtot<-merge(qtot,cfDNAq[,c("SampleID","Conc")], by="SampleID")
bcfDNA_sum_cfDNAq<-qtot$bcfDNA_sum*qtot$Conc*330

pdf(paste0(experiment_name,"bcfDNA_GE_contribution.pdf"))
boxplot(bcfDNA_sum_cfDNAq~qtot$Group)
stripchart(bcfDNA_sum_cfDNAq~qtot$Group, col="blue", vertical=T, method="jitter", add=T, pch=21)
dev.off()