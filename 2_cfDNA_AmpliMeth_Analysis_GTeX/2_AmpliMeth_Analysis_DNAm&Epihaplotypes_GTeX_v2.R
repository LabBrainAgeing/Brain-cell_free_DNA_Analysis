rm(list=ls())

origin <- "/path/to/AMPLIMETH/folder/"
targets <- list.dirs(origin, recursive=F, full.names=F)
trim_mat <- read.table("/path/to/file/TrimmingInput.txt", header=T, sep="\t")

for (t in 1:length(targets)){
  
################# 
amplicon <- targets[t]
print(amplicon)
trim_l <- trim_mat[trim_mat$Target==amplicon,]$TrimLeft
trim_r <- trim_mat[trim_mat$Target==amplicon,]$TrimRight
inpath <- paste0(origin,amplicon,"/")
outpath <- "/path/to/output/folder/"
#################

setwd(outpath)
dir.create(paste0(amplicon,"_filtered"))
setwd(paste0(outpath, amplicon,"_filtered"))
ls_align <- list.files(inpath, ".align$",recursive=T,full.names = T)
ls_out <- list.files(inpath,".out$", recursive=T, full.names = T)
print("Does the folders content correspond?")
print(table(gsub(".align", "", ls_align)==ls_out))

for (i in 1:length(ls_align)) {tryCatch(
{
print(ls_align[i])
sample <- gsub(".out","",strsplit(strsplit(ls_out[i],"/")[[1]][length(strsplit(ls_out[i],"/")[[1]])],"_")[[1]][2])
amplicon <- strsplit(strsplit(ls_out[i],"/")[[1]][length(strsplit(ls_out[i],"/")[[1]])],"_")[[1]][1]
align <- read.table(ls_align[i], header=F, sep="\t")
out <- read.table(ls_out[i], header=F, sep="\t")
print("Do sequencing files share the same number of reads")
print(nrow(out)==(nrow(align)/5))

v1 <- seq(1, nrow(align),by=5)
align_v1 <- align[v1,,drop=F]
align_v1 <- droplevels(align_v1)
align_v1 <- as.data.frame(align_v1)
align_v1 <- as.data.frame(strsplit(as.character(align_v1[,1]), " "))
align_v1 <- as.data.frame(t(align_v1))
colnames(align_v1) <- c("read_ID","experiment_ID","read_length","sample_ID","region_ID")
str(align_v1)

v2 <- seq(2, nrow(align),by=5)
align_v2 <- align[v2,,drop=F]
align_v2 <- droplevels(align_v2)
align_v2 <- as.data.frame(align_v2)
align_v2 <- as.data.frame(strsplit(as.character(align_v2[,1]), " "))
align_v2 <- as.data.frame(t(align_v2))
align_v2 <- align_v2[,c(3,8)]
colnames(align_v2) <- c("fraction_C", "total_C")
align_v2[,1] <- as.numeric(as.character(align_v2[,1]))
align_v2[,2] <- as.numeric(as.character(align_v2[,2]))
str(align_v2)

align <- cbind(align_v1, align_v2)
rownames(align) <- c(1:(nrow(align)))

print(table(align$read_length))
align$read_length <- as.numeric(as.character((substring(align$read_length, 8))))
#align <- align[align$read_length==m,]
#align <- droplevels(align)
dim(align)

print(table(align$fraction_C))
align <- align[align$fraction_C>=0.98,]
align <- droplevels(align)
print(nrow(align))

out <- out[rownames(out) %in% rownames(align),1,drop=F]
out$V1 <- gsub(" ","",out$V1)

out$V1 <- as.factor(out$V1)
epifreq <- data.frame(table(out$V1))
epifreq_filt <- epifreq[epifreq$Freq>=2,,drop=F]
epifreq_filt <- droplevels(epifreq_filt)
epifreq_filt <- epifreq_filt[!grepl("2",epifreq_filt$Var1,fixed = T),]

out_filt <- out[out$V1 %in% levels(epifreq_filt[,1]),1,drop=F]
out_filt <- droplevels(out_filt)
relative_freq <- epifreq_filt$Freq/nrow(out_filt)

epifreq_filt <- data.frame(epifreq_filt, relative_freq)
colnames(epifreq_filt) <- c("Epi_ID", paste0(amplicon,"_",sample, "_occurrence"), paste0(amplicon,"_",sample,"_relative_freq"))

write.table(epifreq_filt, file=paste0(outpath,amplicon,"_filtered/",sample,"_filtered.txt"), row.names=F, sep="\t")
}, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})
}

ls_filt <- list.files()
ss <- read.table("/path/to/GTeX_SampleSheet_Analysis.txt", header=T, sep="\t")
ss <- ss[order(ss$Tissue3),]

data <- read.table(paste0(outpath,amplicon,"_filtered/",ls_filt[1]), header=T, sep="\t", colClasses=c("factor", "numeric", "numeric"), check.names=F)
data[,1] <- substr(data[,1],trim_l+1,nchar(as.character(data[,1]))) #Trimming Left
data[,1] <- substr(data[,1],1,nchar(as.character(data[,1]))-trim_r) #Trimming Right
data[,1] <- as.factor(data[,1])
data_new <- aggregate(data[,2], by=list(Category=data[,1]), FUN=sum)
data_new[,3] <- data_new[,2]/sum(data_new[,2])
colnames(data_new) <- colnames(data)
for (j in 2:length(ls_filt)){
data_temp <- read.table(paste0(outpath,amplicon,"_filtered/",ls_filt[j]), header=T, sep="\t",colClasses=c("factor", "numeric", "numeric"),check.names = F)
if (nrow(data_temp)==0) {next}
data_temp[,1] <- substr(data_temp[,1],trim_l+1,nchar(as.character(data_temp[,1]))) #Trimming Left
data_temp[,1] <- substr(data_temp[,1],1,nchar(as.character(data_temp[,1]))-trim_r) #Trimming Right
data_temp[,1] <- as.factor(data_temp[,1])
data_new_temp <- aggregate(data_temp[,2], by=list(Category=data_temp[,1]), FUN=sum)
data_new_temp[,3] <- data_new_temp[,2]/sum(data_new_temp[,2])
colnames(data_new_temp) <- colnames(data_temp)
data_new <- merge(data_new,data_new_temp, by="Epi_ID", all.x=T, all.y=T)
}
dim(data_new)
data <- data_new

m_data <- as.matrix(data[,seq(3,ncol(data), by=2)])
m_data[is.na(m_data)] <- 0
rownames(m_data) <- data$Epi_ID
m_data <- m_data[order(rownames(m_data)),]

for (n in 1:length(colnames(m_data))){
colnames(m_data)[n] <- strsplit(colnames(m_data),"_")[[n]][2]
}

ss <- ss[c(3:42,1:2,43:141),]
ss <- ss[ss$SampleID %in% colnames(m_data),]
m_data <- m_data[,match(ss$SampleID,colnames(m_data))]
table(ss$SampleID==colnames(m_data))

# library(RColorBrewer)
# sidePalette <- c(heat.colors(13),rep("#8da0cb",38))
# 
# pdf(paste0(outpath,amplicon,"_heatmap_Tissue3.pdf"), height=10, width=20)
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# colSide <- sidePalette[ss$Tissue3]
# #colSide <- colorRampPalette(brewer.pal(11, "BrBG"))(41)[ss$Tissue3]
# heatmap(m_data,scale="column",cexRow = 0.4,cexCol = 0.2, ColSideColors = colSide,col=colMain)
# legend(0.78,0.8, legend=levels(ss$Tissue3),fill=sidePalette, cex=0.75)
# #heatmap.2(m_data,col=mypalette,Rowv=F,Colv=T,dendrogram="column",key=T,symkey=F,density.info="none", margins=c(12,8), srtCol=45, trace="none",scale="none",symm=F)
# dev.off()

n_meth <- c(1:nrow(m_data))
for (i in 1:nrow(m_data)) {
n_meth_temp <- sum(as.numeric(unlist(strsplit(rownames(m_data)[i], ""))))
n_meth[i] <- n_meth_temp
}
n_meth <- as.data.frame(n_meth)
dim(n_meth)
m_data1 <- cbind(m_data,n_meth)
m_data1 <- m_data1[order(m_data1$n_meth, rownames(m_data1)),]
m_data1 <- m_data1[,-ncol(m_data1)]
m_data1 <- as.matrix(m_data1)

print(table(colnames(m_data1)==ss$SampleID))

m_data1_toprint <- m_data1
colnames(m_data1_toprint) <- paste(ss$Tissue2,colnames(m_data1),sep="_")
write.table(m_data1_toprint,file=paste0(outpath,amplicon,"_Supplementary2.txt"),sep="\t",row.names=T)

colnames(m_data1) <- ss$Tissue2

m_data_sum <- cbind(m_data1,apply(m_data1,1,sum))
m_data_sum <- m_data_sum[order(-m_data_sum[,ncol(m_data_sum)]),]
m_data1 <- m_data_sum[1:(ifelse(nrow(m_data_sum)>=10,10,nrow(m_data_sum))),-ncol(m_data_sum)]
n_meth <- c(1:nrow(m_data1))
for (i in 1:nrow(m_data1)) {
  n_meth_temp <- sum(as.numeric(unlist(strsplit(rownames(m_data1)[i], ""))))
  n_meth[i] <- n_meth_temp
}
n_meth <- as.data.frame(n_meth)
dim(n_meth)
m_data2 <- cbind(m_data1,n_meth)
m_data2 <- m_data2[order(m_data2$n_meth, rownames(m_data2)),]

m_data2 <- m_data2[,-ncol(m_data2)]
m_data2 <- as.matrix(m_data2)

ss$Tissue2 <- factor(ss$Tissue2,levels=c(unique(ss$Tissue2)))
levels(ss$Tissue2)

sidePalette <- c(heat.colors(13),rep("#8da0cb",38))


library(RColorBrewer)
library(gplots)
pdf(paste0(outpath,amplicon,"_Figure2.pdf"), height=10, width=20)
colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
colSide <- sidePalette[ss$Tissue2]
#colSide <- colorRampPalette(brewer.pal(11, "BrBG"))(41)[ss$Tissue3]
#heatmap(m_data,scale="column",cexRow = 0.4,cexCol = 0.2, ColSideColors = colSide,col=colMain)
heatmap.2(m_data2,col=colMain,Rowv=F,Colv=T,dendrogram="column",key=T,symkey=F,density.info="none", margins=c(12,12), srtCol=45, trace="none",scale="none",symm=F,cexCol = 0.8,ColSideColors = colSide)
#legend("bottomleft", legend=levels(ss$Tissue2),fill=sidePalette, cex=0.6)
dev.off()
#########

#Clustering Analysis /wo Bootstrap ####
library(dendextend)
eucl_dist <- dist(t(m_data2),method = 'euclidean')
hie_clust <- hclust(eucl_dist,method = 'complete')
plot(hie_clust, cex=.3)
dend <- as.dendrogram(hie_clust)
rule <- as.character(ss$Tissue3[order.dendrogram(dend)])
rule[!grepl("Brain",rule)] <- "#8da0cb"
rule[grep("Brain",rule)] <- "#FFAA00"
labels_colors(dend) <- rule
labels(dend) <- ss$Tissue3[order.dendrogram(dend)]
dend_SIMS <- hang.dendrogram(dend, hang_height = 0.1)
dend_SIMS <- assign_values_to_leaves_nodePar(dend_SIMS, 0.5,"lab.cex")
plot(dend_SIMS)

library(cluster)
set.seed(101)
pamclu=cluster::pam(t(m_data2),k=3)
plot(silhouette(pamclu),main=NULL)

##Optimal number of clusters ####
Ks=sapply(2:7,
          function(i) 
            summary(silhouette(pam(t(m_data2),k=i)))$avg.width)
plot(2:7,Ks,xlab="k",ylab="av. silhouette",type="b",
     pch=19)

##Bootstrap Resampling ####
library(pvclust)
set.seed(42)
result <- pvclust(m_data2, method.dist="euclidean", method.hclust="complete", nboot=2000, parallel=TRUE)


pdf(paste0(outpath,amplicon,"_Dendrogram_Boot2000.pdf"), height=12, width=10)
plot(result, cex = 0.5, cex.pv = 0.5, float = 0.01, main=amplicon)
pvrect(result, alpha=0.95)
dend <- as.dendrogram(result)
rule <- as.character(ss$Tissue3[order.dendrogram(dend)])
rule[!grepl("Brain",rule)] <- "#8da0cb"
rule[grep("Brain",rule)] <- "#FFAA00"
colored_bars(colors = rule, dend = dend, sort_by_labels_order = FALSE,rowLabels = "Clusters", y_shift = -0.3)
dev.off()
#print(result, which="134")

#Parsimony Analysis ####

##Phylogenetic Tree ####
library(phangorn)
#
to_phy_store <- as.data.frame(m_data2)
to_phy <- as.data.frame(m_data2)

for (h in 1:nrow(to_phy)){
  
  brain_temp <- to_phy[h,grep("Brain",colnames(to_phy))]
  other_temp <- to_phy[h,grep("Brain",colnames(to_phy),invert = T)]
  
  thr_brain <- median(as.numeric(brain_temp),na.rm=T)
  perc28 <- quantile(as.numeric(to_phy_store[h,]),probs = c(0.28))
  perc72 <- quantile(as.numeric(to_phy_store[h,]),probs = c(0.72))
  
  if (perc28<thr_brain) {
    to_phy[h,which(as.numeric(to_phy[h,])>=perc72)] <- "Brain"
    to_phy[h,which(as.numeric(to_phy[h,])<perc72)] <- "Other"
  } else {
    to_phy[h,which(as.numeric(to_phy[h,])<=perc28)] <- "Brain"
    to_phy[h,which(as.numeric(to_phy[h,])>perc28)] <- "Other"
  }
}

phy_data <- phyDat(t(to_phy), type = "USER", levels = unique(as.vector(t(to_phy))))
tree <- pratchet(phy_data)
parsimony_tree <- optim.parsimony(tree, phy_data)
##Cluster Analysis ####
eucl_dist <- dist(t(m_data2),method = 'euclidean')
hie_clust <- hclust(eucl_dist,method = 'complete')

###Plots####
library(ape)
library(graphics)
parsimony_tree_phylo <- as.phylo(parsimony_tree)

pdf(paste0(outpath,amplicon,"_Dendrogram_MostParsimony_2.pdf"), height=12, width=10)

#par(mfrow=c(2,1))
par(mar = c(0, 2, 1, 1))

col_phylo <- parsimony_tree_phylo$tip.label
col_phylo[!grepl("Brain",col_phylo)] <- "#8da0cb"
col_phylo[grep("Brain",col_phylo)] <- "#FFAA00"

plot.phylo(parsimony_tree_phylo,type="unrooted", direction = "downwards", cex=.4,label.offset = 0,use.edge.length=F,main=paste0(amplicon,"_Parsimony Dendrogam"),tip.color = col_phylo )

# tip_order <- parsimony_tree_phylo$tip.label[parsimony_tree_phylo$edge[parsimony_tree_phylo$edge[, 2] <= length(parsimony_tree_phylo$tip.label), 2]]
# col_labels <- tip_order
# col_labels[!grepl("Brain",col_labels)] <- "#8da0cb"
# col_labels[grep("Brain",col_labels)] <- "#FFAA00"
# colored_bars(colors = col_labels, dend = dend, sort_by_labels_order = FALSE,rowLabels = "", y_shift = 5,y_scale = 4.1,border = NA)
# axis(2)

par(mar = c(8, 2, 2, 1))  # Bottom, left, top, right
plot(hie_clust, main = paste0(amplicon,"_Cluster Analysis Dendrogram"),cex=.2,xlab="", sub="",hang=-1,labels=F,yaxt='n',ylab=NULL)
text(x = 1:nrow(t(m_data2)), y = rep(min(hie_clust$height) - 0.045, nrow(t(m_data2))),  # Adjust downward position
     labels = labels(hie_clust), srt = 90, adj = 1, cex = 0.4, xpd = TRUE, col=rule)
colored_bars(colors = rule, dend = dend, sort_by_labels_order = FALSE,rowLabels = "", y_shift = -0.005,y_scale = .02)

dev.off()

}
