library(Seurat)
library(ggplot2)
library(dplyr)
library(presto)
library(dplyr)
library(msigdbr)
library(tibble)
library(fgsea)
library(RColorBrewer)
library(GEOquery)
library(WGCNA)
library(tidyverse)

###### cell type specific energy metabolism the scRNA provided by Blish lab########################
setwd("Projects/LLIImmune")
Healthy.integrated.cellType=readRDS("WBscRNA8BlishLab/Integrated/Healthy.integrated.cellType.rds")
metabolism=GSA.read.gmt("EM/MetabolismPathway4Human.gmt") #The metabolism assocaited pathways were manually downloaded from KEGG
fgsea_sets=metabolism$genesets
names(fgsea_sets)=metabolism$geneset.names

DefaultAssay(Healthy.integrated.cellType)="RNA"
Healthy.integratedMarker <- wilcoxauc(Healthy.integrated.cellType, 'cellType')
table(Healthy.integratedMarker$group)
for(cluster in unique(Healthy.integrated.cellType$cellType)){
print (cluster)
clusterCell<- Healthy.integratedMarker %>% dplyr::filter(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0, nPermSimple = 10000)
ranks=na.omit(ranks)
fwrite(fgseaRes, file=paste0("GSEA/EMBetType8BlishLab/",cluster,".txt",sep=""), sep="\t", sep2=c("", " ", ""))
}

data=read.table("GSEA/EMBetType8BlishLabCombine.txt",header=T,sep="\t") #this file were generated via 
data=data[,1:8]
data$sig=ifelse(data$pval<0.05,"Sig","NonSig")
sigPathway=data[data$pval<0.05,"pathway"]
data=data[data$pathway%in%unique(sigPathway),]
data=data[order(data$pval),]
pN=table(data[data$pval<0.05&data$NES>0,"pathway"])
PathwayList=factor(data$pathway,levels=names(sort(pN)))
ClusterList=factor(data$Cluster,levels=c("CD14 Mono ","CD16 Mono ","Doublets","Neutro Mono transitional cell","DC","pDC","Basophil","Eosinophils","Neutrophil","Developing neutrophil","Naive CD4 T","CD4 CTLs","CD4 Treg","MAIT","CD8 T","Platelet","Red blood","NK","CD265 NK","Prolif lymph","Naive B","Memory B","IgG B"))
t=ggplot(data,aes(ClusterList,PathwayList,size=-1*log(pval),colour=NES,shape=sig))+geom_point()+
scale_color_gradient2(low="blue",mid="white",high = "red")+
scale_shape_manual(values=c(3,19))+
theme_bw()+#theme(legend.position="bottom")+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=rel(1.0),colour = "black",angle=90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=rel(1.0),colour = "black"))
pdf("GSEA/EMBetType8BlishLabCombine.pdf",height=6,width=9)
print(t)
dev.off()



###### cell type specific energy metabolism using the scRNA provided by Hashimoto et al ########################
setwd("Projects/LLIImmune")
CENSingleCell=readRDS("Hashimoto/Integrated/CENSingleCell.integrated.cellType.rds")
DefaultAssay(CENSingleCell)="RNA"
Idents(CENSingleCell)=CENSingleCell$cellType

metabolism=GSA.read.gmt("MetabolismPathway4Human.gmt")
fgsea_sets=metabolism$genesets
names(fgsea_sets)=metabolism$geneset.names
DefaultAssay(CENSingleCell)="RNA"

CENSingleCellMarker <- wilcoxauc(CENSingleCell, 'cellType')
table(CENSingleCellMarker$group)
for(cluster in unique(CENSingleCell$cellType)){
print (cluster)
clusterCell<- CENSingleCellMarker %>% dplyr::filter(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0, nPermSimple = 10000)
ranks=na.omit(ranks)
fwrite(fgseaRes, file=paste0("GSEA/EMBetType8Hashimoto/",cluster,".txt",sep=""), sep="\t", sep2=c("", " ", ""))
}


data=read.table("GSEA/EMBetType8Hashimoto.txt",header=T,sep="\t")
data=data[,1:8]
data$sig=ifelse(data$pval<0.05,"Sig","NonSig")
sigPathway=data[data$pval<0.05&data$NES>0,"pathway"]
data=data[data$pathway%in%unique(sigPathway),]
data=data[order(data$pval),]
pN=table(data[data$pval<0.05,"pathway"])
PathwayList=factor(data$pathway,levels=names(sort(pN)))
ClusterList=factor(data$Cluster,levels=c("CD14 Mono","CD16 Mono","DC","MKI67","NK","Megakaryocyte","CD4 T CTLs","CD8 T","Naive CD4 T","B","Red blood cell"))
t=ggplot(data,aes(ClusterList,PathwayList,size=-1*log(pval),colour=NES,shape=sig))+geom_point()+
scale_color_gradient2(low="blue",mid="white",high = "red")+
scale_shape_manual(values=c(3,19))+
theme_bw()+#theme(legend.position="bottom")+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=rel(1.0),colour = "black",angle=90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=rel(1.0),colour = "black"))
pdf("GSEA/EMBetType8Hashimoto.pdf",height=4,width=7)
print(t)
dev.off()


###### cell type specific energy metabolism variation between LLIs and control ########################
setwd("Projects/LLIImmune")
CENSingleCell=readRDS("Hashimoto/Integrated/CENSingleCell.integrated.cellType.rds")
DefaultAssay(CENSingleCell)="RNA"
Idents(CENSingleCell)=CENSingleCell$cellType

cellTypeUnit="CD16 Mono" #"CD14 Mono"
subType=subset(CENSingleCell,cellType%in%cellTypeUnit)
subTypeLLIDEG <- wilcoxauc(subType, 'Group')
clusterCell<- subTypeLLIDEG %>% dplyr::filter(group == "SC") %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0, nPermSimple = 10000)
ranks=na.omit(ranks)
pdf("GSEA/OxphosCD16MonoBetSCandCT.pdf",height=4)
plotEnrichment(fgsea_sets[["Oxidative phosphorylation"]],ranks)
dev.off()

###### WGCNA to identify the variation of OXPHOS in CD14 monocyte  #########
setwd("Projects/LLIImmune")
Sys.setenv(VROOM_CONNECTION_SIZE = 50000000) #Change the size of the connection buffer (131072)
gse=getGEO(filename="MonocyteAging/rawData/GSE56045_series_matrix.txt.gz",getGPL = FALSE) #download the series manually
#Rows: 39982 Columns: 1203
phenotype=gse@phenoData@data
colnames(phenotype)
phenotype=phenotype[,c("title","geo_accession","age:ch1","bcell:ch1","neutro:ch1","nkcell:ch1","tcell:ch1","cac:ch1","plaque:ch1","racegendersite:ch1")]
colnames(phenotype)=c("Title","SampleId","Age","Bcell","Neutro","NKcell","Tcell","Cac","Palque","Racegendersite")
write.table(phenotype,file="sampleInfo.txt",quote=F,sep="\t",row.names=F)

ex <- exprs(gse) #39982  1202
write.table(ex,file="NormalizedExprInfo.txt",quote=F,sep="\t")

#--------WGCNA to identify the co-expressed module--------------------
setwd("Projects/LLIImmune/MonocyteAging/WGCNA")
datExpr=read.table("../NormalizedExprInfo.txt",header=T,row.names=1)
metaInfo=read.table("../sampleInfo.txt",header=T,row.names=1)
all(colnames(datExpr)==(metaInfo$SampleId)) #TRUE
datExprAdjust = limma::removeBatchEffect(datExpr,covariates=model.matrix(~Bcell+Neutro+NKcell+Tcell+Racegendersite,data=metaInfo),design = model.matrix(~Age,data=metaInfo)) # vst transformed and adjusted expression matrix
write.table(datExprAdjust,"../datExprAdjust8Limma.txt",sep="\t",row.names=TRUE,quote=F) 

m.mad <- apply(datExprAdjust,1,mad) 
datExprAdjustVar <- datExprAdjust[which(m.mad >  max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
datExprAdjust <- as.data.frame(t(datExprAdjustVar))

gsg = goodSamplesGenes(datExprAdjust, verbose = 3);
gsg$allOK #TRUE
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExprAdjust), method = "average");
pdf("discreteSample.pdf",width=15)
plot(sampleTree, hang=-1, main = "Sample clustering to detect outliers!", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

clust = cutreeStatic(sampleTree, cutHeight = 100, minSize = 10) 
table(clust)
# clust 1 contains the samples we want to keep.
pdf("discreteSampleCutree.pdf",width=15)
plot(sampleTree, hang=-1, main = "Sample clustering to detect outliers!", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h=100, col = "red");
dev.off()
keepSamples = (clust==1)
datExpr = datExprAdjust[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
write.table(datExpr,"datExprFinalized4WGCNA.txt",sep="\t",row.names=TRUE,quote=F) # vst transformed and adjusted expression matrix


## Powers analysis
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector=powers, networkType="unsigned", verbose=5)
pdf("pickSoftThreshold.pdf",width=12)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red") #need to check
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

power=8
net = blockwiseModules(datExpr, power = power, maxBlockSize = ncol(datExpr),
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.1,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = "pearson", 
                       maxPOutliers=1, loadTOMs=TRUE,
                       saveTOMFileBase = "datExprAdjust8Limma.tom",
                       verbose = 3)
save(net,file="datExprAdjust8Limma.net.RData")

moduleColorsLable=net$colors
moduleColors = labels2colors(moduleColorsLable)
table(moduleColors)
pdf("DendroAndColors.pdf",width=10,height=5)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

metaInfo=read.table("../sampleInfo.txt",header=T,row.names=1)
rownames(metaInfo)=metaInfo$SampleId
metaInfo=metaInfo[rownames(datExpr),]
metaInfo=metaInfo[,c("Bcell","Neutro","NKcell","Tcell","Racegendersite","Age")]

table(net$colors)
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
modTraitCor = cor(MEs_col, metaInfo, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
write.table(modTraitCor,file="Module-trait-Rvalue.txt",sep="\t",quote=F)
write.table(modTraitP,file="Module-trait-Pvalue.txt",sep="\t",quote=F)

#-------------output module genes-----------------------
moduleColorsLable=net$colors
moduleColors = labels2colors(moduleColorsLable)
table(moduleColors)
moduleGene=data.frame(Probe_Id=colnames(datExpr),Module=moduleColors)
GeneInfo=read.csv("../GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt",header=T,sep="\t")
SymbolInfo=GeneInfo[,c("Symbol","Entrez_Gene_ID","Probe_Id")]
moduleGeneInfo=merge(moduleGene,SymbolInfo,by="Probe_Id")
write.table(moduleGeneInfo,file="moduleGeneInfo.txt",sep="\t",quote=F)

#------------visulization of Module-trait-R vauue--------------------
modTraitCor=read.table("Module-trait-Rvalue.txt",header=T,row.names=1)
modTraitCor=modTraitCor[order(modTraitCor$Age),]
modTraitCor$Module=factor(rownames(modTraitCor),levels=rev(rownames(modTraitCor)))
colorList=gsub("ME","",modTraitCor$Module)
p<-ggplot(data=modTraitCor, aes(x=Module, y=Age,fill=Module)) +
   geom_bar(stat="identity")+
   scale_fill_manual(values=rev(colorList))+
   theme_bw()+theme(legend.position="none")+
   theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=rel(1.0),colour = "black",angle=90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=rel(1.0),colour = "black"))
pdf("Module-trait-Rvalue.pdf",width=10,height=3)
print(p)
dev.off()

#------------Metabobolism for each module--------------------
setwd("Projects/LLIImmune/")
Metabolism=read.table("EnergyMetabolism/MetabolismPathway4HumanGene.txt",header=T,sep="\t")
MetabolismList=split(Metabolism$Symbol,Metabolism$PathwayName)
ModuleGene=read.table("MonocyteAging/WGCNA/moduleGeneInfo.txt",header=T,sep="\t")
ModuleGeneList=split(ModuleGene$Symbol,ModuleGene$Module)

Pvalue=matrix(data = NA, nrow = length(MetabolismList)*length(ModuleGeneList), ncol =8)
index=0
for(i in 1:length(MetabolismList)){
  for (j in 1:length(ModuleGeneList)){
    index=index+1
    pathyway=unique(MetabolismList[[i]])
    Module=unique(ModuleGeneList[[j]])
    overlap=intersect(pathyway,Module)
    Pvalue[index,1]=names(MetabolismList[i])
    Pvalue[index,2]=names(ModuleGeneList[j])
    Pvalue[index,3]=1-phyper(length(overlap),length(Module),length(unique(ModuleGene$Symbol))-length(Module),length(pathyway))
    Pvalue[index,4]=length(Module)
    Pvalue[index,5]=length(pathyway)
    Pvalue[index,6]=length(overlap)
    Pvalue[index,7]=length(unique(ModuleGene$Symbol))
    Pvalue[index,8]=toString(overlap)
  }
}
Pvalue=data.frame(Pvalue)
colnames(Pvalue)=c("Pathway","Module","Pvalue","ModuleN","PathwayN","OverlapN","Total","OverlapGene")
write.table(Pvalue,file="MonocyteAging/WGCNA/Pvalue4Metabolism.txt",row.names=F,quote=F,sep="\t")

Oxphs=Pvalue[Pvalue$Pathway%in%c("Oxidative phosphorylation","Citrate cycle (TCA cycle)"),c("Pathway","Module","Pvalue")]
Oxphs$logP=-log10(Oxphs$Pvalue)
OxphsMatrix=reshape2::dcast(Oxphs,Pathway~Module,value.var ="logP")
rownames(OxphsMatrix)=OxphsMatrix$Pathway
OxphsMatrix$Pathway=NULL
OxphsMatrix=OxphsMatrix[,rev(colorList)]
pdf("MonocyteAging/WGCNA/OxphsPvalue.pdf",width=12,height=2)
pheatmap(OxphsMatrix,scale="row",cluster_col=F,cluster_row=F,color = colorRampPalette(c("blue", "white", "Gold"))(50))
dev.off()

