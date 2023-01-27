library(Seurat)
library(ggplot2)
library(pheatmap)
library(viridis)
library(DEGreport)
library(RColorBrewer)

###### variation of the LHA genes during congnitive impairment ########################
setwd("/Projects/LLIImmune/HumanBrainByMathys")
BrainCell.data <- Read10X(data.dir = "FilterData/")
BrainCell <- CreateSeuratObject(counts = BrainCell.data, project = "BrainCell")
Info=read.table("FilterData/Phenotype4Cell.txt",header=TRUE,row.names=1,sep="\t")
BrainCell$CellType=Info$broad.cell.type
BrainCell <- NormalizeData(BrainCell, normalization.method = "LogNormalize", scale.factor = 10000)
BrainCell <- FindVariableFeatures(BrainCell, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(BrainCell)
BrainCell <- ScaleData(BrainCell, features = all.genes)
BrainCell <- RunPCA(BrainCell, features = VariableFeatures(object = BrainCell))
BrainCell<- RunUMAP(BrainCell, reduction = "pca", dims = 1:30)
BrainCell<- FindNeighbors(BrainCell, reduction = "pca", dims = 1:30)
BrainCell<- FindClusters(BrainCell, resolution = 0.5)
pdf("BrainCellCluster.pdf",width=8)
DimPlot(BrainCell, reduction = "umap", label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

pdf("BrainCellType.pdf",width=8)
DimPlot(BrainCell, reduction = "umap", group.by="CellType")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

saveRDS(BrainCell,file="BrainCell.rds") #As different R version or suerat version could lead to slightly different result, so using this dataset to reproduce our graph


#----isolated the immunce cells from brain tissue and identify the T cell, microglia, and macrophage-------------
setwd("/data2/deng/Projects/LLIImmune/HumanBrainByMathys")
#As different R version or suerat version could lead to slightly different result, so using this dataset to reproduce our graph
BrainCell=readRDS("BrainCell.rds") 
Immune=subset(BrainCell,seurat_clusters%in%c(13,21))
Immune <- NormalizeData(Immune, normalization.method = "LogNormalize", scale.factor = 10000)
Immune <- FindVariableFeatures(Immune, selection.method = "vst", nfeatures = 2000)
Immune <- ScaleData(Immune, verbose = FALSE)
Immune <- RunPCA(Immune, features = VariableFeatures(object = Immune))
Immune<- RunTSNE(Immune, reduction = "pca", dims = 1:30)
Immune<- FindNeighbors(Immune, reduction = "pca", dims = 1:30)
Immune<- FindClusters(Immune, resolution = 0.1)
#cell type marker for Mic, Tcell and Mac were derived from https://www.nature.com/articles/s41586-021-04369-3/figures/8
pdf("ImmuneMarker.pdf",height=4,width=6)
FeaturePlot(Immune,reduction="tsne",ncol=3,cols=c("lightgray","red"),features=c("SPP1","P2RY12","C3","CD163","MRC1","CD247"))&NoLegend()&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
Immune$CellType=ifelse(Immune$seurat_clusters==0,"Mic",ifelse(Immune$seurat_clusters==1,"Tcell","Mac"))
tiff("ImmuneCellType.tiff",height=400,width=850)
TSNEPlot(Immune,reduction="tsne",group.by="CellType",cols=c("#00BF7D","#00BFC4","#FF62BC"))&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("ImmuneLabel.tiff",height=800,width=850)
TSNEPlot(Immune,reduction="tsne")&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

#----Merged the immune and other cells -------------
NonImmune=subset(BrainCell,seurat_clusters%in%c(0:12,14:20,22))
NonImmune$CellType=NonImmune$OldCellType
BrainCellCom=merge(NonImmune,Immune,add.cell.ids = c("Non-Immune","Immune"))
BrainCellCom <- NormalizeData(BrainCellCom, normalization.method = "LogNormalize", scale.factor = 10000)
BrainCellCom <- FindVariableFeatures(BrainCellCom, selection.method = "vst", nfeatures = 2000)
BrainCellCom <- ScaleData(BrainCellCom, verbose = FALSE)
BrainCellCom <- RunPCA(BrainCellCom, features = VariableFeatures(object = BrainCellCom))
BrainCellCom<- RunTSNE(BrainCellCom, reduction = "pca", dims = 1:30)

tiff("BrainCellComCellType.tiff",height=800,width=850)
TSNEPlot(BrainCellCom,group.by="CellType",label=T)&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

saveRDS(BrainCellCom,file="BrainTypeContainImmune.rds")


###### identify the macrophage specific expressed LHA genes ########################
setwd("/Projects/LLIImmune/HumanBrainByMathys")
Idents(BrainCellCom)=BrainCellCom$CellType
cellTypeExpr <- AverageExpression(BrainCellCom)[["RNA"]]
LHAGene=read.table("../LHAGene.txt",header=T,row.names=1)
geneList=intersect(rownames(LHAGene),rownames(cellTypeExpr)) 
length(geneList)#358
cellTypeExpr4geneList=cellTypeExpr[geneList,]
t=pheatmap(cellTypeExpr4geneList,scale="row",border=NA,clustering_method="ward.D2")
pdf("LHAGeneExprInBrain.pdf",height=10,width=5)
print(t)
dev.off()
tree=data.frame(cutree(t$tree_row,k=10))
all(rownames(tree[geneList,])==rownames(LHAGene[geneList,]))
anno=data.frame(cbind("Pattern"=LHAGene[geneList,],Group=paste0("G",tree[geneList,],sep="")))
rownames(anno)=geneList
ann_colors = list(
    Pattern = c(DownRegulated ="lightblue", UpRegulated="OrangeRed")
)
pdf("LHAGeneInBrainGraph.pdf",height=3.5,width=8)
pheatmap(t(cellTypeExpr4geneList),annotation_colors = ann_colors,scale="column",border=NA,clustering_method="ward.D2",show_colnames=F,annotation_col=anno,color = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"))
dev.off()
anno$cellType="NA"
anno[anno$Group%in%c("G1"),"cellType"]="Neuron"
anno[anno$Group%in%c("G2"),"cellType"]="Opc"
anno[anno$Group%in%c("G3"),"cellType"]="Per"
anno[anno$Group%in%c("G4"),"cellType"]="Neuron"
anno[anno$Group%in%c("G5"),"cellType"]="Oli"
anno[anno$Group%in%c("G6"),"cellType"]="T cell"
anno[anno$Group%in%c("G7"),"cellType"]="Ast"
anno[anno$Group%in%c("G8"),"cellType"]="Mic"
anno[anno$Group%in%c("G9"),"cellType"]="End"
anno[anno$Group%in%c("G10"),"cellType"]="Mac"
GeneList=t$tree_row
GeneOrder=GeneList$labels[GeneList$order]
anno=anno[GeneOrder,]
cellOrder=c("Neuron","Oli","Opc","Ast","Mic","Mac","T cell","Per","End")
anno$cellType=factor(anno$cellType,levels=cellOrder)
anno=anno[order(anno$cellType),]
cellTypeOrder=c("Ex","In","Oli","Opc","Ast","Mic","Mac","Tcell","Per","End")
cellTypeExpr4geneList=cellTypeExpr4geneList[rownames(anno),cellTypeOrder]
LHAGroup=read.table("../WBscRNA8BlishLab/LHAGene/LHAGeneInfoBetCellType.txt",header=T,row.names=1)
LHAGroup=LHAGroup[rownames(anno),]
anno$RawGroup=LHAGroup
pdf("AgingGeneInBrainGraphOrder.pdf",height=3.5,width=8)
pheatmap(t(cellTypeExpr4geneList),annotation_col=anno[,c("Pattern","RawGroup","cellType")],annotation_colors = ann_colors,scale="column",border=NA,show_colnames=F,cluster_row=F,cluster_col=F,color = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"))
dev.off()

MacUp=anno[anno$cellType=="Mac"&anno$Pattern=="UpRegulated",]
write.table(MacUp,file="MacUpGeneInBrain.txt",sep="\t",quote=F)


###### Exploring the expression of macrophage specific-LHA genes along with congnitive impairment ########################
setwd("Projects/LLIImmune/HumanBrainByMathys")
GeneExpr=read.table("Syn3388564/geneExpr.txt",header=TRUE,row.names=1,sep="\t",check.names=F) #this file was generated from ROSMAP_RNAseq_FPKM_gene.tsv
SampleInfo=read.table("Syn3388564/Phenotype.txt",header=TRUE,row.names=1,sep="\t")
Samples=intersect(colnames(GeneExpr),rownames(SampleInfo))
GeneExpr=GeneExpr[rowSums(GeneExpr)>0,Samples]
SampleInfo=SampleInfo[Samples,]
all(colnames(GeneExpr)==rownames(SampleInfo))
SampleInfo=SampleInfo[,c("msex","apoe_genotype","braaksc","ceradscID","cogdx")]
Female=SampleInfo[SampleInfo$msex=="Female",] #only focused on Female

MicSpeAgingGene=read.table("MacUpGeneInBrain.txt")
table(SampleInfo$cogdx)
Female=SampleInfo[SampleInfo$cogdx%in%c(1,2,4,5),]
Female$cogdx=factor(Female$cogdx,levels=c(1,2,4,5))
sampleRemain=intersect(rownames(Female),colnames(GeneExpr))
AgingProtectGeneExprInROSMAP=GeneExpr[rownames(MicSpeAgingGene),sampleRemain]
dim(AgingProtectGeneExprInROSMAP)
all(rownames(Female)==colnames(AgingProtectGeneExprInROSMAP)) #TRUE
AgingProtectGeneExprInROSMAP=na.omit(AgingProtectGeneExprInROSMAP)
res <- degPatterns(AgingProtectGeneExprInROSMAP, Female, time = "cogdx",minc = 1)
write.table(res$normalized,file="MacUpGeneCogdxPattern.txt",quote=F,sep="\t")
resultTable=read.table("MacUpGeneCogdxPattern.txt",header=T)
resultTable$cogdx=factor(as.character(resultTable$cogdx),levels=c(1,2,4,5))
pdf("MacUpGeneCogdxPattern.pdf",height=5,width=6)
degPlotCluster(resultTable,time="cogdx",color="cogdx",boxes = TRUE,min_genes =  1)+geom_line(aes_string(group="genes"),alpha=0.5)+scale_color_brewer(palette="Set2")+theme_bw()
dev.off()


#####Exploring the M1 and M2 in brain marcophage ########################
Mac=subset(BrainCellCom,CellType%in%c("Mac"))
Mac <- NormalizeData(Mac, normalization.method = "LogNormalize", scale.factor = 10000)
Mac <- FindVariableFeatures(Mac, selection.method = "vst", nfeatures = 2000)
Mac <- ScaleData(Mac, verbose = FALSE)
Mac <- RunPCA(Mac, features = VariableFeatures(object = Mac))
Mac<- RunTSNE(Mac, reduction = "pca", dims = 1:30)
Mac<- FindNeighbors(Mac, reduction = "pca", dims = 1:30)
Mac<- FindClusters(Mac, resolution = 1)
M1M2=read.table("M1M2Signature.txt",skip=1,sep="\t")
colnames(M1M2)=c("Symbol","Group")
M1M2List=split(M1M2$Symbol,M1M2$Group)
Mac <- AddModuleScore(
  object = Mac,
  features = M1M2List,
  name = names(M1M2List),
  ctrl = 100
)
pdf("Mac/M1M2.pdf",height=4)
VlnPlot(Mac, pt.size=0.5,features = c("M1.like.genes1","M2.like.genes2"))&NoLegend()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

new.cluster.ids <- c("Macro-M2-like","Macro-M1-like","Macro-Others")
names(new.cluster.ids) <- levels(Mac)
Mac <- RenameIdents(Mac, new.cluster.ids)
Mac$cellType=factor(Idents(Mac),levels=c("Macro-M1-like","Macro-M2-like","Macro-Others"))
pdf("Mac/M1M2ModuleScore.pdf",height=4)
VlnPlot(Mac, pt.size=0.5,features = c("M1.like.genes1","M2.like.genes2"),group.by="cellType")
dev.off()

BrainMac=subset(Mac,CellType%in%c("Macro-M1-like","Macro-M2-like"))
tmp=data.frame(table(BrainMac$CellType,BrainMac$Statues))
colnames(tmp)=c("cellType","Statues","cellNumber")
g=ggplot(tmp, aes(Statues, cellNumber, fill=cellType)) +
  geom_bar(stat="identity",position="fill") +
  scale_y_continuous(expand=c(0,0))+theme_bw()+
  scale_fill_manual(values=c("Violet","CornflowerBlue"))+
  guides(fill = guide_legend(title = "Pattern", title.position = "top"),col = guide_legend(nrow = 1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("Mac/M1M2DisInAD.pdf",width=5)
print(g)
dev.off()
