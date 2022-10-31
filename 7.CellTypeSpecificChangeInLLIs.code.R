library(Seurat)
library(ggplot2)
library(dplyr)
library(presto)
library(dplyr)
library(msigdbr)
library(tibble)
library(fgsea)
library(RColorBrewer)


setwd("Projects/LLIImmune/Hashimoto")
Hashimoto.data=read.table("rawData/01.UMI.txt",header=T,row.names=1,sep="\t",check.names=F)
BarcodeInfo=read.table("rawData/03.Cell.Barcodes.txt",row.names=1)
all(colnames(Hashimoto.data)==rownames(BarcodeInfo)) #should be TRUE

GeneSymbol=read.table("gencode.v38.annotation.geneType.txt",header=T)#gene annotation file was download from Gencode database
head(GeneSymbol)
HashimotoEnsemble=data.frame(Ensemble=rownames(Hashimoto.data))
HashimotoEnsembleSymbol=merge(HashimotoEnsemble,GeneSymbol,by="Ensemble")
HashimotoEnsembleSymbol=HashimotoEnsembleSymbol[!duplicated(HashimotoEnsembleSymbol$Symbol), ]  #remove the duplicated symbol
rownames(HashimotoEnsembleSymbol)=HashimotoEnsembleSymbol$Ensemble
EnsembleId=intersect(rownames(HashimotoEnsembleSymbol),rownames(Hashimoto.data))
HashimotoSymbol.data=Hashimoto.data[EnsembleId,]
all(rownames(HashimotoSymbol.data)==rownames(HashimotoEnsembleSymbol))#should be True
rownames(HashimotoSymbol.data)=HashimotoEnsembleSymbol$Symbol
dim(HashimotoSymbol.data)
#23011 61202
dim(Hashimoto.data)
#23384 61202
#removed 23384-23011=373 genes
CENSingleCell <- CreateSeuratObject(counts = HashimotoSymbol.data, project = "CENSingleCell", min.cells =0, min.features = 0)
all(rownames(CENSingleCell@meta.data)==rownames(BarcodeInfo))
CENSingleCell$SampleID=BarcodeInfo[,1]
CENSingleCell$Group=BarcodeInfo[,2]
CENSingleCell[["percent.mt"]] <- PercentageFeatureSet(CENSingleCell, pattern = "^MT-")
pdf("FirstClustering/HealthyQC.pdf",width=11)
VlnPlot(CENSingleCell, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
CENSingleCellSubet=subset(CENSingleCell,subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt<20) #58093 samples, removed 61202-58093=3109 cells
table(CENSingleCellSubet$SampleID)

CENSingleCellSubet <- NormalizeData(CENSingleCellSubet, normalization.method = "LogNormalize", scale.factor = 10000)
CENSingleCellSubet <- FindVariableFeatures(CENSingleCellSubet, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(CENSingleCellSubet)
CENSingleCellSubet <- ScaleData(CENSingleCellSubet, features = all.genes)
CENSingleCellSubet <- RunPCA(CENSingleCellSubet, features = VariableFeatures(object = CENSingleCellSubet))
pdf("FirstClustering/Elbow.pdf")
ElbowPlot(CENSingleCellSubet,ndims = 50)
dev.off()
CENSingleCellSubet<- RunTSNE(CENSingleCellSubet, reduction = "pca", dims = 1:50)
CENSingleCellSubet<- FindNeighbors(CENSingleCellSubet, reduction = "pca", dims = 1:50)
CENSingleCellSubet<- FindClusters(CENSingleCellSubet, resolution = 3)
pdf("FirstClustering/CENSingleCellSubetCluster.pdf",width=9)
DimPlot(CENSingleCellSubet,reduction="tsne",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

#found the sample batch effect
pdf("FirstClustering/CENSingleCellSample.pdf",width=9)
DimPlot(CENSingleCellSubet,reduction="tsne",group.by="SampleID")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


######  Integration between samples  ##############
CENSingleCell.list <- SplitObject(CENSingleCellSubet, split.by = "SampleID")
CENSingleCell.list <- lapply(X = CENSingleCell.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = CENSingleCell.list)
anchors <- FindIntegrationAnchors(object.list = CENSingleCell.list, reduction = "rpca",dims = 1:50,k.filter = 100)
CENSingleCell.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
CENSingleCell.integrated <- ScaleData(CENSingleCell.integrated, verbose = FALSE)
CENSingleCell.integrated <- RunPCA(CENSingleCell.integrated, verbose = FALSE)
pdf("ElbowIntegrated.pdf")
ElbowPlot(CENSingleCell.integrated,ndims = 50)
dev.off()
CENSingleCell.integrated<- RunTSNE(CENSingleCell.integrated, reduction = "pca", dims = 1:30)
CENSingleCell.integrated<- FindNeighbors(CENSingleCell.integrated, reduction = "pca", dims = 1:30)
DefaultAssay(CENSingleCell.integrated)="integrated"
CENSingleCell.integrated<- FindClusters(CENSingleCell.integrated, resolution = 0.5)
saveRDS(CENSingleCell.integrated,file="CENSingleCell.integrated.rds")

pdf("Integrated/CENSingleCellIntegratedCluster.pdf",width=9)
DimPlot(CENSingleCell.integrated,reduction="tsne",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Integrated/CENSingleCellIntegratedSampleId.pdf",width=8)
DimPlot(CENSingleCell.integrated,reduction="tsne",group.by="SampleID")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Integrated/CENSingleCellIntegratedGroup.pdf",width=16)
DimPlot(CENSingleCell.integrated,reduction="tsne",split.by="Group")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

DefaultAssay(CENSingleCell.integrated)="RNA"
Marker=c("CD3D","CCR7","CD4","CD8A","KLRF1","GZMB","CD14","FCGR3A","MS4A1","HBA1")
pdf("Integrated/MyMarkerGeneFeature4CENSingle.pdf",width=12,height=5)
FeaturePlot(Healthy.integrated, features = Marker, raster=TRUE, min.cutoff = "q9",cols=c("lightgrey","Red"),ncol=5)&NoLegend()&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

CENSingleCell.integrated.subset=subset(CENSingleCell.integrated,idents=c(0:16))
new.cluster.ids <- c("CD14 Mono","NK","NK","Naive CD4 T","CD4 T CTLs","Naive CD4 T","CD8 T","CD8 T","B","CD16 Mono","NK","Megakaryocyte","B","CD8 T","DC","Red blood cell","MKI67")
names(new.cluster.ids) <- levels(CENSingleCell.integrated.subset)
CENSingleCell.integrated.subset <- RenameIdents(CENSingleCell.integrated.subset, new.cluster.ids)
CENSingleCell.integrated.subset$cellType=factor(Idents(CENSingleCell.integrated.subset),levels=c("CD14 Mono","CD16 Mono","NK","CD4 T CTLs","Naive CD4 T","CD8 T","B","Megakaryocyte","DC","Red blood cell","MKI67"))
colorPalette=brewer.pal(length(unique(CENSingleCell.integrated.subset$cellType)),"Set3")
tiff("Integrated/CENSingleCellintegrated.cellType.tiff",width=600)
TSNEPlot(CENSingleCell.integrated.subset,label=T,group.by="cellType",cols =colorPalette)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

tiff("Integrated/CENSingleCell.integrated.cellTypeSplit.tiff",width=1000)
TSNEPlot(CENSingleCell.integrated.subset,label=T,group.by="cellType",split.by="Group",cols =colorPalette)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

saveRDS(CENSingleCell.integrated.subset,file="Integrated/CENSingleCell.integrated.cellType.rds")



###### Validation of the cell type specific variation using the scRNA between LLIs and control  #####################
CENSingleCell=readRDS("CENSingleCell.integrated.subset")
LHAGene=read.table("LHAGeneCombineInfo.txt",header=T)
LHAGene$Pattern=ifelse(LHAGene$FC>0,"UpInLLI","DnInLLI")
fgsea_sets=split(LHAGene$Gene,LHAGene$Pattern)

cluster="B"
subType=subset(CENSingleCell,cellType%in% cluster)
Idents(subType)=subType$Group
CENSingleCellMarker <- wilcoxauc(subType, 'Group')
clusterCell<- CENSingleCellMarker %>% dplyr::filter(group == "SC") %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0, nPermSimple = 10000)
ranks=na.omit(ranks)
pdf("GSEA/C4InBcell.pdf",height=4)
plotEnrichment(fgsea_sets[["G4"]],ranks) 
dev.off()





