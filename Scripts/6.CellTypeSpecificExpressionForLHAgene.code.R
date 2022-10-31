library(Seurat)
library(presto)
library(dplyr)
library(msigdbr)
library(tibble)
library(fgsea)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)

###### Cell type specific for the LHA genes ######
setwd("Projects/LLIImmune/WBscRNA8BlishLab")
#https://www.covid19cellatlas.org/index.patient.html #Web site to download the whole blood scRNAseq
#https://rupress.org/jem/article/218/8/e20210582/212379/Multi-omic-profiling-reveals-widespread  #Reference paper
BlishLab=readRDS("rawData/blish_awilk_covid_seurat.rds") #the raw data download from https://www.covid19cellatlas.org/index.patient.html
Healthy=subset(BlishLab,Status=="Healthy") #selected the normal samples as reference
pdf("FirstClustering/HealthyLabel.pdf",width=11) #varify the cell type identified by authors
DimPlot(Healthy)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("FirstClustering/HealthyQC.pdf",width=11)
VlnPlot(Healthy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
Healthy=subset(Healthy,subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt<20)
Healthy <- NormalizeData(Healthy, normalization.method = "LogNormalize", scale.factor = 10000)
Healthy <- FindVariableFeatures(Healthy, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Healthy)
Healthy <- ScaleData(Healthy, features = all.genes)
Healthy <- RunPCA(Healthy, features = VariableFeatures(object = Healthy))
pdf("FirstClustering/Elbow.pdf")
ElbowPlot(Healthy,ndims = 50)
dev.off()
Healthy<- RunTSNE(Healthy, reduction = "pca", dims = 1:50)
Healthy<- FindNeighbors(Healthy, reduction = "pca", dims = 1:50)
Healthy<- FindClusters(Healthy, resolution = 2)
pdf("FirstClustering/HealthyCluster.pdf",width=9)
DimPlot(Healthy,reduction="tsne",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
#varified the sample batch effect
pdf("FirstClustering/HealthySample.pdf",width=9)
DimPlot(Healthy,reduction="tsne",group.by="orig.ident")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


###### Integration between samples to remove the batch effect and identify the cell types in the whole blood##############
Healthy.list <- SplitObject(Healthy, split.by = "orig.ident")
Healthy.list <- lapply(X = Healthy.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = Healthy.list)
anchors <- FindIntegrationAnchors(object.list = Healthy.list, reduction = "rpca",dims = 1:50,k.filter = 100)
Healthy.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
Healthy.integrated <- ScaleData(Healthy.integrated, verbose = FALSE)
Healthy.integrated <- RunPCA(Healthy.integrated, verbose = FALSE)
pdf("ElbowIntegrated.pdf")
ElbowPlot(Healthy,ndims = 50)
dev.off()
Healthy.integrated<- RunTSNE(Healthy.integrated, reduction = "pca", dims = 1:50)
Healthy.integrated<- FindNeighbors(Healthy.integrated, reduction = "pca", dims = 1:50)
DefaultAssay(Healthy.integrated)="integrated"
Healthy.integrated<- FindClusters(Healthy.integrated, resolution = 3) #a extreme high resolution was set to identify the small subset of the cell types

pdf("Integrated/HealthyIntegratedSample.pdf",width=9) #sample effect was largely removed
DimPlot(Healthy.integrated,reduction="tsne",group.by="orig.ident")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Integrated/HealthyIntegratedCluster.pdf",width=9)
DimPlot(Healthy.integrated,reduction="tsne",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Integrated/HealthyIntegratedCellType8Author.pdf",width=9) #cell types provided by author were treated as reference
DimPlot(Healthy.integrated,reduction="tsne",group.by="cell.type.coarse")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

#validate the cell type by previous validated markers
MyMarkerGene=c("G0S2","PRF1","GZMB","TNFRSF11A","CD8A","LYZ","CD14","FCGR3A","CD3D","CCR4","IL2RA","CCR7","KLRB1","IL7R","CLC","GATA2","CD79A","TCL1A","IGHG1","IGHA1","HBB")
pdf("Integrated/MyMarkerGeneFeature.pdf",width=12,height=10)
FeaturePlot(Healthy.integrated, features = MyMarkerGene, raster=TRUE,min.cutoff = "q9",cols=c("lightgrey","Red"),ncol=5)&NoLegend()&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

#validate the cell type associated markers
Idents(Healthy.integrated.cellType)=Healthy.integrated.cellType$cellType
Healthy.integrated.celltype.markers <- FindAllMarkers(Healthy.integrated.cellType, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Healthy.integrated.celltype.markers,file="Integrated/Healthy.integrated.cellType.markers.txt",sep="\t",quote=F)
Healthy.integrated.celltype.markers %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = avg_log2FC) -> top3
Healthy.integrated.cellType <- ScaleData(Healthy.integrated.cellType, verbose = FALSE)
pdf("Integrated/cellTypeMarkerDoHeatmap.pdf",width=20,height=15)
DoHeatmap(Healthy.integrated.cellType, features = top3$gene) 
dev.off()

Healthy.integrated.cellType=subset(Healthy.integrated,idents=c(0:47)) #remove the last three clusters with extremely lower cells
new.cluster.ids <- c("Neutrophil","NK","CD8 T","CD14 Mono ","CD14 Mono ","Neutrophil","NK","CD4 CTLs","NK","CD14 Mono ","Naive CD4 T","NK","CD8 T","Neutrophil","CD8 T","Neutrophil","CD14 Mono ","Naive CD4 T","Naive CD4 T","Naive B","CD265 NK","Neutrophil","Neutrophil","Naive B","Neutrophil","Naive B","CD16 Mono ","CD4 Treg","Platelet","CD16 Mono ","Naive B","IgG B","Neutro Mono transitional cell","DC","Platelet","MAIT","Red blood","pDC","Eosinophils","Neutrophil","IgA B","Basophil","Developing neutrophil","Prolif lymph","Naive CD4 T","IgG B","Doublets","Basophil")
names(new.cluster.ids) <- levels(Healthy.integrated.cellType)
Healthy.integrated.cellType <- RenameIdents(Healthy.integrated.cellType, new.cluster.ids)
Healthy.integrated.cellType$cellType=Idents(Healthy.integrated.cellType)
tiff("Integrated/Healthy.integrated.cellTypeType.tiff",width=800)
TSNEPlot(Healthy.integrated.cellType,label=T,group.by="cellType")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

saveRDS(Healthy.integrated.cellType,file="Integrated/Healthy.integrated.cellType.rds")



###### verify the cell type specific expression of the LHA genes and common age associated genes ##############
setwd("Projects/LLIImmune")
Healthy.integrated.cellType=readRDS("WBscRNA8BlishLab/Integrated/Healthy.integrated.cellType.rds")
LHAGene=read.table("LHAGene.txt",header=T)
fgsea_sets=split(LHAGene$Gene,LHAGene$Pattern)

CommonAgingGene=read.table("GTEx/CommonAging8GTEx/CommonAgingGene.txt",header=T)
CommonAgingGene$Group=paste0(CommonAgingGene$Pattern,"_G",CommonAgingGene$cluster,sep="")
AgingGene=split(LHAGene$Gene,LHAGene$Pattern)
fgsea_sets=split(CommonAgingGene$genes,CommonAgingGene$Group)
fgsea_sets$UpInLLI=AgingGene$UpInLLI
fgsea_sets$DnInLLI=AgingGene$DnInLLI

AgingGene8Peter=read.table("Peter/AgeRelatedGene.txt")
fgsea_sets$AgingGene8Peter=AgingGene8Peter[,1]
names(fgsea_sets)

Healthy.integratedMarker <- wilcoxauc(Healthy.integrated.cellType, 'cellType')
table(Healthy.integratedMarker$group)
for(cluster in unique(Healthy.integrated.cellType$cellType)){
print (cluster)
clusterCell<- Healthy.integratedMarker %>% dplyr::filter(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0, nPermSimple = 10000)
ranks=na.omit(ranks)
#fwrite(fgseaRes, file=paste0("GSEA/GSEA4AgingProtectedGeneInCellType/",cluster,".txt",sep=""), sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file=paste0("WBscRNA8BlishLab/Integrated/GSEA/LHAgeneCellTypeSpecificInfo/",cluster,".txt",sep=""), sep="\t", sep2=c("", " ", ""))
}


#-------visualization of the LHA gene specific expression among the cell types--------------
WBscRNA8BlishLab/Integrated/GSEA/GSEResultCombine.pl #perl script to combine the cell type specific expression and generate the LHAgeneCellTypeSpecificInfoCombine.txt
data=read.table("WBscRNA8BlishLab/Integrated/GSEA/LHAgeneCellTypeSpecificInfoCombine.txt",header=T,sep="\t") #
data=data[,1:8]
data$sig=ifelse(data$pval<0.05,"Sig","NonSig")
PathwayList=factor(data$pathway,levels=rev(c("AgingGene8Peter","DnInLLI","Dn_G6","Dn_G7","Dn_G17","Dn_G4","Dn_G31","Dn_G21","Dn_G9","Dn_G14","Dn_G34","UpInLLI","Up_G5","Up_G10","Up_G16","Up_G22","Up_G28","Up_G13")))
ClusterList=factor(data$Cluster,levels=c("CD14 Mono ","CD16 Mono ","Doublets","Neutro Mono transitional cell","DC","pDC","Neutrophil","Developing neutrophil","Basophil","Eosinophils","Naive CD4 T","CD4 CTLs","CD4 Treg","MAIT","CD8 T","Platelet","Red blood","NK","CD265 NK","Prolif lymph","Naive B","Memory B","IgG B"))
t=ggplot(data,aes(ClusterList,PathwayList,size=-1*log(pval),colour=NES,shape=sig))+geom_point()+
scale_color_gradient2(low="white",mid="white",high = "red")+
scale_shape_manual(values=c(3,19))+
theme_bw()+#theme(legend.position="bottom")+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=rel(1.0),colour = "black",angle=90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=rel(1.0),colour = "black"))
pdf("WBscRNA8BlishLab/Integrated/GSEA/LHACellTypeSpecificInfoCombine.pdf",height=5,width=7)
print(t)
dev.off()




#-------visualization of the LHA gene specific expression among the cell types--------------
setwd("Projects/LLIImmune")
Idents(Healthy.integrated.cellType)=Healthy.integrated.cellType$cellType
cellTypeAverageExpr <- AverageExpression(Healthy.integrated.cellType)[["RNA"]]
LHAGene=read.table("LHAGeneCombineInfo.txt",header=T) #Generated by 4.LHAgeneIdentification.code.R
LHAGene$pvalue=-log10(LHAGene$PV)
LHAGene$Pattern=ifelse(LHAGene$FC>0,"UpInLLI","DnInLLI")
LHAGene=LHAGene[,c("pvalue","Pattern")]
table(LHAGene$Pattern)
#DnInLLI UpInLLI
#    211     159
geneList=intersect(rownames(cellTypeAverageExpr),rownames(LHAGene))
length(geneList) #367
LHAGeneExpression=cellTypeAverageExpr[geneList,]
LHAGene=LHAGene[geneList,]
LHAGeneExpression=LHAGeneExpression[ , !colnames(LHAGeneExpression) %in% c("Doublets")]
all(rownames(LHAGeneExpression)==rownames(LHAGene)) #TRUE
ScaledMatrix=pheatmap:::scale_rows(LHAGeneExpression)
heatmap=ComplexHeatmap::Heatmap(ScaledMatrix,row_km = 4,clustering_method_columns="ward.D2",clustering_method_rows="ward.D2",cluster_columns=T,show_row_names=F,show_column_names=T,col = colorRampPalette(c("SeaGreen", "white", "OrangeRed"))(100))
gene_cluster <- row_order(heatmap)
clu_df <- lapply(names(gene_cluster), function(i){
  out <- data.frame(Gene = rownames(ScaledMatrix[gene_cluster[[i]],]),Cluster = paste0("G", i), stringsAsFactors = FALSE)
  return(out)
}) %>% do.call(rbind, .)
write.table(clu_df,file="WBscRNA8BlishLab/LHAGene/LHAGeneInfoBetCellType.txt",sep="\t",quote=F,row.names=F)

rownames(clu_df)=clu_df$Gene
clu_df=clu_df[rownames(LHAGene),]
all(rownames(clu_df)==rownames(LHAGene))  #shoulb be TRUE
GeneInfo=cbind(LHAGene,clu_df)
Pattern_colors=c("Orange","   YellowGreen")
names(Pattern_colors)=unique(GeneInfo$Pattern)
GroupPalette<-brewer.pal(length(unique(GeneInfo$Cluster)),"Set3")
names(GroupPalette)=sort(unique(GeneInfo$Cluster))
geneGroupAnno = rowAnnotation( 
    Pattern = GeneInfo$Pattern,
    GeneCluster = GeneInfo$Cluster,
    col = list(Pattern=Pattern_colors,GeneCluster = GroupPalette),
    show_legend = c(TRUE,TRUE,TRUE)
)
GeneShow=GeneInfo
GeneShow$symbol=rownames(GeneInfo)
GeneShow %>% group_by(Cluster) %>% arrange(Cluster, desc(pvalue)) %>% top_frac(n = 0.1, wt = pvalue) -> top5
Index=which(rownames(GeneShow) %in% top5$symbol)
geneSymbolAnno= rowAnnotation(
    GeneSymbol=anno_mark(at = Index,labels = rownames(GeneShow)[Index])
    )
heatmap2=ComplexHeatmap::Heatmap(ScaledMatrix,row_km = 4,left_annotation=geneGroupAnno,right_annotation=geneSymbolAnno,clustering_method_columns="ward.D2",clustering_method_rows="ward.D2",cluster_columns=T,show_row_names=F,show_column_names=T,col = colorRampPalette(c("SeaGreen", "white", "OrangeRed"))(100))

pdf("WBscRNA8BlishLab/LHAGene/ComplexHeatmapRowTmp.pdf",width=7,height=9) #might different if you repeatedly execute this function, but the genes won't changed in the phagocytes
draw(heatmap2)
dev.off()

