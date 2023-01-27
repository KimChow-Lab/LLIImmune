library(ggpubr)
library(DESeq2)
setwd("Projects/LLIImmune/Hainan")
####### DEGs between LLIs and F1SPs from three populations #######
####### Chengmai dataset published in Genome Research #######
GR=read.csv("rawData/GenomeRes/GR.readcounts.cen.f1sp.csv",header=T,row.names=1)
GRSampleInfo=read.csv("rawData/GenomeRes/GR.SAMPLE.infor.cen.f1sp.csv",header=T,row.names=1)
all(colnames(GR)==rownames(GRSampleInfo))
dds <- DESeqDataSetFromMatrix(countData = GR,
                              colData = GRSampleInfo,
                              design= ~ sex+batch+group)

dds <- DESeq(dds)
res <- results(dds, contrast=c("group","CEN","F1SP"))
resOrdered <- res[order(res$pvalue),]
write.table(resOrdered,file="DEG8DeseqInGRSample.txt",quote=F,sep="\t")
summary(res)
vsd <- vst(dds, blind=FALSE)
assay(vsd) = limma::removeBatchEffect(assay(vsd),covariates=model.matrix(~sex+batch,data=GRSampleInfo),design = model.matrix(~group,data=GRSampleInfo))
write.table(format(as.data.frame(assay(vsd)),digits=3),"GR/GR-vsd.adjusted.txt",sep="\t",row.names=TRUE,quote=F) # vst transformed and adjusted expression matrix

####### Lingao dataset published in Science Advances #######
LG=read.csv("rawData/SciAdv/LG.readcounts.csv",header=T,row.names=1)
LGSampleInfo=read.csv("rawData/SciAdv/LG.SAMPLE.infor.csv",header=T,row.names=1)
all(colnames(LG)==rownames(LGSampleInfo))
dds <- DESeqDataSetFromMatrix(countData = LG,
                              colData = LGSampleInfo,
                              design= ~ group)
dds <- DESeq(dds)
res <- results(dds, contrast=c("group","LLI","F1SP"))
resOrdered <- res[order(res$pvalue),]
write.table(resOrdered,file="DEG8DeseqInLGSample.txt",quote=F,sep="\t")
summary(res)
vsd <- vst(dds, blind=FALSE)
write.table(format(as.data.frame(assay(vsd)),digits=3),"LG/LG-vsd.adjusted.txt",sep="\t",row.names=TRUE,quote=F) # vst transformed and adjusted expression matrix

####### Lingshui dataset published in Science Advances #######
setwd("Projects/LLIImmune/Hainan")
LS=read.csv("rawData/SciAdv/LS.readcounts.csv",header=T,row.names=1)
LSSampleInfo=read.csv("rawData/SciAdv/LS.SAMPLE.infor.csv",header=T,row.names=1)
all(colnames(LS)==rownames(LSSampleInfo))
dds <- DESeqDataSetFromMatrix(countData = LS,
                              colData = LSSampleInfo,
                              design= ~ group+lib_type)
dds <- DESeq(dds)
res <- results(dds, contrast=c("group","LLI","F1SP"))
resOrdered <- res[order(res$pvalue),]
write.table(resOrdered,file="DEG8DeseqInLSSample.txt",quote=F,sep="\t")
summary(res)
vsd <- vst(dds, blind=FALSE)
assay(vsd) = limma::removeBatchEffect(assay(vsd),covariates=model.matrix(~lib_type,data=LSSampleInfo),design = model.matrix(~group,data=LSSampleInfo))
write.table(format(as.data.frame(assay(vsd)),digits=3),"LS/LS-vsd.adjusted.txt",sep="\t",row.names=TRUE,quote=F) # vst transformed and adjusted expression matrix



####### Volcano plots of the DEGs between LLIs and F1SPs #######
#----Figure 1A-C--------
resOrdered=read.table("DEG8DeseqInGRSample.txt",header=T,row.names=1)
library(ggrepel)
hs_data=data.frame(resOrdered)
hs_data$threshold = as.factor(ifelse(hs_data$pvalue>0.05,"NoSig",ifelse(abs(hs_data$log2FoldChange)>log2(1.2),ifelse(hs_data$log2FoldChange>log2(1.2),"SigUp","SigDown"),ifelse(hs_data$log2FoldChange>0,"Up","Down"))))
table(hs_data$threshold)
hs_data$threshold=factor(hs_data$threshold,levels=c("SigDown","Down","NoSig","Up","SigUp"))
hs_data$ID=rownames(hs_data)
t=ggplot(data = hs_data, aes(x = log2FoldChange, y = -log10(pvalue), size=log2(baseMean+1),colour=threshold, label =ID )) +
  geom_point(alpha=0.4) +
  theme_bw() + scale_size(range = c(0, 4))+
  scale_color_manual(values=c("blue","RoyalBlue", "grey","Salmon","red")) +
  xlim(c(-4, 4)) + ylim(0,25)+
  geom_vline(xintercept=c(-log2(1.2),log2(1.2)),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
    geom_text_repel(
    data = subset(hs_data, hs_data$pvalue < 0.01 & abs(hs_data$log2FoldChange) >= log2(1.5)),
    aes(label = ID),
    size = 3,
    max.overlaps=8,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("LS/DEGVolcano.pdf",width=9)
print(t)
dev.off()



####### Overlap of the DEGs between LLIs and F1SPs from the three populations #######
setwd("Projects/LLIImmune/Hainan")
LSDEG=read.table("DEG8DeseqInLSSample.txt",header=T,row.names=1)
LSDEG=LSDEG[LSDEG$pvalue<0.05,]
LSDEGUp=LSDEG[LSDEG$log2FoldChange> 0,]
LSDEGDown=LSDEG[LSDEG$log2FoldChange< 0,]

LGDEG=read.table("DEG8DeseqInLGSample.txt",header=T,row.names=1)
LGDEG=LGDEG[LGDEG$pvalue<0.05,]
LGDEGUp=LGDEG[LGDEG$log2FoldChange> 0,]
LGDEGDown=LGDEG[LGDEG$log2FoldChange< 0,]

GRDEG=read.table("DEG8DeseqInGRSample.txt",header=T,row.names=1)
GRDEG=GRDEG[GRDEG$pvalue<0.05,]
GRDEGUp=GRDEG[GRDEG$log2FoldChange> 0,]
GRDEGDown=GRDEG[GRDEG$log2FoldChange< 0,]

#----Figure 1D--------
library(UpSetR)
DEGList=list("LSUp"=rownames(LSDEGUp),"LGUp"=rownames(LGDEGUp),"GRUp"=rownames(GRDEGUp),"LSDn"=rownames(LSDEGDown),"LGDn"=rownames(LGDEGDown),"GRDn"=rownames(GRDEGDown))
upset(fromList(DEGList),nsets = 6, keep.order = T,
	  sets = c("LGUp","LSUp","GRUp","LGDn","LSDn","GRDn"),
	  main.bar.color="RoyalBlue",
	  matrix.color="RoyalBlue",
	  sets.bar.color="LimeGreen",
	    intersections = list(
	  	list("LGUp"),list("LSUp"),list("GRUp"),list("LGDn"),list("LSDn"),list("GRDn"),list("LGUp","LSUp"),list("LGUp","GRUp"),list("LSUp","GRUp"),list("LGUp","LSUp","GRUp"),list("LGDn","LSDn"),list("LGDn","GRDn"),list("LSDn","GRDn"),list("LGDn","LSDn","GRDn")
	  	)
	  )

UpGene=Reduce(intersect, list(DEGList$LSUp,DEGList$LGUp,DEGList$GRUp))
DnGene=Reduce(intersect, list(DEGList$LSDn,DEGList$LGDn,DEGList$GRDn))
upFrame=data.frame(Gene=UpGene,Pattern="UpRegulated")
downFrame=data.frame(Gene=DnGene,Pattern="DownRegulated")
DEG=rbind(upFrame,downFrame)
write.table(DEG,file="SharedDEGList.txt",sep="\t",quote=F,row.names=F)



####### Permutation analysis to obtain the overlap significance #######
setwd("Projects/LLIImmune/Hainan")
GR=read.table("DEG8DeseqInGRSample.txt",header=T,row.names=1)
GRSigGN=GR[GR$pvalue<0.05,]
GRSigUpN=GRSigGN[GRSigGN$log2FoldChange>0,]
GRSigDnN=GRSigGN[GRSigGN$log2FoldChange<0,]
LG=read.table("DEG8DeseqInLGSample.txt",header=T,row.names=1)
LGSigGN=LG[LG$pvalue<0.05,]
LGSigUpN=LGSigGN[LGSigGN$log2FoldChange>0,]
LGSigDnN=LGSigGN[LGSigGN$log2FoldChange<0,]
LS=read.table("DEG8DeseqInLSSample.txt",header=T,row.names=1)
LSSigGN=LS[LS$pvalue<0.05,]
LSSigUpN=LSSigGN[LSSigGN$log2FoldChange>0,]
LSSigDnN=LSSigGN[LSSigGN$log2FoldChange<0,]

DEGOverlapUp=Reduce(intersect, list(rownames(GRSigUpN),rownames(LGSigUpN),rownames(LSSigUpN)))
length(DEGOverlapUp)#317

UpRandLength=array()
for(i in 1:10000){
  GRUpRand=sample(rownames(GR),dim(GRSigUpN)[1])
  LGUpRand=sample(rownames(LG),dim(LGSigUpN)[1])
  LSUpRand=sample(rownames(LS),dim(LSSigUpN)[1])
  DEGOverlapUpRand=Reduce(intersect, list(GRUpRand,LGUpRand,LSUpRand))
  UpRandLength[i]=length(DEGOverlapUpRand)
  i=i+1;
}
mean(UpRandLength) #27.509
UpRandLength=data.frame("RandN"=UpRandLength)
#-----Figure S1A--------
t=ggplot(UpRandLength,aes(x=RandN)) + geom_density(alpha=.5,fill = "Plum1")+xlim(0,length(DEGOverlapUp)+10)+theme_bw()+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title = element_blank())+
geom_segment(aes(x = length(DEGOverlapUp), y = 0.025, xend = length(DEGOverlapUp), yend = 0),arrow = arrow(length = unit(0.5, "cm")))
pdf("DEG/UpDEGOverlapRand.pdf",width=7,height=2)
print(t)
dev.off()

DEGOverlapDown=Reduce(intersect, list(rownames(GRSigDnN),rownames(LGSigDnN),rownames(LSSigDnN)))
length(DEGOverlapDown) #449
DownRandLength=array()
for(i in 1:10000){
  GRDownRand=sample(rownames(GR),dim(GRSigDnN)[1])
  LGDownRand=sample(rownames(LG),dim(LGSigDnN)[1])
  LSDownRand=sample(rownames(LS),dim(LSSigDnN)[1])
  DEGOverlapDownRand=Reduce(intersect, list(GRDownRand,LGDownRand,LSDownRand))
  DownRandLength[i]=length(DEGOverlapDownRand)
  i=i+1;
}
mean(DownRandLength) #22.9166
DownRandLength=data.frame("RandN"=DownRandLength)
#-----Figure S1A--------
t=ggplot(DownRandLength,aes(x=RandN)) + geom_density(alpha=.5,fill = "Plum1")+xlim(0,length(DEGOverlapDown)+10)+theme_bw()+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title = element_blank())+
geom_segment(aes(x = length(DEGOverlapDown), y = 0.025, xend = length(DEGOverlapDown), yend = 0),arrow = arrow(length = unit(0.5, "cm")))
pdf("DEG/DownDEGOverlapRand.pdf",width=7,height=2)
print(t)
dev.off()



############ specific genes visualization by boxplot among all populations########################
setwd("Projects/LLIImmune/Hainan")
GRSampleInfo=read.csv("rawData/GenomeRes/GR.SAMPLE.infor.cen.f1sp.csv",header=T,row.names=1)
ChMExpr=read.table("GR/GR-vsd.adjusted.txt",header=T,row.names=1)
all(rownames(GRSampleInfo)==colnames(ChMExpr))
LGSampleInfo=read.csv("rawData/SciAdv/LG.SAMPLE.infor.csv",header=T,row.names=1)
LGExpr=read.table("LG/LG-vsd.adjusted.txt",header=T,row.names=1)
all(rownames(LGSampleInfo)==colnames(LGExpr))
LSSampleInfo=read.csv("rawData/SciAdv/LS.SAMPLE.infor.csv",header=T,row.names=1)
LSExpr=read.table("LS/LS-vsd.adjusted.txt",header=T,row.names=1)
all(rownames(LSSampleInfo)==colnames(LSExpr))

GRSampleInfo$group=ifelse(GRSampleInfo$group=="CEN","LLI","F1SP")
GRSampleInfo$Population="ChM"
LGSampleInfo$Population="LG"
LSSampleInfo$Population="LS"
sampleInfo=rbind(GRSampleInfo[,c("group","Population")],LGSampleInfo[,c("group","Population")],LSSampleInfo[,c("group","Population")])

TargetGenes="MSR1"
TargetGenesExpr=cbind(ChMExpr[TargetGenes,],LGExpr[TargetGenes,],LSExpr[TargetGenes,])
all(rownames(sampleInfo)==colnames(TargetGenesExpr))
GraphData=cbind(sampleInfo,t(TargetGenesExpr))

t=ggplot(GraphData, aes(Population, y=log2(`MSR1`+1),color=group))+
geom_boxplot()+labs(title="",x="", y = "")+ theme_bw()+
scale_color_brewer(palette="Set2")+
#theme(legend.position="none")+ 
geom_point(aes(color=group),size=1.5, alpha=0.5,show.legend=TRUE, position = position_jitterdodge(jitter.width =0.15,dodge.width = 0.75))+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0),face="bold"))
pdf("SpecificGene/MSR1.pdf",height=3,width=4)
print(t)
dev.off()



######## integrated analysis in the revised version ########

setwd("D:/Aging/PBMC/HaiNan/Qingpeng")
LS=read.csv("rawData/SciAdv/LS.readcounts.csv",header=T,row.names=1)
LSSampleInfo=read.csv("rawData/SciAdv/LS.SAMPLE.infor.csv",header=T,row.names=1)
all(colnames(LS)==rownames(LSSampleInfo))
LG=read.csv("rawData/SciAdv/LG.readcounts.csv",header=T,row.names=1)
LGSampleInfo=read.csv("rawData/SciAdv/LG.SAMPLE.infor.csv",header=T,row.names=1)
all(colnames(LG)==rownames(LGSampleInfo))
GR=read.csv("rawData/GenomeRes/GR.readcounts.cen.f1sp.csv",header=T,row.names=1)
GRSampleInfo=read.csv("rawData/GenomeRes/GR.SAMPLE.infor.cen.f1sp.csv",header=T,row.names=1)
all(colnames(GR)==rownames(GRSampleInfo))
OverlapGene=Reduce(intersect, list(rownames(LS),rownames(LG),rownames(GR)))
LS=LS[OverlapGene,]
LG=LG[OverlapGene,]
GR=GR[OverlapGene,]
all(rownames(LS)==rownames(LG))
all(rownames(LS)==rownames(GR))
allCount=cbind(GR,LS,LG)
GRSampleInfo$lib_type=GRSampleInfo$type
GRSampleInfo$type=NULL
LGSampleInfo$batch="LG"
LSSampleInfo$batch="LS"
GRSampleInfo$group=ifelse(GRSampleInfo$group=="CEN","LLI","F1SP")
LGSampleInfo$sex="F"
LSSampleInfo$sex="F"
phenotye=c("group","lib_type","batch","sex")
allSampleInfo=rbind(GRSampleInfo[,phenotye],LSSampleInfo[,phenotye],LGSampleInfo[,phenotye])
all(colnames(allCount)==rownames(allSampleInfo))
dds <- DESeqDataSetFromMatrix(countData = allCount,
                              colData = allSampleInfo,
                              design= ~ sex+batch+group)

dds <- DESeq(dds)
res <- results(dds, contrast=c("group","LLI","F1SP"))
resOrdered <- res[order(res$pvalue),]
write.table(resOrdered,file="DEG8DeseqInAllSample.txt",quote=F,sep="\t")
