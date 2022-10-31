###### Obtain the LHA genes by removeing the common age-associated genes ######
setwd("Projects/LLIImmune/")
library(ggrepel)
library(ggplot2)
library(UpSetR)

###### common age-associated genes identifed from GTEx ----------######
SharedDEG=read.table("Hainan/SharedDEGList.txt",header=T)
resultTable=read.table("CommonAging8GTEx/AgingGenePattern.txt",header=T)
filtclusters<-resultTable
upc = c(5,10,13,16,22,28)
downc = c(4,6,7,9,14,17,21,31,34) 
UpCommonPattern=filtclusters[filtclusters$cluster%in%upc,]
DnCommonPattern=filtclusters[filtclusters$cluster%in%downc,]
UpCommonAgingGene=unique(UpCommonPattern$genes)
DnCommonAgingGene=unique(DnCommonPattern$genes)

ClusterList=split(SharedDEG$Gene,SharedDEG$Pattern)
ClusterList$CommonUp=UpCommonAgingGene
ClusterList$CommonDn=DnCommonAgingGene
upset(fromList(ClusterList),nsets = 4,queries = list(
  list(query = intersects, params = list("CommonUp", "UpRegulated"),color="Orange",active=T),
  list(query = intersects, params = list("CommonDn", "DownRegulated"),color="Green",active=T)
  ))
table(SharedDEG$Pattern)
#DownRegulated   UpRegulated 
#          449           317
intersect(ClusterList$CommonUp,ClusterList$DownRegulated)
RemoveCommonAgingGene8GTEx=setdiff(SharedDEG$Gene,c(UpCommonAgingGene,DnCommonAgingGene)) #655 Genes

#-----------further removed the age-associated genes identificed by Peter--------------------------
AgingGene8Peter=read.table("Peter/AgeRelatedGene.txt") #1497 genes
length(intersect(RemoveCommonAgingGene8GTEx,AgingGene8Peter[,1])) #285 overlap
LHAGeneGene=setdiff(RemoveCommonAgingGene8GTEx,AgingGene8Peter[,1]) #370 genes left
LHAGeneGene=data.frame(Gene=LHAGeneGene)
LHAGeneGenePattern=merge(SharedDEG,LHAGeneGene,by="Gene")
LHAGeneGenePattern=LHAGeneGenePattern[order(LHAGeneGenePattern$Pattern),]
table(LHAGeneGenePattern$Pattern)
#DownRegulated   UpRegulated 
#          211           159
write.table(LHAGeneGenePattern,file="LHAGene.txt",sep="\t",quote=F,row.names=F)



###### visualization of the LHA genes and common aging genes  ######
LSDEG=read.table("HaiNan/DEG8DeseqInLSSample.txt",header=T,row.names=1)
LGDEG=read.table("HaiNan/DEG8DeseqInLGSample.txt",header=T,row.names=1)
GRDEG=read.table("HaiNan/DEG8DeseqInGRSample.txt",header=T,row.names=1)
SharedGene=intersect(rownames(LSDEG),rownames(LGDEG))
SharedGene=intersect(SharedGene,rownames(GRDEG))
APGeneInLS=LSDEG[SharedGene,c("baseMean","log2FoldChange","pvalue")]
APGeneInLG=LGDEG[SharedGene,c("baseMean","log2FoldChange","pvalue")]
APGeneInGR=GRDEG[SharedGene,c("baseMean","log2FoldChange","pvalue")]
all(rownames(APGeneInLS)==rownames(APGeneInLG))
all(rownames(APGeneInLS)==rownames(APGeneInGR))
APGeneInLS$Gene=rownames(APGeneInLS)
APGeneInLG$Gene=rownames(APGeneInLG)
APGeneInGR$Gene=rownames(APGeneInGR)
APGeneInLS$Population="LS"
APGeneInLG$Population="LG"
APGeneInGR$Population="ChM"
IntegratedDEG4Order=cbind(APGeneInLS,APGeneInLG,APGeneInGR)
IntegratedDEG4Order$FC=(IntegratedDEG4Order[,2]+IntegratedDEG4Order[,7]+IntegratedDEG4Order[,12])/3
IntegratedDEG4Order$PV=(IntegratedDEG4Order[,3]+IntegratedDEG4Order[,8]+IntegratedDEG4Order[,13])/3
IntegratedDEG4Order=IntegratedDEG4Order[order(IntegratedDEG4Order$FC,IntegratedDEG4Order$PV),]


resultTable=read.table("GTex/CommonAging8GTEx/AgingGenePattern.txt",header=T)
filtclusters<-resultTable
upc = c(5,10,13,16,22,28)
downc = c(4,6,7,9,14,17,21,31,34) 
UpCommonPattern=filtclusters[filtclusters$cluster%in%upc,]
DnCommonPattern=filtclusters[filtclusters$cluster%in%downc,]
UpCommonAgingGene=unique(UpCommonPattern$genes)
DnCommonAgingGene=unique(DnCommonPattern$genes)

CommonAgingGeneDEGInfo=read.table("GTex/CommonAging8GTEx/LTR_Result_DESeq2.txt",header=T,row.names=1)
CommonAgingGeneDEGInfo=CommonAgingGeneDEGInfo[,c("padj","Symbol")]
colnames(CommonAgingGeneDEGInfo)=c("CommonAgingGenePadj","Symbol")
IntegratedDEG4Order$Symbol=rownames(IntegratedDEG4Order)
IntegratedDEG4Order=IntegratedDEG4Order[,c("FC","PV","Symbol")]

AgingGene8Peter=read.table("Peter/AgeRelatedGene.txt") 

library(ggrepel)
hs_data=merge(IntegratedDEG4Order,CommonAgingGeneDEGInfo,by="Symbol")
hs_data$threshold = as.factor(ifelse(hs_data$Symbol %in% UpCommonAgingGene,"UpCommon8GTEx",ifelse(hs_data$Symbol%in%DnCommonAgingGene,"DnCommon8GTEx",ifelse(hs_data$Symbol %in% AgingGene8Peter[,1],"AgingGene8Peter","NotCommonAgingGene"))))
table(hs_data$threshold)
#----Figure 1H----
t=ggplot(data = hs_data, aes(x = FC, y = -log10(PV), colour=threshold, size=-log10(CommonAgingGenePadj),label =Symbol)) +
  geom_point(alpha=0.4) +
  theme_bw() + scale_size(range = c(0,3))+
  scale_color_manual(values=c("orange","blue","lightgrey","red")) +
  #xlim(c(-5, 5)) + ylim(0,12.5)+
  geom_vline(xintercept=c(-log2(1.2),log2(1.2)),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
    geom_text_repel(
    data = subset(hs_data, hs_data$PV < 0.01 & abs(hs_data$FC) >= log2(1.2)),
    aes(label = Symbol),
    size = 3,
    max.overlaps=5,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("Volcano_CommonAgingGeneFromLHAGene.pdf",width=10)
print(t)
dev.off()


####### visualization for LHA genes #######
setwd("Projects/LLIImmune/")
LHAGene=read.table("LHAGene.txt",header=T)#from GTEX.AgeGene.code.R
LSDEG=read.table("Hainan/DEG8DeseqInLSSample.txt",header=T,row.names=1)
LGDEG=read.table("Hainan/DEG8DeseqInLGSample.txt",header=T,row.names=1)
GRDEG=read.table("Hainan/DEG8DeseqInGRSample.txt",header=T,row.names=1)

APGeneInLS=LSDEG[LHAGene$Gene,c("baseMean","log2FoldChange","pvalue")]
APGeneInLG=LGDEG[LHAGene$Gene,c("baseMean","log2FoldChange","pvalue")]
APGeneInGR=GRDEG[LHAGene$Gene,c("baseMean","log2FoldChange","pvalue")]
all(rownames(APGeneInLS)==rownames(APGeneInLG))
all(rownames(APGeneInLS)==rownames(APGeneInGR))
APGeneInLS$Gene=rownames(APGeneInLS)
APGeneInLG$Gene=rownames(APGeneInLG)
APGeneInGR$Gene=rownames(APGeneInGR)
APGeneInLS$Population="LS"
APGeneInLG$Population="LG"
APGeneInGR$Population="ChM"

LHAGeneDEG4Order=rbind(APGeneInLS,APGeneInLG,APGeneInGR)
LHAGeneDEG4Order=LHAGeneDEG4Order[order(LHAGeneDEG4Order$Gene,LHAGeneDEG4Order$Population),]

LHAGeneDEG4Order=cbind(APGeneInLS,APGeneInLG,APGeneInGR)
LHAGeneDEG4Order$FC=(LHAGeneDEG4Order[,2]+LHAGeneDEG4Order[,7]+LHAGeneDEG4Order[,12])/3
LHAGeneDEG4Order$PV=(LHAGeneDEG4Order[,3]+LHAGeneDEG4Order[,8]+LHAGeneDEG4Order[,13])/3
LHAGeneDEG4Order=LHAGeneDEG4Order[order(LHAGeneDEG4Order$FC,LHAGeneDEG4Order$PV),]


hs_data=data.frame(LHAGeneDEG4Order)
hs_data$threshold = as.factor(ifelse(hs_data$PV > 0.05,"NoSig",ifelse(abs(hs_data$FC)>log2(1.2),ifelse(hs_data$FC>log2(1.2),"SigUp","SigDown"),ifelse(hs_data$FC>0,"Up","Down"))))
write.table(hs_data,file="LHAGeneCombineInfo.txt",sep="\t",quote=F) #output the LHAGenes with average p value and FC

table(hs_data$threshold)
#----Figure 1I----
t=ggplot(data = hs_data, aes(x = FC, y = -log10(PV), colour=threshold, size=log2(baseMean+1), label =Gene)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  scale_color_manual(values=c( "RoyalBlue","blue","red","Salmon")) +
  #xlim(c(-5, 5)) + ylim(0,12.5)+
  geom_vline(xintercept=c(-log2(1.2),log2(1.2)),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
    geom_text_repel(
    data = subset(hs_data, hs_data$PV < 0.01 & abs(hs_data$FC) >= log2(1.2)),
    aes(label = Gene),
    size = 3,
    max.overlaps=5,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("Volcano_LHAGene.pdf",width=9)
print(t)
dev.off()
