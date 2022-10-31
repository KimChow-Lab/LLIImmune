####### Defination of the common age-associated genes based on GTEx ######
setwd("Projects/LLIImmune/GTEx")
#######  identify the signficantly changed age-associated genes ######
library(DESeq2)
GTExSampleInfo=read.csv("rawData/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",header=T,sep="\t",check.names=F,row.names=1)
GTExGeneCount=read.table("rawData/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct",header=T,sep="\t",row.names=1,skip=2,check.names=F)
WBSample=GTExSampleInfo[GTExSampleInfo$SMTSD %in% "Whole Blood",]
WBList=intersect(rownames(WBSample),colnames(GTExGeneCount))
WBGeneCount=GTExGeneCount[,c("Description",WBList)]
AgeInfo=read.table("rawData/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",header=T,row.names=1)
AgeInfo$Prefix=rownames(AgeInfo)
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(WBList),'-')))
result=data.frame("SampleID"=WBList,"Prefix"=paste0(TmpInfo[,1],"-",TmpInfo[,2],sep=""))
SampleInfo=merge(result,AgeInfo,by="Prefix")
rownames(SampleInfo)=SampleInfo$SampleID
SampleInfo=SampleInfo[WBList,]
all(rownames(SampleInfo)==colnames(WBGeneCount)[-1]) #TRUE
SampleInfo$Gender=ifelse(SampleInfo$SEX==1,"Male","Female")
SampleInfo$DTHHRDY = as.character(SampleInfo$DTHHRDY)
WBGeneExpr=WBGeneCount
WBGeneExpr$Description=NULL
SampleInfo$Prefix=NULL
SampleInfo$SampleID=NULL
SampleInfo$SEX=NULL
WBGeneExpr=[rowMeans(WBGeneExpr)>0,]
dds <- DESeqDataSetFromMatrix(countData=WBGeneExpr, colData=SampleInfo, design=~Gender+DTHHRDY+AGE) 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds = DESeq(dds,test="LRT",reduced= ~Gender+DTHHRDY) #obtain the age-associted DEGs
res_LRT <- results(dds)
res_LRT<-res_LRT[complete.cases(res_LRT),]
res_LRT$Ensemble=rownames(res_LRT)
Symbol=data.frame(Ensemble=rownames(WBGeneCount),Symbol=WBGeneCount$Description)
res_LRT_Symbol=merge(data.frame(res_LRT),Symbol,by="Ensemble")
sigGene=res_LRT_Symbol[res_LRT_Symbol$padj<0.05,] #definition of the significant changed age-associted DEGs
write.table(res_LRT_Symbol,file="CommonAging8GTEx/LTR_Result_DESeq2.txt",row.names=F,sep="\t",quote=F)

#######  Identify the adjust gene expression matrix ######
assay(vsd) = limma::removeBatchEffect(assay(vsd),covariates=model.matrix(~Gender+DTHHRDY,data=SampleInfo),design = model.matrix(~AGE,data=SampleInfo))
vsd.adjust=format(as.data.frame(assay(vsd)),digits=3)
Symbol=data.frame(Ensemble=rownames(WBGeneCount),Symbol=WBGeneCount$Description)
vsd.adjust$Ensemble=rownames(vsd.adjust)
vsd.adjust.symbol=merge(Symbol,vsd.adjust,by="Ensemble")
rownames(vsd.adjust.symbol)=vsd.adjust.symbol$Ensemble
vsd.adjust.symbol$Ensemble=NULL
write.table(vsd.adjust.symbol,"GTEx/vsd.adjusted.txt",sep="\t",row.names=TRUE,quote=F) # vst transformed and adjusted expression matrix

#######  Identify the common age-associated genes  ######
library(DEGreport)
library(RColorBrewer)
WBSampleInfo=SampleInfo
all(rownames(WBSampleInfo)==colnames(sigGeneExpr)) # forever TRUE lalalala
WBSampleInfo$AGE=factor(WBSampleInfo$AGE,levels=c("20-29","30-39","40-49","50-59","60-69","70-79"))
res <- degPatterns(sigGeneExpr, WBSampleInfo, time = "AGE",minc = 1)

write.table(res$normalized,file="CommonAging8GTEx/AgingGenePattern.txt",quote=F,sep="\t")
resultTable=read.table("CommonAging8GTEx/AgingGenePattern.txt",header=T)
pdf("CommonAging8GTEx/FullAgingGenePattern.pdf",height=20,width=30)
degPlotCluster(resultTable,time="AGE",color="AGE",boxes = FALSE,min_genes =  1)+geom_line(aes_string(group="genes"),alpha=0.5)+scale_color_brewer(palette="YlOrRd")+theme_bw()
dev.off()

filtclusters<-resultTable
upc = c(5,10,13,16,22,28)
downc = c(4,6,7,9,14,17,21,31,34) 
UpCommonPattern=filtclusters[filtclusters$cluster%in%upc,]
DnCommonPattern=filtclusters[filtclusters$cluster%in%downc,]
pdf("CommonAging8GTEx/UpCommonPattern.pdf",height=4,width=7)
degPlotCluster(UpCommonPattern,time="AGE",color="AGE",boxes = FALSE)+geom_line(aes_string(group="genes"),alpha=0.5)+scale_color_brewer(palette="YlOrRd")+theme_bw()
dev.off()
pdf("CommonAging8GTEx/DnCommonPattern.pdf",height=4,width=7)
degPlotCluster(DnCommonPattern,time="AGE",color="AGE",boxes = FALSE)+geom_line(aes_string(group="genes"),alpha=0.5)+scale_color_brewer(palette="YlGn")+theme_bw()
dev.off()

UpCommonAgingGene=unique(UpCommonPattern$genes)
DnCommonAgingGene=unique(DnCommonPattern$genes)
CommonAgingGene=filtclusters[filtclusters$cluster%in% c(upc,downc),c("genes","cluster")]
CommonAgingGene=unique(CommonAgingGene)
CommonAgingGene$Pattern=ifelse(CommonAgingGene$cluster%in%upc,"CommonUp","CommonDn")
write.table(CommonAgingGene,file="CommonAging8GTEx/CommonAgingGene.txt",sep="\t",quote=F,row.names=F)



#######  Visualizatoin the common age-associated genes  ######
WBGeneExprAdj=read.table("GTEx/vsd.adjusted.txt",header=T,row.names=1,check.names=F)
WBSampleInfo=read.table("GTEx/rawData/WholeBloodSampleInfo.txt",header=T,row.names=1)
table(WBSampleInfo$AGE)
sampleList=intersect(rownames(WBSampleInfo),colnames(WBGeneExprAdj))
length(sampleList)
WBGeneExprAdj=WBGeneExprAdj[,c("Symbol",sampleList)]

geneList="MSR1"
ShareListExprInGTEx=WBGeneExprAdj[WBGeneExprAdj$Symbol%in%geneList,]
rownames(ShareListExprInGTEx)=ShareListExprInGTEx$Symbol
ShareListExprInGTEx$Symbol=NULL
all(colnames(ShareListExprInGTEx)==rownames(WBSampleInfo)) #TRUE

ShareListExprWithPhenotype=cbind(WBSampleInfo,t(ShareListExprInGTEx))
t=ggplot(ShareListExprWithPhenotype, aes(AGE, y=log2(`MSR1`+1),color=as.character(AGE)))+
geom_boxplot()+labs(title="",x="", y = "")+ theme_bw()+
scale_color_brewer(palette="Set2")+
stat_compare_means(method = "anova")+
theme(legend.position="none")+ geom_jitter(shape=1, position=position_jitter(0.2))+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90),axis.text.y = element_text(size=rel(1.0),face="bold"))
pdf("SpecificGene/MSR1InGETx.pdf",height=3,width=5)
print(t)
dev.off()
