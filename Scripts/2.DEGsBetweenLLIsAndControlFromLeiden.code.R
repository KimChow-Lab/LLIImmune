setwd("Projects/LLIImmune/Leiden")
#obtain the sample information from GEO
gse <- GEOquery::getGEO("GSE16717", GSEMatrix = TRUE)
phenotype=gse[[1]]@phenoData@data[,c("age:ch1","family:ch1","gender:ch1","group:ch1")]
colnames(phenotype)=c("Age","family","Sex","BioGroup")
write.table(phenotype,file="/sampleInfo.txt",sep="\t",quote=F)

####### obtain the DEGs from Leiden by limma package

data=read.table("rawData/combine.txt",header=T,row.names=1,sep="\t") #the normalized gene expression were downloaded from GSE16717
SampleInfor=read.table("sampleInfo.txt",header=T,row.names=1)
CompairSample=SampleInfor[SampleInfor$BioGroup%in%c("control","longLived"),]
CompairSample=CompairSample[order(CompairSample$BioGroup),]
sampleList=intersect(rownames(CompairSample),colnames(data))
exprSet=data[,sampleList]
all(rownames(CompairSample)==colnames(exprSet))


####### obtain the DEGs from Leiden by limma package
f <- factor(CompairSample$BioGroup, levels=c("control","longLived"))
design <- model.matrix(~0+f+Sex,CompairSample)
colnames(design) <- c("control","longLived","Sex")
rownames(design) = rownames(CompairSample)
fit <- lmFit(exprSet, design)
cont.matrix <- makeContrasts("control-longLived",  levels=design)
fit2  <- contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit2)
top.table <- topTable(fit3, sort.by = "P", n = Inf)
head(top.table, 20)

anno=read.csv("GPL2895.annot",header=T,sep="\t",check.names=F) #gene annotation file
anno=anno[,c("ID","Gene symbol","Gene ID")]
colnames(anno)=c("ID","GeneSymbol","GeneID")
top.table$ID=rownames(top.table)
degInfo=merge(top.table,anno,by="ID")
sigGene=degInfo[degInfo$P.Value<0.01,]
write.table(sigGene,file="sigGene8limma.txt",sep="\t",quote=F) 


####### Output the gene experssion matrix 
probleID=intersect(sigGene$ID,rownames(exprSet))
length(probleID)
sigGeneExpr=exprSet[probleID,]
all(rownames(CompairSample)==colnames(sigGeneExpr))
sigGeneExpr$ID=rownames(sigGeneExpr)
sigGeneExprSymbol=merge(anno,sigGeneExpr,by="ID")
sigGeneExprSymbol=sigGeneExprSymbol[sigGeneExprSymbol$GeneID!="",]
sigGeneExprSymbol$ID=NULL
sigGeneExprSymbol$GeneID=NULL
write.table(sigGeneExprSymbol,file="sigGeneExpr.txt",sep="\t",quote=F,row.names=F)



####### expresion for the focused genes ##############
data=read.table("rawData/combine.txt",header=T,row.names=1,sep="\t")
SampleInfor=read.table("sampleInfo.txt",header=T,row.names=1)
CompairSample=SampleInfor[SampleInfor$BioGroup%in%c("control","longLived"),]
CompairSample=CompairSample[order(CompairSample$BioGroup),]
sampleList=intersect(rownames(CompairSample),colnames(data))
exprSet=data[,sampleList]
CompairSample=CompairSample[sampleList,]
all(rownames(CompairSample)==colnames(exprSet))

anno=read.csv("GPL2895.annot",header=T,sep="\t",check.names=F)
anno=anno[,c("ID","Gene symbol","Gene ID","GenBank Accession")]
colnames(anno)=c("ID","GeneSymbol","GeneID","GBAC")
exprSet$ID=rownames(exprSet)
exprSetSymbol=merge(anno,exprSet,by="ID")
exprSetSymbol$ID=NULL
exprSetSymbol$GeneID=NULL
all(rownames(CompairSample)==colnames(exprSetSymbol)[-c(1:2)])

FocusedGene=c("MSR1")
FocusedGeneExpr=exprSetSymbol[exprSetSymbol$GeneSymbol%in%FocusedGene,]
dim(FocusedGeneExpr)
rownames(FocusedGeneExpr)=paste0(FocusedGeneExpr$GeneSymbol,"_",FocusedGeneExpr$GBAC,sep="")
FocusedGeneExpr$GeneSymbol=NULL
FocusedGeneExpr$GBAC=NULL
all(rownames(CompairSample)==colnames(FocusedGeneExpr))
GraphData=cbind(CompairSample,t(FocusedGeneExpr))
library(reshape2)
GraphData=melt(GraphData,id=c(1:4))
colnames(GraphData)=c("Age","family","Sex","BioGroup","TranscriptID","Expr")
#GraphData=GraphData[GraphData$Expr>1.5,] #only for TRABD2A
#GraphData=GraphData[GraphData$Expr>2,] #only for LMO2
GraphData=na.omit(GraphData)
t=ggplot(GraphData, aes(BioGroup, y=log2(Expr),color=BioGroup))+
geom_boxplot(outlier.shape = NA)+labs(title="",x="", y = "")+ theme_bw()+
scale_color_brewer(palette="Set2")+
geom_point(aes(color=BioGroup),size=1.5, alpha=0.5,show.legend=TRUE, position = position_jitterdodge(jitter.width =0.15,dodge.width = 0.75))+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0),face="bold"))+
facet_wrap(.~TranscriptID,scale="free_y")

pdf("Leiden/specificGene/MSR1InLeiden.pdf",width=5,height=3)
print(t)
dev.off()

