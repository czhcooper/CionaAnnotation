


library(FoldGO)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(topGO)
library(dplyr)
library(reshape2)
library(tidyverse)
library(apeglm)
library(Rmisc)
library(topGO)
library(qqplotr)
library(pheatmap)
#GAFReader("/Users/cooper/bio/ciona_s/genome/aniseed/Functional_Annotations/Cisavi_slimTunicate.gaf",geneid_col=6)->csgaf
 #csgaf@annotation->Csgaf
#csGO2geneID<- getAnnotation(csgaf)
#csgeneID2GO<-inverseList(csGO2geneID)


#read.table("/Volumes/data_cooper/hxn/SL/merge.tsv",header = F)->csSLgeneID
#readRDS("/Volumes/data_cooper/run_ssd_cooper/ciona_i/FRY/SL/ciSLgeneID")->ciSLgeneID




GOenrich.cs<-function(data,category,method="fisher",cutOff=0.05,padjust= TRUE,algorithm="classic") {

  require(topGO)
  require(reshape2)
  require(tidyverse)

  if(is.data.frame(data) ==T ){
    genelist<-factor(as.integer(unique(names(csgeneID2GO) ) %in% data[,1]))
  } else {
    genelist<-factor(as.integer(unique(names(csgeneID2GO) ) %in% data))
  }

  names(genelist)<-unique(names(csgeneID2GO))
  GOdata<-new("topGOdata", ontology=category,
              allGenes=genelist,
              annot=annFUN.gene2GO,
              gene2GO=csgeneID2GO,
              geneSel=selection,
              nodeSize=5
  )

  allDEGs_resultFis<-runTest(GOdata,algorithm = algorithm,statistic = method)

  #summarising the results
  allDEGs_Res<- GenTable(GOdata, classic = allDEGs_resultFis,orderBy = "classic", ranksOf = "classic", topNodes = 100)
  as.numeric(allDEGs_Res$classic)->allDEGs_Res$classic
  GOID2Description<-data.frame(GOID=Cigaf$V5,Description=Cigaf$V18)
  GOID2Description[GOID2Description$GOID %in% allDEGs_Res$GO.ID,]->a
  allDEGs_Res[allDEGs_Res$GO.ID %in% a$GOID,]->allDEGs_Res
  allDEGs_Res$Term <-a[match(allDEGs_Res$GO.ID,a$GOID),2]
  allDEGs_Res %>% arrange(desc(Significant))->allDEGs_Res
  allDEGs_Res$Term<-factor(allDEGs_Res$Term,levels =rev(allDEGs_Res$Term) )
  allDEGs_Res$classic<-as.numeric(allDEGs_Res$classic)

  #adjusting pvalue for multiple test
  if (padjust==TRUE){
    allDEGs_Res[p.adjust(allDEGs_Res$classic,method = "BH")< cutOff,]->allDEGs_Res
  } else { allDEGs_Res ->allDEGs_Res}

  return(allDEGs_Res)
}


#GOenrich.cs(csSLgeneID[csSLgeneID$V2>0,],category = "BP")->csSL.BP
#GOenrich.cs(csSLgeneID[csSLgeneID$V2>0,],category = "MF")->csSL.MF
#GOenrich.cs(csSLgeneID[csSLgeneID$V2>0,],category = "CC")->csSL.CC

#csSL.BP$type<-"Biological Process"
#csSL.MF$type<-"Molecular Function"
#csSL.CC$type<-"Cellular Component"

#rbind(csSL.BP[1:10,],csSL.CC[1:10,],csSL.MF[1:10,])->csSL.GO

#ggplot(csSL.GO,aes(x = Term,y=Significant,fill=type))+geom_col()+
#  theme(axis.text.x=element_text(angle=-40, hjust=0),legend.title = element_blank() )+ facet_grid(.~type, scale="free_x")+
#  ylab("Number of Genes")+xlab("GO Term")



#GOenrich(ciSLgeneID,category = "BP")->ciSL.BP
#GOenrich(ciSLgeneID,category = "MF")->ciSL.MF
#GOenrich(ciSLgeneID,category = "CC")->ciSL.CC

#ciSL.BP$type<-"Biological Process"
#ciSL.MF$type<-"Molecular Function"
#ciSL.CC$type<-"Cellular Component"

#rbind(ciSL.BP[1:10,],ciSL.CC[1:10,],ciSL.MF[1:10,])->ciSL.GO

#ggplot(ciSL.GO,aes(x = Term,y=Significant,fill=type))+geom_col()+
#  theme(axis.text.x=element_text(angle=-40, hjust=0),legend.title = element_blank() )+ facet_grid(.~type, scale="free_x")+
#  ylab("Number of Genes")+xlab("GO Term")

#read.delim("~/bio/ciona_s/expression/featureCount/aniseed/gene_feature.txt",row.names = 1,skip = 1)->cs.rawCount


#library(DESeq2)
#library(ggplot2)
#read.table("~/bio/ciona_s/expression/featureCount/aniseed/sampleData",header = T)->cs.sampleData




 # ch1<-"C48"
#  ch2<-"M48"


 # cs.rawCount[,c(grep(paste(ch1,"_","[0-9]",sep = ""),colnames(cs.rawCount)),grep(paste(ch2,"_","[0-9]",sep = ""),colnames(cs.rawCount)))]->cs.rawCount.ch
#  round(cs.rawCount.ch)->cs.rawCount.ch
 # cs.sampleData[cs.sampleData$condition %in% c(ch1,ch2),]->cs.sampleData.ch
#  cs.deseq2Data.ch<-DESeqDataSetFromMatrix(countData =cs.rawCount.ch,colData =cs.sampleData.ch ,design = ~condition)
 # cs.deseq2Data.ch<-DESeq(cs.deseq2Data.ch)

#  deseq2Results <- results(cs.deseq2Data.ch)
 # summary(deseq2Results)
#  res<-as.data.frame(deseq2Results)
  #res$operon<-ifelse(rownames(res) %in% csSLID,"SL","nonSL")
 # res$condition<-ch2
#  res->res.M48


#  rbind(res.M1,res.M24,res.M48)->res.cs
#  res.cs$up<-ifelse(res.cs$log2FoldChange>1,"up","down")
#  na.omit(res.cs)->aa

#  ggplot(aa[abs(aa$log2FoldChange)>1 & aa$padj<0.01 ,],aes(x=condition,fill=operon))+geom_bar()+xlab("")



#  ggplot(res.cs,aes(x=condition,y=up,fill=up))+geom_col()






