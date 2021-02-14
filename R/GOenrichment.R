# Library required
library(topGO)
library(reshape2)
library(tidyverse)
library(ggplot2)


#preparing data
data("cigeneID2GO","Cigaf","ci_kegg","ciGO2geneID","ciGeneAlias2KHID")

#
GOenrich<-function(data,category,method="fisher",cutOff=0.05,padjust= TRUE,algorithm="classic") {

  require(topGO)
  require(reshape2)
  require(tidyverse)

  genelist<-factor(as.integer(unique(names(cigeneID2GO) ) %in% data[,1]))
  names(genelist)<-unique(names(cigeneID2GO))
  GOdata<-new("topGOdata", ontology=category,
              allGenes=genelist,
              annot=annFUN.gene2GO,
              gene2GO=cigeneID2GO,
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

GOplot<- function(data,lim=20){
  require(ggplot2)

  if(dim(data)[1] > lim ){
    p<-ggplot(data[1:lim,],aes(x=Significant,y=Term,fill=classic))+geom_bar(stat="identity") +
      scale_fill_gradient("pvalue",low = "red",high = "blue")+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))+
      xlab("Count")+ylab("")+guides(colour=guide_legend(title = ""))+
      theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank())+theme()
  } else {
    p<-ggplot(data,aes(x=Significant,y=Term,fill=classic))+geom_bar(stat="identity") +
      scale_fill_gradient("pvalue",low = "red",high = "blue")+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))+
      xlab("Count")+ylab("")+guides(colour=guide_legend(title = ""))+
      theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank())+theme()
  }

  return(p)
}








