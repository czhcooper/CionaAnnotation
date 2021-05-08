##' calculate the expression correlation between  gene pairs
##'
##'
##' @title Expression Similarity
##' @param
##' @param
##' @return Correlation coefficient
##' @author Chen Zaohuang
##' @export



#expressionData<-KH.tpm
#group<-data.frame(samples=colnames(KH.tpm)[-27])
#group$condition<-c(rep("0h",5),rep("120h",4),rep("144h",3),rep("24h",4),rep("48h",6),rep("96h",4))

# HT gene pairs
ci.ht



#operon gene pairs

ci.operon.genePairs<-bind_rows(ci.operon_2gene) %>% select(Praent,GeneID)
colnames(ci.operon.genePairs)[1]<-"typeID"


gene_expressionSim<-function(genePairs,expressionData,group){
    expressionData <-expressionData %>% select(-GeneID) %>% apply( 1, normalized.tpm) %>% t() %>% as.data.frame() %>% add_column(GeneID=expressionData$GeneID)
     ea<-expressionData[match(genePairs$Geneid,expressionData$GeneID),-which(names(expressionData) %in% c("GeneID"))]
      gea<-cbind(genePairs,ea)
      aa<-gea[,names(gea) %in% c(group$samples)]
      aa$typeID<-gea$typeID
    cc<-split(aa,aa$typeID)

    cc<-lapply(cc,function(x){
      t(x[,-which(names(x) %in% c("typeID"))])
    })

   ee<- lapply(cc,function(x){
      aa<- as.data.frame(x)
      aa$group<-group$condition[match(rownames(aa),group$samples)]
      aa<-split(aa[,c(1,2)],aa$group,drop = T)
      dd<-lapply(aa,function(x){
        a.cor<-cor.test(x[,1],x[,2])
        a.cor<-data.frame(cor=a.cor$estimate,p.value=a.cor$p.value,row.names = F)
        })
       bind_rows(dd,.id="group")

    })
  #combine each gene pair's expression correlation
      es.cor<-bind_rows(ee,.id = "typeID")

   return(es.cor)

}








cor.df<-function(x){
  cor(x[,1],x[,2])

}

aa<-lapply(ci.HH,function(x){
  x[c(1,2,3,4,5)]
})




#' @title Z-score Normalization
#' @param X should be a numeric data.frame or matrix
#' @author Chen Zaohuang
#' @export
normalized.tpm<-function(x){
  (x-mean(x))/sd(x)
}






