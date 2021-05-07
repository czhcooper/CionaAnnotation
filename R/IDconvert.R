##' convert gene ID among KH,KY,gene name, symbol,
##'
##'
##' @title Gff2GeneTable
##' @param gffFile GFF file
##' @param compress compress file or not
##' @return file save.
##' @export
##' @author Chen Zaohuang



convertgeneID<-function(gene,fromType="KYgeneID",toType="KHgeneID",Db=KY_geneID2Symbol){

  #Library required
  library(tidyverse)

  if(all(c(fromType,toType) %in% c("KYgeneID","KHgeneID","EnGeneID","ALIAS","ENTREZID","geneName","hs_Symbol","HGNC","ci_Symbol") )){
    if(class(gene)=="character"){

      gene<-sub("(KY2019:|KH2012:|KH2013:)","",gene)   #remove prefix of gene id

      id_res<-lapply(gene, function(x){
        unique(KY_2_KH_Orthology[,toType][KY_2_KH_Orthology[,fromType] %in% x])

      })

      names(id_res)<-gene
      return(id_res)

    } else if (class(gene)=="data.frame"){
      # This function requires complement in future.

    } else (cat("Input of \'gene\' should be a class of character or data.frame"))


  } else { cat("ID type should be terms among (\"KYgeneID\",\"KHgeneID\",\"ENSEMBL\",\"HGNC\",\"geneName\",\"ci_Symbol\",\"hs_Symbol\"). ")}



}










##' convert gene ID among KH,KY,gene name, symbol,
##'
##'
##' @title Gff2GeneTable
##' @param gffFile GFF file
##' @param compress compress file or not
##' @return file save.
##' @export
##' @author Chen Zaohuang



get_KY_GO<-function(gene){
  require(clusterProfiler)

  gene<-sub("KY2019:|KH2012:|KH2013","",gene)
  KY.GO<-lapply(gene, function(x){
    if(x %in% KY_2_KH_Orthology[,"KYgeneID"]){
      KH.id<-unique(KY_2_KH_Orthology[,"KHgeneID"][KY_2_KH_Orthology[,"KYgeneID"] %in% x])
      unique( unlist(cigeneID2GO[sub("KH2012:","",names(cigeneID2GO)) %in% KH.id],use.names = F))

    } else {
      gene.name<-unique(KY_2_Symbol[,"hs_Symbol"][KY_2_Symbol[,"Parent"] %in% gene])
      gene.GO<- bitr(gene.name,fromType = "SYMBOL",toType = "GO",OrgDb = "org.Hs.eg.db")
      return(gene.GO$GO)
    }
  })
  names(KY.GO)<-gene
  return(KY.GO)
}



















