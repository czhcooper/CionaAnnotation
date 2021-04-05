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




  if(class(gene)=="character"){

    gene<-sub("(KY2019:|KH2012:|KH2013:)","",gene)   #remove prefix of gene id

    id_res<-lapply(gene, function(x){
     KY_2_KH_geneID[,toType][KY_2_KH_geneID[,fromType] %in% x]

   })

   names(id_res)<-gene
   return(id_res)

  } else if (class(gene)=="data.frame"){


  } else (cat("Input of \'gene\' should be a class of character or data.frame"))

}




get_KY_GO<-function(gene){

  KY.GO<-lapply(gene, function(x){
    KH.id<-KY_2_KH_geneID[,toType][KY_2_KH_geneID[,fromType] %in% x]
   unique( unlist(cigeneID2GO[sub("KH2012:","",names(cigeneID2GO)) %in% KH.id],use.names = F))

  })
  names(KY.GO)<-gene

  return(KY.GO)


}
