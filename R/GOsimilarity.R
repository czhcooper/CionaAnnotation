

length(ci.HH)

ci.HH.bothGO<-ci.HH[unlist(lapply(ci.HH,function(x){
  all(x[[1]] %in% names(cigeneID2GO))
}),use.names = F)]

ci.HH.GO.sim.BP<-lapply(ci.HH.bothGO,function(x){

  if( all(x[[1]] %in% names(cigeneID2GO))) {
    cal_geneSim(x[[1]][1],x[[1]][2],method = "max",simFun = "CoutoResnik",category = "BP")
  } else {return(NA)}

})

ci.operon.bothGO<-ci.operon_2gene[unlist(lapply(ci.operon_2gene,function(x){
  all(x[["GeneID"]] %in% names(cigeneID2GO))
}),use.names = F)]


ci.Operon.GO.sim<-lapply(ci.operon.bothGO,function(x){

  if( all(x[["GeneID"]] %in% names(cigeneID2GO))) {
    cal_geneSim(x[["GeneID"]][1],x[["GeneID"]][2],method = "max",simFun = "Lin",category = "BP")
  } else {return(NA)}

})

for ( i in 1:10000) {
  gene<-sample(names(cigeneID2GO),2,replace = F)
  ci.random.GO.sim[i]<-cal_geneSim(gene[1],gene[2],method = "max",simFun = "Lin",category = "BP")
}

ci.HT.GO.sim<-lapply(ci.HT,function(x){

  if( all(x[[1]] %in% names(cigeneID2GO))) {
    cal_geneSim(x[[1]][1],x[[1]][2],method = "max",simFun = "CoutoResnik",category = "BP")
  } else {return(NA)}

})

ci.TT.GO.sim<-lapply(ci.TT,function(x){

  if( all(x[[1]] %in% names(cigeneID2GO))) {
    cal_geneSim(x[[1]][1],x[[1]][2],method = "max",simFun = "CoutoResnik",category = "BP")
  } else {return(NA)}

})
