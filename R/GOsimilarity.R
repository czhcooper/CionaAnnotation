
#libaray requaried
library(GOSim)



# This function refers to the "GOsim" R package
precomputeTermSims<- fucntion (x, y, similarityTerm="JiangConrath", verbose=FALSE){
  if(verbose)
    message("precomputing term similarities ...")
    gotermsx<-as.vector(unique(unlist(cigeneID2GO[x])))
    gotermsy<-as.vector(unique(unlist(cigeneID2GO[y])))

      STerm<-matrix(0, nrow=length(gotermsx), ncol=length(gotermsy))
      rownames(STerm)=gotermsx
      colnames(STerm)=gotermsy
      for(i in 1:length(gotermsx)){
        for(j in 1:length(gotermsy)){
          STerm[i,j]<-upper.tri(getTermSim(c(gotermsx[i],gotermsy[j]),method = "relevance"))
          }
      }

}
