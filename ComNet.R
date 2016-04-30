library(igraph)

CommonModules <- function(G1,G2,method,iter=1000){
  
  if(method == "fgr"){
    MoA <- cluster_fast_greedy(G1)
    MoB <- cluster_fast_greedy(G2)
  }else if(method == "ebet"){
    MoA <- cluster_edge_betweenness(G1,directed = F)
    MoB <- cluster_edge_betweenness(G2,directed = F)
  }else if(method == "leig"){
    MoA <- cluster_leading_eigen(G1,options = list(maxiter = iter))
    MoB <- cluster_leading_eigen(G2,options = list(maxiter = iter))
  }else if(method == "walk"){
    MoA <- cluster_walktrap(G1)
    MoB <- cluster_walktrap(G2)
  }
  
  
  if(length(MoA) >= length(MoB)){
    graph1 <- MoA
    graph2 <- MoB
    bigGraph <- G1
    smallGraph <- G2
  }else{
    graph1 <- MoB
    graph2 <- MoA
    bigGraph <- G2
    smallGraph <- G1
  }
  
  iqual <- vector(mode = 'list')
  counter <- 1
  
  for(i in 1:length(graph1)){
    subA <- induced.subgraph(bigGraph,vids = as.vector(unlist(graph1[i])))
    
    for(j in 1:length(graph2)){
      subB <- induced.subgraph(smallGraph,vids = as.vector(unlist(graph2[j])))
      
      if(ecount(subA) > 0 && ecount(subB) > 0){
        if(length(names(subA[2])) == length(names(subB[2]))){
          dif <- graph.difference(subA,subB)
          if(names(subA[2]) == names(subB[2]) && ecount(dif) == 0){
            iqual[[counter]] <- names(dif[2])
            counter <- counter + 1
          }
        }
      }
    }
  }
  return(iqual)
}

SubModules <- function(G1,G2,method,iter=1000){
  common <- function(G1,G2,ModMa,ModMe){
    subA <- induced.subgraph(G1,vids = as.vector(ModMa))
    subB <- induced.subgraph(G2,vids = as.vector(ModMe))
    
    if(ecount(graph.intersection(subA,subB,keep.all.vertices = F)) == ecount(subA) && 
       vcount(graph.intersection(subA,subB,keep.all.vertices = F)) == vcount(subA)){
      l <- list(G1=ModMa,G2=ModMe)
    }else if(ecount(graph.intersection(subA,subB,keep.all.vertices = F)) == ecount(subB) && 
             vcount(graph.intersection(subA,subB,keep.all.vertices = F)) == vcount(subB)){
      l <- list(G1=ModMa,G2=ModMe)
    }
    return(l)
  }
  
  if(method == "fgr"){
    ModA <- cluster_fast_greedy(G1)
    ModB <- cluster_fast_greedy(G2)
  }else if(method == "ebet"){
    ModA <- cluster_edge_betweenness(G1,directed = F)
    ModB <- cluster_edge_betweenness(G2,directed = F)
  }else if(method == "leig"){
    ModA <- cluster_leading_eigen(G1,options = list(maxiter = iter))
    ModB <- cluster_leading_eigen(G2,options = list(maxiter = iter))
  }else if(method == "walk"){
    ModA <- cluster_walktrap(G1)
    ModB <- cluster_walktrap(G2)
  }
  
  counter <- 1
  total <- list()
  
  for(i in 1:length(ModA)){
    ModMa <- as.vector(unlist(ModA[i]))
    for(j in 1:length(ModB)){
      ModMe <- as.vector(unlist(ModB[j]))
      if(length(ModMa) != length(ModMe)){
        if(is.na(table(ModMa %in% ModMe)[2]) != TRUE && table(ModMa %in% ModMe)[2] == length(ModMe)){
          t <- common(G1,G2,ModMa,ModMe)
          total[[counter]] <- t
          counter <- counter + 1
        }else if(is.na(table(ModMe %in% ModMa)[2]) != TRUE && table(ModMe %in% ModMa)[2] == length(ModMa)){
          t <- common(G1,G2,ModMa,ModMe)
          total[[counter]] <- t
          counter <- counter + 1
        }
      }
    }
  }
  return(total)
}

#################################################################

ADN <- read.graph("ALZ.txt",format = "ncol")
PDN <- read.graph("PRK.txt",format = "ncol")
MSN <- read.graph("ESMU.txt",format = "ncol")

SubAvsSubP <- SubModules(ADN,PDN,method = "fgr")
SubAvsSubM <- SubModules(ADN,MSN,method = "fgr")
SubPvsSubM <- SubModules(PDN,MSN,method = "fgr")

SubAvsSubP2 <- SubModules(ADN,PDN,method = "walk")
SubAvsSubM2 <- SubModules(ADN,MSN,method = "walk")
SubPvsSubM2 <- SubModules(PDN,MSN,method = "walk")

SubAvsSubP3 <- SubModules(ADN,PDN,method = "leig",iter = 50000)
SubAvsSubM3 <- SubModules(ADN,MSN,method = "leig",iter = 50000)
SubPvsSubM3 <- SubModules(PDN,MSN,method = "leig",iter = 50000)

################################################################

AvsP <- CommonModules(ADN,PDN,method = "fgr")
PvsM <- CommonModules(PDN,MSN,method = "fgr")
AvsM <- CommonModules(ADN,MSN,method = "fgr")

AvsP2 <- CommonModules(ADN,PDN,method = "walk")
PvsM2 <- CommonModules(PDN,MSN,method = "walk")
AvsM2 <- CommonModules(ADN,MSN,method = "walk")

AvsP3 <- CommonModules(ADN,PDN,method = "leig",iter = 50000)
PvsM3 <- CommonModules(PDN,MSN,method = "leig",iter = 50000)
AvsM3 <- CommonModules(ADN,MSN,method = "leig",iter = 50000)

AvsP4 <- CommonModules(ADN,PDN,method = "ebet")
PvsM4 <- CommonModules(PDN,MSN,method = "ebet")
AvsM4 <- CommonModules(ADN,MSN,method = "ebet")

DAN <- read.graph("AD.txt",format = "ncol")
PAN <- read.graph("PD.txt",format = "ncol")
MAN <- read.graph("Merge.txt",format = "ncol")

PDAD <- CommonModules(DAN,PAN,method = "fgr")
MEAD <- CommonModules(DAN,MAN,method = "fgr")
MEPD <- CommonModules(DAN,MAN,method = "fgr")

####################################################

SubC <- function(Ma,Md,Me){
  if((is.na(table(as.vector(unlist(Ma)) %in% as.vector(unlist(Me)))[2]) != TRUE &&
      table(as.vector(unlist(Ma)) %in% as.vector(unlist(Me)))[2] == length(as.vector(Me)))
     &&
     (is.na(table(as.vector(unlist(Ma)) %in% as.vector(unlist(Md)))[2]) != TRUE &&
      table(as.vector(unlist(Ma)) %in% as.vector(unlist(Me)))[2] == length(as.vector(Md)))
     &&
     (is.na(table(as.vector(unlist(Md)) %in% as.vector(unlist(Me)))[2]) != TRUE &&
      table(as.vector(unlist(Md)) %in% as.vector(unlist(Me)))[2] == length(as.vector(Me)))
     ){
    return("existe")
  }
}

AD <- cluster_walktrap(ADN)
PD <- cluster_walktrap(PDN)
MS <- cluster_walktrap(MSN)

for(i in 1:length(AD)){
  for(j in 1:length(PD)){
    for(k in 1:length(MS)){
      if((length(as.vector(unlist(AD[i]))) < length(as.vector(unlist(PD[j]))) &&
         length(as.vector(unlist(AD[i]))) < length(as.vector(unlist(MS[k])))) &&
         length(as.vector(unlist(MS[k]))) < length(as.vector(unlist(PD[j])))){
        t <- SubC(PD[j],MS[k],AD[i])
        if(is.null(t) != TRUE){
          a <- AD[i]
          b <- PD[j]
          c <- MS[k]
          stop("aparecio")
        }
      }else if((length(as.vector(unlist(AD[i]))) < length(as.vector(unlist(PD[j]))) &&
                length(as.vector(unlist(AD[i]))) < length(as.vector(unlist(MS[k])))) &&
               length(as.vector(unlist(PD[j]))) < length(as.vector(unlist(MS[k])))){
        t <- SubC(MS[k],PD[j],AD[i])
        if(is.null(t) != TRUE){
          a <- AD[i]
          b <- PD[j]
          c <- MS[k]
          stop("aparecio")
        }
      }else if((length(as.vector(unlist(PD[j]))) < length(as.vector(unlist(AD[i]))) &&
                length(as.vector(unlist(PD[j]))) < length(as.vector(unlist(MS[k])))) &&
               length(as.vector(unlist(AD[i]))) < length(as.vector(unlist(MS[k])))){
        t <- SubC(MS[k],AD[i],PD[j])
        if(is.null(t) != TRUE){
          a <- AD[i]
          b <- PD[j]
          c <- MS[k]
          stop("aparecio")
        }
      }else if((length(as.vector(unlist(PD[j]))) < length(as.vector(unlist(AD[i]))) &&
                length(as.vector(unlist(PD[j]))) < length(as.vector(unlist(MS[k])))) &&
               length(as.vector(unlist(MS[k]))) < length(as.vector(unlist(AD[i])))){
        t <- SubC(AD[i],MS[k],PD[j])
        if(is.null(t) != TRUE){
          a <- AD[i]
          b <- PD[j]
          c <- MS[k]
          stop("aparecio")
        }
      }else if((length(as.vector(unlist(MS[k]))) < length(as.vector(unlist(AD[i]))) &&
                length(as.vector(unlist(MS[k]))) < length(as.vector(unlist(PD[j])))) &&
               length(as.vector(unlist(PD[j]))) < length(as.vector(unlist(AD[i])))){
        t <- SubC(AD[i],PD[j],MS[k])
        if(is.null(t) != TRUE){
          a <- AD[i]
          b <- PD[j]
          c <- MS[k]
          stop("aparecio")
        }
      }else if((length(as.vector(unlist(MS[k]))) < length(as.vector(unlist(AD[i]))) &&
                length(as.vector(unlist(MS[k]))) < length(as.vector(unlist(PD[j])))) &&
               length(as.vector(unlist(AD[i]))) < length(as.vector(unlist(PD[j])))){
        t <- SubC(PD[j],AD[i],MS[k])
        if(is.null(t) != TRUE){
          a <- AD[i]
          b <- PD[j]
          c <- MS[k]
          stop("aparecio")
        }
      }
    }
  }
}

####################################################
common <- function(G1,G2,ModMa,ModMe){
  subA <- induced.subgraph(G1,vids = as.vector(ModMa))
  subB <- induced.subgraph(G2,vids = as.vector(ModMe))
  
  if(ecount(graph.intersection(subA,subB,keep.all.vertices = F)) == ecount(subA) && 
     vcount(graph.intersection(subA,subB,keep.all.vertices = F)) == vcount(subA)){
    l <- list(G1=ModMa,G2=ModMe)
  }else if(ecount(graph.intersection(subA,subB,keep.all.vertices = F)) == ecount(subB) && 
           vcount(graph.intersection(subA,subB,keep.all.vertices = F)) == vcount(subB)){
    l <- list(G1=ModMa,G2=ModMe)
  }
  return(l)
}

G1 <- ADN
G2 <- PDN

ModA <- cluster_fast_greedy(G1)
ModB <- cluster_fast_greedy(G2)

counter <- 1
total <- list()

for(i in 1:length(ModA)){
  ModMa <- as.vector(unlist(ModA[i]))
  for(j in 1:length(ModB)){
    ModMe <- as.vector(unlist(ModB[j]))
    if(length(ModMa) != length(ModMe)){
      if(is.na(table(ModMa %in% ModMe)[2]) != TRUE && table(ModMa %in% ModMe)[2] == length(ModMe)){
        t <- common(G1,G2,ModMa,ModMe)
        total[[counter]] <- t
        counter <- counter + 1
      }else if(is.na(table(ModMe %in% ModMa)[2]) != TRUE && table(ModMe %in% ModMa)[2] == length(ModMa)){
        t <- common(G1,G2,ModMa,ModMe)
        total[[counter]] <- t
        counter <- counter + 1
      }
    }
  }
}

###############################################################################