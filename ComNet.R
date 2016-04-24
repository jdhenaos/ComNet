library(igraph)

CommonModules <- function(G1,G2){
  MoA <- cluster_fast_greedy(G1)
  MoB <- cluster_fast_greedy(G2)
  
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
      
      if(length(names(subA[2])) > length(names(subB[2]))){
        dif <- graph.difference(subA,subB)
      }else if(length(names(subA[2])) < length(names(subB[2]))){
        dif <- graph.difference(subB,subA)
      }else{
        dif <- graph.difference(subA,subB)
        if(names(subA[2]) == names(subB[2]) && ecount(dif) == 0){
          iqual[[counter]] <- names(dif[2])
          counter <- counter + 1
        }
      }
    }
  }
  return(iqual)
}

ADN <- read.graph("ALZ.txt",format = "ncol")
PDN <- read.graph("PRK.txt",format = "ncol")
MSN <- read.graph("ESMU.txt",format = "ncol")

AvsP <- CommonModules(ADN,PDN)
PvsM <- CommonModules(PDN,MSN)
AvsM <- CommonModules(ADN,MSN)

####################################################

MoA <- cluster_fast_greedy(ADN)
MoB <- cluster_fast_greedy(PDN)

if(length(MoAD) >= length(MoA)){
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
    
    if(length(names(subA[2])) > length(names(subB[2]))){
      print("subA1 > subP1")
      dif <- graph.difference(subA,subB)
    }else if(length(names(subA[2])) < length(names(subB[2]))){
      print("subP1 > subA1")
      dif <- graph.difference(subB,subA)
    }else{
      print("subA1 = subP1")
      dif <- graph.difference(subA,subB)
      if(names(subA[2]) == names(subB[2]) && ecount(dif) == 0){
        iqual[[counter]] <- names(dif[2])
        counter <- counter + 1
      }
    }
  }
}

####################################################