library(igraph)

ADN <- read.graph("ALZ.txt",format = "ncol")
PDN <- read.graph("PRK.txt",format = "ncol")
MSN <- read.graph("ESMU.txt",format = "ncol")

MoAD <- cluster_fast_greedy(ADN)
MoPD <- cluster_fast_greedy(PDN)
MoMS <- cluster_fast_greedy(MSN)

if(length(MoAD) >= length(MoPD)){
  graph1 <- MoAD
  graph2 <- MoPD
  bigGraph <- ADN
  smallGraph <- PDN
}else{
  graph1 <- MoPD
  graph2 <- MoAD
  bigGraph <- PDN
  smallGraph <- ADN
}

iqual <- c()
  
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
        if(length(iqual) == 0){
          iqual <- c(as.vector(names(dif[2]),mode = "character"))
        }else{
          iqual <- c(iqual,as.vector(names(dif[2]),mode = "character"))
        }
      }
    }
  }
}

####################################################

iqual <- data.frame()

subA <- induced.subgraph(ADN,vids = as.vector(unlist(MoAD[27])))
subB <- induced.subgraph(ADN,vids = as.vector(unlist(MoAD[27])))

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
      if(length(iqual) == 0){
        iqual <- rbind(names(dif[2]))
      }else{
        iqual <- rbind(iqual,names(dif[2]))
      }
    }
}

MoAD <- cluster_edge_betweenness(ADN,directed = F)
MoPD <- cluster_edge_betweenness(PDN,directed = F)
MoMS <- cluster_edge_betweenness(MSN,directed = F)