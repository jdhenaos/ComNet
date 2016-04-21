library(igraph)

ADN <- read.graph("./ALZ.txt",format = "ncol")
PDN <- read.graph("./PRK.txt",format = "ncol")
MSN <- read.graph("./ESMU.txt",format = "ncol")

MoAD <- cluster_fast_greedy(ADN)
MoPD <- cluster_fast_greedy(PDN)
MoMS <- cluster_fast_greedy(MSN)

subA1 <- induced.subgraph(ADN,vids = as.vector(unlist(MoAD[27])))
subP1 <- induced.subgraph(PDN,vids = as.vector(unlist(MoPD[3])))

if(length(unlist(MoAD[27])) > length(unlist(MoPD[3]))){
  print("subA1 > subP1")
  dif <- graph.difference(subA1,subP1)
}else if(length(unlist(MoAD[27])) < length(unlist(MoPD[3]))){
  print("subP1 > subA1")
  dif <- graph.difference(subP1,subA1)
}else{
  print("subA1 = subP1")
  dif <- graph.difference(subA1,subP1)
    if(ecount(dif) == 0 && vcount(dif) == 0){
      break
    }
}

MoAD <- cluster_edge_betweenness(ADN,directed = F)
MoPD <- cluster_edge_betweenness(PDN,directed = F)
MoMS <- cluster_edge_betweenness(MSN,directed = F)