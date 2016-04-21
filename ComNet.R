library(igraph)

ADN <- read.graph("./ALZ.txt",format = "ncol")
PDN <- read.graph("./PRK.txt",format = "ncol")
MSN <- read.graph("./ESMU.txt",format = "ncol")

if(method=="fg"){
  MoAD <- cluster_fast_greedy(ADN)
  MoPD <- cluster_fast_greedy(PDN)
  MoMS <- cluster_fast_greedy(MSN)
}else if(method=="eb"){
  MoAD <- cluster_edge_betweenness(ADN,directed = F)
  MoPD <- cluster_edge_betweenness(PDN,directed = F)
  MoMS <- cluster_edge_betweenness(MSN,directed = F)
}

MoAD <- cluster_fast_greedy(ADN)
MoPD <- cluster_fast_greedy(PDN)
MoMS <- cluster_fast_greedy(MSN)

MoAD <- cluster_edge_betweenness(ADN,directed = F)
MoPD <- cluster_edge_betweenness(PDN,directed = F)
MoMS <- cluster_edge_betweenness(MSN,directed = F)