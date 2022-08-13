#in R

wdImport <- c("E:/working/SCI/Soil Micro/SCI/Figures/16S/QIIME2-DAD2/FourStages-DAD2-F295+R205-Filtered1/Import")
wdImport2 <- c("E:/working/SCI/Soil Micro/SCI/Figures/16S/QIIME2-DAD2/FourStages-DAD2-F295+R205-Filtered1/Network/FastSpar-Genus/IP")
wdOutput <- c("E:/working/SCI/Soil Micro/SCI/Figures/16S/QIIME2-DAD2/FourStages-DAD2-F295+R205-Filtered1/Network/FastSpar-Genus/IP/1000times/0.77_p0.01")
library(stringr)#字符串处理.string manipulation
library(tidyr)
library(dplyr)

setwd(wdImport2)
data.set.name2 <- "_FastSpar_FourStages_Filtered1_IP_0.77_p0.01"
cor <- read.delim('GenusCountIPFastSparWithoutAllZero_FourStages_Filtered1_median_correlation.tsv', row.names = 1, sep = '\t', check.names = FALSE)
pvals <- read.delim('GenusCountIPFourStages_Filtered1_pvalues.tsv', row.names = 1, sep = '\t', check.names = FALSE)

#filtering
cor[abs(cor) <= 0.77] <- 0

pvals[pvals>=0.01] <- -1
pvals[pvals<0.01 & pvals>=0] <- 1
pvals[pvals==-1] <- 0

#adjacent matrix after filtering
network_adj<- as.matrix(cor) * as.matrix(pvals)
sum(network_adj)
diag(network_adj) <- 0    
network_adj[network_adj!=0]
View(network_adj)
setwd(wdOutput)
write.table(data.frame(network_adj, check.names = FALSE), paste('network_adj',data.set.name2,'.txt',sep = ''), col.names = NA, sep = '\t', quote = FALSE)

library(igraph)
neetwork_adj <- read.delim(paste('network_adj',data.set.name2,'.txt',sep = ''), row.names = 1, sep = '\t', check.names = FALSE)
head(neetwork_adj)[1:6]    

ig <- graph_from_adjacency_matrix(as.matrix(neetwork_adj), mode = 'undirected', weighted = TRUE, diag = FALSE)
ig    #igraph adjacent list

E(ig)$sparcc <- E(ig)$weight
E(ig)$weight <- abs(E(ig)$weight)

adj_matrix <- as.matrix(get.adjacency(ig, attr = 'sparcc'))
setwd(wdOutput)
write.table(data.frame(adj_matrix, check.names = FALSE), paste('network.adj_matrix',data.set.name2,'.csv',sep = ''), col.names = NA, sep = '\t', quote = FALSE)

#graphml for gephi
setwd(wdOutput)
write.graph(ig, paste('network',data.set.name2,'.graphml',sep = ""), format = 'graphml')

#gml for cytoscape
setwd(wdOutput)
write.graph(ig, paste('network',data.set.name2,'.gml',sep = ""), format = 'gml')

#edges list
edge <- data.frame(as_edgelist(ig))

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(ig)$weight,
  sparcc = E(ig)$sparcc
)
head(edge_list)
CorType <- edge_list$sparcc
CorType[CorType>0] <-"positive"
CorType[CorType<0] <-"negative"

unique(CorType)
CorType
edge_list$CorType <- CorType
unique(union(edge_list$source,edge_list$target))
setwd(wdOutput)
write.table(edge_list, paste('network.edge_list',data.set.name2,'.csv',sep = ""), sep = '\t', row.names = FALSE, quote = FALSE)

#edge_list_NetShift
edge_list_NetShift<-edge_list[,-c(4,5)]
write.table(edge_list_NetShift, paste('network.edge_list',data.set.name2,'_NetShift','.csv',sep = ""), sep = '\t', row.names = FALSE, col.names=F,quote = FALSE)

节#Nodes list
node_list <- data.frame(
  Id = V(ig)$name,    #nodes name
  degree = degree(ig)    #nodes degree
)
head(node_list)
node_list <- node_list%>% filter(degree!=0)
row.names(node_list)
#check
all(unique(sort(union(edge_list$source,edge_list$target)))==row.names(node_list))

#Nodes coloring
setwd(wdImport)
ASVCount <- read.table("feature_table_FourStages_DAD2_F295_R205_Filtered1.txt",header=T,row.names = 1,na.strings = c("NA"))
SampleData <- read.table("FourStages_SampleData.txt",header=T,row.names = 1,na.strings = c("NA"))
taxonomy_FourStages_Filtered1 <- read.table("Sup_taxonomy_FourStages_DAD2_F295_R205_Filtered1.txt",sep=",",header=T,row.names = 1,fill=TRUE, na.strings = "")
dim(ASVCount)
rownames(ASVCount) <- row.names(taxonomy_FourStages_Filtered1)
colnames(ASVCount) <- row.names(SampleData)
ASVCountIP <- ASVCount%>%select(starts_with("IP"))

ASVCountIPTaxonmy <- cbind(ASVCountIP,taxonomy_FourStages_Filtered1)
GenusCountIP <- as.data.frame(ASVCountIPTaxonmy %>% group_by(Genus,Domain,Phylum,Class,Order,Family)%>%
                                summarise_at(vars(IP1_1:IP4_3),funs(sum)))
GenusCountIP$Genus<-substring(text=GenusCountIP$Genus,first=4,last = 1000000L)#从第4位开始取值
GenusCountIP$Phylum<-substring(text=GenusCountIP$Phylum,first=4,last = 1000000L)
rownames(GenusCountIP) <- GenusCountIP$Genus
GenusIPNode <- GenusCountIP[node_list$Id,]
GenusIPNodeTax <- GenusIPNode[,c(2:6)]
dim(GenusIPNodeTax)
NodePhylum <- GenusIPNodeTax$Phylum
NodePhylum
#Top 10 phyla
PhylumTop10 <- c("Proteobacteria","Actinobacteriota","Chloroflexi","Acidobacteriota","Bacteroidota",
                 "Patescibacteria","Gemmatimonadota","Firmicutes","Myxococcota","Verrucomicrobiota")

NodePhylumOther <- setdiff(unique(NodePhylum),PhylumTop10)
NodePhylumOther
# Naming others
NodePhylum[pmatch(NodePhylumOther,NodePhylum,duplicates.ok=T)]<-"Others"
NodePhylum[pmatch(NodePhylumOther,NodePhylum,duplicates.ok=T)]<-"Others"
NodePhylum[pmatch(NodePhylumOther,NodePhylum,duplicates.ok=T)]<-"Others"
NodePhylum[pmatch(NodePhylumOther,NodePhylum,duplicates.ok=T)]<-"Others"
NodePhylum[pmatch(NodePhylumOther,NodePhylum,duplicates.ok=T)]<-"Others"
NodePhylum[pmatch(NodePhylumOther,NodePhylum,duplicates.ok=T)]<-"Others"
NodePhylum[pmatch(NodePhylumOther,NodePhylum,duplicates.ok=T)]<-"Others"
NodePhylum[pmatch(NodePhylumOther,NodePhylum,duplicates.ok=T)]<-"Others"
NodePhylum[pmatch(NodePhylumOther,NodePhylum,duplicates.ok=T)]<-"Others"
NodePhylum[pmatch(NodePhylumOther,NodePhylum,duplicates.ok=T)]<-"Others"
NodePhylum[pmatch(NodePhylumOther,NodePhylum,duplicates.ok=T)]<-"Others"
NodePhylum[pmatch(NodePhylumOther,NodePhylum,duplicates.ok=T)]<-"Others"
setdiff(unique(NodePhylum),PhylumTop10)
unique(NodePhylum)
NodePhylum
#
PhylumTop10ColorValue <- c("#e31a1c","#1f78b4","#b2df8a","#33a02c","#fb9a99","#a6cee3","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a")
NodeColor <- NodePhylum
i=1
for(i in 1:10){
  NodeColor <-gsub(pattern = PhylumTop10[i],replacement = PhylumTop10ColorValue[i],x=NodeColor)
}
NodeColor <- gsub(pattern = c("Others"),replacement = c("#c0c0c0"),x=NodeColor)
unique(NodeColor)
GenusIPNodeTax$color <- NodeColor
node_listTax <- cbind(node_list,GenusIPNodeTax)
setwd(wdOutput)
write.table(node_listTax, paste('network.node_list',data.set.name2,'.csv',sep = ""), sep = '\t', row.names = FALSE, quote = FALSE)

