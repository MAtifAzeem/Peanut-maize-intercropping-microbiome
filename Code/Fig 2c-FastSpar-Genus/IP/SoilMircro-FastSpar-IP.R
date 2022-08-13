#in R

library(stringr)#字符串处理.string manIPulation
library(tidyr)
library(dplyr)
wdImport <- c("E:/Study/SCI/Soil Micro/SCI/Figures/16S/QIIME2-DAD2/FourStages-DAD2-F295+R205-Filtered1/Import")
wdOutput <- c("E:/Study/SCI/Soil Micro/SCI/Figures/16S/QIIME2-DAD2/FourStages-DAD2-F295+R205-Filtered1/Network/FastSpar-Genus/IP")

data.set.name = '_FourStages_Filtered1' 
setwd(wdImport)
ASVCount <- read.table("feature_table_FourStages_DAD2_F295_R205_Filtered1.txt",header=T,row.names = 1,na.strings = c("NA"))
SampleData <- read.table("FourStages_SampleData.txt",header=T,row.names = 1,na.strings = c("NA"))
taxonomy_FourStages_Filtered1 <- read.table("Sup_taxonomy_FourStages_DAD2_F295_R205_Filtered1.txt",header=T,row.names = 1,sep=",", na.strings = "")
dim(ASVCount)
dim(taxonomy_FourStages_Filtered1)
rownames(ASVCount) <- row.names(taxonomy_FourStages_Filtered1)
colnames(ASVCount) <- row.names(SampleData)
ASVCountIP <- ASVCount%>%select(starts_with("IP"))
ASVCountIPTaxonmy <- cbind(ASVCountIP,taxonomy_FourStages_Filtered1)
GenusCountIPFastSpar <- as.data.frame(ASVCountIPTaxonmy %>% group_by(Genus)%>%
                                         summarise_at(vars(IP1_1:IP4_3),funs(sum)))
GenusCountIPFastSpar$Genus<-substring(text=GenusCountIPFastSpar$Genus,first=4,last = 1000000L)
dim(GenusCountIPFastSpar)

rownames(GenusCountIPFastSpar)<-GenusCountIPFastSpar$Genus
GenusCountIPFastSpar <- subset(GenusCountIPFastSpar,select=-c(Genus))
GenusCountIPFastSparWithoutAllZero<- GenusCountIPFastSpar[apply(GenusCountIPFastSpar,1,function(x) mean(x)>0),]
GenusCountIPFastSparWithoutAllZero$Genus <-row.names(GenusCountIPFastSparWithoutAllZero)
GenusCountIPFastSparWithoutAllZero<-GenusCountIPFastSparWithoutAllZero[,c(13,seq(1:12))]
dim(GenusCountIPFastSparWithoutAllZero)
names(GenusCountIPFastSparWithoutAllZero)<-c("OTU ID",rep("IP",12))
setwd(wdOutput)
getwd()
write.table(GenusCountIPFastSparWithoutAllZero, paste("GenusCountIPFastSparWithoutAllZero",data.set.name,".tsv",sep=""), row.names = F, sep = '\t', quote = FALSE)
