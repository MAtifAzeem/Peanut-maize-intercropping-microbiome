# in R

library(stringr)#字符串处理.string manmpulation
library(tidyr)
library(dplyr)
wdImport <- c("E:/Study/SCI/Soil Micro/SCI/Figures/16S/QIIME2-DAD2/FourStages-DAD2-F295+R205-Filtered1/Import")
wdOutput <- c("E:/Study/SCI/Soil Micro/SCI/Figures/16S/QIIME2-DAD2/FourStages-DAD2-F295+R205-Filtered1/Network/FastSpar-Genus/MP")

data.set.name = '_FourStages_Filtered1' 
setwd(wdImport)
ASVCount <- read.table("feature_table_FourStages_DAD2_F295_R205_Filtered1.txt",header=T,row.names = 1,na.strings = c("NA"))
SampleData <- read.table("FourStages_SampleData.txt",header=T,row.names = 1,na.strings = c("NA"))
taxonomy_FourStages_Filtered1 <- read.table("Sup_taxonomy_FourStages_DAD2_F295_R205_Filtered1.txt",header=T,row.names = 1,sep=",", na.strings = "")
dim(ASVCount)
dim(taxonomy_FourStages_Filtered1)
rownames(ASVCount) <- row.names(taxonomy_FourStages_Filtered1)
colnames(ASVCount) <- row.names(SampleData)
ASVCountMP <- ASVCount%>%select(starts_with("MP"))
ASVCountMPTaxonmy <- cbind(ASVCountMP,taxonomy_FourStages_Filtered1)
GenusCountMPFastSpar <- as.data.frame(ASVCountMPTaxonmy %>% group_by(Genus)%>%
                                         summarise_at(vars(MP1_1:MP4_3),funs(sum)))
GenusCountMPFastSpar$Genus<-substring(text=GenusCountMPFastSpar$Genus,first=4,last = 1000000L)
dim(GenusCountMPFastSpar)
#把全为0的行去掉
rownames(GenusCountMPFastSpar)<-GenusCountMPFastSpar$Genus
GenusCountMPFastSpar <- subset(GenusCountMPFastSpar,select=-c(Genus))
GenusCountMPFastSparWithoutAllZero<- GenusCountMPFastSpar[apply(GenusCountMPFastSpar,1,function(x) mean(x)>0),]
GenusCountMPFastSparWithoutAllZero$Genus <-row.names(GenusCountMPFastSparWithoutAllZero)
GenusCountMPFastSparWithoutAllZero<-GenusCountMPFastSparWithoutAllZero[,c(13,seq(1:12))]
dim(GenusCountMPFastSparWithoutAllZero)
names(GenusCountMPFastSparWithoutAllZero)<-c("OTU ID",rep("MP",12))
setwd(wdOutput)
getwd()
write.table(GenusCountMPFastSparWithoutAllZero, paste("GenusCountMPFastSparWithoutAllZero",data.set.name,".tsv",sep=""), row.names = F, sep = '\t', quote = FALSE)

