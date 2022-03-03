#### 1. loading required libraries####
library(ggplot2)#作图 plot
library(dplyr)#数据清洗，Data cleaning
library(plyr)#数据清洗，Data cleaning
library(readxl)#读入 excel, read excel
library(data.table)
library(tidyr)
wdImport<-("E:/Study/SCI/Soil Micro/SCI/Figures/Data/Data for submit")
data.set.name = '_FourStages_Filtered1' 
setwd(wdImport)
ASVCount <- read.table("feature_table_FourStages_DAD2_F295_R205_Filtered1.txt",header=T,row.names = 1,na.strings = c("NA"))
SampleData <- read.table("FourStages_SampleData.txt",header=T,row.names = 1,na.strings = c("NA"))
Tax <- read.table("Sup_taxonomy_FourStages_DAD2_F295_R205_Filtered1.txt",header=T,row.names = 1,na.strings = c("NA"),sep=",")
Tax <- subset(Tax,select=-c(ASV,Species))
colnames(Tax)[1] <- 'Kingdom'
Tax$Kingdom<-substring(text=Tax$Kingdom,first=4,last = 1000000L)#从第4位开始取值
Tax$Phylum<-substring(text=Tax$Phylum,first=4,last = 1000000L)#从第4位开始取值
Tax$Class<-substring(text=Tax$Class,first=4,last = 1000000L)#从第4位开始取值
Tax$Order<-substring(text=Tax$Order,first=4,last = 1000000L)#从第4位开始取值
Tax$Family<-substring(text=Tax$Family,first=4,last = 1000000L)#从第4位开始取值
Tax$Genus<-substring(text=Tax$Genus,first=4,last = 1000000L)#从第4位开始取值
library(amplicon)
####Peanut####
ASVCountPeanut <- ASVCount[,c(31:42,55:66)]
SampleDataPeanut <-SampleData%>%filter(Species=="Peanut"&Rhizocompartments=="Rhizosphere")
wdOutputPeanut <- c("E:/Study/SCI/Soil Micro/SCI/Figures/16S/QIIME2-DAD2/FourStages-DAD2-F295+R205-Filtered1/Lefse/Peanut")
setwd(wdOutputPeanut)
format2lefse2(
  ASVCountPeanut,
  Tax,
  SampleDataPeanut,
  thre = 0.0,
  groupID = "SystemSpecies",
  output = "LEfSePeanutRhizo_0.txt"
)

####Maize####
ASVCountMaize <- ASVCount[,c(19:30,43:54)]
SampleDataMaize <-SampleData%>%filter(Species=="Maize"&Rhizocompartments=="Rhizosphere")
wdOutputMaize <- c("E:/Study/SCI/Soil Micro/SCI/Figures/16S/QIIME2-DAD2/FourStages-DAD2-F295+R205-Filtered1/Lefse/Maize")
setwd(wdOutputMaize)
format2lefse2(
  ASVCountMaize,
  Tax,
  SampleDataMaize,
  thre = 0.0,
  groupID = "SystemSpecies",
  output = "LEfSeMaizeRhizo_0.txt"
)


####MPvsMM####
ASVCount_MPvsMM <- ASVCount[,c(43:66)]
SampleData_MPvsMM <-SampleData[c(43:66),]
wdOutputMaize <- c("E:/Study/SCI/Soil Micro/SCI/Figures/16S/QIIME2-DAD2/FourStages-DAD2-F295+R205-Filtered1/Lefse/MPvsMM")
setwd(wdOutputMaize)
format2lefse2(
  ASVCount_MPvsMM,
  Tax,
  SampleData_MPvsMM,
  thre = 0.0,
  groupID = "SystemSpecies",
  output = "LEfSeRhizo_MPvsMM_0.txt"
)

