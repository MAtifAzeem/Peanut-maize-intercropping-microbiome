#### 1. loading required libraries ####
library(ggplot2)#作图 plot
library(ggpubr)#添加显著性标记, Add the significance marker
library(ggsignif)#添加显著性标记, Add the significance marker
library(dplyr)#数据清洗，Data cleaning
library(plyr)#数据清洗，Data cleaning
library(ggthemes)#ggplot所用主题，Themes for ggplot2
library(readxl)#读入 excel, read excel
library(ggsci)#配色，color scheme
library(showtext)#字体设置, font setting
library(extrafont)#使用系统字体，Using the system fonts
library(sysfonts)#加载系统字体，loading the system fonts
library(Cairo)#抗锯齿,anti-aliasing
library(ape)
library(phyloseq)#microbiome analysis
library(vegan)#Adonis analysis
#### 2. Setting themes and working dictionary path ####
loadfonts()
Sys.setenv(R_GSCMD = "C:/Program Files (x86)/gs/gs9.50/bin/gswin32c.exe")

mytheme1 <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#4d4d4d",size=0.2),
                              panel.background = element_rect(colour = "#4d4d4d",size=0.2),
                              text = element_text(family = "Arial"),
                              strip.text = element_text(size = 6,hjust = 0.5),
                              plot.title = element_text(size = 6,hjust = 0.5),
                              axis.text=element_text(size=6,color = "#4D4D4D"),
                              axis.title=element_text(size = 6),
                              legend.text = element_text(size = 6),
                              legend.title = element_text(size = 6),
                              legend.background = element_blank(),
                              panel.border = element_rect(colour = NA),
                              axis.line = element_line(color = "#4D4D4D",size=0.2),
                              axis.ticks.length = unit(0.8, "mm"))#移除整体的边???

mytheme_bigfonts <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#4d4d4d",size=0.2),
                              panel.background = element_rect(colour = "#4d4d4d",size=0.2),
                              text = element_text(family = "Arial"),
                              strip.text = element_text(size = 9,hjust = 0.5),
                              plot.title = element_text(size = 9,hjust = 0.5),
                              axis.text=element_text(size=9,color = "#4D4D4D"),
                              axis.title=element_text(size = 9),
                              legend.text = element_text(size = 9),
                              legend.title = element_text(size = 9),
                              legend.background = element_blank(),
                              panel.border = element_rect(colour = NA),
                              axis.line = element_line(color = "#4D4D4D",size=0.2),
                              axis.ticks.length = unit(0.8, "mm"))#移除整体的边???

wdImport<-("E:/Study/SCI/Soil Micro/SCI/Figures/Data/Data for submit")
wdOutput<- ("E:/Study/SCI/Soil Micro/SCI/Figures/Figures from R/Figure2")
data.set.name = '_FourStages_DAD2_F295_R205_Filtered1' 


#### 3 Import and process data ####
setwd(wdImport)
ASVCount <- read.table("feature_table_FourStages_DAD2_F295_R205_Filtered1.txt",header=T,row.names = 1,na.strings = c("NA"))
SampleData <- read.table("FourStages_SampleData.txt",header=T,row.names = 1,na.strings = c("NA"))
colnames(ASVCount)<- row.names(SampleData)
RootedTree_FourStages_Filtered1<- read.tree("rooted_tree_FourStages_DAD2_F295_R205_Filtered1.nwk")
Df <- phyloseq(otu_table(ASVCount, taxa_are_rows = T),sample_data(SampleData),phy_tree(RootedTree_FourStages_Filtered1))
Df
Dfr <-transform_sample_counts(Df, function(x) x / sum(x) )
unifrac_PCoA <- ordinate(Dfr, "PCoA", "unifrac")
unifrac_PCoA_Vector<-as.data.frame(unifrac_PCoA$vectors)
setwd(wdOutput)
write.table(unifrac_PCoA_Vector, paste("unifrac_vector",data.set.name,".txt",sep=""), col.names = NA, sep = '\t', quote = FALSE)
sample_data(Dfr)$RhizocompartmentsSystemSpecies <- factor(sample_data(Dfr)$RhizocompartmentsSystemSpecies,level=c("Bulk","MP","IP","IM","MM"))
sample_data(Dfr)$Species <- factor(sample_data(Dfr)$Species,level=c("Bulk","Peanut","Maize"))
sample_data(Dfr)$Days<-factor(sample_data(Dfr)$Days)
Dfr_rhizo<-subset_samples(Dfr,Rhizocompartments=="Rhizosphere")

#### 4. Days ####
#### 4.1 46d ####
Dfr_46 <- subset_samples(Dfr_rhizo,Dfr_rhizo@sam_data$Days=="46")
Dfr_46
#adonis
unifrac_dis_46 <- distance(Dfr_46, method = 'unifrac')
SampleData_rhizo_46 <-SampleData[c(19:21,31:33,43:45,55:57),]
adonis2(formula=unifrac_dis_46 ~ SystemSpecies, 
        data=SampleData_rhizo_46,
        permutations=9999)
#
ASVCount_rhizo_46<-ASVCount[,c(19:21,31:33,43:45,55:57)]
tree=RootedTree_FourStages_Filtered1
amplicon::MicroTest(otu=ASVCount_rhizo_46,map=SampleData_rhizo_46,group="RhizocompartmentsSystemSpecies",
                    Micromet="adonis",dist="unifrac")

#plots#
Dfr_rhizo <- subset_samples(Dfr,Dfr@sam_data$Rhizocompartments=="Rhizosphere")
Dfr_rhizo_46 <- subset_samples(Dfr_rhizo,Dfr_rhizo@sam_data$Days=="46")
Dfr_rhizo_46
sample_data(Dfr_rhizo_46)$SystemSpecies <- factor(sample_data(Dfr_rhizo_46)$SystemSpecies,level=c("MP","IP","IM","MM"))

unifrac_PCoA_rhizo_46 <- ordinate(Dfr_rhizo_46, method="PCoA", distance="unifrac")
unifrac_PCoAPoints_rhizo_46 <- plot_ordination(Dfr_rhizo_46, unifrac_PCoA_rhizo_46,type = "samples", color = "SystemSpecies" )+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point()+
  scale_color_manual(values = c("#4D78B2","#BE1A17","#E19896","#B3CDE2"))+
  scale_x_continuous(limits=c(-0.3,0.3),breaks=c(-0.2,0.0,0.2))+
  scale_y_continuous(limits=c(-0.3,0.45),breaks=c(-0.2,0.0,0.2))+
  geom_polygon(aes(colour = SystemSpecies),fill=NA)+
  mytheme1
unifrac_PCoAPoints_rhizo_46
setwd(wdOutput)
getwd()
ggsave(paste("unifrac_PCoAPoints_rhizo_46",data.set.name,".pdf",sep=""),
       unifrac_PCoAPoints_rhizo_46,
       device=cairo_pdf,width=80,height=60,dpi = 600,units = "mm")

#### 4.2 53d ####
Dfr_53 <- subset_samples(Dfr_rhizo,Dfr_rhizo@sam_data$Days=="53")
Dfr_53
#adonis
unifrac_dis_53 <- distance(Dfr_53, method = 'unifrac')
SampleData_rhizo_53 <-SampleData[c(22:24,34:36,46:48,58:60),]
adonis2(formula=unifrac_dis_53 ~ SystemSpecies, 
        data=SampleData_rhizo_53,
        permutations=9999)
#
ASVCount_rhizo_53<-ASVCount[,c(22:24,34:36,46:48,58:60)]
tree=RootedTree_FourStages_Filtered1
amplicon::MicroTest(otu=ASVCount_rhizo_53,map=SampleData_rhizo_53,group="RhizocompartmentsSystemSpecies",
                    Micromet="adonis",dist="unifrac")
ASVCount_rhizo_73_Peanut<-ASVCount[,c(34:36,58:60)]
tree=RootedTree_FourStages_Filtered1



#plots#
Dfr_rhizo <- subset_samples(Dfr,Dfr@sam_data$Rhizocompartments=="Rhizosphere")
Dfr_rhizo_53 <- subset_samples(Dfr_rhizo,Dfr_rhizo@sam_data$Days=="53")
Dfr_rhizo_53
sample_data(Dfr_rhizo_53)$SystemSpecies <- factor(sample_data(Dfr_rhizo_53)$SystemSpecies,level=c("MP","IP","IM","MM"))

unifrac_PCoA_rhizo_53 <- ordinate(Dfr_rhizo_53, method="PCoA", distance="unifrac")
unifrac_PCoAPoints_rhizo_53 <- plot_ordination(Dfr_rhizo_53, unifrac_PCoA_rhizo_53,type = "samples", color = "SystemSpecies" )+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point(size=1)+
  scale_color_manual(values = c("#4D78B2","#BE1A17","#E19896","#B3CDE2"))+
  scale_x_continuous(limits=c(-0.3,0.3),breaks=c(-0.2,0.0,0.2))+
  scale_y_continuous(limits=c(-0.3,0.45),breaks=c(-0.2,0.0,0.2))+
  geom_polygon(aes(colour = SystemSpecies),fill=NA)+
  mytheme1
unifrac_PCoAPoints_rhizo_53
setwd(wdOutput)
getwd()
ggsave(paste("unifrac_PCoAPoints_rhizo_53",data.set.name,".pdf",sep=""),
       unifrac_PCoAPoints_rhizo_53,
       device=cairo_pdf,width=80,height=60,dpi = 600,units = "mm")

#### 4.3 63d ####
#adonis
Dfr_63 <- subset_samples(Dfr_rhizo,Dfr_rhizo@sam_data$Days=="63")
Dfr_63
unifrac_dis_63 <- distance(Dfr_63, method = 'unifrac')
SampleData_rhizo_63 <-SampleData[c(25:27,37:39,49:51,61:63),]
adonis2(formula=unifrac_dis_63 ~ SystemSpecies, 
        data=SampleData_rhizo_73,
        permutations=9999)
#
ASVCount_rhizo_63<-ASVCount[,c(25:27,37:39,49:51,61:63)]
tree=RootedTree_FourStages_Filtered1
amplicon::MicroTest(otu=ASVCount_rhizo_63,map=SampleData_rhizo_63,group="RhizocompartmentsSystemSpecies",
                    Micromet="adonis",dist="unifrac")

#plots#
Dfr_rhizo <- subset_samples(Dfr,Dfr@sam_data$Rhizocompartments=="Rhizosphere")
Dfr_rhizo_63 <- subset_samples(Dfr_rhizo,Dfr_rhizo@sam_data$Days=="63")
Dfr_rhizo_63
sample_data(Dfr_rhizo_63)$SystemSpecies <- factor(sample_data(Dfr_rhizo_63)$SystemSpecies,level=c("MP","IP","IM","MM"))

unifrac_PCoA_rhizo_63 <- ordinate(Dfr_rhizo_63, method="PCoA", distance="unifrac")
unifrac_PCoAPoints_rhizo_63 <- plot_ordination(Dfr_rhizo_63, unifrac_PCoA_rhizo_63,type = "samples", color = "SystemSpecies" )+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point(size=1)+
  scale_color_manual(values = c("#4D78B2","#BE1A17","#E19896","#B3CDE2"))+
  scale_x_continuous(limits=c(-0.3,0.3),breaks=c(-0.2,0.0,0.2))+
  scale_y_continuous(limits=c(-0.3,0.45),breaks=c(-0.2,0.0,0.2))+
  geom_polygon(aes(colour = SystemSpecies),fill=NA)+
  mytheme1
unifrac_PCoAPoints_rhizo_63
setwd(wdOutput)
getwd()
ggsave(paste("unifrac_PCoAPoints_rhizo_63",data.set.name,".pdf",sep=""),
       unifrac_PCoAPoints_rhizo_63,
       device=cairo_pdf,width=80,height=60,dpi = 600,units = "mm")

#### 4.4 73d ####
#adonis
Dfr_73 <- subset_samples(Dfr_rhizo,Dfr_rhizo@sam_data$Days=="73")
Dfr_73
unifrac_dis_73 <- distance(Dfr_73, method = 'unifrac')
SampleData_rhizo_73 <-SampleData[c(28:30,40:42,52:54,64:66),]
adonis2(formula=unifrac_dis_73 ~ SystemSpecies, 
        data=SampleData_rhizo_73,
        permutations=9999)
#
ASVCount_rhizo_73<-ASVCount[,c(28:30,40:42,52:54,64:66)]
tree=RootedTree_FourStages_Filtered1
amplicon::MicroTest(otu=ASVCount_rhizo_73,map=SampleData_rhizo_73,group="RhizocompartmentsSystemSpecies",
                    Micromet="adonis",dist="unifrac")

#plots#
Dfr_rhizo_73 <- subset_samples(Dfr_rhizo,Dfr_rhizo@sam_data$Days=="73")
Dfr_rhizo_73
sample_data(Dfr_rhizo_73)$SystemSpecies <- factor(sample_data(Dfr_rhizo_73)$SystemSpecies,level=c("MP","IP","IM","MM"))

unifrac_PCoA_rhizo_73 <- ordinate(Dfr_rhizo_73, method="PCoA", distance="unifrac")
unifrac_PCoAPoints_rhizo_73 <- plot_ordination(Dfr_rhizo_73, unifrac_PCoA_rhizo_73,type = "samples", color = "SystemSpecies" )+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point(size=1)+
  scale_color_manual(values = c("#4D78B2","#BE1A17","#E19896","#B3CDE2"))+
  scale_x_continuous(limits=c(-0.3,0.3),breaks=c(-0.2,0.0,0.2))+
  scale_y_continuous(limits=c(-0.3,0.45),breaks=c(-0.2,0.0,0.2))+
  geom_polygon(aes(colour = SystemSpecies),fill=NA)+
  mytheme1
unifrac_PCoAPoints_rhizo_73
setwd(wdOutput)
getwd()
ggsave(paste("unifrac_PCoAPoints_rhizo_73",data.set.name,".pdf",sep=""),
       unifrac_PCoAPoints_rhizo_73,
       device=cairo_pdf,width=80,height=60,dpi = 600,units = "mm")

#### 7. Rhizo ####
#adonis
Dfr_rhizo <- subset_samples(Dfr,Dfr@sam_data$Rhizocompartments=="Rhizosphere")
Dfr_rhizo@sam_data$Species
unifrac_dis_rhizo <- distance(Dfr_rhizo, method = 'unifrac')
SampleData_rhizo <-SampleData[-c(1:18),]
adonis2(formula=unifrac_dis_rhizo ~ RhizocompartmentsSystemSpecies+Species, 
        data=SampleData_rhizo,
        permutations=9999)
ASVCount_rhizo<-ASVCount[,-c(1:18)]
tree=RootedTree_FourStages_Filtered1
amplicon::MicroTest(otu=ASVCount_rhizo,map=SampleData_rhizo,group="RhizocompartmentsSystemSpecies",
          Micromet="anosim",dist="unifrac")
#plots#
Dfr_rhizo <- subset_samples(Dfr,Rhizocompartments=="Rhizosphere")
sample_data(Dfr_rhizo)$Days<-factor(sample_data(Dfr_rhizo)$Days)
unifrac_PCoA <- ordinate(Dfr_rhizo, method="NMDS", distance="unifrac")
sample_data(Dfr_rhizo)$Days<-factor(sample_data(Dfr_rhizo)$Days)
unifrac_PCoAPoints_rhizo <- plot_ordination(Dfr_rhizo, unifrac_PCoA,shape = "Days", color="RhizocompartmentsSystemSpecies" )+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point(size=1)+
  scale_x_continuous(limits=c(-0.3,0.22),breaks=c(-0.2,0.0,0.2))+
  scale_y_continuous(limits=c(-0.22,0.4),breaks=c(-0.2,0.0,0.2,0.3))+
  scale_color_manual(values = c("#4D78B2","#BE1A17","#E19896","#B3CDE2"))+
  scale_shape_few()+
  mytheme1
unifrac_PCoAPoints_rhizo
setwd(wdOutput)
getwd()
ggsave(paste("unifrac_PCoAPoints_rhizo",data.set.name,".pdf",sep=""),unifrac_PCoAPoints_rhizo,
       device=cairo_pdf,width=80,height=70,dpi = 600,units = "mm")

#### 4.unifrac_similarity ####
### 46 days ###
## unifrac_distance_46days ##
setwd(wdImport)
RootedTree_FourStages_Filtered1 <- read.tree("rooted_tree_FourStages_DAD2_F295_R205_Filtered1.nwk")
Df <- phyloseq(otu_table(ASVCount, taxa_are_rows = T),sample_data(SampleData),phy_tree(RootedTree_FourStages_Filtered1))
Df
Dfr <-transform_sample_counts(Df, function(x) x / sum(x) )
Dfr_rhizo<-subset_samples(Dfr,Rhizocompartments=="Rhizosphere")
Unifrac_dis <- distance(Dfr_rhizo, method = 'unifrac') 
Unifrac_disA <- phyloseq::UniFrac(Dfr_rhizo,weighted=F)#和上面的结果一样的
all(Unifrac_dis==Unifrac_disA)
### 46days ###
## unifrac_similarity_46days ##
Dfr_rhizo_46<-subset_samples(Dfr_rhizo,Days==46)
Dfr_rhizo_46
unifrac_dis_rhizo_46 <- phyloseq::UniFrac(Dfr_rhizo_46,weighted=F)
unifrac_sim_rhizo_46 <- as.matrix(1 - unifrac_dis_rhizo_46)
## IPvsIM46 ##
IPvsIM46_unifrac_sim<-as.numeric(unifrac_sim_rhizo_46[c("IP1_1","IP1_2","IP1_3"),c("IM1_1","IM1_2","IM1_3")])
## MPvsMM46 ##
MPvsMM46_unifrac_sim<-as.numeric(unifrac_sim_rhizo_46[c("MP1_1","MP1_2","MP1_3"),c("MM1_1","MM1_2","MM1_3")])

### 53days ###
## unifrac_similarity_53days ##
Dfr_rhizo_53<-subset_samples(Dfr_rhizo,Days==53)
Dfr_rhizo_53
unifrac_dis_rhizo_53 <- phyloseq::UniFrac(Dfr_rhizo_53,weighted=F)
unifrac_sim_rhizo_53 <- as.matrix(1 - unifrac_dis_rhizo_53)
## IPvsIM53 ##
IPvsIM53_unifrac_sim<-as.numeric(unifrac_sim_rhizo_53[c("IP2_1","IP2_3","IP2_4"),c("IM2_1","IM2_2","IM2_3")])
## MPvsMM53 ##
MPvsMM53_unifrac_sim<-as.numeric(unifrac_sim_rhizo_53[c("MP2_1","MP2_2","MP2_3"),c("MM2_2","MM2_3","MM2_4")])

### 63days ###
## unifrac_similarity_63days ##
Dfr_rhizo_63<-subset_samples(Dfr_rhizo,Days==63)
Dfr_rhizo_63
unifrac_dis_rhizo_63 <- phyloseq::UniFrac(Dfr_rhizo_63,weighted=F)
unifrac_sim_rhizo_63 <- as.matrix(1 - unifrac_dis_rhizo_63)
## IPvsIM63 ##
IPvsIM63_unifrac_sim<-as.numeric(unifrac_sim_rhizo_63[c("IP3_1","IP3_3","IP3_4"),c("IM3_1","IM3_2","IM3_3")])
## MPvsMM63 ##
MPvsMM63_unifrac_sim<-as.numeric(unifrac_sim_rhizo_63[c("MP3_1","MP3_2","MP3_4"),c("MM3_2","MM3_3","MM3_4")])

### 73days ###
## unifrac_similarity_73days ##
Dfr_rhizo_73<-subset_samples(Dfr_rhizo,Days==73)
Dfr_rhizo_73
unifrac_dis_rhizo_73 <- phyloseq::UniFrac(Dfr_rhizo_73,weighted=F)
unifrac_sim_rhizo_73 <- as.matrix(1 - unifrac_dis_rhizo_73)
## IPvsIM73 ##
IPvsIM73_unifrac_sim<-as.numeric(unifrac_sim_rhizo_73[c("IP4_1","IP4_2","IP4_3"),c("IM4_1","IM4_2","IM4_3")])
## MPvsMM73 ##
MPvsMM73_unifrac_sim<-as.numeric(unifrac_sim_rhizo_73[c("MP4_1","MP4_2","MP4_3"),c("MM4_1","MM4_2","MM4_3")])

Group<-factor(c(rep("IP vs IM46",9),rep("IP vs IM53",9),rep("IP vs IM63",9),rep("IP vs IM73",9)))
VSIn<-factor(c(rep("IP vs IM",36)))
Days<-as.numeric(c(rep(46,9),rep(53,9),rep(63,9),rep(73,9)))
Class<-c(rep("All taxa",36))
unifrac_sim_IPvsIM<-c(IPvsIM46_unifrac_sim,IPvsIM53_unifrac_sim,IPvsIM63_unifrac_sim,IPvsIM73_unifrac_sim)
IPvsIMAllTaxaunifrac_sim <-data.frame(Group,VS=VSIn,Days,Class=Class,Similarity=unifrac_sim_IPvsIM)
sapply(IPvsIMAllTaxaunifrac_sim,class)
library(tidyverse)
library(broom)#将统计模型结果整理成数据框形式
##线性回归模型的统计量
#All taxa
lm(Similarity~Days,IPvsIMAllTaxaunifrac_sim)
IPvsIMunifracAlllmfit <- lm(Similarity~Days,IPvsIMAllTaxaunifrac_sim)
summary(IPvsIMunifracAlllmfit)
tidy(IPvsIMunifracAlllmfit)
augment(IPvsIMunifracAlllmfit)
glance(IPvsIMunifracAlllmfit)
cor(IPvsIMAllTaxaunifrac_sim$Days,IPvsIMAllTaxaunifrac_sim$Similarity)
#### IM vsMM ####
##All taxa
Group<-factor(c(rep("MP vs MM46",9),rep("IM vsMM53",9),rep("IM vsMM63",9),rep("IM vsMM73",9)))
VSMo<-factor(c(rep("MP vs MM",36)))
Days<-as.numeric(c(rep(46,9),rep(53,9),rep(63,9),rep(73,9)))
Class<-c(rep("All taxa",36))
unifrac_sim_MPvsMM<-c(MPvsMM46_unifrac_sim,MPvsMM53_unifrac_sim,MPvsMM63_unifrac_sim,MPvsMM73_unifrac_sim)
MPvsMMAllTaxaunifrac_sim<-data.frame(Group,VS=VSMo,Days,Class=Class,Similarity=unifrac_sim_MPvsMM)
sapply(MPvsMMAllTaxaunifrac_sim,class)
library(tidyverse)
library(broom)#将统计模型结果整理成数据框形式
##线性回归模型的统计量
#All taxa
lm(Similarity~Days,MPvsMMAllTaxaunifrac_sim)
MPvsMMunifracAlllmfit <- lm(Similarity~Days,MPvsMMAllTaxaunifrac_sim)
summary(MPvsMMAllTaxaunifrac_sim)
tidy(MPvsMMunifracAlllmfit)
augment(MPvsMMunifracAlllmfit)
glance(MPvsMMunifracAlllmfit)
cor(MPvsMMAllTaxaunifrac_sim$Days,MPvsMMAllTaxaunifrac_sim$Similarity)
#### Plots ####
IPvsIM_MPvsMM_unifrac_sim<-rbind(IPvsIMAllTaxaunifrac_sim,MPvsMMAllTaxaunifrac_sim)
IPvsIM_MPvsMM_unifrac_sim$VS <- factor(IPvsIM_MPvsMM_unifrac_sim$VS,levels = c("IP vs IM","MP vs MM"))

IPvsIM_MPvsMM_unifrac_simPointLM <- ggplot(IPvsIM_MPvsMM_unifrac_sim,aes(factor(Days),Similarity,group=VS))+
  geom_point(aes(color=VS),size=0.1,position=position_dodge(0.15),stat="identity")+
  stat_smooth(aes(color=VS,fill=VS),method = loess,
              size=0.5,position=position_dodge(0.15),alpha=0.1,span=2)+
  mytheme1+
  theme(legend.title = element_blank(),legend.key.width=unit(3,'mm'))+
  scale_y_continuous(limits = c(0.4,0.75),breaks=c(0.4,0.5,0.6,0.7))+
  scale_color_npg()+scale_fill_npg()+
  xlab("Days")+ylab("unifrac-Curits Similarity")+
  guides(color=FALSE)
IPvsIM_MPvsMM_unifrac_simPointLM
setwd(wdOutput)
getwd()
ggsave(paste("IPvsIM_MPvsMM_unifrac_simPointLM_FourStages_Filtered1",data.set.name,".pdf",sep=""),
       IPvsIM_MPvsMM_unifrac_simPointLM,device=cairo_pdf,width=80,height=40,dpi = 600,units = "mm")

###statistcal analysis
library(rcompanion)
IPvsIM_MPvsMM_unifrac_sim$VS <- factor(IPvsIM_MPvsMM_unifrac_sim$VS)
IPvsIM_MPvsMM_unifrac_sim$Days <- factor(IPvsIM_MPvsMM_unifrac_sim$Days)
str(IPvsIM_MPvsMM_unifrac_sim)
head(IPvsIM_MPvsMM_unifrac_sim)
scheirerRayHare(Similarity~VS*Days+Days:VS, data = IPvsIM_MPvsMM_unifrac_sim)

library(HH)
interaction2wt(Similarity~Days*VS, data = IPvsIM_MPvsMM_unifrac_sim)##saved as 5*12 IPvsIM_MPvsMM_unifracDistance2wt

# 46 days #
IPvsIM_MPvsMM_unifrac_sim_46<-filter(IPvsIM_MPvsMM_unifrac_sim,Days=="46")
compare_means(data=IPvsIM_MPvsMM_unifrac_sim_46,Similarity~VS,method="wilcox.test")

# 53 days #
IPvsIM_MPvsMM_unifrac_sim_53<-filter(IPvsIM_MPvsMM_unifrac_sim,Days=="53")
compare_means(data=IPvsIM_MPvsMM_unifrac_sim_53,Similarity~VS,method="wilcox.test")

# 63 days #
IPvsIM_MPvsMM_unifrac_sim_63<-filter(IPvsIM_MPvsMM_unifrac_sim,Days=="63")
compare_means(data=IPvsIM_MPvsMM_unifrac_sim_63,Similarity~VS,method="wilcox.test")


# 73 days #
IPvsIM_MPvsMM_unifrac_sim_73<-filter(IPvsIM_MPvsMM_unifrac_sim,Days=="73")
compare_means(data=IPvsIM_MPvsMM_unifrac_sim_73,Similarity~VS,method="wilcox.test")


