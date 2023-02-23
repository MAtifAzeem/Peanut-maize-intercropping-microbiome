#### 1. loading required libraries ####
library(ggplot2)#作图 plot
library(ggpubr)#添加显著性标记, Add the significance marker
library(ggsignif)#添加显著性标记, Add the significance marker
library(dplyr)#数据清洗，Data cleaning
library(plyr)#数据清洗，Data cleaning
library(reshape2)#数据清洗，Data cleaning
library(ggthemes)#ggplot所用主题，Themes for ggplot2
library(ggsci)#TOP期刊配色方案, color plates from top journal
library(readxl)#读入 excel, read excel
library(showtext)#字体设置, font setting
library(extrafont)#使用系统字体，Using the system fonts
library(sysfonts)#加载系统字体，loading the system fonts
library(Cairo)#抗锯齿,anti-aliasing
library(agricolae)#多重比较，Multiple comparisons.
library(rcompanion)
library(car)

#### 2. Setting themes and working dictionary path ####
loadfonts()
Sys.setenv(R_GSCMD = "C:/Program Files (x86)/gs/gs9.50/bin/gswin32c.exe")

mytheme <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#4d4d4d"),
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
                             axis.ticks.length = unit(0.8, "mm"))

mytheme1 <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#4d4d4d",linewidth=0.2),
                              panel.background = element_rect(colour = "#4d4d4d",linewidth=0.2),
                              text = element_text(family = "Arial"),
                              strip.text = element_text(size = 6,hjust = 0.5),
                              plot.title = element_text(size = 6,hjust = 0.5),
                              axis.text=element_text(size=6,color = "#4D4D4D"),
                              axis.title=element_text(size = 6),
                              legend.text = element_text(size = 6),
                              legend.title = element_text(size = 6),
                              legend.background = element_blank(),
                              axis.line = element_line(color = "#4D4D4D",size=0.2),
                              axis.ticks.length = unit(0.8, "mm"))

wdImport<-("E:/working/SCI/Soil Micro/SCI/Figures/Data/Data for submit")
wdOutput <- ("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Supplemental materials/αdiversity_similarity")

#### 4. Fig. s3a αdiversity ####
### 4.1 Import and process data ###
setwd(wdImport)
alpha_diversity <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                    sheet = "fig s3 α diversity")

alpha_diversity$Species <- factor(alpha_diversity$Species,levels = c("Peanut","Maize"))
alpha_diversity$System <- factor(alpha_diversity$System,levels = c("Intercropping","Monocropping"))
alpha_diversity$SystemSpecies <- factor(alpha_diversity$SystemSpecies,levels = c("IP","MP","IM","MM"))

### 4.1 chao1_Peanut ###
## 4.1.1 Statistic analysis ##
alpha_diversity_Peanut <- alpha_diversity%>%filter(Species=="Peanut")
scheirerRayHare(data=alpha_diversity_Peanut, chao1~System*Days)

## 4.1.2 Plots ##
chao1_Peanut_line<-ggplot(alpha_diversity_Peanut, aes(x=Days, y=chao1, fill=SystemSpecies,group=SystemSpecies)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line",size=0.5,aes(col=SystemSpecies)) +
  stat_summary(fun="mean", geom="point",size=0.5,aes(col=SystemSpecies)) +
  mytheme1+
  scale_x_continuous(breaks=c(46,53,63,73),limits=c(38,80))+
  scale_y_continuous(limits=c(0,3500))+
  labs(x="dps",
       y="chao1",parse =T)+
  scale_color_npg()+scale_fill_npg()+
  theme(legend.position = "top")+
  guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))
chao1_Peanut_line
setwd(wdOutput)
getwd()
ggsave("chao1_Peanut_line.pdf",
       chao1_Peanut_line,
       device=cairo_pdf,width=40,height=50,dpi = 300,units = "mm")

### 4.2 chao1_Peanut ###
## 4.2.1 Statistic analysis ##
alpha_diversity_Peanut <- alpha_diversity%>%filter(Species=="Peanut")
scheirerRayHare(data=alpha_diversity_Peanut, shannon~System*Days)

## 4.2.2 Plots ##
shannon_Peanut_line<-ggplot(alpha_diversity_Peanut, aes(x=Days, y=shannon, fill=SystemSpecies,group=SystemSpecies)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line",size=0.5,aes(col=SystemSpecies)) +
  stat_summary(fun="mean", geom="point",size=0.5,aes(col=SystemSpecies)) +
  mytheme1+
  scale_x_continuous(breaks=c(46,53,63,73),limits=c(38,80))+
  scale_y_continuous(limits=c(0,11))+
  labs(x="dps",
       y="chao1",parse =T)+
  scale_color_npg()+scale_fill_npg()+
  theme(legend.position = "top")+
  guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))
shannon_Peanut_line
setwd(wdOutput)
getwd()
ggsave("shannon_Peanut_line.pdf",
       shannon_Peanut_line,
       device=cairo_pdf,width=40,height=53,dpi = 300,units = "mm")

### 4.3 chao1_Maize ###
## 4.3.1 Statistic analysis ##
alpha_diversity_Maize <- alpha_diversity%>%filter(Species=="Maize")
scheirerRayHare(data=alpha_diversity_Maize, chao1~System*Days)

## 4.3.2 Plots ##
chao1_Maize_line<-ggplot(alpha_diversity_Maize, aes(x=Days, y=chao1, fill=SystemSpecies,group=SystemSpecies)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line",size=0.5,aes(col=SystemSpecies)) +
  stat_summary(fun="mean", geom="point",size=0.5,aes(col=SystemSpecies)) +
  mytheme1+
  scale_x_continuous(breaks=c(46,53,63,73),limits=c(38,80))+
  scale_y_continuous(limits=c(0,3500))+
  labs(x="dps",
       y="chao1",parse =T)+
  scale_color_npg()+scale_fill_npg()+
  theme(legend.position = "top")+
  guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))
chao1_Maize_line
setwd(wdOutput)
getwd()
ggsave("chao1_Maize_line.pdf",
       chao1_Maize_line,
       device=cairo_pdf,width=40,height=50,dpi = 300,units = "mm")

### 4.4 chao1_Maize ###
## 4.4.1 Statistic analysis ##
alpha_diversity_Maize <- alpha_diversity%>%filter(Species=="Maize")
scheirerRayHare(data=alpha_diversity_Maize, shannon~System*Days)

## 4.4.2 Plots ##
shannon_Maize_line<-ggplot(alpha_diversity_Maize, aes(x=Days, y=shannon, fill=SystemSpecies,group=SystemSpecies)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line",size=0.5,aes(col=SystemSpecies)) +
  stat_summary(fun="mean", geom="point",size=0.5,aes(col=SystemSpecies)) +
  mytheme1+
  scale_x_continuous(breaks=c(46,53,63,73),limits=c(38,80))+
  scale_y_continuous(limits=c(0,11))+
  labs(x="dps",
       y="chao1",parse =T)+
  scale_color_npg()+scale_fill_npg()+
  theme(legend.position = "top")+
  guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))
shannon_Maize_line
setwd(wdOutput)
getwd()
ggsave("shannon_Maize_line.pdf",
       shannon_Maize_line,
       device=cairo_pdf,width=40,height=53,dpi = 300,units = "mm")

data.set.name = '_FourStages_rhizosphere_DAD2_F295_R205_Filtered1' 

#### 5. Fig. s3b similarity ####
data.set.name = '_FourStages_rhizosphere_DAD2_F295_R205_Filtered1'
#### 5.1.unifrac_similarity ####
### 46 days ###
## unifrac_distance_46days ##
setwd(wdImport)
RootedTree_FourStages_rhizosphere_Filtered1 <- read.tree("rooted_tree_FourStages_rhizosphere_DAD2_F295_R205_Filtered1.nwk")
Df <- phyloseq(otu_table(ASVCount, taxa_are_rows = T),sample_data(SampleData),phy_tree(RootedTree_FourStages_rhizosphere_Filtered1))
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
  stat_smooth(aes(color=VS,fill=VS),method =  stats::loess,
              size=0.5,position=position_dodge(0.15),alpha=0.1,span=2)+
  mytheme1+
  theme(legend.title = element_blank(),legend.key.width=unit(3,'mm'))+
  scale_y_continuous(limits = c(0,0.62))+
  scale_color_npg()+scale_fill_npg()+
  xlab("Days")+ylab("unifrac-Curits Similarity")+
  guides(color="none")
IPvsIM_MPvsMM_unifrac_simPointLM
setwd(wdOutput)
getwd()
ggsave(paste("IPvsIM_MPvsMM_unifrac_simPointLM_FourStages_rhizosphere_Filtered1",data.set.name,".pdf",sep=""),
       IPvsIM_MPvsMM_unifrac_simPointLM,device=cairo_pdf,width=160,height=70,dpi = 600,units = "mm")

###statistcal analysis
library(rcompanion)
IPvsIM_MPvsMM_unifrac_sim$VS <- factor(IPvsIM_MPvsMM_unifrac_sim$VS)
IPvsIM_MPvsMM_unifrac_sim$Days <- factor(IPvsIM_MPvsMM_unifrac_sim$Days)
str(IPvsIM_MPvsMM_unifrac_sim)
head(IPvsIM_MPvsMM_unifrac_sim)
scheirerRayHare(Similarity~VS*Days+Days:VS, data = IPvsIM_MPvsMM_unifrac_sim)
