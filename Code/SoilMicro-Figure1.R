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
library(PMCMRplus)
library(car)
library(rcompanion)
library(forecast)#box-cox数据变换
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
                             axis.line = element_line(color = "#4D4D4D",linewidth=0.2),
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
wdOutput_Figure1 <- ("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Figure1")
#### 4. Fig. 1b-Field-SPAD and Acitve Fe ####
setwd(wdImport)
IPvsMP_Field_SPAD<- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                               sheet = "fig 1b SPAD")

leveneTest(YL_SPAD ~ Treatment, data = IPvsMP_Field_SPAD)#p>0.05，则满足方差齐性
shapiro.test(IPvsMP_Field_SPAD$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
IPvsMP_Field_SPAD<-IPvsMP_Field_SPAD%>%mutate(boxcox_YL_SPAD=BoxCox(IPvsMP_Field_SPAD$YL_SPAD,lambda = "auto"))
leveneTest(boxcox_YL_SPAD ~ Treatment, data = IPvsMP_Field_SPAD)#p>0.05，则满足方差齐性
shapiro.test(IPvsMP_Field_SPAD$boxcox_YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
compare_means(YL_SPAD~Treatment,IPvsMP_Field_SPAD,method="wilcox.test")

IPvsMP_Field_ActiveFe<- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                   sheet = "fig 1b ActiveFe")

leveneTest(YL_ActiveFe ~ Treatment, data = IPvsMP_Field_ActiveFe)#p>0.05，则满足方差齐性
shapiro.test(IPvsMP_Field_ActiveFe$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
IPvsMP_Field_ActiveFe<-IPvsMP_Field_ActiveFe%>%mutate(boxcox_YL_ActiveFe=BoxCox(IPvsMP_Field_ActiveFe$YL_ActiveFe,lambda="auto"))
leveneTest(boxcox_YL_ActiveFe ~ Treatment, data = IPvsMP_Field_ActiveFe)#p>0.05，则满足方差齐性
shapiro.test(IPvsMP_Field_ActiveFe$boxcox_YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
compare_means(YL_ActiveFe~Treatment,IPvsMP_Field_ActiveFe,method="wilcox.test")

#### 5. Fig. 1c-Field-LER-NetEffect-SamePositon####
### 5.1 Import and process data ###
setwd(wdImport)
Filed_LER_NE_SamePosition <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                        sheet = "fig 1c same farm")

sapply(Filed_LER_NE_SamePosition,class)
Filed_LER_NE_SamePosition[(6:10),"Value"]<-0-Filed_LER_NE_SamePosition[(6:10),"Value"]
Filed_LER_NE_SamePosition$Year<-factor(Filed_LER_NE_SamePosition$Year, levels = c(2015,2014,2013,2012,2011))
### 5.1 Plot ###
Filed_LER_NE_SamePosition_Barplot <- ggplot(Filed_LER_NE_SamePosition,
                                            aes(x=Year, y=Value, fill=Indicator)) + 
  geom_bar(stat="identity", position="identity",width = 0.9,color="black",size=0.1)+
  mytheme1+coord_flip()+guides(fill="none")+
  scale_fill_manual(values = c("#bebada","#8dd3c7"))+
  theme(panel.background = element_rect(fill="#fbb4ae",size = 0.1))
Filed_LER_NE_SamePosition_Barplot
setwd(wdOutput_figure1)
getwd()
ggsave("Filed_LER_NE_SamePosition_Barplot.pdf",
       Filed_LER_NE_SamePosition_Barplot,
       device=cairo_pdf,width=29,height=27,dpi = 300,units = "mm")

#### 6. fig. 1c-Field-LER-NetEffect-same year####
### 6.1 Import and process data ###
setwd(wdImport)
Filed_LER_NE_SameYear <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                    sheet = "fig 1c same year")

sapply(Filed_LER_NE_SameYear,class)
Filed_LER_NE_SameYear[(6:10),"Value"]<-0-Filed_LER_NE_SameYear[(6:10),"Value"]
Filed_LER_NE_SameYear$Farm<-factor(Filed_LER_NE_SameYear$Farm, levels = c("#5","#4","#3","#2","#1"))
### 6.2 Plot ###
Filed_LER_NE_SameYear_Barplot <- ggplot(Filed_LER_NE_SameYear,
                                        aes(x=Farm, y=Value, fill=Indicator)) + 
  geom_bar(stat="identity", position="identity",width = 0.9,color="black",size=0.1)+
  mytheme1+coord_flip()+guides(fill="none")+
  scale_fill_manual(values = c("#bebada","#8dd3c7"))+
  theme(panel.background = element_rect(fill="#ffffb3"))
Filed_LER_NE_SameYear_Barplot
setwd(wdOutput_Figure1)
getwd()
ggsave("Filed_LER_NE_SameYear_Barplot.pdf",
       Filed_LER_NE_SameYear_Barplot,
       device=cairo_pdf,width=25,height=25,dpi = 300,units = "mm")

#### 7. Fig. 1d-Peanut-Young leaves-SPAD value in young leaves####
### 7.1 Import and process data ###
setwd(wdImport)
NormalPotPeanutYLSPAD <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                    sheet = "fig 1d")
### 7.2 Summary by grouping ###
NormalPotPeanutYLSPAD_summary <- NormalPotPeanutYLSPAD%>%
  group_by(System,Days,Species,SystemSpecies,SpeciesDays,SystemSpeciesDays)%>%
  summarise_at("YL_SPAD",funs(mean,sd))
NormalPotPeanutYLSPAD_summary

### 7.3 Statistic analysis ###
#Statistic analysis by days post sowing(dps)
#53 dps
NormalPotPeanutYLSPAD_53 <-NormalPotPeanutYLSPAD%>%filter(Days=="53")
leveneTest(YL_SPAD ~ SystemSpeciesDays, data = NormalPotPeanutYLSPAD_53)#p>0.05，则满足方差齐性
shapiro.test(NormalPotPeanutYLSPAD_53$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
compare_means(YL_SPAD~SystemSpeciesDays,NormalPotPeanutYLSPAD_53,method="t.test")

#63 dps
NormalPotPeanutYLSPAD_63 <-NormalPotPeanutYLSPAD%>%filter(Days=="63")
leveneTest(YL_SPAD ~ SystemSpeciesDays, data = NormalPotPeanutYLSPAD_63)#p>0.05，则满足方差齐性
shapiro.test(NormalPotPeanutYLSPAD_63$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
NormalPotPeanutYLSPAD_63<-NormalPotPeanutYLSPAD_63%>%mutate(boxcox_YL_SPAD=BoxCox(NormalPotPeanutYLSPAD_63$YL_SPAD,lambda="auto"))
leveneTest(boxcox_YL_SPAD ~ SystemSpeciesDays, data = NormalPotPeanutYLSPAD_63)#p>0.05，则满足方差齐性
shapiro.test(NormalPotPeanutYLSPAD_63$boxcox_YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
compare_means(YL_SPAD~SystemSpeciesDays,NormalPotPeanutYLSPAD_63,method="wilcox.test")

#73 dps
NormalPotPeanutYLSPAD_73 <-NormalPotPeanutYLSPAD%>%filter(Days=="73")
leveneTest(YL_SPAD ~ SystemSpeciesDays, data = NormalPotPeanutYLSPAD_73)#p>0.05，则满足方差齐性
shapiro.test(NormalPotPeanutYLSPAD_73$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
compare_means(YL_SPAD~SystemSpeciesDays,NormalPotPeanutYLSPAD_73,method="t.test")

#### 8. Fig. 1e Peanut-Intercropping-Pot-Fe status ####
### 8.1 Import and process data ###
setwd(wdImport)
NormalPotPeanutFe <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                sheet = "fig 1e")
NormalPotPeanutFe$Species <- factor(NormalPotPeanutFe$Species,levels = c("Peanut","Maize"))
NormalPotPeanutFe$System <- factor(NormalPotPeanutFe$System,levels = c("Intercropping","Monocropping"))
NormalPotPeanutFe$SystemSpecies <- factor(NormalPotPeanutFe$SystemSpecies,levels = c("IP","MP"))

### 8.2 Fig.1E Peanut-Intercropping-Pot-Active Fe concentration in young leaves ###
## 8.2.1 Statistic analysis ##
leveneTest(YL_ActiveFe ~ SystemSpeciesDays, data = NormalPotPeanutFe)#p>0.05，则满足方差齐性
shapiro.test(NormalPotPeanutFe$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
bb<-compare_means(data=NormalPotPeanutFe,YL_ActiveFe~SystemSpeciesDays,method = "t.test")
bb

summary(aov(data=NormalPotPeanutFe,YL_ActiveFe~System+Days+Days*System))
summary.aov(lm(data=NormalPotPeanutFe, YL_ActiveFe~System*Days))

leveneTest(AvailableFe ~ SystemSpeciesDays, data = NormalPotPeanutFe)#p>0.05，则满足方差齐性
shapiro.test(NormalPotPeanutFe$AvailableFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
summary(aov(data=NormalPotPeanutFe,AvailableFe~System+Days+Days*System))
## 8.2.2 Plots ##
Pot_NormalPotPeanut_AvailableFe_line<-ggplot(NormalPotPeanutFe, aes(x=Days, y=AvailableFe, fill=SystemSpecies,group=SystemSpecies)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line",size=0.5,aes(col=SystemSpecies)) +
  stat_summary(fun="mean", geom="point",size=0.5,aes(col=SystemSpecies)) +
  mytheme1+
  scale_x_continuous(breaks=c(46,53,63,73),limits=c(38,80))+
  scale_y_continuous(limits=c(0,6))+
  labs(x="dps",
       y="Concentration",parse =T)+
  scale_color_npg()+scale_fill_npg()+
  theme(legend.position = "top")+
  guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))
Pot_NormalPotPeanut_AvailableFe_line
setwd(wdOutput_Figure1)
getwd()
ggsave("Pot_NormalPotPeanut_AvailableFe_line.pdf",
       Pot_NormalPotPeanut_AvailableFe_line,
       device=cairo_pdf,width=44,height=44,dpi = 300,units = "mm")

Pot_NormalPotPeanut_YL_ActiveFe_line<-ggplot(NormalPotPeanutFe, aes(x=Days, y=YL_ActiveFe, fill=SystemSpecies,group=SystemSpecies)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line",size=0.5,aes(col=SystemSpecies)) +
  stat_summary(fun="mean", geom="point",size=0.5,aes(col=SystemSpecies)) +
  mytheme1+
  scale_x_continuous(breaks=c(46,53,63,73),limits=c(38,80))+
  scale_y_continuous(limits=c(0,14))+
  labs(x="dps",
       y="Concentration",parse =T)+
  scale_color_npg()+scale_fill_npg()+
  theme(legend.position = "top")+
  guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))
Pot_NormalPotPeanut_YL_ActiveFe_line
setwd(wdOutput_Figure1)
getwd()
ggsave("Pot_NormalPotPeanut_YL_ActiveFe_line.pdf",
       Pot_NormalPotPeanut_YL_ActiveFe_line,
       device=cairo_pdf,width=44,height=44,dpi = 300,units = "mm")

#### 9. Fig. 1e-Peanut-sterilization-YL-SPAD####
### 9.1 Import and process data ###
setwd(wdImport)
PotSterilization<- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                              sheet = "fig 1g")
PotSterilizationSPAD<-PotSterilization[,-c(8:9)]
PotSterilizationSPAD<-na.omit(PotSterilizationSPAD)
PotSterilizationSPAD$SoilTreatment <-factor(PotSterilizationSPAD$SoilTreatment,level=c("Normal","Sterilization"))
PotSterilizationSPAD$System <-factor(PotSterilizationSPAD$System,level=c("MP","IP"))

### 9.2 Statistic analysis ###
leveneTest(YL_SPAD ~ Treatment, data = PotSterilizationSPAD)#p>0.05，则满足方差齐性
shapiro.test(PotSterilizationSPAD$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
summary(aov(data=PotSterilizationSPAD, YL_SPAD~System+SoilTreatment+System*SoilTreatment))
aov_model_PotSterilization_SPAD<-aov(data=PotSterilization,YL_SPAD~Treatment)
LSD_model_PotSterilization_SPAD<- LSD.test(aov_model_PotSterilization_SPAD,"Treatment",p.adj = "BH")
LSD_model_PotSterilization_SPAD

### 9.3 Plots ###
Pot_Sterilizatiion_YL_SPAD_Bar<-ggplot(PotSterilizationSPAD,aes(x=SoilTreatment,y=YL_SPAD,group=System))+
  stat_summary(fun=mean, geom='bar',aes(fill=System),color="black",width=.25,position = position_dodge(0.6),size=0.15)+
  stat_summary(fun.data = mean_sdl,geom='errorbar',width=.1,position = position_dodge(0.6),size=0.2)+
  geom_point(aes(fill=System),shape=21,position = position_jitterdodge(0.18),size=0.8,stroke = 0.2,alpha=0.6)+
  mytheme+theme(legend.position = "top",legend.direction="horizontal")+
  scale_fill_manual(values = c("#4DBBD5","#E64B35"))+
  scale_y_continuous(limits=c(0,48))
Pot_Sterilizatiion_YL_SPAD_Bar
setwd(wdOutput_Figure1)
getwd()
ggsave(paste("Pot_Sterilizatiion_YL_SPAD_Bar",".pdf",sep=""),
       Pot_Sterilizatiion_YL_SPAD_Bar,device=cairo_pdf,width=45,height=40,dpi = 300,units = "mm")

#### 10. Fig. 1e-Peanut-sterilization-YL-ActiveFe ####
### 10.1 Import and process data ###
setwd(wdImport)
PotSterilization<- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                              sheet = "fig 1g")
PotSterilizationActiveFe<-PotSterilization[,-c(7,9)]
PotSterilizationActiveFe<-na.omit(PotSterilizationActiveFe)
PotSterilizationActiveFe$SoilTreatment <-factor(PotSterilizationActiveFe$SoilTreatment,level=c("Normal","Sterilization"))
PotSterilizationActiveFe$System <-factor(PotSterilizationActiveFe$System,level=c("MP","IP"))
PotSterilizationActiveFe<-PotSterilizationActiveFe
dim(PotSterilizationActiveFe)
### 10.2 Statistic analysis ###
leveneTest(YL_ActiveFe ~ Treatment, data = PotSterilizationActiveFe)#p>0.05，则满足方差齐性
shapiro.test(PotSterilizationActiveFe$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
#BoxCox.lambda(PotSterilizationActiveFe$YL_ActiveFe,method="loglik")#a number indicating the Box-Cox transformation parameter
PotSterilizationActiveFe<-PotSterilizationActiveFe%>%mutate(boxcox_YL_ActiveFe=BoxCox(PotSterilizationActiveFe$YL_ActiveFe,lambda="auto"))
shapiro.test(PotSterilizationActiveFe$boxcox_YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
summary(aov(data=PotSterilizationActiveFe,boxcox_YL_ActiveFe~System+SoilTreatment+System*SoilTreatment))
aov_model_PotSterilization_boxcox_ActiveFe<-aov(data=PotSterilizationActiveFe,boxcox_YL_ActiveFe~Treatment)
LSD_model_PotSterilization_boxcox_ActiveFe<- LSD.test(aov_model_PotSterilization_boxcox_ActiveFe,"Treatment",p.adj = "BH")
LSD_model_PotSterilization_boxcox_ActiveFe

### 10.3 Plots ###
Pot_Sterilizatiion_YL_ActiveFe_Bar<-ggplot(PotSterilizationActiveFe,aes(x=SoilTreatment,y=YL_ActiveFe,group=System))+
  stat_summary(fun=mean, geom='bar',aes(fill=System),color="black",width=.25,position = position_dodge(0.6),size=0.15)+
  stat_summary(fun.data = mean_sdl,geom='errorbar',width=.1,position = position_dodge(0.6),size=0.2)+
  geom_point(aes(fill=System),shape=21,position = position_jitterdodge(0.18),size=0.8,stroke = 0.2)+
  mytheme+theme(legend.position = "top",legend.direction="horizontal")+
  scale_fill_manual(values = c("#4DBBD5","#E64B35"))+
  scale_y_continuous(limits=c(0,15))
Pot_Sterilizatiion_YL_ActiveFe_Bar

setwd(wdOutput_Figure1)
getwd()
ggsave(paste("Pot_Sterilizatiion_YL_ActiveFe_Bar",".pdf",sep=""),
       Pot_Sterilizatiion_YL_ActiveFe_Bar,device=cairo_pdf,width=45,height=40,dpi = 300,units = "mm")

#### 11. Fig. 1e-Peanut-sterilization-YL-AvailableFe ####
### 11.1 Import and process data ###
setwd(wdImport)
PotSterilization<- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                              sheet = "fig 1g")
PotSterilizationAvailableFe<-PotSterilization[,-c(7:8)]
PotSterilizationAvailableFe<-na.omit(PotSterilizationAvailableFe)
PotSterilizationAvailableFe$SoilTreatment <-factor(PotSterilizationAvailableFe$SoilTreatment,level=c("Normal","Sterilization"))
PotSterilizationAvailableFe$System <-factor(PotSterilizationAvailableFe$System,level=c("MP","IP"))

### 11.2 Statistic analysis ###
leveneTest(AvailableFe ~ Treatment, data = PotSterilizationAvailableFe)#p>0.05，则满足方差齐性
shapiro.test(PotSterilizationAvailableFe$AvailableFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
summary(aov(data=PotSterilizationAvailableFe, AvailableFe~System+SoilTreatment+System*SoilTreatment))
aov_model_PotSterilization_AvailableFe<-aov(data=PotSterilizationAvailableFe,AvailableFe~Treatment)
summary(aov_model_PotSterilization_AvailableFe)
LSD_model_PotSterilization_AvailableFe<- LSD.test(aov_model_PotSterilization_AvailableFe,"Treatment",p.adj = "BH")
LSD_model_PotSterilization_AvailableFe

### 11.3 Plots ###
Pot_Sterilizatiion_AvailableFe_Bar<-ggplot(PotSterilizationAvailableFe,aes(x=SoilTreatment,y=AvailableFe,group=System))+
  stat_summary(fun=mean, geom='bar',aes(fill=System),color="black",width=.25,position = position_dodge(0.6),size=0.15)+
  stat_summary(fun.data = mean_sdl,geom='errorbar',width=.1,position = position_dodge(0.6),size=0.2)+
  geom_point(aes(fill=System),shape=21,position = position_jitterdodge(0.18),size=0.8,stroke = 0.2)+
  mytheme+theme(legend.position = "top",legend.direction="horizontal")+
  scale_fill_manual(values = c("#4DBBD5","#E64B35"))+
  scale_y_continuous(limits=c(0,15))
Pot_Sterilizatiion_AvailableFe_Bar

setwd(wdOutput_Figure1)
getwd()
ggsave(paste("Pot_Sterilizatiion_AvailableFe_Bar",".pdf",sep=""),
       Pot_Sterilizatiion_AvailableFe_Bar,device=cairo_pdf,width=45,height=40,dpi = 300,units = "mm")

