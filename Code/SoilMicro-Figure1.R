#### 1. loading required libraries ####
library(ggplot2)#作图 plot
library(ggpubr)#添加显著性标记, Add the significance marker
library(ggsignif)#添加显著性标记, Add the significance marker
library(dplyr)#数据清洗，Data cleaning
library(plyr)#数据清洗，Data cleaning
library(reshape2)#数据清洗，Data cleaning
library(ggthemes)#ggplot所用主题，Themes for ggplot2
library(readxl)#读入 excel, read excel
library(showtext)#字体设置, font setting
library(extrafont)#使用系统字体，Using the system fonts
library(sysfonts)#加载系统字体，loading the system fonts
library(Cairo)#抗锯齿,anti-aliasing

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
                              axis.line = element_line(color = "#4D4D4D",size=0.2),
                              axis.ticks.length = unit(0.8, "mm"))#移除整体的边???

wdImport<-("E:/Study/SCI/Soil Micro/SCI/Figures/Data/Data for submit")
wdOutput_Figure1 <- ("E:/Study/SCI/Soil Micro/SCI/Figures/Figures from R/Figure1")


#### 4. Fig. 1C-Field-LER-NetEffect-SamePositon####
#### 4.1 Import and process data ####
setwd(wdImport)
Filed_LER_NE_SamePosition <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                    sheet = "Fig 1C same position")

sapply(Filed_LER_NE_SamePosition,class)
Filed_LER_NE_SamePosition[(6:10),"Value"]<-0-Filed_LER_NE_SamePosition[(6:10),"Value"]
Filed_LER_NE_SamePosition$Year<-factor(Filed_LER_NE_SamePosition$Year, levels = c(2015,2014,2013,2012,2011))
#### 4.1 Plot ####
Filed_LER_NE_SamePosition_Barplot <- ggplot(Filed_LER_NE_SamePosition,
                                             aes(x=Year, y=Value, fill=Indicator)) + 
  geom_bar(stat="identity", position="identity",width = 0.9,color="black",size=0.1)+
  mytheme1+coord_flip()+guides(fill="none")+
  scale_fill_manual(values = c("#bebada","#8dd3c7"))+
  theme(panel.background = element_rect(fill="#fbb4ae",size = 0.1))
Filed_LER_NE_SamePosition_Barplot
setwd(wdOutput_Figure1)
getwd()
ggsave("Filed_LER_NE_SamePosition_Barplot.pdf",
       Filed_LER_NE_SamePosition_Barplot,
       device=cairo_pdf,width=29,height=27,dpi = 300,units = "mm")

#### 5. Fig. 1C-Field-LER-NetEffect-same year####
#### 5.1 Import and process data ####
setwd(wdImport)
Filed_LER_NE_SameYear <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                         sheet = "Fig 1C same year")

sapply(Filed_LER_NE_SameYear,class)
Filed_LER_NE_SameYear[(6:10),"Value"]<-0-Filed_LER_NE_SameYear[(6:10),"Value"]
Filed_LER_NE_SameYear$Farm<-factor(Filed_LER_NE_SameYear$Farm, levels = c("#5","#4","#3","#2","#1"))
#### 5.2 Plot ####
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

#### 6. Fig. 1D-Peanut-Young leaves-SPAD value in young leaves####
#### 6.1 Import and process data ####
setwd(wdImport)
NormalPotPeanutYLSPAD <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                       sheet = "Fig 1D")
#### 6.2 Summary by grouping ####
NormalPotPeanutYLSPADStat <- NormalPotPeanutYLSPAD%>%
  group_by(System,Days,Species,SystemSpecies,SpeciesDays,SystemSpeciesDays)%>%
  summarise_at("YL_SPAD",funs(mean,sd))
PotPeanutYLSPAD2015Stat

#### 6.3 Student's t test ####
#46 days post sowing(dps)
Peanut46YLSPAD <- subset(NormalPotPeanutYLSPAD,SpeciesDays=="Peanut46")
compare_means(data=Peanut46YLSPAD,YL_SPAD~System,method = "t.test")

#53 dps
Peanut53YLSPAD <- subset(NormalPotPeanutYLSPAD,SpeciesDays=="Peanut53")
compare_means(data=Peanut53YLSPAD,YL_SPAD~System,method = "t.test")

#63 dps
Peanut63YLSPAD <- subset(NormalPotPeanutYLSPAD,SpeciesDays=="Peanut63")
compare_means(data=Peanut63YLSPAD,YL_SPAD~System,method = "t.test")

#73 dps
Peanut73YLSPAD <- subset(NormalPotPeanutYLSPAD,SpeciesDays=="Peanut73")
compare_means(data=Peanut73YLSPAD,YL_SPAD~System,method = "t.test")

#80 dps
Peanut80YLSPAD <- subset(NormalPotPeanutYLSPAD,SpeciesDays=="Peanut80")
compare_means(data=Peanut80YLSPAD,YL_SPAD~System,method = "t.test")

#### 7. Fig. 1E Peanut-Intercropping-Pot-Fe status ####
#### 7.1 Import and process data ####
setwd(wdImport)
NormalPotPeanutFe <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                  sheet = "Fig 1E")

NormalPotPeanutFe$Species <- factor(NormalPotPeanutFe$Species,levels = c("Peanut","Maize"))
NormalPotPeanutFe$System <- factor(NormalPotPeanutFe$System,levels = c("Monocropping","Intercropping"))
NormalPotPeanutFe$SystemSpecies <- factor(NormalPotPeanutFe$SystemSpecies,levels = c("MP","IP"))

#### 7.2 Fig.1E Peanut-Intercropping-Pot-Available Fe concentration in rhizosphere####
#### 7.2.1 Statistic analysis ####
ActiveFe_aov<-aov(data=NormalPotPeanutFe,YL_ActiveFe~System+Days+Days*System)
summary(ActiveFe_aov)

#### 7.2.2 Plots ####
Pot_NormalPotPeanut_ActiveFe_line<-ggplot(NormalPotPeanutFe, aes(x=Days, y=YL_ActiveFe, fill=SystemSpecies,group=SystemSpecies)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line",size=0.5,aes(col=SystemSpecies)) +
  stat_summary(fun="mean", geom="point",size=0.5,aes(col=SystemSpecies)) +
  FacetTheme+
  scale_x_continuous(breaks=c(46,53,63,73,80),limits=c(40,84))+
  labs(x="dps",
       y="Concentration",parse =T)+
  scale_color_npg()+scale_fill_npg()+
  theme(legend.position = "top")+
  guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))
Pot_NormalPotPeanut_ActiveFe_line
setwd(wdOutput_Figure1)
getwd()
ggsave("Pot_NormalPotPeanut_ActiveFe_line.pdf",
       Pot_NormalPotPeanut_ActiveFe_line,
       device=cairo_pdf,width=52,height=44,dpi = 300,units = "mm")

#### 7.3 Fig.1E Peanut-Intercropping-Pot-Active Fe concentration in young leaves####
#### 7.3.1 Statistic analysis ####
AvailableFe_aov<-aov(data=NormalPotPeanutFe,AvailableFe~System+Days+Days*System)
summary(AvailableFe_aov)

#### 7.3.2 Plots ####
Pot_NormalPotPeanut_AvailableFe_line<-ggplot(NormalPotPeanutFe, aes(x=Days, y=AvailableFe, fill=SystemSpecies,group=SystemSpecies)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line",size=0.5,aes(col=SystemSpecies)) +
  stat_summary(fun="mean", geom="point",size=0.5,aes(col=SystemSpecies)) +
  FacetTheme+
  scale_x_continuous(breaks=c(46,53,63,73,80),limits=c(40,84))+
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
       device=cairo_pdf,width=52,height=44,dpi = 300,units = "mm")

#### 8. Fig. 1F-Peanut-sterilization-biomass####
#### 8.1 Import and process data ####
setwd(wdImport)
PotSterilization<- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                    sheet = "Fig 1F")
PotSterilizationBiomass<-PotSterilization[,-c(7,8)]
PotSterilizationBiomass<-na.omit(PotSterilizationBiomass)
PotSterilizationBiomass$SoilTreatment <-factor(PotSterilizationBiomass$SoilTreatment,level=c("Sterilization","Normal"))
PotSterilizationBiomass$System <-factor(PotSterilizationBiomass$System,level=c("MP","IP"))

#### 8.2 Statistic analysis####
PotSterilizationBiomass_IP<-filter(PotSterilizationBiomass,System=="IP")
compare_means(data=PotSterilizationBiomass_IP,Biomass~SoilTreatment,method = "t.test")
PotSterilizationBiomass_MP<-filter(PotSterilizationBiomass,System=="MP")
compare_means(data=PotSterilizationBiomass_MP,Biomass~SoilTreatment,method = "t.test")
PotSterilizationBiomass_Sterilization<-filter(PotSterilizationBiomass,SoilTreatment=="Sterilization")
compare_means(data=PotSterilizationBiomass_Sterilization,Biomass~System,method = "t.test")
PotSterilizationBiomass_Normal<-filter(PotSterilizationBiomass,SoilTreatment=="Normal")
compare_means(data=PotSterilizationBiomass_Normal,Biomass~System,method = "t.test")

#### 8.3 Plots ####
Pot_Sterilizatiion_Bimass_Bar<-ggplot(PotSterilizationBiomass,aes(x=SoilTreatment,y=Biomass,color=System,group=System))+
  stat_summary(fun=mean, geom='bar',fill="white",width=.5,position = position_dodge(0.6),size=0.12)+
  stat_summary(fun.data = mean_sdl,geom='errorbar',width=.1,position = position_dodge(0.6),size=0.12)+
  geom_point(shape=1,position = position_jitterdodge(0.1),size=0.1)+
  mytheme+theme(legend.position = "top",legend.direction="horizontal")+
  scale_color_manual(values = c("#000000","#DD7843"))+
  scale_y_continuous(limits=c(0,16))
Pot_Sterilizatiion_Bimass_Bar

setwd(wdOutput_Figure1)
getwd()
ggsave(paste("Pot_Sterilizatiion_Bimass_Bar",".pdf",sep=""),
       Pot_Sterilizatiion_Bimass_Bar,device=cairo_pdf,width=45,height=40,dpi = 300,units = "mm")

#### 9. Fig. 1F-Peanut-sterilization-YL-SPAD####
#### 9.1 Import and process data ####
setwd(wdImport)
PotSterilization<- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                              sheet = "Fig 1F")
PotSterilizationSPAD<-PotSterilization[,-c(8,9)]
PotSterilizationSPAD<-na.omit(PotSterilizationSPAD)
PotSterilizationSPAD$SoilTreatment <-factor(PotSterilizationSPAD$SoilTreatment,level=c("Sterilization","Normal"))
PotSterilizationSPAD$System <-factor(PotSterilizationSPAD$System,level=c("MP","IP"))

#### 9.2 Statistic analysis####
PotSterilizationSPAD_IP<-filter(PotSterilizationSPAD,System=="IP")
compare_means(data=PotSterilizationSPAD_IP,YL_SPAD~SoilTreatment,method = "t.test")
PotSterilizationSPAD_MP<-filter(PotSterilizationSPAD,System=="MP")
compare_means(data=PotSterilizationSPAD_MP,YL_SPAD~SoilTreatment,method = "t.test")
PotSterilizationSPAD_Sterilization<-filter(PotSterilizationSPAD,SoilTreatment=="Sterilization")
compare_means(data=PotSterilizationSPAD_Sterilization,YL_SPAD~System,method = "t.test")
PotSterilizationSPAD_Normal<-filter(PotSterilizationSPAD,SoilTreatment=="Normal")
compare_means(data=PotSterilizationSPAD_Normal,YL_SPAD~System,method = "t.test")

#### 9.3 Plots ####
Pot_Sterilizatiion_YL_SPAD_Bar<-ggplot(PotSterilizationSPAD,aes(x=SoilTreatment,y=YL_SPAD,color=System,group=System))+
  stat_summary(fun=mean, geom='bar',fill="white",width=.5,position = position_dodge(0.6),size=0.12)+
  stat_summary(fun.data = mean_sdl,geom='errorbar',width=.1,position = position_dodge(0.6),size=0.12)+
  geom_point(shape=1,position = position_jitterdodge(0.18),size=0.1)+
  mytheme+theme(legend.position = "top",legend.direction="horizontal")+
  scale_color_manual(values = c("#000000","#DD7843"))+
  scale_y_continuous(limits=c(0,48))
Pot_Sterilizatiion_YL_SPAD_Bar

setwd(wdOutput_Figure1)
getwd()
ggsave(paste("Pot_Sterilizatiion_YL_SPAD_Bar",".pdf",sep=""),
       Pot_Sterilizatiion_YL_SPAD_Bar,device=cairo_pdf,width=45,height=40,dpi = 300,units = "mm")

#### 10. Fig. 1F-Peanut-sterilization-YL-ActiveFe####
#### 10.1 Import and process data ####
setwd(wdImport)
PotSterilization<- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                              sheet = "Fig 1F")
PotSterilizationActiveFe<-PotSterilization[,-c(7,9)]
PotSterilizationActiveFe<-na.omit(PotSterilizationActiveFe)
PotSterilizationActiveFe$SoilTreatment <-factor(PotSterilizationActiveFe$SoilTreatment,level=c("Sterilization","Normal"))
PotSterilizationActiveFe$System <-factor(PotSterilizationActiveFe$System,level=c("MP","IP"))

#### 10.2 Statistic analysis####
PotSterilizationActiveFe_IP<-filter(PotSterilizationActiveFe,System=="IP")
compare_means(data=PotSterilizationActiveFe_IP,YL_ActiveFe~SoilTreatment,method = "t.test")
PotSterilizationActiveFe_MP<-filter(PotSterilizationActiveFe,System=="MP")
compare_means(data=PotSterilizationActiveFe_MP,YL_ActiveFe~SoilTreatment,method = "t.test")
PotSterilizationActiveFe_Sterilization<-filter(PotSterilizationActiveFe,SoilTreatment=="Sterilization")
compare_means(data=PotSterilizationActiveFe_Sterilization,YL_ActiveFe~System,method = "t.test")
PotSterilizationActiveFe_Normal<-filter(PotSterilizationActiveFe,SoilTreatment=="Normal")
compare_means(data=PotSterilizationActiveFe_Normal,YL_ActiveFe~System,method = "t.test")

#### 10.3 Plots ####
Pot_Sterilizatiion_YL_ActiveFe_Bar<-ggplot(PotSterilizationActiveFe,aes(x=SoilTreatment,y=YL_ActiveFe,color=System,group=System))+
  stat_summary(fun=mean, geom='bar',fill="white",width=.5,position = position_dodge(0.6),size=0.12)+
  stat_summary(fun.data = mean_sdl,geom='errorbar',width=.1,position = position_dodge(0.6),size=0.12)+
  geom_point(shape=1,position = position_jitterdodge(0.18),size=0.1)+
  mytheme+theme(legend.position = "top",legend.direction="horizontal")+
  scale_color_manual(values = c("#000000","#DD7843"))+
  scale_y_continuous(limits=c(0,15))
Pot_Sterilizatiion_YL_ActiveFe_Bar

setwd(wdOutput_Figure1)
getwd()
ggsave(paste("Pot_Sterilizatiion_YL_ActiveFe_Bar",".pdf",sep=""),
       Pot_Sterilizatiion_YL_ActiveFe_Bar,device=cairo_pdf,width=45,height=40,dpi = 300,units = "mm")

