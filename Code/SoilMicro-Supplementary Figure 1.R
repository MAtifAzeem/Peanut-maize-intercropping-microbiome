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
wdOutput <- ("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Supplemental materials/IMvsMM_iron")

#### 4. Fig. 1d-Maize-Young leaves-SPAD value in young leaves####
### 4.1 Import and process data ###
setwd(wdImport)
NormalPotMaizeYLSPAD <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                    sheet = "fig s1 SPAD")
NormalPotMaizeYLSPAD$System <- factor(NormalPotMaizeYLSPAD$System,levels = c("Intercropping","Monocropping"))
NormalPotMaizeYLSPAD$SystemSpecies <- factor(NormalPotMaizeYLSPAD$SystemSpecies,levels = c("IM","MM"))

### 4.2 Summary by grouping ###
NormalPotMaizeYLSPAD_summary <- NormalPotMaizeYLSPAD%>%
  group_by(System,Days,Species,SystemSpecies,SpeciesDays,SystemSpeciesDays)%>%
  summarise_at("YL_SPAD",funs(mean,sd))
NormalPotMaizeYLSPAD_summary
## 4.2.1Statistic analysis ###
leveneTest(YL_SPAD ~ RhizocompartmentsSystemSpeciesDays, data = NormalPotMaizeYLSPAD)#p>0.05，则满足方差齐性
shapiro.test(NormalPotMaizeYLSPAD$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
NormalPotMaizeYLSPAD_aov<-aov(data=NormalPotMaizeYLSPAD,YL_SPAD~System+Days+Days*System)
summary(NormalPotMaizeYLSPAD_aov)
## 4.3 Plots ##
Pot_NormalPotMaize_SPAD_line<-ggplot(NormalPotMaizeYLSPAD, aes(x=Days, y=YL_SPAD, fill=SystemSpecies,group=SystemSpecies)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line",size=0.5,aes(col=SystemSpecies)) +
  stat_summary(fun="mean", geom="point",size=0.5,aes(col=SystemSpecies)) +
  mytheme1+
  scale_x_continuous(breaks=c(46,53,63,73),limits=c(38,80))+
  scale_y_continuous(limits=c(limits=c(0,45)))+
  labs(x="dps",
       y="Concentration",parse =T)+
  scale_color_npg()+scale_fill_npg()+
  theme(legend.position = "top")+
  guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))
Pot_NormalPotMaize_SPAD_line
setwd(wdOutput)
getwd()
ggsave("Pot_NormalPotMaize_SPAD_line.pdf",
       Pot_NormalPotMaize_SPAD_line,
       device=cairo_pdf,width=44,height=44,dpi = 300,units = "mm")

#### 5. Fig. 1d-Maize-Young leaves-active Fe in young leaves####
### 5.1 Import and process data ###
setwd(wdImport)
NormalPotMaize_iron <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                   sheet = "fig s1 iron")

NormalPotMaize_iron$System <- factor(NormalPotMaize_iron$System,levels = c("Intercropping","Monocropping"))
NormalPotMaize_iron$SystemSpecies <- factor(NormalPotMaize_iron$SystemSpecies,levels = c("IM","MM"))

### 5.2 Summary by grouping ###
NormalPotMaize_iron_summary <- NormalPotMaize_iron%>%
  group_by(System,Days,Species,SystemSpecies,SpeciesDays,SystemSpeciesDays)%>%
  summarise_at(c("YL_ActiveFe","AvailableFe"),funs(mean,sd))
NormalPotMaize_iron_summary
## 5.2.1 Statistic analysis ###
leveneTest(YL_ActiveFe ~ RhizocompartmentsSystemSpeciesDays, data = NormalPotMaize_iron)#p>0.05，则满足方差齐性
shapiro.test(NormalPotMaize_iron$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
NormalPotMaize_activeIron_aov<-aov(data=NormalPotMaize_iron,YL_ActiveFe~System+Days+Days*System)
summary(NormalPotMaize_activeIron_aov)
NormalPotMaize_activeIron_lm.fit<-lm(data=NormalPotMaize_iron, YL_ActiveFe~System*Days)
summary.aov(NormalPotMaize_activeIron_lm.fit)

## 5.3 Plots ##
Pot_NormalPotMaize_activeIron_line<-ggplot(NormalPotMaize_iron, aes(x=Days, y=YL_ActiveFe, fill=SystemSpecies,group=SystemSpecies)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line",size=0.5,aes(col=SystemSpecies)) +
  stat_summary(fun="mean", geom="point",size=0.5,aes(col=SystemSpecies)) +
  mytheme1+
  scale_x_continuous(breaks=c(46,53,63,73),limits=c(38,80))+
  scale_y_continuous(limits=c(0,16))+
  labs(x="dps",
       y="Concentration",parse =T)+
  scale_color_npg()+scale_fill_npg()+
  theme(legend.position = "top")+
  guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))
Pot_NormalPotMaize_activeIron_line
setwd(wdOutput)
getwd()
ggsave("Pot_NormalPotMaize_activeIron_line.pdf",
       Pot_NormalPotMaize_activeIron_line,
       device=cairo_pdf,width=44,height=44,dpi = 300,units = "mm")

#### 6. Fig. 1d-Maize-Young leaves-available Fe####
### 6.1 Import and process data ###
setwd(wdImport)
NormalPotMaize_iron <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                  sheet = "fig s1 iron")

NormalPotMaize_iron$System <- factor(NormalPotMaize_iron$System,levels = c("Intercropping","Monocropping"))
NormalPotMaize_iron$SystemSpecies <- factor(NormalPotMaize_iron$SystemSpecies,levels = c("IM","MM"))

### 6.2 Summary by grouping ###
NormalPotMaize_iron_summary <- NormalPotMaize_iron%>%
  group_by(System,Days,Species,SystemSpecies,SpeciesDays,SystemSpeciesDays)%>%
  summarise_at(c("YL_ActiveFe","AvailableFe"),funs(mean,sd))
NormalPotMaize_iron_summary

## 6.2.1 Statistic analysis ###
leveneTest(AvailableFe ~ RhizocompartmentsSystemSpeciesDays, data = NormalPotMaize_iron)#p>0.05，则满足方差齐性
shapiro.test(NormalPotMaize_iron$AvailableFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
NormalPotMaize_AvailableFe_aov<-aov(data=NormalPotMaize_iron,AvailableFe~System+Days+Days*System)
summary(NormalPotMaize_AvailableFe_aov)

## 6.3 Plots ##
Pot_NormalPotMaize_AvailableFe_line<-ggplot(NormalPotMaize_iron, aes(x=Days, y=AvailableFe, fill=SystemSpecies,group=SystemSpecies)) +
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
Pot_NormalPotMaize_AvailableFe_line
setwd(wdOutput)
getwd()
ggsave("Pot_NormalPotMaize_AvailableFe_line.pdf",
       Pot_NormalPotMaize_AvailableFe_line,
       device=cairo_pdf,width=44,height=44,dpi = 300,units = "mm")
