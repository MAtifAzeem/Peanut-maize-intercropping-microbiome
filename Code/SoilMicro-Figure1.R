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
IPvsMP_Field_SPAD<-IPvsMP_Field_SPAD%>%mutate(boxcox_YL_SPAD=BoxCox(IPvsMP_Field_SPAD$YL_SPAD,BoxCox.lambda(IPvsMP_Field_SPAD$YL_SPAD)))
leveneTest(boxcox_YL_SPAD ~ Treatment, data = IPvsMP_Field_SPAD)#p>0.05，则满足方差齐性
shapiro.test(IPvsMP_Field_SPAD$boxcox_YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
compare_means(YL_SPAD~Treatment,IPvsMP_Field_SPAD,method="wilcox.test")

IPvsMP_Field_ActiveFe<- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                   sheet = "fig 1b ActiveFe")

leveneTest(YL_ActiveFe ~ Treatment, data = IPvsMP_Field_ActiveFe)#p>0.05，则满足方差齐性
shapiro.test(IPvsMP_Field_ActiveFe$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
IPvsMP_Field_ActiveFe<-IPvsMP_Field_ActiveFe%>%mutate(boxcox_YL_ActiveFe=BoxCox(IPvsMP_Field_ActiveFe$YL_ActiveFe,BoxCox.lambda(IPvsMP_Field_ActiveFe$YL_ActiveFe)))
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