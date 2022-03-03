#### 1. loading required libraries ####
library(ggplot2)#作图 plot
library(ggpubr)#添加显著性标记, Add the significance marker
library(ggsignif)#添加显著性标记, Add the significance marker
library(dplyr)#数据清洗，Data cleaning
library(plyr)#数据清洗，Data cleaning
library(reshape2)#数据清洗，Data cleaning
library(ggthemes)#ggplot所用主题，Themes for ggplot2
library(grid)#分面和嵌合图，facet and Mosaic graph
library(agricolae)#多重比较，Multiple comparisons.
library(readxl)#读入 excel, read excel
library(ggsci)#配色，color scheme
library(showtext)#字体设置, font setting
library(car)#方差齐性检验，homogeneity test of variance, levene test
library(extrafont)#使用系统字体，Using the system fonts
library(sysfonts)#加载系统字体，loading the system fonts
library(Cairo)#抗锯齿,anti-aliasing
library(psych)
loadfonts()
Sys.setenv(R_GSCMD = "C:/Program Files (x86)/gs/gs9.50/bin/gswin32c.exe")
#### 2. Setting themes and working dictionary path ####
mytheme <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#000000"),
                             text = element_text(family = "Arial"),
                             strip.text = element_text(size = 8,hjust = 0.5),
                             plot.title = element_text(size = 8,hjust = 0.5),
                             axis.text=element_text(size=8,color = "#626262"),
                             axis.title=element_text(size = 8,color = "black"),
                             legend.text = element_text(size = 8),
                             legend.title = element_text(size = 8),
                             legend.background = element_blank(),
                             panel.border = element_rect(colour = NA),
                             axis.line = element_line(color = "black",size=1))#移除整体的边???

FacetTheme <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#000000"),
                                text = element_text(family = "Arial"),
                                strip.text = element_text(size = 8,hjust = 0.5),
                                plot.title = element_text(size = 8,hjust = 0.5),
                                axis.text =element_text(size=8,color = "#626262"),
                                axis.title =element_text(size = 8,color = "black"),
                                legend.text = element_text(size = 8),
                                legend.title = element_text(size = 8),
                                legend.background = element_blank(),
                                axis.line = element_line(color = "black",size=0.4))#移除整体的边???
wdImport<-("E:/Study/SCI/Soil Micro/SCI/Figures/Data/Data for submit")
wdOutput<- ("E:/Study/SCI/Soil Micro/SCI/Figures/Figures from R/Supplemental materials/Pot2015_available_active Fe")
####Peanut-YL-SPAD-Pot2020####
setwd(wdImport)
Iron_Peanut <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                       sheet = "fig s5")
Iron_Peanut

Peanut_2015_Avaliable_Active_Fe_cor<-corr.test(Iron_Peanut$YL_ActiveFe,Iron_Peanut$AvailableFe,
                            method="pearson",minlength=5)
Peanut_2015_Avaliable_Active_Fe_cor
#Plot#
Peanut_2015_Avaliable_Active_Fe <- ggscatter(Iron_Peanut, x = "AvailableFe", y = "YL_ActiveFe",color = "#E64B35",
                             cor.method = "pearson",add.params = list(color = "#E64B35",fill="lightgrey",size=1),size = 1,conf.int = T,
                             add = "reg.line")+
  labs(x=expression('Available Fe (μg '*g^{-1}*')'),parse =T,
       y=expression('Active Fe (μg '*g^{-1}*')'),parse =T,
       title="")+
  scale_y_continuous(limits=c(4,15))+
  FacetTheme+
  stat_cor(aes(
    label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    size=2.5,method = "pearson",cor.coef.name="r",
    label.x.npc ="left",label.y.npc = "top")
Peanut_2015_Avaliable_Active_Fe
setwd(wdOutput)
ggsave("Peanut_2015_Avaliable_Active_Fe.pdf",device=cairo_pdf,width=60,height=60,dpi = 300,units = "mm")
getwd()

