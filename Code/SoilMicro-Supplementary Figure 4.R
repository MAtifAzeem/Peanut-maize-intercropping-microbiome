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
library(ComplexHeatmap)#heatmapping
library(pheatmap)#heatmapping
library(cluster)#hierarchical clustering
library(psych)#correlationship analysis#
library(VennDiagram)#Venn plot

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
                              strip.text = element_text(size = 5,hjust = 0.5),
                              plot.title = element_text(size = 5,hjust = 0.5),
                              axis.text=element_text(size=5,color = "#4D4D4D"),
                              axis.title=element_text(size = 5),
                              legend.text = element_text(size = 5),
                              legend.title = element_text(size = 5),
                              legend.background = element_blank(),
                              axis.line = element_line(color = "#4D4D4D",size=0.2),
                              axis.ticks.length = unit(0.8, "mm"))#移除整体的边???

wdImport<-("E:/working/SCI/Soil Micro/SCI/Figures/Data/Data for submit")
wdOutput_Figure3 <- ("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Figure3")

#### 3. Fig. 3a-LefSe_MPvsIP_results####
### 3.1 Import and process data ###
setwd(wdImport)
LefSe_MPvsIP <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                    sheet = "fig s4")
LefSe_MPvsIP$ID<-factor(LefSe_MPvsIP$ID,levels=LefSe_MPvsIP$ID)
### 3.2 Plots ###
LefSe_MPvsIP_LDA <- ggplot(LefSe_MPvsIP,aes(x=ID, y=LDA_plus_minus, fill=group)) + 
  geom_bar(stat="identity", position="identity",width = 0.8,color="black",size=0.1)+
  scale_fill_manual(values=c("#E31A1C","#1F78B4"))+
  scale_y_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3,4))+
  mytheme1+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
LefSe_MPvsIP_LDA
setwd(wdOutput_Figure3)
getwd()
ggsave(paste("LefSe_MPvsIP_LDA_3.0_p0.05",".pdf",sep=""),
       LefSe_MPvsIP_LDA,device=cairo_pdf,width=140,height=90,dpi = 600,units = "mm")
