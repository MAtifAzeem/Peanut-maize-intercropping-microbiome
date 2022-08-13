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
                              axis.ticks.length = unit(0.8, "mm"))
wdImport<-("E:/working/SCI/Soil Micro/SCI/Figures/Data/Data for submit")
wdOutput <- ("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Supplemental materials/divers_IPvsMP")

#### 3 Import and process data ####
setwd(wdImport)
getwd()
Iron <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                             sheet = "fig s5")

peanut_activeFe_available_cor<-corr.test(Iron$YL_ActiveFe,Iron$AvailableFe,
                                              method="pearson",adjust="BH",minlength=5)
peanut_activeFe_available_cor
peanut_activeFe_available_cor$r
peanut_activeFe_available_cor$p

##plots
peanut_activeFe_available_cor_line <- ggplot(Iron, aes(x= AvailableFe, y= YL_ActiveFe))+
  geom_point(color="#E64B35")+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(label.x = 3, label.y = 14)+
  mytheme
peanut_activeFe_available_cor_line

formula <- y ~ x
peanut_activeFe_available_cor_line <- ggplot(Iron, aes(x= AvailableFe, y= YL_ActiveFe))+
  geom_point(color="#E64B35")+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula = formula,label.x = 3,label.y = 14) +
  stat_cor (method = "pearson",
            label.x = 3)+
  mytheme
peanut_activeFe_available_cor_line
setwd(wdOutput)
getwd()
ggsave(paste("peanut_activeFe_available_cor_line",data.set.name,".pdf",sep=""),
       peanut_activeFe_available_cor_line,
       device=cairo_pdf,width=80,height=60,dpi = 600,units = "mm")

