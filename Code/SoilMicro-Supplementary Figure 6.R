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

wdImport<-("E:/Study/SCI/Soil Micro/SCI/Figures/Data/Data for submit")
wdOutput <- ("E:/Study/SCI/Soil Micro/SCI/Figures/Figures from R/Supplemental materials/Top10drivers_MPvsIP")


#### 3. Heatmap_Top10_MPvsIP####
### 3.1 Import and process data ###
setwd(wdImport)
getwd()
Top10drivers_MPvsIP <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                           sheet = "fig s4")
Top10drivers_MPvsIP_mean<-Top10drivers_MPvsIP%>%group_by(SystemSpeciesDays,SystemSpecies,System,Species,Days)%>%
  summarise_at(vars(Pseudomonas:AKYG1722),funs(mean))
Top10drivers_MPvsIP_mean<-Top10drivers_MPvsIP_mean[c(13:16,5:8,9:12,1:4),]
Top10drivers_MPvsIP_RelativeAbundance <- Top10drivers_MPvsIP_mean[,-c(1:5)]
Top10drivers_MPvsIP_RelativeAbundance_MP_IP_MM<-Top10drivers_MPvsIP_RelativeAbundance[c(1:12),]
Top10drivers_MPvsIP_RelativeAbundance_MP_IP_MM_Normalization<-t(scale(Top10drivers_MPvsIP_RelativeAbundance_MP_IP_MM,center = T))
colnames(Top10drivers_MPvsIP_RelativeAbundance_MP_IP_MM_Normalization)<-Top10drivers_MPvsIP_mean$SystemSpeciesDays[c(1:12)]
View(Top10drivers_MPvsIP_RelativeAbundance_MP_IP_MM_Normalization)
### 4.2 Complexheatmap ###
## 4.2.1 color plate ##
max(Top10drivers_MPvsIP_RelativeAbundance_MP_IP_MM_Normalization)
min(Top10drivers_MPvsIP_RelativeAbundance_MP_IP_MM_Normalization)
col_fun = circlize::colorRamp2(c(-2, 0, 3), c("#2166ac", "white", "#b2182b"))
col_fun(seq(-3, 3))
## 4.2.2 split ##
column_split_order<-c(rep("MP",4),rep("IP",4),rep("MM",4))
column_split_order<-factor(column_split_order, levels = c("IP","MP","MM"))
## 4.2.3 annotation ##
Top10drivers_MPvsIP_IP <-filter(Top10drivers_MPvsIP,SystemSpecies=="IP")
Top10drivers_MPvsIP_IP_value <-Top10drivers_MPvsIP_IP[,-c(1:13)]
Top10drivers_MPvsIP_IP_value<-t(Top10drivers_MPvsIP_IP_value)
ha = rowAnnotation(RA_mean = anno_boxplot(Top10drivers_MPvsIP_IP_value,height = unit(4, "cm"),
                                          box_width = 0.5))
ha
## 4.2.4 plots ##
Heatmap(Top10drivers_MPvsIP_RelativeAbundance_MP_IP_MM_Normalization,col = col_fun,
           cluster_columns = F,column_split = column_split_order,
           name="Relative abundance",row_km = 4,row_gap = unit(c(2,2,2), "mm"),
           column_title_gp = gpar(fill = c("#E31A1C", "#1F78B4", "#A6CEE3")),
           right_annotation = ha, row_dend_reorder = TRUE,
           clustering_distance_rows = "spearman")
setwd(wdOutput)

#save as width*height7*4

#### 4. fig 3b Correlationship with active Fe and available Fe in peanuts ####
### 4.1 Import and process data ###
setwd(wdImport)
getwd()
Top10drivers_MPvsIP <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                              sheet = "fig s6")
Top10drivers_MPvsIP_peanut<-filter(Top10drivers_MPvsIP,Species=="Peanut")

## 5.2 correlationship anlaysis ##
Top10drivers_MPvsIP_peanut_activeFe_cor<-corr.test(Top10drivers_MPvsIP_peanut$YLActiveFe,Top10drivers_MPvsIP_peanut[,14:length(Top10drivers_MPvsIP_peanut[1,])],
                              method="spearman",adjust="BH",minlength=5)
Top10drivers_MPvsIP_peanut_activeFe_cor$r
Top10drivers_MPvsIP_peanut_activeFe_cor$p
Top10drivers_MPvsIP_peanut_AvailableFe_cor<-corr.test(Top10drivers_MPvsIP_peanut$AvailableFe,Top10drivers_MPvsIP_peanut[,14:length(Top10drivers_MPvsIP_peanut[1,])],
                                                  method="spearman",adjust="BH",minlength=5)
Top10drivers_MPvsIP_peanut_AvailableFe_cor$r
Top10drivers_MPvsIP_peanut_AvailableFe_cor$p

Top10drivers_MPvsIP_peanut_cor_results<-as.data.frame(t(rbind(Top10drivers_MPvsIP_peanut_activeFe_cor$r,
                                                          Top10drivers_MPvsIP_peanut_activeFe_cor$p,
                                                          Top10drivers_MPvsIP_peanut_AvailableFe_cor$r,
                                                          Top10drivers_MPvsIP_peanut_AvailableFe_cor$p)))
colnames(Top10drivers_MPvsIP_peanut_cor_results)<-c("activeFe_r","activeFe_p_value","avaliableFe_r","availableFe_p_value")
Top10drivers_MPvsIP_peanut_cor_results$Genus<-row.names(Top10drivers_MPvsIP_peanut_cor_results)
Top10drivers_MPvsIP_peanut_cor_r<-as.data.frame(t(rbind(Top10drivers_MPvsIP_peanut_activeFe_cor$r,
                                                      Top10drivers_MPvsIP_peanut_AvailableFe_cor$r)))
colnames(Top10drivers_MPvsIP_peanut_cor_r)<-c("activeFe_r","avaliableFe_r")
### 5.3. Complexheatmap ###
## 5.3.1 color plate ##
max(Top10drivers_MPvsIP_peanut_cor_r)
min(Top10drivers_MPvsIP_peanut_cor_r)
col_cor = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#4d9221", "white", "#c51b7d"))
## 5.3.2 split ##
row_order_setting<-c("Pseudomonas","Aeromicrobium","Solirubrobacter","Iamia","unclassified_f__Nocardioidaceae","Subgroup_7",
                     "Pontibacter","Microbacterium","Noviherbaspirillum","AKYG1722")
row_order_setting
Top10drivers_MPvsIP_peanut_cor_r
length(row_order_setting)

# 5.3.3 plots #
Heatmap(Top10drivers_MPvsIP_peanut_cor_r,col = col_cor,
        cluster_columns = F,cluster_rows = F,row_order = row_order_setting)
setwd(wdOutput_Figure3)

#saving as width*hight 7*4
