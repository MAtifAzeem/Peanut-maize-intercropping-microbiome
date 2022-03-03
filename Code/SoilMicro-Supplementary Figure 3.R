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
wdOutput <- ("E:/Study/SCI/Soil Micro/SCI/Figures/Figures from R/Supplemental materials/Biomarkers_PhylumToFamily_MPvsIP")
#### 3. Biomarks_PhylumToFamily ####
### 3.1. Import and process data ###
setwd(wdImport)
getwd()
Biomarks_MPvsIP_PhylumToFamilyToFamily <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                           sheet = "fig s3",col_types = c(rep("text",11),rep("numeric",27)))
apply(Biomarks_MPvsIP_PhylumToFamilyToFamily,2,class)
Biomarks_MPvsIP_PhylumToFamily_mean<-Biomarks_MPvsIP_PhylumToFamilyToFamily%>%group_by(SystemSpeciesDays,SystemSpecies,System,Species,Days)%>%
  summarise_at(vars(p__Gemmatimonadota:f__Geitlerinemataceae),funs(mean))
Biomarks_MPvsIP_PhylumToFamily_mean<-Biomarks_MPvsIP_PhylumToFamily_mean[c(13:16,5:8,9:12,1:4),]
Biomarks_MPvsIP_PhylumToFamily_RelativeAbundance <- Biomarks_MPvsIP_PhylumToFamily_mean[,-c(1:5)]
Biomarks_MPvsIP_PhylumToFamily_RelativeAbundance_MP_IP_MM<-Biomarks_MPvsIP_PhylumToFamily_RelativeAbundance[c(1:12),]
Biomarks_MPvsIP_PhylumToFamily_RelativeAbundance_MP_IP_MM_Normalization<-t(scale(Biomarks_MPvsIP_PhylumToFamily_RelativeAbundance_MP_IP_MM,center = T))
colnames(Biomarks_MPvsIP_PhylumToFamily_RelativeAbundance_MP_IP_MM_Normalization)<-Biomarks_MPvsIP_PhylumToFamily_mean$SystemSpeciesDays[c(1:12)]
View(Biomarks_MPvsIP_PhylumToFamily_RelativeAbundance_MP_IP_MM_Normalization)
### 4.2 Complexheatmap ###
## 4.2.1 color plate ##
max(Biomarks_MPvsIP_PhylumToFamily_RelativeAbundance_MP_IP_MM_Normalization)
min(Biomarks_MPvsIP_PhylumToFamily_RelativeAbundance_MP_IP_MM_Normalization)
col_fun = circlize::colorRamp2(c(-2, 0, 3), c("#2166ac", "white", "#b2182b"))
col_fun(seq(-3, 3))
## 4.2.2 split ##
column_split_order<-c(rep("MP",4),rep("IP",4),rep("MM",4))
column_split_order<-factor(column_split_order, levels = c("IP","MP","MM"))
## 4.2.3 annotation ##
Biomarks_MPvsIP_PhylumToFamily_IP <-filter(Biomarks_MPvsIP_PhylumToFamilyToFamily,SystemSpecies=="IP")
Biomarks_MPvsIP_PhylumToFamily_IP_value <-Biomarks_MPvsIP_PhylumToFamily_IP[,-c(1:13)]
Biomarks_MPvsIP_PhylumToFamily_IP_value_t<-t(Biomarks_MPvsIP_PhylumToFamily_IP_value)
ha = rowAnnotation(RA_mean = anno_boxplot(Biomarks_MPvsIP_PhylumToFamily_IP_value_t,height = unit(4, "cm"),
                                          box_width = 0.5))
## 4.2.4 plots ##
Heatmap(Biomarks_MPvsIP_PhylumToFamily_RelativeAbundance_MP_IP_MM_Normalization,col = col_fun,
           cluster_columns = F,column_split = column_split_order,
           name="Relative abundance",row_km = 2,row_gap = unit(c(2), "mm"),
           column_title_gp = gpar(fill = c("#E31A1C", "#1F78B4", "#A6CEE3")),
           right_annotation = ha, row_dend_reorder = TRUE,
           clustering_distance_rows = "spearman")
setwd(wdOutput)

#save as width*height8*5

#### 4. Supplementary fig 3 Correlationship with active Fe and available Fe in peanuts ####
### 4.1 Import and process data ###
setwd(wdImport)
getwd()
Biomarks_MPvsIP_PhylumToFamilyToFamily <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                              sheet = "fig s4")
apply(Biomarks_MPvsIP_PhylumToFamilyToFamily,2,class)
Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut<-filter(Biomarks_MPvsIP_PhylumToFamilyToFamily,Species=="Peanut")
AA<-Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut[,14:length(Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut[1,])]

## 4.2 correlationship anlaysis ##
Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_activeFe_cor<-corr.test(Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut$YLActiveFe,Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut[,14:length(Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut[1,])],
                              method="spearman",adjust="BH",minlength=5)
Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_activeFe_cor$r
Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_activeFe_cor$p
Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_AvailableFe_cor<-corr.test(Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut$AvailableFe,Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut[,14:length(Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut[1,])],
                                                  method="spearman",adjust="BH",minlength=5)
Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_AvailableFe_cor$r
Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_AvailableFe_cor$p

Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_cor_results<-as.data.frame(t(rbind(Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_activeFe_cor$r,
                                                          Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_activeFe_cor$p,
                                                          Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_AvailableFe_cor$r,
                                                          Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_AvailableFe_cor$p)))
colnames(Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_cor_results)<-c("activeFe_r","activeFe_p_value","avaliableFe_r","availableFe_p_value")
Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_cor_results$Genus<-row.names(Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_cor_results)
Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_cor_r<-as.data.frame(t(rbind(Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_activeFe_cor$r,
                                                      Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_AvailableFe_cor$r)))
colnames(Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_cor_r)<-c("activeFe_r","avaliableFe_r")
### 4.3. Complexheatmap ###
## 4.3.1 color plate ##
max(Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_cor_r)
min(Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_cor_r)
col_cor = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#4d9221", "white", "#c51b7d"))
## 4.3.2 split ##
row_order_setting<-c("p__Cyanobacteria","f__Pseudomonadaceae","o__Pseudomonadales","c__Cyanobacteriia",
                     "o__Geitlerinematales","f__Geitlerinemataceae","f__Phormidiaceae","o__Cyanobacteriales",
                     "f__Devosiaceae","f__Rhizobiaceae","f__Comamonadaceae","p__Verrucomicrobiota","c__Verrucomicrobiae","o__Verrucomicrobiales",
                     "f__Rubritaleaceae","f__Saccharimonadales","f__Geodermatophilaceae","o__Frankiales","o__IMCC26256",
                     "f__IMCC26256","f__Gemmatimonadaceae","c__Gemmatimonadetes","o__Gemmatimonadales",
                     "p__Gemmatimonadota","f__Hymenobacteraceae")
row_order_setting
Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_cor_r
length(row_order_setting)

# 4.3.3 plots #
Heatmap(Biomarks_MPvsIP_PhylumToFamilyToFamily_peanut_cor_r,col = col_cor,
        cluster_columns = F,cluster_rows = F,row_order = row_order_setting)
setwd(wdOutput_Figure3)

#saving as width*hight 7*4
