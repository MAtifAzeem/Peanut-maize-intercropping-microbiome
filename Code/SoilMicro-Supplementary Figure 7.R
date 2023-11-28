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
wdOutput <- ("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Supplemental materials/Biomarkers_MMvsIM")

#### 3. Fig. 3a-LefSe_MMvsIM_results####
### 3.1 Import and process data ###
setwd(wdImport)
LefSe_MMvsIM <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                    sheet = "fig s7a")
LefSe_MMvsIM$ID<-factor(LefSe_MMvsIM$ID,levels=LefSe_MMvsIM$ID)
### 3.2 Plots ###
LefSe_MMvsIM_LDA <- ggplot(LefSe_MMvsIM,aes(x=ID, y=LDA, fill=group)) + 
  geom_bar(stat="identity", position="identity",width = 0.8,color="black",size=0.1)+
  scale_fill_manual(values=c("#E19896","#A6CEE3"))+
  scale_y_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3,4))+
  mytheme1+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
LefSe_MMvsIM_LDA
setwd(wdOutput)
getwd()
ggsave(paste("LefSe_MMvsIM_LDA_3.0_p0.05",".pdf",sep=""),
       LefSe_MMvsIM_LDA,device=cairo_pdf,width=140,height=90,dpi = 600,units = "mm")

#### 4. Fig. 3b-Heatmap_genus_biomarker_MMvsIM####
### 4.1 Import and process data ###
setwd(wdImport)
getwd()
genus_biomarker <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                              sheet = "fig s7b")
genus_biomarker_mean<-genus_biomarker%>%group_by(SystemSpeciesDays,SystemSpecies,System,Species,Days)%>%
  summarise_at(vars(p__Firmicutes:g__Bryobacter),funs(mean))
genus_biomarker_mean<-genus_biomarker_mean[c(9:12,1:8,13:16),]
genus_biomarker_RelativeAbundance <- genus_biomarker_mean[,-c(1:5)]
genus_biomarker_RelativeAbundance_MM_IM_MP<-genus_biomarker_RelativeAbundance[c(1:12),]
genus_biomarker_RelativeAbundance_MM_IM_MP_Normalization<-t(scale(genus_biomarker_RelativeAbundance_MM_IM_MP,center = T))
genus_biomarker_RelativeAbundance_MM_IM_MP_Normalization
colnames(genus_biomarker_RelativeAbundance_MM_IM_MP_Normalization)<-genus_biomarker_mean$SystemSpeciesDays[c(1:12)]
View(genus_biomarker_RelativeAbundance_MM_IM_MP_Normalization)
### 4.2 Complexheatmap ###
## 4.2.1 color plate ##
max(genus_biomarker_RelativeAbundance_MM_IM_MP_Normalization)
min(genus_biomarker_RelativeAbundance_MM_IM_MP_Normalization)
col_fun = circlize::colorRamp2(c(-2, 0, 3), c("#2166ac", "white", "#b2182b"))
col_fun(seq(-3, 3))
## 4.2.2 split ##
column_split_order<-c(rep("MM",4),rep("IM",4),rep("MP",4))
column_split_order<-factor(column_split_order, levels = c("IM","MM","MP"))
## 4.2.3 annotation ##
genus_biomarker_IM <-filter(genus_biomarker,SystemSpecies=="IM")
genus_biomarker_IM_value <-genus_biomarker_IM[,-c(1:13)]
genus_biomarker_IM_value<-t(genus_biomarker_IM_value)
ha = rowAnnotation(RA_mean = anno_boxplot(genus_biomarker_IM_value,height = unit(4, "cm"),
                                          box_width = 0.5))
ha
View(genus_biomarker_IM_value)
## 4.2.4 plots ##
Heatmap(genus_biomarker_RelativeAbundance_MM_IM_MP_Normalization,col = col_fun,
           cluster_columns = F,column_split = column_split_order,
           name="Relative abundance",row_km = 4,row_gap = unit(c(2,4,2), "mm"),
           column_title_gp = gpar(fill = c("#CC9694", "#BBCCE1", "#5F78B0")), row_dend_reorder = TRUE,
           clustering_method_rows = "average",right_annotation = ha)
setwd(wdOutput)

#save as height*width8*6
