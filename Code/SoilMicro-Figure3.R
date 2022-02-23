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
wdOutput_Figure3 <- ("E:/Study/SCI/Soil Micro/SCI/Figures/Figures from R/Figure3")

#### 3. Fig. 3a-LefSe_MPvsIP_results####
### 3.1 Import and process data ###
setwd(wdImport)
LefSe_MPvsIP <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                    sheet = "fig 3a")
LefSe_MPvsIP$ID<-factor(LefSe_MPvsIP$ID,levels=LefSe_MPvsIP$ID)
### 3.2 Plots ###
LefSe_MPvsIP_LDA <- ggplot(LefSe_MPvsIP,aes(x=ID, y=LDA, fill=group)) + 
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

#### 4. Fig. 3b-Heatmap_genus_biomarker_MPvsIP####
### 4.1 Import and process data ###
setwd(wdImport)
getwd()
genus_biomarker <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                           sheet = "fig 3b")
genus_biomarker_mean<-genus_biomarker%>%group_by(SystemSpeciesDays,SystemSpecies,System,Species,Days)%>%
  summarise_at(vars(unclassified_f__Comamonadaceae:IMCC26256),funs(mean))
genus_biomarker_mean<-genus_biomarker_mean[c(13:16,5:8,9:12,1:4),]
genus_biomarker_RelativeAbundance <- genus_biomarker_mean[,-c(1:5)]
genus_biomarker_RelativeAbundance_MP_IP_MM<-genus_biomarker_RelativeAbundance[c(1:12),]
genus_biomarker_RelativeAbundance_MP_IP_MM_Normalization<-t(scale(genus_biomarker_RelativeAbundance_MP_IP_MM,center = T))
colnames(genus_biomarker_RelativeAbundance_MP_IP_MM_Normalization)<-genus_biomarker_mean$SystemSpeciesDays[c(1:12)]
View(genus_biomarker_RelativeAbundance_MP_IP_MM_Normalization)
### 4.2 Complexheatmap ###
## 4.2.1 color plate ##
max(genus_biomarker_RelativeAbundance_MP_IP_MM_Normalization)
min(genus_biomarker_RelativeAbundance_MP_IP_MM_Normalization)
col_fun = circlize::colorRamp2(c(-2, 0, 3), c("#2166ac", "white", "#b2182b"))
col_fun(seq(-3, 3))
## 4.2.2 split ##
column_split_order<-c(rep("MP",4),rep("IP",4),rep("MM",4))
column_split_order<-factor(column_split_order, levels = c("IP","MP","MM"))
## 4.2.3 annotation ##
genus_biomarker_IP <-filter(genus_biomarker,SystemSpecies=="IP")
genus_biomarker_IP_value <-genus_biomarker_IP[,-c(1:13)]
genus_biomarker_IP_value<-t(genus_biomarker_IP_value)
ha = rowAnnotation(RA_mean = anno_boxplot(genus_biomarker_IP_value,height = unit(4, "cm"),
                                          box_width = 0.5))
ha
## 4.2.4 plots ##
Heatmap(genus_biomarker_RelativeAbundance_MP_IP_MM_Normalization,col = col_fun,
           cluster_columns = F,column_split = column_split_order,
           name="Relative abundance",row_km = 4,row_gap = unit(c(2,4,2), "mm"),
           column_title_gp = gpar(fill = c("#E31A1C", "#1F78B4", "#A6CEE3")),
           right_annotation = ha, row_dend_reorder = TRUE,
           clustering_distance_rows = "spearman")
setwd(wdOutput_Figure3)

#save as width*height7*4

#### 5. fig 3b Correlationship with active Fe and available Fe in peanuts ####
### 5.1 Import and process data ###
setwd(wdImport)
getwd()
genus_biomarker <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                              sheet = "fig 3b")
genus_biomarker_peanut<-filter(genus_biomarker,Species=="Peanut")

## 5.2 correlationship anlaysis ##
genus_biomarker_peanut_activeFe_cor<-corr.test(genus_biomarker_peanut$YLActiveFe,genus_biomarker_peanut[,14:length(genus_biomarker_peanut[1,])],
                              method="spearman",adjust="BH",minlength=5)
genus_biomarker_peanut_activeFe_cor$r
genus_biomarker_peanut_activeFe_cor$p
genus_biomarker_peanut_AvailableFe_cor<-corr.test(genus_biomarker_peanut$AvailableFe,genus_biomarker_peanut[,14:length(genus_biomarker_peanut[1,])],
                                                  method="spearman",adjust="BH",minlength=5)
genus_biomarker_peanut_AvailableFe_cor$r
genus_biomarker_peanut_AvailableFe_cor$p

genus_biomarker_peanut_cor_results<-as.data.frame(t(rbind(genus_biomarker_peanut_activeFe_cor$r,
                                                          genus_biomarker_peanut_activeFe_cor$p,
                                                          genus_biomarker_peanut_AvailableFe_cor$r,
                                                          genus_biomarker_peanut_AvailableFe_cor$p)))
colnames(genus_biomarker_peanut_cor_results)<-c("activeFe_r","activeFe_p_value","avaliableFe_r","availableFe_p_value")
genus_biomarker_peanut_cor_results$Genus<-row.names(genus_biomarker_peanut_cor_results)
genus_biomarker_peanut_cor_r<-as.data.frame(t(rbind(genus_biomarker_peanut_activeFe_cor$r,
                                                      genus_biomarker_peanut_AvailableFe_cor$r)))
colnames(genus_biomarker_peanut_cor_r)<-c("activeFe_r","avaliableFe_r")
### 5.3. Complexheatmap ###
## 5.3.1 color plate ##
max(genus_biomarker_peanut_cor_r)
min(genus_biomarker_peanut_cor_r)
col_cor = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#4d9221", "white", "#c51b7d"))
## 5.3.2 split ##

row_order_setting<-c("Massilia","Sphingomonas","Saccharimonadales","IMCC26256","Ellin6067","uncultured_f__Gemmatimonadaceae",
                     "Ohtaekwangia","unclassified_f__Pseudomonadaceae","Pseudomonas","Geitlerinema","unclassified_f__Comamonadaceae",
                     "Allorhizobium_Neorhizobium_Pararhizobium_Rhizobium","Sphingobium","Devosia","Pseudoxanthomonas","Luteolibacter")

length(row_order_setting)

# 5.3.3 plots #
Heatmap(genus_biomarker_peanut_cor_r,col = col_cor,
        cluster_columns = F,cluster_rows = F,
        row_order = row_order_setting
        )
setwd(wdOutput_Figure3)

#saving as width*hight 7*4

#### 6. Fig 3d ####
### 6.1 Import and process data ###
setwd(wdImport)
getwd()
NetShift <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                              sheet = "fig 3d")
NetShift_IP<-filter(NetShift,DelBet>0)
NetShift_IP<-NetShift_IP[order(NetShift_IP$NESH_score,decreasing = T),]
NetShift_IP_top10<-NetShift_IP[c(1:10),]
NetShift_IP_top10<-NetShift_IP_top10[order(NetShift_IP_top10$NESH_score),]
NetShift_IP_top10$Genus<-factor(NetShift_IP_top10$Genus,levels=NetShift_IP_top10$Genus)
### 6.2 plots ###
NetShift_IP_top10_lollipop<-ggplot(NetShift_IP_top10, aes(Genus, NESH_score)) +
  geom_segment( aes(x=Genus, xend=Genus, y=0, yend=NESH_score), color="black",size=0.3) +
  geom_point( color="black", size=1) +
  mytheme1+
  coord_flip()
NetShift_IP_top10_lollipop
setwd(wdOutput_Figure3)
getwd()
ggsave(paste("NetShift_IP_top10_lollipop",".pdf",sep=""),
       NetShift_IP_top10_lollipop,device=cairo_pdf,width=60,height=35,dpi = 600,units = "mm")

#### 7. Fig. 3e####
### 6.1 Import and process data ###
setwd(wdImport)
getwd()
Venn <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                       sheet = "fig 3e")
Venn_data<-list(A=na.omit(Venn$match),B=na.omit(Venn$Cor_active),
                C=na.omit(Venn$Cor_available),D=na.omit(Venn$Keystone))
Venn_data

if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
Venn_biomarker_IPvsMP<-ggvenn(
  Venn_data, 
  fill_color = c("#CD534CFF", "#EFC000FF", "#868686FF", "#0073C2FF"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F,
)
Venn_biomarker_IPvsMP
setwd(wdOutput_Figure3)
getwd()
ggsave(paste("Venn_biomarker_IPvsMP_FourStages_filtered1",".pdf",sep=""),
       Venn_biomarker_IPvsMP,device=cairo_pdf,width=60,height=30,dpi = 600,units = "mm")



