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
wdOutput <- ("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Supplemental materials/divers_IPvsMP")

#### 3 Import and process data ####
setwd(wdImport)
getwd()
NetShift <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                             sheet = "fig s8a")

drivers_IPvsMP <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                           sheet = "fig s8b")

#### 4 ####
### 4.1 NetShift ###
### 4.1.1 processing data ####
NetShift<-NetShift[order(-NetShift$DelBet,-NetShift$NESH_score),]
NetShift_IP<-filter(NetShift,DelBet>0)
NetShift_IP<-NetShift_IP[order(NetShift_IP$NESH_score,decreasing = T),]
NetShift_IP$Genus<-factor(NetShift_IP$Genus,levels=NetShift_IP$Genus)
write.table(NetShift_IP, paste("NetShift_IP",".csv",sep=""), row.names=T,sep = '\t', quote = FALSE)
### 4.1.1 plots ###
NetShift_IP_lollipop<-ggplot(NetShift_IP, aes(Genus, NESH_score)) +
  geom_segment( aes(x=Genus, xend=Genus, y=0, yend=NESH_score), color="black",size=0.3) +
  geom_point( color="black", size=1)+
  mytheme1+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
NetShift_IP_lollipop
setwd(wdOutput)
getwd()
ggsave(paste("NetShift_IP_lollipop",".pdf",sep=""),
       NetShift_IP_lollipop,device=cairo_pdf,width=110,height=60,dpi = 600,units = "mm")


drivers_IPvsMP_mean<-drivers_IPvsMP%>%group_by(SystemSpeciesDays,SystemSpecies,System,Species,Days)%>%
  summarise_at(vars(unclassified_f__Nocardioidaceae:uncultured_p__Armatimonadota),funs(mean))
drivers_IPvsMP_mean<-drivers_IPvsMP_mean[c(13:16,5:8,9:12,1:4),]
drivers_IPvsMP_RelativeAbundance <- drivers_IPvsMP_mean[,-c(1:5)]
drivers_IPvsMP_RelativeAbundance_MP_IP_MM<-drivers_IPvsMP_RelativeAbundance[c(1:12),]
drivers_IPvsMP_RelativeAbundance_MP_IP_MM_Normalization<-t(scale(drivers_IPvsMP_RelativeAbundance_MP_IP_MM,center = T))
colnames(drivers_IPvsMP_RelativeAbundance_MP_IP_MM_Normalization)<-drivers_IPvsMP_mean$SystemSpeciesDays[c(1:12)]
View(drivers_IPvsMP_RelativeAbundance_MP_IP_MM_Normalization)
### 4.2 Complexheatmap ###
## 4.2.1 color plate ##
max(drivers_IPvsMP_RelativeAbundance_MP_IP_MM_Normalization)
min(drivers_IPvsMP_RelativeAbundance_MP_IP_MM_Normalization)
col_fun = circlize::colorRamp2(c(-2, 0, 3), c("#2166ac", "white", "#b2182b"))
col_fun(seq(-3, 3))
## 4.2.2 split ##
column_split_order<-c(rep("MP",4),rep("IP",4),rep("MM",4))
column_split_order<-factor(column_split_order, levels = c("IP","MP","MM"))
## 4.2.3 annotation ##
drivers_IPvsMP_IP <-filter(drivers_IPvsMP,SystemSpecies=="IP")
drivers_IPvsMP_IP_value <-drivers_IPvsMP_IP[,-c(1:13)]
drivers_IPvsMP_IP_value<-t(drivers_IPvsMP_IP_value)
ha = rowAnnotation(RA_mean = anno_boxplot(drivers_IPvsMP_IP_value,height = unit(4, "cm"),
                                          box_width = 0.5))
ha
## 4.2.4 plots ##
Heatmap(drivers_IPvsMP_RelativeAbundance_MP_IP_MM_Normalization,col = col_fun,
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
drivers_IPvsMP <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                              sheet = "fig s6b")
drivers_IPvsMP_peanut<-filter(drivers_IPvsMP,Species=="Peanut")

## 5.2 correlationship anlaysis ##
drivers_IPvsMP_peanut_activeFe_cor<-corr.test(drivers_IPvsMP_peanut$YLActiveFe,drivers_IPvsMP_peanut[,14:length(drivers_IPvsMP_peanut[1,])],
                              method="spearman",adjust="BH",minlength=5)
drivers_IPvsMP_peanut_activeFe_cor$r
drivers_IPvsMP_peanut_activeFe_cor$p
drivers_IPvsMP_peanut_AvailableFe_cor<-corr.test(drivers_IPvsMP_peanut$AvailableFe,drivers_IPvsMP_peanut[,14:length(drivers_IPvsMP_peanut[1,])],
                                                  method="spearman",adjust="BH",minlength=5)
drivers_IPvsMP_peanut_AvailableFe_cor$r
drivers_IPvsMP_peanut_AvailableFe_cor$p

drivers_IPvsMP_peanut_cor_results<-as.data.frame(t(rbind(drivers_IPvsMP_peanut_activeFe_cor$r,
                                                          drivers_IPvsMP_peanut_activeFe_cor$p,
                                                          drivers_IPvsMP_peanut_AvailableFe_cor$r,
                                                          drivers_IPvsMP_peanut_AvailableFe_cor$p)))
colnames(drivers_IPvsMP_peanut_cor_results)<-c("activeFe_r","activeFe_p_value","avaliableFe_r","availableFe_p_value")
drivers_IPvsMP_peanut_cor_results$Genus<-row.names(drivers_IPvsMP_peanut_cor_results)
drivers_IPvsMP_peanut_cor_r<-as.data.frame(t(rbind(drivers_IPvsMP_peanut_activeFe_cor$r,
                                                      drivers_IPvsMP_peanut_AvailableFe_cor$r)))
colnames(drivers_IPvsMP_peanut_cor_r)<-c("activeFe_r","avaliableFe_r")
### 5.3. Complexheatmap ###
## 5.3.1 color plate ##
max(drivers_IPvsMP_peanut_cor_r)
min(drivers_IPvsMP_peanut_cor_r)
col_cor = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#4d9221", "white", "#c51b7d"))
## 5.3.2 split ##
row_order_setting<-c("Adhaeribacter","Geodermatophilus","67-14","uncultured_o__Ardenticatenales","AKYG1722","Sphingomonas","Nocardioides","Aeromicrobium",
                     "uncultured_f__Blastocatellaceae","unclassified_f__Blastocatellaceae","Luteitalea","uncultured_p__Armatimonadota","Iamia","IMCC26256",
                     "RB41","Subgroup_7","unclassified_f__Nocardioidaceae","Microbacterium","unclassified_o__Gaiellales","CCD24","uncultured_f__Gemmatimonadaceae",
                     "bacteriap25","MND1","Noviherbaspirillum","unclassified_f__Xanthomonadaceae","unclassified_f__Planococcaceae","Paenibacillus","Paenisporosarcina",
                     "Pseudomonas")
row_order_setting
drivers_IPvsMP_peanut_cor_r
length(row_order_setting)

# 5.3.3 plots #
Heatmap(drivers_IPvsMP_peanut_cor_r,col = col_cor,
        cluster_columns = F,cluster_rows = F,row_order = row_order_setting)
setwd(wdOutput)


#saving as width*hight 7*4
