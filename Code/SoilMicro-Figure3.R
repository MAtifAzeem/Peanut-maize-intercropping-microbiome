#### 1. loading required libraries ####
library(ggplot2)#作图 plot
library(ggpubr)#添加显著性标记, Add the significance marker
library(ggsignif)#添加显著性标记, Add the significance marker
library(dplyr)#数据清洗，Data cleaning
library(plyr)#数据清洗，Data cleaning
library(ggthemes)#ggplot所用主题，Themes for ggplot2
library(readxl)#读入 excel, read excel
library(ggsci)#配色，color scheme
library(showtext)#字体设置, font setting
library(extrafont)#使用系统字体，Using the system fonts
library(sysfonts)#加载系统字体，loading the system fonts
library(Cairo)#抗锯齿,anti-aliasing
library(ape)
library(phyloseq)#microbiome analysis
library(vegan)#Adonis analysis
library(ComplexHeatmap)#heatmapping
library(cluster)#hierarchical clustering
library(psych)#correlationship analysis#
#### 2. Setting themes and working dictionary path ####
loadfonts()
Sys.setenv(R_GSCMD = "C:/Program Files (x86)/gs/gs9.50/bin/gswin32c.exe")

mytheme1 <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#4d4d4d",linewidth=0.2),
                              panel.background = element_rect(colour = "#4d4d4d",size=0.2),
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
                              axis.ticks.length = unit(0.8, "mm"))#移除整体的边???

mytheme_bigfonts <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#4d4d4d",linewidth=0.2),
                              panel.background = element_rect(colour = "#4d4d4d",size=0.2),
                              text = element_text(family = "Arial"),
                              strip.text = element_text(size = 9,hjust = 0.5),
                              plot.title = element_text(size = 9,hjust = 0.5),
                              axis.text=element_text(size=9,color = "#4D4D4D"),
                              axis.title=element_text(size = 9),
                              legend.text = element_text(size = 9),
                              legend.title = element_text(size = 9),
                              legend.background = element_blank(),
                              panel.border = element_rect(colour = NA),
                              axis.line = element_line(color = "#4D4D4D",linewidth=0.2),
                              axis.ticks.length = unit(0.8, "mm"))#移除整体的边???

wdImport<-("E:/working/SCI/Soil Micro/SCI/Figures/Data/Data for submit")
wdOutput<- ("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Figure3")
data.set.name = '_FourStages_rhizosphere_DAD2_F295_R205_Filtered1' 

#### 3 Import and process data ####
setwd(wdImport)
ASVCount <- read.table("feature_table_FourStages_rhizosphere_DAD2_F295_R205_Filtered1.txt",header=T,row.names = 1,na.strings = c("NA"))
SampleData <- read.table("FourStages_rhizosphere_SampleData.txt",header=T,row.names = 1,na.strings = c("NA"))
colnames(ASVCount)<- row.names(SampleData)
RootedTree_FourStages_rhizosphere_Filtered1<- read.tree("rooted_tree_FourStages_rhizosphere_DAD2_F295_R205_Filtered1.nwk")
Df <- phyloseq(otu_table(ASVCount, taxa_are_rows = T),sample_data(SampleData),phy_tree(RootedTree_FourStages_rhizosphere_Filtered1))
Df
Dfr <-transform_sample_counts(Df, function(x) x / sum(x) )
unifrac_PCoA <- ordinate(Dfr, "PCoA", "unifrac")
unifrac_PCoA_Vector<-as.data.frame(unifrac_PCoA$vectors)
setwd(wdOutput)
write.table(unifrac_PCoA_Vector, paste("unifrac_vector",data.set.name,".txt",sep=""), col.names = NA, sep = '\t', quote = FALSE)
sample_data(Dfr)$RhizocompartmentsSystemSpecies <- factor(sample_data(Dfr)$RhizocompartmentsSystemSpecies,level=c("MP","IP","IM","MM"))
sample_data(Dfr)$Species <- factor(sample_data(Dfr)$Species,level=c("Peanut","Maize"))
sample_data(Dfr)$Days<-factor(sample_data(Dfr)$Days)
Dfr_rhizo<-subset_samples(Dfr,Rhizocompartments=="Rhizosphere")

#### 6. Days ####
### 6.1 46d ###
Dfr_46 <- subset_samples(Dfr_rhizo,Dfr_rhizo@sam_data$Days=="46")
Dfr_46
SampleData_rhizo_46<-filter(SampleData,Days==46)
#adonis
ASVCount_rhizo_46<-ASVCount[,c(1:3,13:15,25:27,37:39)]
setwd(wdImport)
getwd()
tree=RootedTree_FourStages_rhizosphere_Filtered1
amplicon::MicroTest(otu=ASVCount_rhizo_46,map=SampleData_rhizo_46,group="RhizocompartmentsSystemSpecies",
                    Micromet="adonis",dist="unifrac")

#plots#
Dfr_rhizo <- subset_samples(Dfr,Dfr@sam_data$Rhizocompartments=="Rhizosphere")
Dfr_rhizo_46 <- subset_samples(Dfr_rhizo,Dfr_rhizo@sam_data$Days=="46")
Dfr_rhizo_46
sample_data(Dfr_rhizo_46)$SystemSpecies <- factor(sample_data(Dfr_rhizo_46)$SystemSpecies,level=c("MP","IP","IM","MM"))

unifrac_PCoA_rhizo_46 <- ordinate(Dfr_rhizo_46, method="PCoA", distance="unifrac")
unifrac_PCoAPoints_rhizo_46 <- plot_ordination(Dfr_rhizo_46, unifrac_PCoA_rhizo_46,type = "samples", color = "SystemSpecies" )+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point()+
  scale_color_manual(values = c("#4D78B2","#BE1A17","#E19896","#B3CDE2"))+
  scale_x_continuous(limits=c(-0.3,0.3),breaks=c(-0.2,0.0,0.2))+
  scale_y_continuous(limits=c(-0.45,0.45),breaks=c(-0.2,0.0,0.2))+
  geom_polygon(aes(colour = SystemSpecies),fill=NA)+
  mytheme1
unifrac_PCoAPoints_rhizo_46
setwd(wdOutput)
getwd()
ggsave(paste("unifrac_PCoAPoints_rhizo_46",data.set.name,".pdf",sep=""),
       unifrac_PCoAPoints_rhizo_46,
       device=cairo_pdf,width=80,height=60,dpi = 600,units = "mm")

#### 6.2 53d ####
Dfr_53 <- subset_samples(Dfr_rhizo,Dfr_rhizo@sam_data$Days=="53")
Dfr_53
#adonis
SampleData_rhizo_53 <-filter(SampleData,Days==53)
ASVCount_rhizo_53<-ASVCount[,c(4:6,16:18,28:30,40:42)]
tree=RootedTree_FourStages_rhizosphere_Filtered1
amplicon::MicroTest(otu=ASVCount_rhizo_53,map=SampleData_rhizo_53,group="RhizocompartmentsSystemSpecies",
                    Micromet="adonis",dist="unifrac")
#plots#
Dfr_rhizo <- subset_samples(Dfr,Dfr@sam_data$Rhizocompartments=="Rhizosphere")
Dfr_rhizo_53 <- subset_samples(Dfr_rhizo,Dfr_rhizo@sam_data$Days=="53")
Dfr_rhizo_53
sample_data(Dfr_rhizo_53)$SystemSpecies <- factor(sample_data(Dfr_rhizo_53)$SystemSpecies,level=c("MP","IP","IM","MM"))

unifrac_PCoA_rhizo_53 <- ordinate(Dfr_rhizo_53, method="PCoA", distance="unifrac")
unifrac_PCoAPoints_rhizo_53 <- plot_ordination(Dfr_rhizo_53, unifrac_PCoA_rhizo_53,type = "samples", color = "SystemSpecies" )+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point(size=1)+
  scale_color_manual(values = c("#4D78B2","#BE1A17","#E19896","#B3CDE2"))+
  scale_x_continuous(limits=c(-0.3,0.3),breaks=c(-0.2,0.0,0.2))+
  scale_y_continuous(limits=c(-0.45,0.45),breaks=c(-0.2,0.0,0.2))+
  geom_polygon(aes(colour = SystemSpecies),fill=NA)+
  mytheme1
unifrac_PCoAPoints_rhizo_53
setwd(wdOutput)
getwd()
ggsave(paste("unifrac_PCoAPoints_rhizo_53",data.set.name,".pdf",sep=""),
       unifrac_PCoAPoints_rhizo_53,
       device=cairo_pdf,width=80,height=60,dpi = 600,units = "mm")

#### 6.3 63d ####
#adonis
SampleData_rhizo_63 <-filter(SampleData,Days==63)
ASVCount_rhizo_63<-ASVCount[,c(7:9,19:21,31:33,43:45)]
tree=RootedTree_FourStages_rhizosphere_Filtered1
amplicon::MicroTest(otu=ASVCount_rhizo_63,map=SampleData_rhizo_63,group="RhizocompartmentsSystemSpecies",
                    Micromet="adonis",dist="unifrac")

#plots#
Dfr_rhizo <- subset_samples(Dfr,Dfr@sam_data$Rhizocompartments=="Rhizosphere")
Dfr_rhizo_63 <- subset_samples(Dfr_rhizo,Dfr_rhizo@sam_data$Days=="63")
Dfr_rhizo_63
sample_data(Dfr_rhizo_63)$SystemSpecies <- factor(sample_data(Dfr_rhizo_63)$SystemSpecies,level=c("MP","IP","IM","MM"))

unifrac_PCoA_rhizo_63 <- ordinate(Dfr_rhizo_63, method="PCoA", distance="unifrac")
unifrac_PCoAPoints_rhizo_63 <- plot_ordination(Dfr_rhizo_63, unifrac_PCoA_rhizo_63,type = "samples", color = "SystemSpecies" )+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point(size=1)+
  scale_color_manual(values = c("#4D78B2","#BE1A17","#E19896","#B3CDE2"))+
  scale_x_continuous(limits=c(-0.3,0.3),breaks=c(-0.2,0.0,0.2))+
  scale_y_continuous(limits=c(-0.45,0.45),breaks=c(-0.2,0.0,0.2))+
  geom_polygon(aes(colour = SystemSpecies),fill=NA)+
  mytheme1
unifrac_PCoAPoints_rhizo_63
setwd(wdOutput)
getwd()
ggsave(paste("unifrac_PCoAPoints_rhizo_63",data.set.name,".pdf",sep=""),
       unifrac_PCoAPoints_rhizo_63,
       device=cairo_pdf,width=80,height=60,dpi = 600,units = "mm")

#### 6.4 73d ####
#adonis
SampleData_rhizo_73 <-filter(SampleData,Days==73)
ASVCount_rhizo_73<-ASVCount[,c(10:12,22:24,34:36,46:48)]
tree=RootedTree_FourStages_rhizosphere_Filtered1
amplicon::MicroTest(otu=ASVCount_rhizo_73,map=SampleData_rhizo_73,group="RhizocompartmentsSystemSpecies",
                    Micromet="adonis",dist="unifrac")

#plots#
Dfr_rhizo_73 <- subset_samples(Dfr_rhizo,Dfr_rhizo@sam_data$Days=="73")
Dfr_rhizo_73
sample_data(Dfr_rhizo_73)$SystemSpecies <- factor(sample_data(Dfr_rhizo_73)$SystemSpecies,level=c("MP","IP","IM","MM"))

unifrac_PCoA_rhizo_73 <- ordinate(Dfr_rhizo_73, method="PCoA", distance="unifrac")
unifrac_PCoAPoints_rhizo_73 <- plot_ordination(Dfr_rhizo_73, unifrac_PCoA_rhizo_73,type = "samples", color = "SystemSpecies" )+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point(size=1)+
  scale_color_manual(values = c("#4D78B2","#BE1A17","#E19896","#B3CDE2"))+
  scale_x_continuous(limits=c(-0.3,0.3),breaks=c(-0.2,0.0,0.2))+
  scale_y_continuous(limits=c(-0.45,0.45),breaks=c(-0.2,0.0,0.2))+
  geom_polygon(aes(colour = SystemSpecies),fill=NA)+
  mytheme1
unifrac_PCoAPoints_rhizo_73
setwd(wdOutput)
getwd()
ggsave(paste("unifrac_PCoAPoints_rhizo_73",data.set.name,".pdf",sep=""),
       unifrac_PCoAPoints_rhizo_73,
       device=cairo_pdf,width=80,height=60,dpi = 600,units = "mm")

#### 7. Fig. 2b-Heatmap_genus_biomarker_MPvsIP####
### 7.1 Import and process data ###
setwd(wdImport)
getwd()
genus_biomarker <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                              sheet = "fig 3b")
genus_biomarker_mean<-genus_biomarker%>%group_by(SystemSpeciesDays,SystemSpecies,System,Species,Days)%>%
  summarise_at(vars(IMCC26256:Luteolibacter),funs(mean))
genus_biomarker_mean<-genus_biomarker_mean[c(13:16,5:8,9:12,1:4),]
genus_biomarker_RelativeAbundance <- genus_biomarker_mean[,-c(1:5)]
genus_biomarker_RelativeAbundance_MP_IP_MM<-genus_biomarker_RelativeAbundance[c(1:12),]
genus_biomarker_RelativeAbundance_MP_IP_MM_Normalization<-t(scale(genus_biomarker_RelativeAbundance_MP_IP_MM,center = T))
colnames(genus_biomarker_RelativeAbundance_MP_IP_MM_Normalization)<-genus_biomarker_mean$SystemSpeciesDays[c(1:12)]
View(genus_biomarker_RelativeAbundance_MP_IP_MM_Normalization)
### 7.2 Complexheatmap ###
## 7.2.1 color plate ##
max(genus_biomarker_RelativeAbundance_MP_IP_MM_Normalization)
min(genus_biomarker_RelativeAbundance_MP_IP_MM_Normalization)
col_fun = circlize::colorRamp2(c(-2, 0, 3), c("#2166ac", "white", "#b2182b"))
col_fun(seq(-3, 3))
## 7.2.2 split ##
column_split_order<-c(rep("MP",4),rep("IP",4),rep("MM",4))
column_split_order<-factor(column_split_order, levels = c("IP","MP","MM"))
## 7.2.3 annotation ##
genus_biomarker_IP <-filter(genus_biomarker,SystemSpecies=="IP")
genus_biomarker_IP_value <-genus_biomarker_IP[,-c(1:13)]
genus_biomarker_IP_value<-t(genus_biomarker_IP_value)
ha = rowAnnotation(RA_mean = anno_boxplot(genus_biomarker_IP_value,height = unit(4, "cm"),
                                          box_width = 0.5))
ha
## 7.2.4 plots ##
Heatmap(genus_biomarker_RelativeAbundance_MP_IP_MM_Normalization,col = col_fun,
        cluster_columns = F,column_split = column_split_order,
        name="Relative abundance",row_km = 4,row_gap = unit(c(2,4,2), "mm"),
        column_title_gp = gpar(fill = c("#E31A1C", "#1F78B4", "#A6CEE3")),
        right_annotation = ha, row_dend_reorder = TRUE,
        clustering_distance_rows = "spearman")
setwd(wdOutput)

#save as width*height7*4

#### 8. fig 3b Correlationship with active Fe and available Fe in peanuts ####
### 8.1 Import and process data ###
setwd(wdImport)
getwd()
genus_biomarker <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                              sheet = "fig 3b")
genus_biomarker_peanut<-filter(genus_biomarker,Species=="Peanut")

## 8.2 correlationship anlaysis ##
genus_biomarker_peanut_activeFe_cor<-corr.test(genus_biomarker_peanut$YLActiveFe,genus_biomarker_peanut[,14:length(genus_biomarker_peanut[1,])],
                                               method="spearman",adjust="BH",minlength=5)
genus_biomarker_peanut_activeFe_cor$r
genus_biomarker_peanut_activeFe_cor$p
genus_biomarker_peanut_activeFe_cor$p.adj
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

### 8.3. Complexheatmap ###
## 8.3.1 color plate ##
max(genus_biomarker_peanut_cor_r)
min(genus_biomarker_peanut_cor_r)
col_cor = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#4d9221", "white", "#d55e00"))
## 8.3.2 split ##
row_order_setting<-c("Sphingomonas","Massilia","Saccharimonadales","IMCC26256",
                     "uncultured_f__Gemmatimonadaceae","MND1","Ellin6067","MB-A2-108",
                     "Ohtaekwangia","Geitlerinema","unclassified_f__Pseudomonadaceae","Pseudomonas","unclassified_f__Comamonadaceae",
                     "Pseudoxanthomonas","Luteolibacter","Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium","Sphingobium","Devosia")
length(row_order_setting)
row_order_setting

# 8.3.3 plots #
Heatmap(genus_biomarker_peanut_cor_r,col = col_cor,
        cluster_columns = F,cluster_rows = F,
        row_order = row_order_setting
)
setwd(wdOutput)
#saving as width*hight 7*4

