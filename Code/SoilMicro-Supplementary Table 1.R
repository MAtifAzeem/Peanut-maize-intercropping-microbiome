#### 1.Loading packet####
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
library(stringr)#字符串处理.string manipulation
library(graphics)#坐标轴表达式，expression for axis
library(vegan)
library(data.table)

####2 set theme and filepath ####
loadfonts()
Sys.setenv(R_GSCMD = "C:/Program Files (x86)/gs/gs9.50/bin/gswin32c.exe")

mytheme <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#000000"),
                             text = element_text(family = "Arial"),
                             strip.text = element_text(size=8,hjust = 0.5),
                             plot.title = element_text(size=8,hjust = 0.5),
                             axis.text=element_text(size=8,color = "black"),
                             axis.title=element_text(size=8),
                             legend.text = element_text(size=8),
                             legend.title = element_text(size=8),
                             legend.background = element_blank(),
                             panel.border = element_rect(colour = NA),
                             axis.line = element_line(color = "black",size=0.4))#移除整体的边???

FacetTheme <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#000000"),
                                text = element_text(family = "Arial"),
                                strip.text = element_text(size=8,hjust = 0.5),
                                plot.title = element_text(size=8,hjust = 0.5),
                                axis.text =element_text(size=8,color = "black"),
                                axis.title =element_text(size=8,color = "black"),
                                legend.text = element_text(size=8,color = "black"),
                                legend.title = element_text(size=8,color = "black"),
                                legend.background = element_blank(),
                                axis.line = element_line(color = "black",size=0.4))#移除整体的边???

wdImport<- c("E:/working/SCI/Soil Micro/SCI/Figures/Data/Data for submit")
wdOutput <- c("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Supplemental materials")

#### 3. Loading and progressing data####
#### 3.1 Genus Iron correlation####
setwd(wdImport)
Iron_All_Genus_Peanut <- read_excel("./Intercropping-microbiome-Data for submit.xlsx",
                       sheet = "table s1")
Iron_All_Genus_Peanut$Species <- factor(Iron_All_Genus_Peanut$Species,levels = c("Peanut","Maize"))
Iron_All_Genus_Peanut$System <- factor(Iron_All_Genus_Peanut$System,levels = c("Monocropping","Intercropping"))
Iron_All_Genus_Peanut$SystemSpecies <- factor(Iron_All_Genus_Peanut$SystemSpecies,levels = c("MP","IP","IM","MM"))
library(psych)
ActiveFe_Genus_cor<-corr.test(Iron_All_Genus_Peanut[8],Iron_All_Genus_Peanut[10:655],
                                      method="spearman",adjust="BH",minlength=5)
ActiveFe_Genus_cor
ActiveFe_Genus_cor$p
ActiveFe_Genus_cor$p.adj
ActiveFe_Genus_cor_Tab<-as.data.frame(t(rbind(ActiveFe_Genus_cor$r,ActiveFe_Genus_cor$p,ActiveFe_Genus_cor$p.adj)))
colnames(ActiveFe_Genus_cor_Tab)<- c("Spearman_ρ_YL_Active_Fe","p-value_YL_Active_Fe","adjust_p_value_YL_Active_Fe")
setwd(wdOutput)
write.table(ActiveFe_Genus_cor_Tab, paste("ActiveFe_Genus_cor_Tab",".txt",sep=""), row.names=F,sep = '\t', quote = FALSE)


AvailableFe_Genus_cor<-corr.test(Iron_All_Genus_Peanut[9],Iron_All_Genus_Peanut[10:655],
                              method="spearman",adjust="BH",minlength=5)
AvailableFe_Genus_cor$r
AvailableFe_Genus_cor$p
AvailableFe_Genus_cor$p.adj
AvailableFe_Genus_cor_Tab<-as.data.frame(t(rbind(AvailableFe_Genus_cor$r,AvailableFe_Genus_cor$p,AvailableFe_Genus_cor$p.adj)))
colnames(AvailableFe_Genus_cor_Tab)<- c("Spearman_ρ_Available_Fe","p-value_Available_Fe","adjust_p_value_Available_Fe")
setwd(wdOutput)
write.table(AvailableFe_Genus_cor_Tab, paste("AvailableFe_Genus_cor_Tab",".txt",sep=""), row.names=F,sep = '\t', quote = FALSE)

#### 3.2 Genus Iron regression####
Pseudomonas_YlActiveFe_reg <- lm(Pseudomonas~YLActiveFe,data=Iron_All_Genus_Peanut)
Pseudomonas_YlActiveFe_reg
Pseudomonas_YlActiveFe_reg_sum<-summary(Pseudomonas_YlActiveFe_reg)
Pseudomonas_YlActiveFe_reg_sum
Pseudomonas_YlActiveFe_reg_sum$coefficients[2,4]
class(Pseudomonas_YlActiveFe_reg_sum$coefficients)

Pseudomonas_AvailableFe_reg <- lm(Pseudomonas~AvailableFe,data=Iron_All_Genus_Peanut)
Pseudomonas_AvailableFe_reg

### 3.2.1 wide to long format
Iron_All_Genus_Peanut_melt <- melt(Iron_All_Genus_Peanut,id.vars =c("ID"), 
                                   measure.vars=colnames(Iron_All_Genus_Peanut)[10:655],
                                   variable.name="genus",
                                   value.name = "relative_abundance")
Iron_All_Genus_Peanut_melt$YLActiveFe<-rep(Iron_All_Genus_Peanut$YLActiveFe,length(colnames(Iron_All_Genus_Peanut)[10:655]))

Iron_All_Genus_Peanut_melt$AvailableFe<-rep(Iron_All_Genus_Peanut$AvailableFe,length(colnames(Iron_All_Genus_Peanut)[10:655]))

### 3.2.2 split the dataframe based on genus name
Iron_All_Genus_Peanut_melt_list <- split(Iron_All_Genus_Peanut_melt,Iron_All_Genus_Peanut_melt$genus)
View(Iron_All_Genus_Peanut_melt_list)
## 3.2.3 YLActiveFe ##
res_YLActiveFe <- lapply(Iron_All_Genus_Peanut_melt_list, function(x) lm(x$relative_abundance~x$YLActiveFe))
res_YLActiveFe_summary <- lapply(res_YLActiveFe, function(x) summary(x))
res_YLActiveFe_summary_R2<-lapply(res_YLActiveFe_summary, function(x) data.frame("lm.R2"=x$r.squared,"lm.adj.R2"=x$adj.r.squared,"p.value"=x$coefficients[2,4]))
res_YLActiveFe_summary_R2_detail <- do.call(rbind,res_YLActiveFe_summary_R2)
coef_res_YLActiveFe <- lapply(res_YLActiveFe, coef)
coef_res_YLActiveFe_detail <- do.call(rbind,coef_res_YLActiveFe)
res_YLActiveFe_final_detail <- cbind(res_YLActiveFe_summary_R2_detail,coef_res_YLActiveFe_detail)
colnames(res_YLActiveFe_final_detail) <- c("lm.R2_YL_Active_Fe","lm.adj.R2_YL_Active_Fe","lm.p.value_YL_Active_Fe","intercept_YL_Active_Fe","slope_YL_Active_Fe")
View(res_YLActiveFe_final_detail)
## 3.2.4 AvailableFe ##
res_AvailableFe <- lapply(Iron_All_Genus_Peanut_melt_list, function(x) lm(x$relative_abundance~x$AvailableFe))
res_AvailableFe_summary <- lapply(res_AvailableFe, function(x) summary(x))
res_AvailableFe_summary_R2<-lapply(res_AvailableFe_summary, function(x) data.frame("lm.R2"=x$r.squared,"lm.adj.R2"=x$adj.r.squared,"p.value"=x$coefficients[2,4]))
res_AvailableFe_summary_R2_detail <- do.call(rbind,res_AvailableFe_summary_R2)
coef_res_AvailableFe <- lapply(res_AvailableFe, coef)
coef_res_AvailableFe_detail <- do.call(rbind,coef_res_AvailableFe)
res_AvailableFe_final_detail <- cbind(res_AvailableFe_summary_R2_detail,coef_res_AvailableFe_detail)
colnames(res_AvailableFe_final_detail) <- c("lm.R2_Available_Fe","lm.adj.R2_Available_Fe","lm.p.value_Available_Fe","intercept_Available_Fe","slope_Available_Fe")
View(res_AvailableFe_final_detail)
FourStages_Genus_cor_lm_ActiveFe_AvailableFe<-cbind(ActiveFe_Genus_cor_Tab,AvailableFe_Genus_cor_Tab,res_YLActiveFe_final_detail,res_AvailableFe_final_detail)
setwd(wdOutput)
getwd()
write.table(FourStages_Genus_cor_lm_ActiveFe_AvailableFe, paste("FourStages_Genus_cor_lm_ActiveFe_AvailableFe",".csv",sep=""), row.names=T,sep = '\t', quote = FALSE)
