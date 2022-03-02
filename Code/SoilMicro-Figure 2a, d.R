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

mytheme1 <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#4d4d4d",size=0.1),
                              panel.background = element_rect(colour = "#4d4d4d",size=0.1),
                              text = element_text(family = "Arial"),
                              strip.text = element_text(size = 6,hjust = 0.5),
                              plot.title = element_text(size = 6,hjust = 0.5),
                              axis.text=element_text(size=6,color = "#4D4D4D"),
                              axis.title=element_text(size = 6),
                              legend.text = element_text(size = 6),
                              legend.title = element_text(size = 6),
                              legend.background = element_blank(),
                              axis.line = element_line(color = "#4D4D4D",size=0.1),
                              axis.ticks.length = unit(0.8, "mm"))#移除整体的边???

wdImport<-("E:/Study/SCI/Soil Micro/SCI/Figures/Data/Data for submit")
wdOutput_Figure2 <- ("E:/Study/SCI/Soil Micro/SCI/Figures/Figures from R/Figure2")

#### 3. Fig. 2a-αdiversity####
### 3.1 Import and process data ###
setwd(wdImport)
αDiversity <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                    sheet = "fig 2a")

sapply(αDiversity,class)
αDiversityRhizosphere<-filter(αDiversity,Rhizocompartments=="Rhizosphere")
### 3.2 Summary by grouping ###
αDiversityRhizosphere_summary <- αDiversityRhizosphere%>%
  group_by(RhizocompartmentsSystemSpeciesDays,System,Days,Species,SystemSpecies,SpeciesDays,SystemSpeciesDays)%>%
  summarise_at(c("chao1","shannon_entropy"),funs(mean,sd))
αDiversityRhizosphere_summary_peanut<-filter(αDiversityRhizosphere_summary,Species=="Peanut")
αDiversityRhizosphere_summary_maize<-filter(αDiversityRhizosphere_summary,Species=="Maize")
### 3.3 Statistic analysis ###
## 3.3.1 Peanut-by days-Statistic analysis ###
αDiversityRhizosphere_peanut <- filter(αDiversityRhizosphere,Species=="Peanut")
View(αDiversityRhizosphere_peanut)
#chao1
chao1_peanut_stat_byDays<-compare_means(data=αDiversityRhizosphere_peanut, chao1~SystemSpeciesDays,method ="t.test")
View(chao1_peanut_stat_byDays)
#shannon
shannon_peanut_stat_byDays<-compare_means(data=αDiversityRhizosphere_peanut, shannon_entropy~SystemSpeciesDays,method ="t.test")
View(shannon_peanut_stat_byDays)
## 3.3.2 maize-by days-Statistic analysis ##
αDiversityRhizosphere_maize <- filter(αDiversityRhizosphere,Species=="Maize")
#chao1
chao1_maize_stat_byDays<-compare_means(data=αDiversityRhizosphere_maize, chao1~SystemSpeciesDays,method ="t.test")
View(chao1_maize_stat_byDays)
#shannon
shannon_maize_stat_byDays<-compare_means(data=αDiversityRhizosphere_maize, shannon_entropy~SystemSpeciesDays,method ="t.test")
View(shannon_maize_stat_byDays)

### 3.4 Plots ###
## 3.4.1chao1_peanut ##
View(αDiversityRhizosphere_peanut)
chao1_peanut_plot <- ggplot(αDiversityRhizosphere_peanut, aes(Days,chao1,color=System,group=System))+
  geom_point(position = position_dodge(8),shape=1,size=0.3)+
  stat_summary(fun.data = mean_sdl,geom='errorbar',position = position_dodge(3),size=0.3,width=2)+
  stat_summary(fun=mean, geom='point',position = position_dodge(3),size=0.3)+
  scale_x_continuous(breaks = c(46,53,63,73,80))+
  scale_y_continuous(limits = c(680,3200))+
  scale_color_npg()+scale_fill_npg()+
  mytheme1
chao1_peanut_plot
setwd(wdOutput_Figure2)
getwd()
ggsave(paste("chao1_peanut_plot",".pdf",sep=""),
       chao1_peanut_plot,device=cairo_pdf,width=63,height=32,dpi = 600,units = "mm")

## 3.4.2shannon_peanut ##
shannon_peanut_plot <- ggplot(αDiversityRhizosphere_peanut, aes(Days,shannon_entropy,color=System,group=System))+
  geom_point(position = position_dodge(8),shape=1,size=0.3)+
  stat_summary(fun.data = mean_sdl,geom='errorbar',position = position_dodge(3),size=0.3,width=2)+
  stat_summary(fun=mean, geom='point',position = position_dodge(3),size=0.3)+
  scale_x_continuous(breaks = c(46,53,63,73,80))+
  scale_y_continuous(limits = c(7.8,11))+
  scale_color_npg()+scale_fill_npg()+
  mytheme1
shannon_peanut_plot
setwd(wdOutput_Figure2)
getwd()
ggsave(paste("shannon_peanut_plot",".pdf",sep=""),
       shannon_peanut_plot,device=cairo_pdf,width=63,height=32,dpi = 600,units = "mm")

## 3.4.3chao1_maize ##
chao1_maize_plot <- ggplot(αDiversityRhizosphere_maize, aes(Days,chao1,color=System,group=System))+
  geom_point(position = position_dodge(8),shape=1,size=0.3)+
  stat_summary(fun.data = mean_sdl,geom='errorbar',position = position_dodge(3),size=0.3,width=2)+
  stat_summary(fun=mean, geom='point',position = position_dodge(3),size=0.3)+
  scale_x_continuous(breaks = c(46,53,63,73,80))+
  scale_y_continuous(limits = c(680,3200))+
  scale_color_npg()+scale_fill_npg()+
  mytheme1
chao1_maize_plot
setwd(wdOutput_Figure2)
getwd()
ggsave(paste("chao1_maize_plot",".pdf",sep=""),
       chao1_maize_plot,device=cairo_pdf,width=63,height=32,dpi = 600,units = "mm")

## 3.4.4shannon_maize ##
shannon_maize_plot <- ggplot(αDiversityRhizosphere_maize, aes(Days,shannon_entropy,color=System,group=System))+
  geom_point(position = position_dodge(8),shape=1,size=0.3)+
  stat_summary(fun.data = mean_sdl,geom='errorbar',position = position_dodge(3),size=0.3,width=2)+
  stat_summary(fun=mean, geom='point',position = position_dodge(3),size=0.3)+
  scale_x_continuous(breaks = c(46,53,63,73,80))+
  scale_y_continuous(limits = c(7.8,11))+
  scale_color_npg()+scale_fill_npg()+
  mytheme1
shannon_maize_plot
setwd(wdOutput_Figure2)
getwd()
ggsave(paste("shannon_maize_plot",".pdf",sep=""),
       shannon_maize_plot,device=cairo_pdf,width=63,height=32,dpi = 600,units = "mm")

#### 4. Fig.2d ####
### 4.1 Import and process data ###
setwd(wdImport)
NST_results <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                      sheet = "fig 2e")

NST_results$group<-factor(NST_results$group,levels = c("Bulk","IP","IM","MP","MM"))

### 4.2 ###
compare_means(data=NST_results,NST_score~group,method = "kruskal.test")
compare_means(data=NST_results,NST_score~group,method = "wilcox.test")


### 4.3 plots ###
NST_results_Pots<-ggplot(NST_results,aes(x=group,y=NST_score))+
  geom_violin(aes(fill=group),width=1) +
  geom_jitter(aes(fill=group),width = 0.05,shape=21,color="black")+
  geom_boxplot(width = 0.35,outlier.alpha = 0.5,alpha=0.5,notch=T)+
  scale_color_manual(values = c("#c0c0c0","#e31a1c","#E19896","#4D78B2","#B3CDE2"))+
  scale_fill_manual(values = c("#c0c0c0","#e31a1c","#E19896","#4D78B2","#B3CDE2"))+
  mytheme1
NST_results_Pots
setwd(wdOutput_Figure2)
getwd()
ggsave(paste("NST_results_Pots",".pdf",sep=""),
       NST_results_Pots,device=cairo_pdf,width=140,height=80,dpi = 600,units = "mm")

