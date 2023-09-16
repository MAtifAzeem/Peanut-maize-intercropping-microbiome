#### 1. Loading packet####
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
library(rcompanion)
library(forecast)#box-cox数据变换
library(PMCMRplus)
#### 2. setting theme and filepath ####
loadfonts()
Sys.setenv(R_GSCMD = "C:/Program Files (x86)/gs/gs9.50/bin/gswin32c.exe")

mytheme <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#000000"),
                             text = element_text(family = "Arial"),
                             strip.text = element_text(size=7,hjust = 0.5),
                             plot.title = element_text(size=7,hjust = 0.5),
                             axis.text=element_text(size=7,color = "#808080"),
                             axis.title=element_text(size=7),
                             legend.text = element_text(size=7),
                             legend.title = element_text(size=7),
                             legend.background = element_blank(),
                             panel.border = element_rect(colour = NA),
                             axis.line = element_line(color="black",size=0.2),
                             axis.ticks = element_line(color="black",size=0.2,lineend = 0.1),
                             axis.ticks.length = unit(0.8, "mm"))#移除整体的边???

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
wdOutput_Figure5 <- c("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Figure5")

#### 3. PAO1 ####
#### 3.1 Import and process data ####
setwd(wdImport)
PAO1_SPAD <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                        sheet = "fig 5b SPAD")

PAO1_SPAD$treatment<-factor(PAO1_SPAD$treatment,levels=c("CK","mutant","wt","1502IPR-01"))

PAO1_ActiveFe <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                            sheet = "fig 5b activeFe")
PAO1_ActiveFe$treatment<-factor(PAO1_ActiveFe$treatment,levels=c("CK","mutant","wt","1502IPR-01"))

PAO1_AvailableFe <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                               sheet = "fig 5b availableFe")
PAO1_AvailableFe$treatment<-factor(PAO1_AvailableFe$treatment,levels=c("CK","mutant","wt","1502IPR-01"))

PAO1_Biomass <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                           sheet = "fig 5b biomass")
PAO1_Biomass$treatment<-factor(PAO1_Biomass$treatment,levels=c("CK","mutant","wt","1502IPR-01"))

PAO1_Reduction <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                             sheet = "fig 5b reduction")
PAO1_Reduction$treatment<-factor(PAO1_Reduction$treatment,levels=c("CK","mutant","wt","1502IPR-01"))

#### 4.1 PAO1 SPAD ####
# 4.1.1 statistical analysis #
leveneTest(YL_SPAD ~ treatment, data = PAO1_SPAD)#p>0.05，则满足方差齐性
shapiro.test(PAO1_SPAD$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_PAO1_SPAD<-aov(data=PAO1_SPAD,YL_SPAD~treatment)
summary(aov_model_PAO1_SPAD)
LSD_model_PAO1_SPAD<-LSD.test(aov_model_PAO1_SPAD,"treatment",p.adj = "BH")
LSD_model_PAO1_SPAD

# 4.2 Plots #
PAO1_SPAD_Bar<-ggplot(PAO1_SPAD,aes(treatment,YL_SPAD))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=treatment),width=0.5,size=0.15)+
  geom_point(aes(fill=treatment),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2,alpha=0.6)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  scale_y_continuous(limits = c(0,45))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9","#00A087"))+
  mytheme+
  guides(fill="none")
PAO1_SPAD_Bar
setwd(wdOutput_Figure5)
getwd()
ggsave(paste("PAO1_SPAD_Bar",".pdf",sep=""),
       PAO1_SPAD_Bar,device=cairo_pdf,width=35,height=36,dpi = 300,units = "mm")

## 4.2 PAO1-Acitve Fe ###
# 4.2.1 statistical analysis #
leveneTest(YL_ActiveFe ~ treatment, data = PAO1_ActiveFe)#p>0.05，则满足方差齐性
shapiro.test(PAO1_ActiveFe$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_PAO1_ActiveFe<-aov(data=PAO1_ActiveFe,YL_ActiveFe~treatment)
summary(aov_model_PAO1_ActiveFe)
LSD_model_PAO1_ActiveFe<-LSD.test(aov_model_PAO1_ActiveFe,"treatment",p.adj = "BH")
LSD_model_PAO1_ActiveFe

# 4.2.2 statistical analysis #
PAO1_YL_ActiveFe_Bar<-ggplot(PAO1_ActiveFe,aes(treatment,YL_ActiveFe))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=treatment),width=0.5,size=0.15)+
  geom_point(aes(fill=treatment),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",
       y=expression('Active Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,16))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9","#00A087"))+
  mytheme+
  guides(fill="none")
PAO1_YL_ActiveFe_Bar
setwd(wdOutput_Figure5)
getwd()
ggsave(paste("PAO1_YL_ActiveFe_Bar",".pdf",sep=""),
       PAO1_YL_ActiveFe_Bar,device=cairo_pdf,width=35,height=36,dpi = 300,units = "mm")

## 4.3 PAO1-AvailableFe ###
# 4.3.1 statistical analysis #
leveneTest(availableFe ~ treatment, data = PAO1_AvailableFe)#p>0.05，则满足方差齐性
shapiro.test(PAO1_AvailableFe$availableFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_PAO1_AvailableFe<-aov(data=PAO1_AvailableFe,availableFe~treatment)
summary(aov_model_PAO1_AvailableFe)
LSD_model_PAO1_AvailableFe<-LSD.test(aov_model_PAO1_AvailableFe,"treatment",p.adj = "BH")
LSD_model_PAO1_AvailableFe

# 4.3.2 Plots #
PAO1_AvailableFe_Bar<-ggplot(PAO1_AvailableFe,aes(treatment,availableFe))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=treatment),width=0.5,size=0.15)+
  geom_point(aes(fill=treatment),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",
       y=expression('Available Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,8))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9","#00A087"))+
  mytheme+
  guides(fill="none")
PAO1_AvailableFe_Bar
setwd(wdOutput_Figure5)
getwd()
ggsave(paste("PAO1_AvailableFe_Bar",".pdf",sep=""),
       PAO1_AvailableFe_Bar,device=cairo_pdf,width=35,height=36,dpi = 300,units = "mm")

## 4.4. PAO1-Biomass ###
# 4.4.1 statistical analysis #
leveneTest(biomass ~ treatment, data = PAO1_Biomass)#p>0.05，则满足方差齐性
shapiro.test(PAO1_Biomass$biomass)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_PAO1_Biomass<-aov(data=PAO1_Biomass,biomass~treatment)
LSD_model_PAO1_Biomass<-LSD.test(aov_model_PAO1_Biomass,"treatment",p.adj = "BH")
LSD_model_PAO1_Biomass

# 4.4.2 statistical analysis #
PAO1_Biomass_Bar<-ggplot(PAO1_Biomass,aes(treatment,biomass))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=treatment),width=0.5,size=0.15)+
  geom_point(aes(fill=treatment),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",
       y=expression('Biomass(g)'),parse =T)+
  scale_y_continuous(limits = c(0,6.5))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9","#00A087"))+
  mytheme+
  guides(fill="none")
PAO1_Biomass_Bar
setwd(wdOutput_Figure5)
getwd()
ggsave(paste("PAO1_Biomass_Bar",".pdf",sep=""),
       PAO1_Biomass_Bar,device=cairo_pdf,width=35,height=36,dpi = 300,units = "mm")

## 4.5. PAO1-Reduction ###
# 4.5.1 statistical analysis #
leveneTest(reduction ~ treatment, data = PAO1_Reduction)#p>0.05，则满足方差齐性
shapiro.test(PAO1_Reduction$reduction)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
PAO1_Reduction<-PAO1_Reduction%>%mutate(boxcox_reduction =BoxCox(PAO1_Reduction$reduction,lambda="auto"))
leveneTest(boxcox_reduction ~ treatment, data = PAO1_Reduction)#p>0.05，则满足方差齐性
shapiro.test(PAO1_Reduction$boxcox_reduction)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
#Non-parameter test
kruskal.test(reduction~treatment, data = PAO1_Reduction)
aov_model_NMP_Reduction<-aov(data=PAO1_Reduction,reduction~treatment)
dunnettT3Test(aov_model_NMP_Reduction,p.adjust.method = "BH")


# 4.5.2 statistical analysis #
PAO1_Reduction_Bar<-ggplot(PAO1_Reduction,aes(treatment,reduction))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=treatment),width=0.5,size=0.15)+
  geom_point(aes(fill=treatment),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",
       y=expression('Active Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,0.6))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9","#00A087"))+
  mytheme+
  guides(fill="none")
PAO1_Reduction_Bar
setwd(wdOutput_Figure5)
getwd()
ggsave(paste("PAO1_Reduction_Bar",".pdf",sep=""),
       PAO1_Reduction_Bar,device=cairo_pdf,width=35,height=36,dpi = 300,units = "mm")
