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
                                axis.line = element_line(color = "black",linewidth=0.4))#移除整体的边???
wdImport<- c("E:/working/SCI/Soil Micro/SCI/Figures/Data/Data for submit")
wdOutput <- ("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Supplemental materials/treatments_reduction_DMA")

#### 3. Peanut-Pot and field ####
#### 3.1 Import and process data ####
setwd(wdImport)
SterilePotreduction<- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                             sheet = "fig s25 reduction")
SterilePotreduction$Treatment3<-factor(SterilePotreduction$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

SterilePot_MArelease <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                  sheet = "fig s25 DMA")
SterilePot_MArelease$Treatment3<-factor(SterilePot_MArelease$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

#### 4 SIP ####
## 4.1 SIP-reduction ###
SIP_reduction<-SterilePotreduction%>%filter(Treatment2=="SIP")
# 4.1.1 statistical analysis #
leveneTest(ReductionValue ~ Treatment3, data = SIP_reduction)#p>0.05，则满足方差齐性
shapiro.test(SIP_reduction$ReductionValue)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_SIP_reduction<-aov(data=SIP_reduction,ReductionValue~Treatment3)
summary(aov_model_SIP_reduction)
LSD_model_SIP_reduction<- LSD.test(aov_model_SIP_reduction,"Treatment3",p.adj = "BH")
LSD_model_SIP_reduction
# 4.1.2 Plots #
SIP_reduction_Bar<-ggplot(SIP_reduction,aes(Treatment3,ReductionValue))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.5,size=0.15)+
  geom_point(aes(fill=Treatment3),shape=21,position = position_jitterdodge(1),size=1,stroke = 0.2,alpha=0.7)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",y='reduction')+
  scale_y_continuous(limits = c(0,120))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
SIP_reduction_Bar
setwd(wdOutput)
getwd()
ggsave(paste("SIP_reduction_Bar",".pdf",sep=""),
       SIP_reduction_Bar,device=cairo_pdf,width=40,height=42,dpi = 300,units = "mm")

#### 5 NIP ####
## 5.1 NIP-reduction ###
NIP_reduction<-SterilePotreduction%>%filter(Treatment2=="NIP")
# 5.1.1 statistical analysis #
leveneTest(ReductionValue ~ Treatment3, data = NIP_reduction)#p>0.05，则满足方差齐性
shapiro.test(NIP_reduction$ReductionValue)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_NIP_reduction<-aov(data=NIP_reduction,ReductionValue~Treatment3)
summary(aov_model_NIP_reduction)
LSD_model_NIP_reduction<- LSD.test(aov_model_NIP_reduction,"Treatment3",p.adj = "BH")
LSD_model_NIP_reduction
# 5.1.2 Plots #
NIP_reduction_Bar<-ggplot(NIP_reduction,aes(Treatment3,ReductionValue))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.5,size=0.15)+
  geom_point(aes(fill=Treatment3),shape=21,position = position_jitterdodge(1),size=1,stroke = 0.2,alpha=0.7)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",y='reduction')+
  scale_y_continuous(limits = c(0,200))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
NIP_reduction_Bar
setwd(wdOutput)
getwd()
ggsave(paste("NIP_reduction_Bar",".pdf",sep=""),
       NIP_reduction_Bar,device=cairo_pdf,width=40,height=42,dpi = 300,units = "mm")

#### 6 SMP ####
## 6.1 SMP-reduction ###
SMP_reduction<-SterilePotreduction%>%filter(Treatment2=="SMP")
# 6.1.1 statistical analysis #
leveneTest(ReductionValue ~ Treatment3, data = SMP_reduction)#p>0.05，则满足方差齐性
shapiro.test(SMP_reduction$ReductionValue)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_SMP_reduction<-aov(data=SMP_reduction,ReductionValue~Treatment3)
summary(aov_model_SMP_reduction)
LSD_model_SMP_reduction<- LSD.test(aov_model_SMP_reduction,"Treatment3",p.adj = "BH")
LSD_model_SMP_reduction
# 6.1.2 Plots #
SMP_reduction_Bar<-ggplot(SMP_reduction,aes(Treatment3,ReductionValue))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.5,size=0.15)+
  geom_point(aes(fill=Treatment3),shape=21,position = position_jitterdodge(1),size=1,stroke = 0.2,alpha=0.7)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",y='reduction')+
  scale_y_continuous(limits = c(0,250))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
SMP_reduction_Bar
setwd(wdOutput)
getwd()
ggsave(paste("SMP_reduction_Bar",".pdf",sep=""),
       SMP_reduction_Bar,device=cairo_pdf,width=40,height=42,dpi = 300,units = "mm")

#### 7 NMP ####
## 7.1 NMP-reduction ###
NMP_reduction<-SterilePotreduction%>%filter(Treatment2=="NMP")
# 7.1.1 statistical analysis #
leveneTest(ReductionValue ~ Treatment3, data = NMP_reduction)#p>0.05，则满足方差齐性
shapiro.test(NMP_reduction$ReductionValue)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
NMP_reduction<-NMP_reduction%>%mutate(boxcox_ReductionValue=BoxCox(NMP_reduction$ReductionValue,lambda = "auto"))
leveneTest(boxcox_ReductionValue ~ Treatment3, data = NMP_reduction)#p>0.05，则满足方差齐性
shapiro.test(NMP_reduction$boxcox_ReductionValue)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_NMP_reduction<-aov(data=NMP_reduction,boxcox_ReductionValue~Treatment3)
summary(aov_model_NMP_reduction)
LSD_model_NMP_reduction<- LSD.test(aov_model_SMP_reduction,"Treatment3",p.adj = "BH")
LSD_model_NMP_reduction
# 7.1.2 Plots #
NMP_reduction_Bar<-ggplot(NMP_reduction,aes(Treatment3,ReductionValue))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.5,size=0.15)+
  geom_point(aes(fill=Treatment3),shape=21,position = position_jitterdodge(1),size=1,stroke = 0.2,alpha=0.7)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",y='reduction')+
  scale_y_continuous(limits = c(0,300))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
NMP_reduction_Bar
setwd(wdOutput)
getwd()
ggsave(paste("NMP_reduction_Bar",".pdf",sep=""),
       NMP_reduction_Bar,device=cairo_pdf,width=40,height=42,dpi = 300,units = "mm")

## 8. SIM-MArelease ###
SIM_MArelease<-SterilePot_MArelease%>%filter(Treatment2=="SIM")
SIM_MArelease
# 8.1 statistical analysis #
leveneTest(MAReleaseValue ~ Treatment3, data = SIM_MArelease)#p>0.05，则满足方差齐性
shapiro.test(SIM_MArelease$MAReleaseValue)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_SIM_MArelease<-aov(data=SIM_MArelease,MAReleaseValue~Treatment3)
summary(aov_model_SIM_MArelease)
LSD_model_SIM_MArelease<- LSD.test(aov_model_SIM_MArelease,"Treatment3",p.adj = "BH")
LSD_model_SIM_MArelease

# 8.2 statistical analysis #
SIM_MArelease_Bar<-ggplot(SIM_MArelease,aes(Treatment3,MAReleaseValue))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.5,size=0.15)+
  geom_point(aes(fill=Treatment3),shape=21,position = position_jitterdodge(1),size=1,stroke = 0.2,alpha=0.7)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",y='MArelease')+
  scale_y_continuous(limits = c(0,600))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
SIM_MArelease_Bar
setwd(wdOutput)
getwd()
ggsave(paste("SIM_MArelease_Bar",".pdf",sep=""),
       SIM_MArelease_Bar,device=cairo_pdf,width=40,height=42,dpi = 300,units = "mm")

## 5.4. NIM-MArelease ###
NIM_MArelease<-SterilePot_MArelease%>%filter(Treatment2=="NIM")
NIM_MArelease
# 5.4.1 statistical analysis #
leveneTest(MAReleaseValue ~ Treatment3, data = NIM_MArelease)#p>0.05，则满足方差齐性
shapiro.test(NIM_MArelease$MAReleaseValue)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_NIM_MArelease<-aov(data=NIM_MArelease,MAReleaseValue~Treatment3)
summary(aov_model_NIM_MArelease)
LSD_model_NIM_MArelease<- LSD.test(aov_model_NIM_MArelease,"Treatment3",p.adj = "BH")
LSD_model_NIM_MArelease

# 5.4.2 statistical analysis #
NIM_MArelease_Bar<-ggplot(NIM_MArelease,aes(Treatment3,MAReleaseValue))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.5,size=0.15)+
  geom_point(aes(fill=Treatment3),shape=21,position = position_jitterdodge(1),size=1,stroke = 0.2,alpha=0.7)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",y='MArelease')+
  scale_y_continuous(limits = c(0,700))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
NIM_MArelease_Bar
setwd(wdOutput)
getwd()
ggsave(paste("NIM_MArelease_Bar",".pdf",sep=""),
       NIM_MArelease_Bar,device=cairo_pdf,width=40,height=42,dpi = 300,units = "mm")

