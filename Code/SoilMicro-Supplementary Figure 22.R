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
wdOutput <- ("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Supplemental materials/treatmens_maize")

#### 3. Peanut-Pot and field ####
#### 3.1 Import and process data ####
setwd(wdImport)
SterilePotSPAD <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                             sheet = "fig s22 SPAD")
SterilePotSPAD$Treatment3<-factor(SterilePotSPAD$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

SterilePot_ActiveFe <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                  sheet = "fig s22 activeFe")
SterilePot_ActiveFe$Treatment3<-factor(SterilePot_ActiveFe$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

SterilePot_Biomass <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                 sheet = "fig s22 biomass")
SterilePot_Biomass$Treatment3<-factor(SterilePot_Biomass$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

#SterilePot_AvailableFe <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
#                                 sheet = "fig 4b availableFe in pot")
#SterilePot_AvailableFe$Treatment3<-factor(SterilePot_AvailableFe$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

#### 4 SIM ####
## 4.1 SIM-SPAD ###
SIM_SPAD<-SterilePotSPAD%>%filter(Treatment2=="SIM")
# 4.1.1 statistical analysis #
leveneTest(YL_SPAD ~ Treatment3, data = SIM_SPAD)#p>0.05，则满足方差齐性
shapiro.test(SIM_SPAD$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
SIM_SPAD<-SIM_SPAD%>%mutate(boxcox_YL_SPAD=BoxCox(SIM_SPAD$YL_SPAD,lambda="auto"))
shapiro.test(SIM_SPAD$boxcox_YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_SIM_SPAD<-aov(data=SIM_SPAD,YL_SPAD~Treatment3)
summary(aov_model_SIM_SPAD)
#Nonparametric tests
kruskal.test(YL_SPAD~Treatment3, data = SIM_SPAD)
aov_model_SIM_SPAD<-aov(data=SIM_SPAD,YL_SPAD~Treatment3)
dunnettT3Test(aov_model_SIM_SPAD,p.adjust.method = "BH")

# 4.1.2 Plots #
SIM_SPAD_Bar<-ggplot(SIM_SPAD,aes(Treatment3,YL_SPAD))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.5,size=0.15)+
  geom_point(aes(fill=Treatment3),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2,alpha=0.7)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",y='SPAD')+
  scale_y_continuous(limits = c(0,40))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
SIM_SPAD_Bar
setwd(wdOutput)
getwd()
ggsave(paste("SIM_SPAD_Bar",".pdf",sep=""),
       SIM_SPAD_Bar,device=cairo_pdf,width=35,height=36,dpi = 300,units = "mm")

## 4.2. SIM-Acitve Fe ###
SIM_ActiveFe<-SterilePot_ActiveFe%>%filter(Treatment2=="SIM")
SIM_ActiveFe
# 4.2.1 statistical analysis #
leveneTest(YL_ActiveFe ~ Treatment3, data = SIM_ActiveFe)#p>0.05，则满足方差齐性
shapiro.test(SIM_ActiveFe$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
SIM_ActiveFe<-SIM_ActiveFe%>%mutate(boxcox_YL_ActiveFe=BoxCox(SIM_ActiveFe$YL_ActiveFe,BoxCox.lambda(SIM_ActiveFe$YL_ActiveFe)))
leveneTest(boxcox_YL_ActiveFe ~ Treatment3, data = SIM_ActiveFe)#p>0.05，则满足方差齐性
shapiro.test(SIM_ActiveFe$boxcox_YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_SIM_ActiveFe<-aov(data=SIM_ActiveFe,boxcox_YL_ActiveFe~Treatment3)
summary(aov_model_SIM_ActiveFe)
LSD.test(aov_model_SIM_ActiveFe,"Treatment3",p.adj = "BH",console = T)
LSD.test(aov_model_SIM_ActiveFe,"Treatment3",p.adj = "BH",console = T,group = F)

# 4.2.2 statistical analysis #
SIM_YL_ActiveFe_Bar<-ggplot(SIM_ActiveFe,aes(Treatment3,YL_ActiveFe))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.5,size=0.15)+
  geom_point(aes(fill=Treatment3),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",
       y=expression('Active Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,15))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
SIM_YL_ActiveFe_Bar
setwd(wdOutput)
getwd()
ggsave(paste("SIM_YL_ActiveFe_Bar",".pdf",sep=""),
       SIM_YL_ActiveFe_Bar,device=cairo_pdf,width=35,height=36,dpi = 300,units = "mm")

# ## 4.3 SIM-AvailableFe ###
# SIM_AvailableFe<-SterilePot_AvailableFe%>%filter(Treatment2=="SIM")
# # 4.3.1 statistical analysis #
# leveneTest(availableFe ~ Treatment3, data = SIM_AvailableFe)#p>0.05，则满足方差齐性
# shapiro.test(SIM_AvailableFe$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
# SIM_AvailableFe<-SIM_AvailableFe%>%mutate(boxcox_YL_ActiveFe=BoxCox(SIM_ActiveFe$YL_ActiveFe,lambda="auto"))
# shapiro.test(SIM_ActiveFe$boxcox_YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
# aov_model_SIM_ActiveFe<-aov(data=SIM_ActiveFe,boxcox_YL_ActiveFe~Treatment3)
# compare_means(data=SIM_AvailableFe,availableFe~Treatment3,method = "t.test")
# aov_model_SIM_AvailableFe<-aov(data=SIM_AvailableFe,availableFe~Treatment3)
# summary(aov_model_SIM_AvailableFe)
# duncan_model_SIM_AvailableFe<- duncan.test(aov_model_SIM_AvailableFe,"Treatment3")
# duncan_model_SIM_AvailableFe
# 
# # 4.3.2 Plots #
# SIM_AvailableFe_Bar<-ggplot(SIM_AvailableFe,aes(Treatment3,availableFe))+
#   geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.5,size=0.15)+
#   geom_point(aes(fill=Treatment3),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
#   stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
#                geom='errorbar', width=0.2,size=0.2)+
#   labs(x="",
#        y=expression('Available Fe (μg '*g^{-1}*')'),parse =T)+
#   scale_y_continuous(limits = c(0,8))+
#   scale_color_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
#   mytheme+
#   guides(color="none")
# SIM_AvailableFe_Bar
# setwd(wdOutput)
# getwd()
# ggsave(paste("SIM_AvailableFe_Bar",".pdf",sep=""),
#        SIM_AvailableFe_Bar,device=cairo_pdf,width=50,height=59,dpi = 300,units = "mm")

## 4.5. SIM-Biomass ###
SIM_Biomass<-SterilePot_Biomass%>%filter(Treatment2=="SIM")
SIM_Biomass
# 4.5.1 statistical analysis #
leveneTest(Total ~ Treatment3, data = SIM_Biomass)#p>0.05，则满足方差齐性
shapiro.test(SIM_Biomass$Total)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_SIM_Biomass<-aov(data=SIM_Biomass,Total~Treatment3)
summary(aov_model_SIM_Biomass)
LSD.test(aov_model_SIM_Biomass,"Treatment3",p.adj = "BH",console = T)
LSD.test(aov_model_SIM_Biomass,"Treatment3",p.adj = "BH",console = T,group = F)

# 4.5.2 plots #
SIM_Biomass_Bar<-ggplot(SIM_Biomass,aes(Treatment3,Total))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.5,size=0.15)+
  geom_point(aes(fill=Treatment3),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",
       y=expression('Biomass(g)'),parse =T)+
  scale_y_continuous(limits = c(0,25))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
SIM_Biomass_Bar
setwd(wdOutput)
getwd()
ggsave(paste("SIM_Biomass_Bar",".pdf",sep=""),
       SIM_Biomass_Bar,device=cairo_pdf,width=35,height=36,dpi = 300,units = "mm")

#### 5 NIM ####
## 5.1 NIM-SPAD ###
NIM_SPAD<-SterilePotSPAD%>%filter(Treatment2=="NIM")
# 5.1.1 statistical analysis #
leveneTest(YL_SPAD ~ Treatment3, data = NIM_SPAD)#p>0.05，则满足方差齐性
shapiro.test(NIM_SPAD$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
NIM_SPAD<-NIM_SPAD%>%mutate(boxcox_YL_SPAD=BoxCox(NIM_SPAD$YL_SPAD,BoxCox.lambda(NIM_SPAD$YL_SPAD)))
leveneTest(boxcox_YL_SPAD ~ Treatment3, data = NIM_SPAD)#p>0.05，则满足方差齐性
shapiro.test(NIM_SPAD$boxcox_YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_NIM_SPAD<-aov(data=NIM_SPAD,YL_SPAD~Treatment3)
summary(aov_model_NIM_SPAD)
#Nonparametric tests
kruskal.test(YL_SPAD~Treatment3, data = NIM_SPAD)
aov_model_NIM_SPAD<-aov(data=NIM_SPAD,YL_SPAD~Treatment3)
dunnettT3Test(aov_model_NIM_SPAD,p.adjust.method = "BH")

# 5.1.2 Plots #
NIM_SPAD_Bar<-ggplot(NIM_SPAD,aes(Treatment3,YL_SPAD))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.5,size=0.15)+
  geom_point(aes(fill=Treatment3),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2,alpha=0.7)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",y='SPAD')+
  scale_y_continuous(limits = c(0,40))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
NIM_SPAD_Bar
setwd(wdOutput)
getwd()
ggsave(paste("NIM_SPAD_Bar",".pdf",sep=""),
       NIM_SPAD_Bar,device=cairo_pdf,width=35,height=36,dpi = 300,units = "mm")

## 5.2. NIM-Acitve Fe ###
NIM_ActiveFe<-SterilePot_ActiveFe%>%filter(Treatment2=="NIM")
NIM_ActiveFe
# 5.2.1 statistical analysis #
leveneTest(YL_ActiveFe ~ Treatment3, data = NIM_ActiveFe)#p>0.05，则满足方差齐性
shapiro.test(NIM_ActiveFe$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
NIM_ActiveFe<-NIM_ActiveFe%>%mutate(boxcox_YL_ActiveFe=BoxCox(NIM_ActiveFe$YL_ActiveFe,BoxCox.lambda(NIM_ActiveFe$YL_ActiveFe)))
leveneTest(boxcox_YL_ActiveFe ~ Treatment3, data = NIM_ActiveFe)#p>0.05，则满足方差齐性
shapiro.test(NIM_ActiveFe$boxcox_YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_NIM_ActiveFe<-aov(data=NIM_ActiveFe,YL_ActiveFe~Treatment3)
summary(aov_model_NIM_ActiveFe)
LSD.test(aov_model_NIM_ActiveFe,"Treatment3",p.adj = "BH",console = T)
LSD.test(aov_model_NIM_ActiveFe,"Treatment3",p.adj = "BH",console = T,group = F)

# 5.2.2 statistical analysis #
NIM_YL_ActiveFe_Bar<-ggplot(NIM_ActiveFe,aes(Treatment3,YL_ActiveFe))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.5,size=0.15)+
  geom_point(aes(fill=Treatment3),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",
       y=expression('Active Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,17))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
NIM_YL_ActiveFe_Bar
setwd(wdOutput)
getwd()
ggsave(paste("NIM_YL_ActiveFe_Bar",".pdf",sep=""),
       NIM_YL_ActiveFe_Bar,device=cairo_pdf,width=35,height=36,dpi = 300,units = "mm")

# ## 5.3 NIM-AvailableFe ###
# NIM_AvailableFe<-SterilePot_AvailableFe%>%filter(Treatment2=="NIM")
# # 5.3.1 statistical analysis #
# compare_means(data=NIM_AvailableFe,availableFe~Treatment3,method = "t.test")
# aov_model_NIM_AvailableFe<-aov(data=NIM_AvailableFe,availableFe~Treatment3)
# summary(aov_model_NIM_AvailableFe)
# duncan_model_NIM_AvailableFe<- duncan.test(aov_model_NIM_AvailableFe,"Treatment3")
# duncan_model_NIM_AvailableFe
# 
# # 5.3.2 Plots #
# NIM_AvailableFe_Bar<-ggplot(NIM_AvailableFe,aes(Treatment3,availableFe))+
#   geom_bar(stat = "summary", fun = "mean",color="black",fill="white",width=0.65,size=0.2)+
#   stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
#                geom='errorbar', width=0.15,size=0.15)+
#   geom_jitter(aes(Treatment3,color=Treatment3),width = 0.2, height = 0.2,size=0.25)+
#   labs(x="",
#        y=expression('Available Fe (μg '*g^{-1}*')'),parse =T)+
#   scale_y_continuous(limits = c(0,8))+
#   scale_color_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
#   mytheme+
#   guides(color=F)
# NIM_AvailableFe_Bar
# setwd(wdOutput)
# getwd()
# ggsave(paste("NIM_AvailableFe_Bar",".pdf",sep=""),
#        NIM_AvailableFe_Bar,device=cairo_pdf,width=50,height=59,dpi = 300,units = "mm")

## 5.5. NIM-Biomass ###
NIM_Biomass<-SterilePot_Biomass%>%filter(Treatment2=="NIM")
NIM_Biomass
# 5.5.1 statistical analysis #
leveneTest(Total ~ Treatment3, data = NIM_Biomass)#p>0.05，则满足方差齐性
shapiro.test(NIM_Biomass$Total)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_NIM_Biomass<-aov(data=NIM_Biomass,Total~Treatment3)
summary(aov_model_NIM_Biomass)
LSD.test(aov_model_NIM_Biomass,"Treatment3",p.adj = "BH",console = T)
LSD.test(aov_model_NIM_Biomass,"Treatment3",p.adj = "BH",console = T,group = F)

# 5.5.2 statistical analysis #
NIM_Biomass_Bar<-ggplot(NIM_Biomass,aes(Treatment3,Total))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.5,size=0.15)+
  geom_point(aes(fill=Treatment3),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(x="",
       y=expression('Biomass(g)'),parse =T)+
  scale_y_continuous(limits = c(0,25))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
NIM_Biomass_Bar
setwd(wdOutput)
getwd()
ggsave(paste("NIM_Biomass_Bar",".pdf",sep=""),
       NIM_Biomass_Bar,device=cairo_pdf,width=35,height=36,dpi = 300,units = "mm")

