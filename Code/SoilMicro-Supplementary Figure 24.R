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
                             axis.line = element_line(color = "black",linewidth=0.4))#移除整体的边???

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
wdOutput <- ("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Supplemental materials/NIP")

#### 3. Peanut-Pot and field ####
#### 3.1 Import and process data ####
setwd(wdImport)
NIP_SPAD <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                             sheet = "fig s24 SPAD of NIP")
NIP_SPAD$Treatment3<-factor(NIP_SPAD$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

NIP_Iron <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                  sheet = "fig s24 iron of NIP")
NIP_Iron$Treatment3<-factor(NIP_Iron$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

NIP_Biomass <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                       sheet = "fig s24 biomass of NIP")
NIP_Biomass$Treatment3<-factor(NIP_Biomass$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

#### 4 NIP ####
## 4.1 NIP-SPAD ###
# 4.1.1 statistical analysis #
leveneTest(YL_SPAD ~ Treatment3, data = NIP_SPAD)#p>0.05，则满足方差齐性
shapiro.test(NIP_SPAD$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
compare_means(data=NIP_SPAD,YL_SPAD~Treatment3,method = "t.test")
aov_model_NIP_SPAD<-aov(data=NIP_SPAD,YL_SPAD~Treatment3)
summary(aov_model_NIP_SPAD)
LSD_model_NIP_SPAD<- LSD.test(aov_model_NIP_SPAD,"Treatment",p.adj = "BH")
LSD_model_NIP_SPAD

# 4.1.2 Plots #
NIP_SPAD_Bar<-ggplot(NIP_SPAD,aes(Treatment3,YL_SPAD))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",y='SPAD')+
  scale_y_continuous(limits = c(0,50))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
NIP_SPAD_Bar
setwd(wdOutput)
getwd()
ggsave(paste("NIP_SPAD_Bar",".pdf",sep=""),
       NIP_SPAD_Bar,device=cairo_pdf,width=40,height=72,dpi = 300,units = "mm")

## 4.2. NIP-Acitve Fe ###
# 4.2.1 statistical analysis #
leveneTest(YL_ActiveFe ~ Treatment3, data = NIP_Iron)#p>0.05，则满足方差齐性
shapiro.test(NIP_Iron$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
NIP_Iron<-NIP_Iron%>%mutate(boxcox_YL_ActiveFe=BoxCox(NIP_Iron$YL_ActiveFe,lambda="auto"))
shapiro.test(NIP_Iron$boxcox_YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_NIP_ActiveFe<-aov(data=NIP_Iron,boxcox_YL_ActiveFe~Treatment3)
summary(aov_model_NIP_ActiveFe)
LSD_model_NIP_ActiveFe<- LSD.test(aov_model_NIP_ActiveFe,"Treatment3",p.adj = "BH")
LSD_model_NIP_ActiveFe

# 4.2.2 statistical analysis #
NIP_YL_ActiveFe_Bar<-ggplot(NIP_Iron,aes(Treatment3,YL_ActiveFe))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",
       y=expression('Active Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,15))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
NIP_YL_ActiveFe_Bar
setwd(wdOutput)
getwd()
ggsave(paste("NIP_YL_ActiveFe_Bar",".pdf",sep=""),
       NIP_YL_ActiveFe_Bar,device=cairo_pdf,width=40,height=60,dpi = 300,units = "mm")

## 4.3 NIP-AvailableFe ###
# 4.3.1 statistical analysis #
leveneTest(availableFe ~ Treatment3, data = NIP_Iron)#p>0.05，则满足方差齐性
shapiro.test(NIP_Iron$availableFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_NIP_AvailableFe<-aov(data=NIP_Iron,availableFe~Treatment3)
summary(aov_model_NIP_AvailableFe)
LSD_model_NIP_AvailableFe<- LSD.test(aov_model_NIP_AvailableFe,"Treatment3",p.adj = "BH")
LSD_model_NIP_AvailableFe

# 4.3.2 Plots #
NIP_AvailableFe_Bar<-ggplot(NIP_Iron,aes(Treatment3,availableFe))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",
       y=expression('Available Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,8))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
NIP_AvailableFe_Bar
setwd(wdOutput)
getwd()
ggsave(paste("NIP_AvailableFe_Bar",".pdf",sep=""),
       NIP_AvailableFe_Bar,device=cairo_pdf,width=50,height=60,dpi = 300,units = "mm")

## 4.4 NIP-Biomass ###
# 4.4.1 statistical analysis #
leveneTest(Total ~ Treatment3, data = NIP_Biomass)#p>0.05，则满足方差齐性
shapiro.test(NIP_Biomass$Total)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_NIP_Biomass<-aov(data=NIP_Biomass,Total~Treatment3)
summary(aov_model_NIP_Biomass)
LSD_model_NIP_Biomass<- LSD.test(aov_model_NIP_Biomass,"Treatment3",p.adj = "BH")
LSD_model_NIP_Biomass

# 4.4.2 Plots #
NIP_Biomass_Bar<-ggplot(NIP_Biomass,aes(Treatment3,Total))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",
       y=expression('Available Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,17))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
NIP_Biomass_Bar
setwd(wdOutput)
getwd()
ggsave(paste("NIP_Biomass_Bar",".pdf",sep=""),
       NIP_Biomass_Bar,device=cairo_pdf,width=50,height=60,dpi = 300,units = "mm")
