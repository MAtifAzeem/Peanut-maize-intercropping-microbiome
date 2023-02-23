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
wdOutput_Figure4 <- c("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Figure4")

#### 3. Peanut-Pot and field ####
#### 3.1 Import and process data ####
setwd(wdImport)
SterilePotSPAD <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                             sheet = "fig 4b SPAD in pot")
SterilePotSPAD$Treatment3<-factor(SterilePotSPAD$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

SterilePot_ActiveFe <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                  sheet = "fig 4b activeFe in pot")
SterilePot_ActiveFe$Treatment3<-factor(SterilePot_ActiveFe$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

SterilePot_Biomass<-read_excel("Intercropping-microbiome-Data for submit.xlsx",
                               sheet = "fig 4b biomass in pot")
SterilePot_Biomass$Treatment3<-factor(SterilePot_Biomass$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

SterilePot_AvailableFe <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                  sheet = "fig 4b availableFe in pot")
SterilePot_AvailableFe$Treatment3<-factor(SterilePot_AvailableFe$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

Field_SPAD <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                         sheet = "fig 4d SPAD in field") 
Field_SPAD$Treatment<-factor(Field_SPAD$Treatment,levels=c("CK","1502IPR-01","Pyoverdine","EDTA-Fe"))


Field_Iron <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                         sheet = "fig 4d iron in field")
Field_Iron$Treatment<-factor(Field_Iron$Treatment,levels=c("CK","1502IPR-01","Pyoverdine","EDTA-Fe"))

Field_Yield <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                          sheet = "fig 4d yield in field")
Field_Yield$Treatment<-factor(Field_Yield$Treatment,levels=c("CK","1502IPR-01","Pyoverdine","EDTA-Fe"))

#### 4 SIP ####
## 4.1 SIP-SPAD ###
SIP_SPAD<-SterilePotSPAD%>%filter(Treatment2=="SIP")
# 4.1.1 statistical analysis #
leveneTest(YL_SPAD ~ Treatment3, data = SIP_SPAD)#p>0.05，则满足方差齐性
shapiro.test(SIP_SPAD$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
#BoxCox.lambda(PotSterilizationActiveFe$YL_ActiveFe,method="loglik")#a number indicating the Box-Cox transformation parameter
SIP_SPAD<-SIP_SPAD%>%mutate(boxcox_YL_SPAD =BoxCox(SIP_SPAD$YL_SPAD,lambda="auto"))
leveneTest(boxcox_YL_SPAD ~ Treatment3, data = SIP_SPAD)#p>0.05，则满足方差齐性
shapiro.test(SIP_SPAD$boxcox_YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_SIP_SPAD<-aov(data=SIP_SPAD,boxcox_YL_SPAD~Treatment3)
summary(aov_model_SIP_SPAD)
#Non-parameter test
kruskal.test(YL_SPAD~Treatment3, data = SIP_SPAD)
aov_model_SIP_SPAD<-aov(data=SIP_SPAD,YL_SPAD~Treatment3)
dunnettT3Test(aov_model_SIP_SPAD,p.adjust.method = "BH")
# 4.1.2 Plots #
SIP_SPAD_Bar<-ggplot(SIP_SPAD,aes(Treatment3,YL_SPAD))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.65,size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",y='SPAD')+
  scale_y_continuous(limits = c(0,40))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
SIP_SPAD_Bar
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("SIP_SPAD_Bar",".pdf",sep=""),
       SIP_SPAD_Bar,device=cairo_pdf,width=50,height=72,dpi = 300,units = "mm")

## 4.2. SIP-Acitve Fe ###
SIP_ActiveFe<-SterilePot_ActiveFe%>%filter(Treatment2=="SIP")
SIP_ActiveFe
# 4.2.1 statistical analysis #
leveneTest(YL_ActiveFe ~ Treatment3, data = SIP_ActiveFe)#p>0.05，则满足方差齐性
shapiro.test(SIP_ActiveFe$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_SIP_ActiveFe<-aov(data=SIP_ActiveFe,YL_ActiveFe~Treatment3)
summary(aov_model_SIP_ActiveFe)
LSD_model_SIP_ActiveFe<-LSD.test(aov_model_SIP_ActiveFe,"Treatment3",p.adj = "BH")
LSD_model_SIP_ActiveFe

# 4.2.2 statistical analysis #
SIP_YL_ActiveFe_Bar<-ggplot(SIP_ActiveFe,aes(Treatment3,YL_ActiveFe))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.65,size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",
       y=expression('Active Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,7))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
SIP_YL_ActiveFe_Bar
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("SIP_YL_ActiveFe_Bar",".pdf",sep=""),
       SIP_YL_ActiveFe_Bar,device=cairo_pdf,width=50,height=72,dpi = 300,units = "mm")

## 4.3 SIP-biomass ###
SIP_Biomass<-SterilePot_Biomass%>%filter(Treatment2=="SIP")
# 4.3.1 statistical analysis #
leveneTest(Total ~ Treatment3, data = SIP_Biomass)#p>0.05，则满足方差齐性
shapiro.test(SIP_Biomass$Total)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
SIP_Biomass<-SIP_Biomass%>%mutate(boxcox_Total =BoxCox(SIP_Biomass$Total,lambda="auto"))
leveneTest(boxcox_Total ~ Treatment3, data = SIP_Biomass)#p>0.05，则满足方差齐性
shapiro.test(SIP_Biomass$boxcox_Total)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_SIP_Biomass<-aov(data=SIP_Biomass,boxcox_Total~Treatment3)
summary(aov_model_SIP_Biomass)
LSD_model_SIP_Biomass<-LSD.test(aov_model_SIP_Biomass,"Treatment3",p.adj = "BH")
LSD_model_SIP_Biomass

# 4.3.2 Plots #
SIP_Biomass_Bar<-ggplot(SIP_Biomass,aes(Treatment3,Total))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.65,size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",
       y="Biomas (g/plant)",parse =T)+
  scale_y_continuous(limits = c(0,15))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
SIP_Biomass_Bar
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("SIP_Biomass_Bar",".pdf",sep=""),
       SIP_Biomass_Bar,device=cairo_pdf,width=50,height=59,dpi = 300,units = "mm")

## 4.4 SIP-AvailableFe ###
SIP_AvailableFe<-SterilePot_AvailableFe%>%filter(Treatment2=="SIP")
# 4.4.1 statistical analysis #
leveneTest(availableFe ~ Treatment3, data = SIP_AvailableFe)#p>0.05，则满足方差齐性
shapiro.test(SIP_AvailableFe$availableFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_SIP_AvailableFe<-aov(data=SIP_AvailableFe,availableFe~Treatment3)
summary(aov_model_SIP_AvailableFe)
LSD_model_SIP_AvailableFe<-LSD.test(aov_model_SIP_AvailableFe,"Treatment3",p.adj = "BH")
LSD_model_SIP_AvailableFe

# 4.4.2 Plots #
SIP_AvailableFe_Bar<-ggplot(SIP_AvailableFe,aes(Treatment3,availableFe))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.65,size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",
       y=expression('Available Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,8))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
SIP_AvailableFe_Bar
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("SIP_AvailableFe_Bar",".pdf",sep=""),
       SIP_AvailableFe_Bar,device=cairo_pdf,width=50,height=59,dpi = 300,units = "mm")

#### 5 SMP ####
## 5.2.1 SMP-SPAD ###
SMP_SPAD<-SterilePotSPAD%>%filter(Treatment2=="SMP")
# 5.2.1.1 statistical analysis #
leveneTest(YL_SPAD ~ Treatment3, data = SMP_SPAD)#p>0.05，则满足方差齐性
shapiro.test(SMP_SPAD$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
SMP_SPAD<-SMP_SPAD%>%mutate(boxcox_YL_SPAD =BoxCox(SMP_SPAD$YL_SPAD,lambda="auto"))
leveneTest(boxcox_YL_SPAD ~ Treatment3, data = SMP_SPAD)#p>0.05，则满足方差齐性
shapiro.test(SMP_SPAD$boxcox_YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_SMP_SPAD<-aov(data=SMP_SPAD,boxcox_YL_SPAD~Treatment3)
summary(aov_model_SMP_SPAD)
LSD_model_SMP_SPAD<-LSD.test(aov_model_SMP_SPAD,"Treatment3",p.adj = "BH")
LSD_model_SMP_SPAD

# 5.2.1.2 Plots #
SMP_SPAD_Bar<-ggplot(SMP_SPAD,aes(Treatment3,YL_SPAD))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.65,size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",y='SPAD')+
  scale_y_continuous(limits = c(0,40))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
SMP_SPAD_Bar
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("SMP_SPAD_Bar",".pdf",sep=""),
       SMP_SPAD_Bar,device=cairo_pdf,width=50,height=72,dpi = 300,units = "mm")

## 5.2.2 SMP-Acitve Fe ###
SMP_ActiveFe<-SterilePot_ActiveFe%>%filter(Treatment2=="SMP")
# 5.2.2.1 statistical analysis #
leveneTest(YL_ActiveFe ~ Treatment3, data = SMP_ActiveFe)#p>0.05，则满足方差齐性
shapiro.test(SMP_ActiveFe$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_SMP_ActiveFe<-aov(data=SMP_ActiveFe,YL_ActiveFe~Treatment3)
summary(aov_model_SMP_ActiveFe)
LSD_model_SMP_ActiveFe<-LSD.test(aov_model_SMP_ActiveFe,"Treatment3",p.adj = "BH")
LSD_model_SMP_ActiveFe

# 5.2.2.2 statistical analysis #
SMP_YL_ActiveFe_Bar<-ggplot(SMP_ActiveFe,aes(Treatment3,YL_ActiveFe))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.65,size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",
       y=expression('Active Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,7))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
SMP_YL_ActiveFe_Bar
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("SMP_YL_ActiveFe_Bar",".pdf",sep=""),
       SMP_YL_ActiveFe_Bar,device=cairo_pdf,width=50,height=72,dpi = 300,units = "mm")

## 5.3 SMP-biomass ###
SMP_Biomass<-SterilePot_Biomass%>%filter(Treatment2=="SMP")
# 5.3.1 statistical analysis #
leveneTest(Total ~ Treatment3, data = SMP_Biomass)#p>0.05，则满足方差齐性
shapiro.test(SMP_Biomass$Total)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_SMP_Biomass<-aov(data=SMP_Biomass,Total~Treatment3)
summary(aov_model_SMP_Biomass)
LSD_model_SMP_Biomass<-LSD.test(aov_model_SMP_Biomass,"Treatment3",p.adj = "BH")
LSD_model_SMP_Biomass

# 5.3.2 Plots #
SMP_Biomass_Bar<-ggplot(SMP_Biomass,aes(Treatment3,Total))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.65,size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",
       y="Biomas (g/plant)",parse =T)+
  scale_y_continuous(limits = c(0,15))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
SMP_Biomass_Bar
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("SMP_Biomass_Bar",".pdf",sep=""),
       SMP_Biomass_Bar,device=cairo_pdf,width=50,height=59,dpi = 300,units = "mm")

## 5.4 SMP-AvailableFe ###
SMP_AvailableFe<-SterilePot_AvailableFe%>%filter(Treatment2=="SMP")
# 5.4.1 statistical analysis #
leveneTest(availableFe ~ Treatment3, data = SMP_AvailableFe)#p>0.05，则满足方差齐性
shapiro.test(SMP_AvailableFe$availableFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_SMP_AvailableFe<-aov(data=SMP_AvailableFe,availableFe~Treatment3)
summary(aov_model_SMP_AvailableFe)
LSD_model_SMP_AvailableFe<-LSD.test(aov_model_SMP_AvailableFe,"Treatment3",p.adj = "BH")
LSD_model_SMP_AvailableFe

# 5.4.2 Plots #
SMP_AvailableFe_Bar<-ggplot(SMP_AvailableFe,aes(Treatment3,availableFe))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.65,size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",
       y=expression('Available Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,8))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
SMP_AvailableFe_Bar
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("SMP_AvailableFe_Bar",".pdf",sep=""),
       SMP_AvailableFe_Bar,device=cairo_pdf,width=50,height=62,dpi = 300,units = "mm")

#### 6 NMP ####
## 6.1 NMP-SPAD ###
NMP_SPAD<-SterilePotSPAD%>%filter(Treatment2=="NMP")
# 6.1.1 statistical analysis #
leveneTest(YL_SPAD ~ Treatment3, data = NMP_SPAD)#p>0.05，则满足方差齐性
shapiro.test(NMP_SPAD$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
NMP_SPAD<-NMP_SPAD%>%mutate(boxcox_YL_SPAD =BoxCox(NMP_SPAD$YL_SPAD,lambda="auto"))
leveneTest(boxcox_YL_SPAD ~ Treatment3, data = NMP_SPAD)#p>0.05，则满足方差齐性
shapiro.test(NMP_SPAD$boxcox_YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_NMP_SPAD<-aov(data=NMP_SPAD,boxcox_YL_SPAD~Treatment3)
summary(aov_model_NMP_SPAD)
#Non-parameter test
kruskal.test(YL_SPAD~Treatment3, data = NMP_SPAD)
aov_model_NMP_SPAD<-aov(data=NMP_SPAD,YL_SPAD~Treatment3)
dunnettT3Test(aov_model_NMP_SPAD,p.adjust.method = "BH")

# 6.1.2 Plots #
NMP_SPAD_Bar<-ggplot(NMP_SPAD,aes(Treatment3,YL_SPAD))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.65,size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",y='SPAD')+
  scale_y_continuous(limits = c(0,40))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
NMP_SPAD_Bar
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("NMP_SPAD_Bar",".pdf",sep=""),
       NMP_SPAD_Bar,device=cairo_pdf,width=50,height=72,dpi = 300,units = "mm")

## 6.2 NMP-Acitve Fe ###
NMP_ActiveFe<-SterilePot_ActiveFe%>%filter(Treatment2=="NMP")
# 6.2.1 statistical analysis #
leveneTest(YL_ActiveFe ~ Treatment3, data = NMP_ActiveFe)#p>0.05，则满足方差齐性
shapiro.test(NMP_ActiveFe$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
NMP_ActiveFe<-NMP_ActiveFe%>%mutate(boxcox_YL_ActiveFe =BoxCox(NMP_ActiveFe$YL_ActiveFe,lambda="auto"))
leveneTest(boxcox_YL_ActiveFe ~ Treatment3, data = NMP_ActiveFe)#p>0.05，则满足方差齐性
shapiro.test(NMP_ActiveFe$boxcox_YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_NMP_ActiveFe<-aov(data=NMP_ActiveFe,YL_ActiveFe~Treatment3)
summary(aov_model_NMP_ActiveFe)
LSD_model_NMP_ActiveFe<-LSD.test(aov_model_NMP_ActiveFe,"Treatment3",p.adj = "BH")
LSD_model_NMP_ActiveFe

#Non-parameter test
kruskal.test(YL_ActiveFe~Treatment3, data = NMP_ActiveFe)
aov_model_NMP_ActiveFe<-aov(data=NMP_ActiveFe,boxcox_YL_ActiveFe~Treatment3)
dunnettT3Test(aov_model_NMP_ActiveFe,p.adjust.method = "BH")

# 6.2.2 Plots #
NMP_YL_ActiveFe_Bar<-ggplot(NMP_ActiveFe,aes(Treatment3,YL_ActiveFe))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.65,size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",
       y=expression('Active Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,7))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
NMP_YL_ActiveFe_Bar
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("NMP_YL_ActiveFe_Bar",".pdf",sep=""),
       NMP_YL_ActiveFe_Bar,device=cairo_pdf,width=50,height=72,dpi = 300,units = "mm")

## 6.3 NMP-biomass ###
NMP_Biomass<-SterilePot_Biomass%>%filter(Treatment2=="NMP")
# 6.3.1 statistical analysis #
leveneTest(Total ~ Treatment3, data = NMP_Biomass)#p>0.05，则满足方差齐性
shapiro.test(NMP_Biomass$Total)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_NMP_Biomass<-aov(data=NMP_Biomass,Total~Treatment3)
summary(aov_model_NMP_Biomass)
LSD_model_NMP_Biomass<-LSD.test(aov_model_NMP_Biomass,"Treatment3",p.adj = "BH")
LSD_model_NMP_Biomass

# 6.3.2 Plots #
NMP_Biomass_Bar<-ggplot(NMP_Biomass,aes(Treatment3,Total))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.65,size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",
       y="Biomas (g/plant)",parse =T)+
  scale_y_continuous(limits = c(0,15))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
NMP_Biomass_Bar
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("NMP_Biomass_Bar",".pdf",sep=""),
       NMP_Biomass_Bar,device=cairo_pdf,width=50,height=59,dpi = 300,units = "mm")

## 6.4 NMP-AvailableFe ###
NMP_AvailableFe<-SterilePot_AvailableFe%>%filter(Treatment2=="NMP")
# 6.4.1 statistical analysis #
leveneTest(availableFe ~ Treatment3, data = NMP_AvailableFe)#p>0.05，则满足方差齐性
shapiro.test(NMP_AvailableFe$availableFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_NMP_AvailableFe<-aov(data=NMP_AvailableFe,availableFe~Treatment3)
summary(aov_model_NMP_AvailableFe)
LSD_model_NMP_ActiveFe<-LSD.test(aov_model_NMP_AvailableFe,"Treatment3",p.adj = "none")
LSD_model_NMP_ActiveFe

# 6.4.2 Plots #
NMP_AvailableFe_Bar<-ggplot(NMP_AvailableFe,aes(Treatment3,availableFe))+
  geom_bar(stat = "summary", fun = "mean",color="black",aes(fill=Treatment3),width=0.65,size=0.2)+
  stat_summary(fun=mean, geom='point',size=1)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  labs(x="",
       y=expression('Available Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,8))+
  scale_fill_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(fill="none")
NMP_AvailableFe_Bar
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("NMP_AvailableFe_Bar",".pdf",sep=""),
       NMP_AvailableFe_Bar,device=cairo_pdf,width=50,height=62,dpi = 300,units = "mm")

#### 7 Peanut-YL SPAD value-Field####
### 7.1 Beijing###
BJ_Field_SPAD<-Field_SPAD%>%filter(Position=="Beijing")
## 7.1.1 statistical analysis##
leveneTest(YL_SPAD ~ Treatment, data = BJ_Field_SPAD)#p>0.05，则满足方差齐性
shapiro.test(BJ_Field_SPAD$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_BJ_Field_SPAD<-aov(data=BJ_Field_SPAD,YL_SPAD~Treatment)
summary(aov_model_BJ_Field_SPAD)
LSD_model_BJ_Field_SPAD<-LSD.test(aov_model_BJ_Field_SPAD,"Treatment",p.adj = "BH")
LSD_model_BJ_Field_SPAD

# 7.1.2 Plot #
BJ_Field_SPAD_box_jitter<-ggplot(BJ_Field_SPAD,aes(Treatment,YL_SPAD))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,notch=T,size=0.2)+
  labs(x="",
       y="SPAD",parse =T)+
  scale_fill_manual(values=c("#BFBF4D","#F99F98","#4DC8F9","#45CB9D"))+
  scale_y_continuous(limits = c(6,45))+
  mytheme+
  guides(color=F,fill=F)
BJ_Field_SPAD_box_jitter
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("BJ_Field_SPAD_box_jitter",".pdf",sep=""),
       BJ_Field_SPAD_box_jitter,device=cairo_pdf,width=40,height=60,dpi = 300,units = "mm")

### 7.2 Puyang ###
PY_Field_SPAD<-Field_SPAD%>%filter(Position=="Puyang")
## 6.2.1 statistical analysis##
leveneTest(YL_SPAD ~ Treatment, data = PY_Field_SPAD)#p>0.05，则满足方差齐性
shapiro.test(PY_Field_SPAD$YL_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_PY_Field_SPAD<-aov(data=PY_Field_SPAD,YL_SPAD~Treatment)
summary(aov_model_PY_Field_SPAD)
LSD_model_PY_Field_SPAD<-LSD.test(aov_model_PY_Field_SPAD,"Treatment",p.adj = "BH")
LSD_model_PY_Field_SPAD
## 6.2.2 Plot##
PY_Field_SPAD_box_jitter<-ggplot(PY_Field_SPAD,aes(Treatment,YL_SPAD))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,notch=T,size=0.2)+
  labs(x="",
       y="SPAD",parse =T)+
  scale_fill_manual(values=c("#BFBF4D","#F99F98","#4DC8F9","#45CB9D"))+
  mytheme+
  scale_y_continuous(limits = c(6,45))+
  guides(color=F,fill=F)
PY_Field_SPAD_box_jitter
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("PY_Field_SPAD_box_jitter",".pdf",sep=""),
       PY_Field_SPAD_box_jitter,device=cairo_pdf,width=40,height=60,dpi = 300,units = "mm")

#### 8 Peanut-YL_ActiveFe-Field####
### 8.1 Beijing ###
BJ_Field_Iron<-Field_Iron%>%filter(Position=="Beijing")
## 8.1.1 statistical analysis##
leveneTest(YL_ActiveFe ~ Treatment, data = BJ_Field_Iron)#p>0.05，则满足方差齐性
shapiro.test(BJ_Field_Iron$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_BJ_Field_ActiveFe<-aov(data=BJ_Field_Iron,YL_ActiveFe~Treatment)
summary(aov_model_BJ_Field_ActiveFe)
LSD_model_BJ_Field_ActiveFe<-LSD.test(aov_model_BJ_Field_ActiveFe,"Treatment",p.adj = "BH")
LSD_model_BJ_Field_ActiveFe
## 8.1.2 Plot ##
BJ_Field_YL_ActiveFe_box_jitter<-ggplot(BJ_Field_Iron,aes(Treatment,YL_ActiveFe))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,notch=F,size=0.2)+
  labs(x="",
       y=expression('Active Fe (μg '*g^{-1}*')'),parse =T)+
  scale_fill_manual(values=c("#BFBF4D","#F99F98","#4DC8F9","#45CB9D"))+
  scale_y_continuous(limits = c(1,9),breaks = c(3,6,9))+
  mytheme+
  guides(color=F,fill=F)
BJ_Field_YL_ActiveFe_box_jitter
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("BJ_Field_YL_ActiveFe_box_jitter",".pdf",sep=""),
       BJ_Field_YL_ActiveFe_box_jitter,device=cairo_pdf,width=40,height=60,dpi = 300,units = "mm")

### 8.2 Puyang ###
PY_Field_Iron<-Field_Iron%>%filter(Position=="Puyang")
## 8.2.1 statistical analysis##
leveneTest(YL_ActiveFe ~ Treatment, data = PY_Field_Iron)#p>0.05，则满足方差齐性
shapiro.test(PY_Field_Iron$YL_ActiveFe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_PY_Field_ActiveFe<-aov(data=PY_Field_Iron,YL_ActiveFe~Treatment)
summary(aov_model_PY_Field_ActiveFe)
LSD_model_PY_Field_ActiveFe<-LSD.test(aov_model_PY_Field_ActiveFe,"Treatment",p.adj = "BH")
LSD_model_PY_Field_ActiveFe
## 8.2.2 Plot ##
PY_Field_YL_ActiveFe_box_jitter<-ggplot(PY_Field_Iron,aes(Treatment,YL_ActiveFe))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,notch=F,size=0.2)+
  labs(x="",
       y=expression('Available Fe (μg '*g^{-1}*')'),parse =T)+
  scale_fill_manual(values=c("#BFBF4D","#F99F98","#4DC8F9","#45CB9D"))+
  mytheme+
  scale_y_continuous(limits = c(1,9),breaks = c(3,6,9))+
  guides(color=F,fill=F)
PY_Field_YL_ActiveFe_box_jitter
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("PY_Field_YL_ActiveFe_box_jitter",".pdf",sep=""),
       PY_Field_YL_ActiveFe_box_jitter,device=cairo_pdf,width=40,height=60,dpi = 300,units = "mm")

#### 9 Peanut-AvailableFe-Field2020####
### 9.1 Beijing ###
BJ_Field_Iron<-Field_Iron[c(1:18),]
## 9.1.1 statistical analysis##
leveneTest(Available_Fe ~ Treatment, data = BJ_Field_Iron)#p>0.05，则满足方差齐性
shapiro.test(PY_Field_Iron$Available_Fe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_BJ_Field_AvailableFe<-aov(data=BJ_Field_Iron,Available_Fe~Treatment)
summary(aov_model_BJ_Field_AvailableFe)
LSD_model_PY_Field_ActiveFe<-LSD.test(aov_model_PY_Field_ActiveFe,"Treatment",p.adj = "BH")
LSD_model_PY_Field_ActiveFe
## 9.1.2 Plot ##
BJ_Field_Available_Fe_box_jitter<-ggplot(BJ_Field_Iron,aes(Treatment,Available_Fe))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,notch=F,size=0.2)+
  labs(x="",
       y=expression('Available Fe (μg '*g^{-1}*')'),parse =T)+
  scale_fill_manual(values=c("#BFBF4D","#F99F98","#4DC8F9","#45CB9D"))+
  mytheme+
  scale_y_continuous(limits = c(2,22),breaks = c(5,10,15))+
  guides(color=F,fill=F)
BJ_Field_Available_Fe_box_jitter
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("BJ_Field_Available_Fe_box_jitter",".pdf",sep=""),
       BJ_Field_Available_Fe_box_jitter,device=cairo_pdf,width=35,height=60,dpi = 300,units = "mm")

### 9.2 Puyang###
PY_Field_AvailableFe<-Field_Iron[c(25:48),]
## 9.2.1 statistical analysis ##
leveneTest(Available_Fe ~ Treatment, data = PY_Field_Iron)#p>0.05，则满足方差齐性
shapiro.test(PY_Field_Iron$Available_Fe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_PY_Field_AvailableFe<-aov(data=PY_Field_AvailableFe,Available_Fe~Treatment)
summary(aov_model_PY_Field_AvailableFe)
LSD_model_PY_Field_AvailableFe<-LSD.test(aov_model_PY_Field_AvailableFe,"Treatment",p.adj = "BH")
LSD_model_PY_Field_AvailableFe
## 9.2.2 Plot ##
PY_Field_Available_Fe_box_jitter<-ggplot(PY_Field_AvailableFe,aes(Treatment,Available_Fe))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,notch=F,size=0.2)+
  labs(x="",
       y=expression('Active Fe (μg '*g^{-1}*')'),parse =T)+
  scale_fill_manual(values=c("#BFBF4D","#F99F98","#4DC8F9","#45CB9D"))+
  mytheme+
  scale_y_continuous(limits = c(2,22),breaks = c(5,10,15))+
  guides(color=F,fill=F)
PY_Field_Available_Fe_box_jitter
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("PY_Field_Available_Fe_box_jitter",".pdf",sep=""),
       PY_Field_Available_Fe_box_jitter,device=cairo_pdf,width=35,height=60,dpi = 300,units = "mm")

#### 10 Peanut-Yield-Field####
### 10.1 Beijing ###
BJ_Field_Yield<-Field_Yield%>%filter(Position=="Beijing")
## 10.1.1 statistical analysis##
leveneTest(YieldPerHa ~ Treatment, data = BJ_Field_Yield)#p>0.05，则满足方差齐性
shapiro.test(BJ_Field_Yield$YieldPerHa)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
BJ_Field_Yield<-BJ_Field_Yield%>%mutate(boxcox_YieldPerHa =BoxCox(BJ_Field_Yield$YieldPerHa,lambda="auto"))
leveneTest(boxcox_YieldPerHa ~ Treatment, data = BJ_Field_Yield)#p>0.05，则满足方差齐性
shapiro.test(BJ_Field_Yield$boxcox_YieldPerHa)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_BJ_Field_Yield<-aov(data=BJ_Field_Yield,YieldPerHa~Treatment)
summary(aov_model_BJ_Field_Yield)
#Non-parameter test
kruskal.test(YieldPerHa~Treatment, data = BJ_Field_Yield)
aov_model_BJ_Field_Yield<-aov(data=BJ_Field_Yield,YieldPerHa~Treatment)
dunnettT3Test(aov_model_BJ_Field_Yield,p.adjust.method = "BH")

## 10.1.1 Plot ##
BJ_Field_Yield_box_jitter<-ggplot(BJ_Field_Yield,aes(Treatment,YieldPerHa))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,notch=F,size=0.2)+
  labs(x="",
       y=expression('The yield (t '*ha^{-1}*')'),parse =T)+
  scale_fill_manual(values=c("#BFBF4D","#F99F98","#4DC8F9","#45CB9D"))+
  mytheme+
  scale_y_continuous(limits = c(1.5,7),breaks = c(2,4,6))+
  guides(color=F,fill=F)
BJ_Field_Yield_box_jitter
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("BJ_Field_Yield_box_jitter",".pdf",sep=""),
       BJ_Field_Yield_box_jitter,device=cairo_pdf,width=40,height=60,dpi = 300,units = "mm")

### 10.2Puyang ###
PY_Field_Yield<-Field_Yield%>%filter(Position=="Puyang")
## 10.2.1 statistical analysis##
leveneTest(YieldPerHa ~ Treatment, data = PY_Field_Yield)#p>0.05，则满足方差齐性
shapiro.test(PY_Field_Yield$YieldPerHa)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_PY_Field_Yield<-aov(data=PY_Field_Yield,YieldPerHa~Treatment)
summary(aov_model_PY_Field_Yield)
LSD_model_PY_Field_Yield<-LSD.test(aov_model_PY_Field_Yield,"Treatment",p.adj = "BH")
LSD_model_PY_Field_Yield
## 10.2.2 Plot ##
PY_Field_Yield_box_jitter<-ggplot(PY_Field_Yield,aes(Treatment,YieldPerHa))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,notch=F,size=0.2)+
  labs(x="",
       y=expression('The yield (t '*ha^{-1}*')'),parse =T)+
  scale_fill_manual(values=c("#BFBF4D","#F99F98","#4DC8F9","#45CB9D"))+
  mytheme+
  scale_y_continuous(limits = c(1.5,5),breaks = c(2,4,6))+
  guides(color=F,fill=F)
PY_Field_Yield_box_jitter
setwd(wdOutput_Figure4)
getwd()
ggsave(paste("PY_Field_Yield_box_jitter",".pdf",sep=""),
       PY_Field_Yield_box_jitter,device=cairo_pdf,width=40,height=60,dpi = 300,units = "mm")
