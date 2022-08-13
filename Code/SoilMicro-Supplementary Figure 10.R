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
wdOutput_Figure5 <- c("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Figure5")

#### 3. Peanut-Pot and field ####
#### 3.1 Import and process data ####
setwd(wdImport)
NIP_SPAD <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                             sheet = "fig s10 SPAD of NIP")
NIP_SPAD$Treatment3<-factor(NIP_SPAD$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

NIP_Iron <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                                  sheet = "fig s10 iron of NIP")
NIP_Iron$Treatment3<-factor(NIP_Iron$Treatment3,levels=c("CK","1502IPR-01","Pyoverdine"))

#### 4 NIP ####
## 4.1 NIP-SPAD ###
# 4.1.1 statistical analysis #
compare_means(data=NIP_SPAD,YL_SPAD~Treatment3,method = "t.test")
aov_model_NIP_SPAD<-aov(data=NIP_SPAD,YL_SPAD~Treatment3)
summary(aov_model_NIP_SPAD)
duncan_model_NIP_SPAD<- duncan.test(aov_model_NIP_SPAD,"Treatment3")
duncan_model_NIP_SPAD

# 4.1.2 Plots #
NIP_SPAD_Bar<-ggplot(NIP_SPAD,aes(Treatment3,YL_SPAD))+
  geom_bar(stat = "summary", fun = "mean",color="black",fill="white",width=0.65,size=0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  geom_jitter(aes(Treatment3,color=Treatment3),width = 0.2, height = 0.2,size=0.25)+
  geom_signif(comparisons = list(c("CK", "1502IPR-01"), 
                                 c("1502IPR-01","Pyoverdine"), 
                                 c("CK","Pyoverdine")),
              textsize=2.5,step_increase = 0.15,size = 0.35,
              map_signif_level = T, #修改参数map_signif_level=TRUE
              test = t.test)+
  labs(x="",y='SPAD')+
  scale_y_continuous(limits = c(0,50))+
  scale_color_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(color=F)
NIP_SPAD_Bar
setwd(wdOutput_Figure5)
getwd()
ggsave(paste("NIP_SPAD_Bar",".pdf",sep=""),
       NIP_SPAD_Bar,device=cairo_pdf,width=40,height=72,dpi = 300,units = "mm")

## 4.2. NIP-Acitve Fe ###
# 4.2.1 statistical analysis #
compare_means(data=NIP_Iron,YL_ActiveFe~Treatment3,method = "t.test")
aov_model_NIP_ActiveFe<-aov(data=NIP_Iron,YL_ActiveFe~Treatment3)
summary(aov_model_NIP_ActiveFe)
duncan_result_NIP_ActiveFe<- duncan.test(aov_model_NIP_ActiveFe,"Treatment3")
duncan_result_NIP_ActiveFe

# 4.2.2 statistical analysis #
NIP_YL_ActiveFe_Bar<-ggplot(NIP_Iron,aes(Treatment3,YL_ActiveFe))+
  geom_bar(stat = "summary", fun = "mean",color="black",fill="white",width=0.65,size=0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  geom_jitter(aes(Treatment3,color=Treatment3),width = 0.2, height = 0.2,size=0.25)+
  geom_signif(comparisons = list(c("CK", "1502IPR-01"), 
                                 c("1502IPR-01","Pyoverdine"), 
                                 c("CK","Pyoverdine")),
              textsize=2.5,step_increase = 0.15,size = 0.35,
              map_signif_level = T, #修改参数map_signif_level=TRUE
              test = t.test)+
  labs(x="",
       y=expression('Active Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,15))+
  scale_color_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(color=F)
NIP_YL_ActiveFe_Bar
setwd(wdOutput_Figure5)
getwd()
ggsave(paste("NIP_YL_ActiveFe_Bar",".pdf",sep=""),
       NIP_YL_ActiveFe_Bar,device=cairo_pdf,width=40,height=60,dpi = 300,units = "mm")

## 4.3 NIP-AvailableFe ###
# 4.3.1 statistical analysis #
compare_means(data=NIP_Iron,availableFe~Treatment3,method = "t.test")
aov_model_NIP_AvailableFe<-aov(data=NIP_Iron,availableFe~Treatment3)
summary(aov_model_NIP_AvailableFe)
duncan_model_NIP_AvailableFe<- duncan.test(aov_model_NIP_AvailableFe,"Treatment3")
duncan_model_NIP_AvailableFe

# 4.3.2 Plots #
NIP_AvailableFe_Bar<-ggplot(NIP_Iron,aes(Treatment3,availableFe))+
  geom_bar(stat = "summary", fun = "mean",color="black",fill="white",width=0.65,size=0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.15,size=0.15)+
  geom_jitter(aes(Treatment3,color=Treatment3),width = 0.2, height = 0.2,size=0.25)+
  geom_signif(comparisons = list(c("CK", "1502IPR-01"), 
                                 c("1502IPR-01","Pyoverdine"), 
                                 c("CK","Pyoverdine")),
              textsize=2.5,step_increase = 0.15,size = 0.35,
              map_signif_level = T, #修改参数map_signif_level=TRUE
              test = t.test)+
  labs(x="",
       y=expression('Available Fe (μg '*g^{-1}*')'),parse =T)+
  scale_y_continuous(limits = c(0,8))+
  scale_color_manual(values=c("#BFBF4D", "#F99F98", "#4DC8F9"))+
  mytheme+
  guides(color=F)
NIP_AvailableFe_Bar
setwd(wdOutput_Figure5)
getwd()
ggsave(paste("NIP_AvailableFe_Bar",".pdf",sep=""),
       NIP_AvailableFe_Bar,device=cairo_pdf,width=50,height=60,dpi = 300,units = "mm")
