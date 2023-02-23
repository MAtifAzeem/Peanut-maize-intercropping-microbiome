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
library(ggtree)
library(rcompanion)#nonparametric two ways test
library(psych)#correlation analysis
library(ape)

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

mytheme1 <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#4d4d4d",linewidth=0.2),
                              panel.background = element_rect(colour = "#4d4d4d",linewidth=0.2),
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
wdOutput_Figure3 <- ("E:/working/SCI/Soil Micro/SCI/Figures/Figures from R/Figure3")

#### 3. Fig. 3b-phylogenetic tree####
## 3.1 Import and process data ##
setwd(wdImport)
Rooted_tree_46Strains <- read.tree("Rooted_tree_46Strains_16s_sequence.nwk")
View(Rooted_tree_46Strains)
Rooted_tree_46Strains_Tib <-as_tibble(Rooted_tree_46Strains)
class(Rooted_tree_46Strains)
StrainScreeningName <- read_excel("Intercropping-microbiome-Data for submit.xlsx",sheet="fig 3b name")

## 3.2 getting Node parameters ##
Rooted_tree_46Strains_Node<-ggtree(Rooted_tree_46Strains,branch.length="none")+geom_text2(aes(subset=!isTip,label=node),hjust=-0.3)+
  geom_tiplab()
Rooted_tree_46Strains_Node
viewClade(Rooted_tree_46Strains_Node,node=66)##node66是Pseudomonas和其他分开的节点


## 3.3 Visualizing Phylogenetic Tree ##
Pse<-StrainScreeningName%>%filter(Genus=="Pseudomonas")
Bac<-StrainScreeningName%>%filter(Genus=="Bacillus")
Ach <-StrainScreeningName%>%filter(Genus=="Achromobacter")
Rhi <-StrainScreeningName%>%filter(Genus=="Rhizobium")
Ent <-StrainScreeningName%>%filter(Genus=="Enterobacter")
Ste <-StrainScreeningName%>%filter(Genus=="Stenotrophomonas")

Genus_Group <- list(Pseudomonas=Pse$Name,Bacillus=Bac$Name,
                    Achromobacter=Ach$Name,Rhizobium=Rhi$Name,
                    Enterobacter=Ent$Name,Stenotrophomonas=Ste$Name)
Genus_Group
Clade_Group<-as_tibble(Rooted_tree_46Strains) %>% groupOTU(Genus_Group)
Clade_Group
Tree_46Strains_LineGroup<-ggtree(Rooted_tree_46Strains, branch.length="none")%>%
  groupOTU(Genus_Group,'Genus')+
  aes(color=Genus)+ylim(NA,7)+
  geom_tiplab(as_ylab=TRUE,align=TRUE)+
  theme(legend.position=c(0.2,0.85))+
  scale_color_manual(values=c("black","#00A087","#3C5488","#4DBBD5","#E64B35","#F39B7F","#7E6148"))
Tree_46Strains_LineGroup
setwd(wdOutput_Figure4)
ggsave("Tree_46Strains_LineGroup.pdf",Tree_46Strains_LineGroup,device=cairo_pdf,width=140,height=140,dpi = 300,units = "mm")
Tree_46Strains<-ggtree(Rooted_tree_46Strains, branch.length="none")+ geom_tiplab(align=TRUE)
Tree_46Strains
View(Rooted_tree_46Strains)
## 3.4 Mapping data to the tree structure ##
setwd(wdImport)
SiderophoreProduction <- read_excel("Intercropping-microbiome-Data for submit.xlsx",sheet="fig 3b SiderophoreProduction")
SiderophoreProductionStat <-SiderophoreProduction%>%
  group_by(GroupID,label,Genus)%>%
  summarise_at("SiderophoreProductionValue",funs(mean,sd))
SiderophoreProductionStat
Sid_Stat<-data.frame(label=SiderophoreProductionStat$label,
                     value=SiderophoreProductionStat$mean)
Sid_Stat
class(Sid_Stat$label)
library(ggstance)
SiderophoreProduction_Bar <- ggplot(SiderophoreProductionStat, aes(mean,label))+
  scale_x_continuous(expand = c(0,0),breaks=c(0,20,40))+theme_tree2()+
  geom_errorbar(aes(xmin=mean-sd,
                    xmax=mean+sd),width=0.2)+geom_colh(fill="#CBCBCB")
SiderophoreProduction_Bar
## 3.5 Composite plots ##
library(aplot)
Rooted_tree_46Strains_Tib
Tree_46Strains_LineGroup
SiderophoreProduction_Bar
Tree_46Strains
Composite_plots_46Strains<-SiderophoreProduction_Bar %>% insert_left(Tree_46Strains_LineGroup,width=6)
Composite_plots_46Strains
setwd(wdOutput_Figure3)
ggsave("Composite_plots_46Strains.pdf",Composite_plots_46Strains,device=cairo_pdf,width=140,height=140,dpi = 300,units = "mm")
wdOutput

####4. Maize-ASV487-Percent-Pot2020####
## 4.1 Import and process data ##
setwd(wdImport)
ASV487_Percent <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                             sheet = "fig 3c")
ASV487_Percent_Maize<-subset(ASV487_Percent,Species=="Maize")
ASV487_Percent_Maize$System<-factor(ASV487_Percent_Maize$System,levels=c("Intercropping","Monocropping"))
## 4.2 two way Scheirer-Ray-Hare test ##

scheirerRayHare(Percent~System+Days+Days:System, data = ASV487_Percent_Maize)

## 4.3 ASV487 Percent plots ##
ASV487_Percent_Maize_Pots<-ggplot(ASV487_Percent_Maize, aes(x=Days, y=Percent, fill=System,group=System)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",size=0,
               alpha=I(.1))+
  stat_summary(fun.data="mean_cl_boot", geom="line",size=1,aes(color=System))+
  geom_point(aes(x=Days,y=Percent,color=System,group=System),size=0.25)+
  labs(x="",
       y="relative abundance(%)",parse =T)+
  scale_y_continuous(limits = c(0,0.8))+
  scale_x_continuous(breaks = c(46,53,63,73))+
  scale_color_manual(values = c("#E64B35","#4DBBD5"))+
  scale_fill_manual(values = c("#E64B35","#4DBBD5"))+
  mytheme1

ASV487_Percent_Maize_Pots
setwd(wdOutput_Figure3)
ggsave(paste("ASV487_Percent_Maize_Pots",".pdf",sep=""),ASV487_Percent_Maize_Pots,device=cairo_pdf,width=80,height=40,dpi = 300,units = "mm")

#### 5. Peanut-ASV487-Percent-Pot2020 ####
## 5.1 Import and process data ##
setwd(wdImport)
getwd()
ASV487_Percent <- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                             sheet = "fig 3c")
ASV487_Percent_Peanut<-subset(ASV487_Percent,Species=="Peanut")
ASV487_Percent_Peanut$System<-factor(ASV487_Percent_Peanut$System,levels=c("Intercropping","Monocropping"))

## 5.2 two way Scheirer-Ray-Hare test ##
scheirerRayHare(Percent~System+Days+Days:System, data = ASV487_Percent_Peanut)

## 5.3 ASV487 Percent Peanut Plots ##

ASV487_Percent_Peanut_Pots<-ggplot(ASV487_Percent_Peanut, aes(x=Days, y=Percent, fill=System,group=System)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",size=0,
               alpha=I(.1))+
  stat_summary(fun.data="mean_cl_boot", geom="line",size=1,aes(color=System))+
  geom_point(aes(x=Days,y=Percent,color=System,group=System),size=0.25)+
  labs(x="",
       y="relative abundance(%)",parse =T)+
  scale_y_continuous(limits = c(0,0.15))+
  scale_x_continuous(breaks = c(46,53,63,73))+
  scale_color_manual(values = c("#E64B35","#4DBBD5"))+
  scale_fill_manual(values = c("#E64B35","#4DBBD5"))+
  mytheme1

ASV487_Percent_Peanut_Pots
setwd(wdOutput_Figure3)
ggsave(paste("ASV487_Percent_Peanut_Pots",".pdf",sep=""),ASV487_Percent_Peanut_Pots,device=cairo_pdf,width=80,height=40,dpi = 300,units = "mm")

#### 6. Pyovdine solublize Fe ####
## 6.1 Import and process data ##
setwd(wdImport)
getwd()
Solubilize_Fe<- read_excel("Intercropping-microbiome-Data for submit.xlsx",
                             sheet = "fig 3f")


## 6.2 correlationship anlaysis ##
Solubilize_Fe_pyoverdine_Fe<-corr.test(Solubilize_Fe$Additon_sid_nmol,Solubilize_Fe$Iron_content_nmol,
                                               method="pearson",adjust="BH",minlength=5)
Solubilize_Fe_pyoverdine_Fe
Solubilize_Fe_pyoverdine_Fe$r
Solubilize_Fe_pyoverdine_Fe$p
Solubilize_Fe_pyoverdine_Fe$p.adj
## 6.3 Plots ##
Solubilize_Fe_line <- ggplot(Solubilize_Fe,aes(x=Additon_sid_nmol,y=Iron_content_nmol))+
  geom_point(color="#33a02c",size=0.5,shape=21)+
  stat_smooth(method = lm,color="#6a3d9a")+
  mytheme1
Solubilize_Fe_line
setwd(wdOutput_Figure3)
getwd()
ggsave(paste("Solubilize_Fe_line",".pdf",sep=""),
       Solubilize_Fe_line,device=cairo_pdf,width=120,height=40,dpi = 300,units = "mm")

