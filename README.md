---
title: "Microbiome analysis of peanut/maize intercropping"
author:
- Nanqi Wang
---

## Contents

[Overview](#overview)

[System Requirements](#requirements)

[Installation guide](#installation)

## Overview {#overview}

All software dependencies and operating systems using in 'Microbiome convergence enables siderophore-secreting-rhizobacteria to improve iron nutrition and yield of peanut intercropped with maize'

## System Requirements {#requirements}

### Hardware requirements

This analysis require a standard computer with enough RAM to support the in-memory operations.

### Software requirements

#### OS Requirements

This analysis is supported for *windows* and *Linux*, tested on the following systems:

Windows: 10

Linux: Ubuntu 20.10

#### R Package

```{bash}
 phyloseq=1.34.0 #microbiome analysis
 stats=4.2.2 #stats package
 car=3.1-1 #companion to applied regression
 forecast=8.2.0 ##box-cox数据变换 Data transformation
 ggpubr=0.4.0 #"ggplot2" based publication ready plots
 agricolae=1.3-5 #statistical procedures for agricultural research
 PMCMRplus=1.9.6 #calculate pairwise mutiple comparisons 
 rcompanion=2.4.21 #fucntions to support extension education proguram evaluaion
 ggplot2=3.4.0 #plotting
 vegan=2.5-7 # community ecology package
 psych=2.2.5 # correlation analysis
 ggplot2=3.4.0 #作图 plot
 dplyr=2.3.0 #数据清洗，Data cleaning
 plyr=1.8.8 #数据清洗，Data cleaning
 reshaped2=1.4.4 #数据清洗，Data cleaning
 ggthemes=4.2.4 #ggplot所用主题，Themes for ggplot2
 readxl=1.4.1 #读入 excel, read excel files
 ggsci=2.9 #配色，Themed color paletttes for "ggplto2"
 showtext=0.9-5 #字体设置, font setting
 extrafont=0.18 #Tools for using fonts
 sysfonts=0.8.8 #加载系统字体，loading the system fonts into R
 Cairo=1.6-0 #抗锯齿,anti-aliasing
 ggtree=3.4.4 #visualization of tree and annotation data
 stringr=1.5.0 #字符串处理string operations
 amplicon=1.14.2 # statistic and visualization for microbiome data
```

#### Linux Package

```{bash}
 QIIME2=2021.8 #microbiome analysis
 TBLASTN+=2.10.0 #Phylogenetic tree construction 
 MAFFT=7.487 #Multiple sequence alignment
 IQ-TREE2=2.1.4-beta #phylogenetic inference
```

# Installation Guide: {#installation}

### R Package installation

installation time: 1\~1.5 h

```{bash}
#phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE))
BiocManager::install("phyloseq") 
install.packages("devtools") #devtools
devtools::install_version("car", version = "3.1-1") # car
devtools::install_version("forecast", version = "8.2.0") # forecast
devtools::install_version("ggpubr", version = "0.4.0") # ggpubr
devtools::install_version("agricolae", version = "1.3-5") # agricolae
devtools::install_version("PMCMRplus", version = "1.9.6") # PMCMRplus
devtools::install_version("rcompanion", version = "2.4.21") #rcompanion
devtools::install_version("vegan", version = "2.5-7") #vegan
devtools::install_version("psych", version = "2.2.5") #psych
devtools::install_version("ggplot2", version = "3.4.0") #ggplot2
devtools::install_version("dplyr", version = "2.3.0") #dplyr
devtools::install_version("plyr", version = "1.8.8") #dplyr
devtools::install_version("reshaped2", version = "1.4.4") #dplyr
devtools::install_version("ggthemes", version = "4.2.4") #dplyr
devtools::install_version("readxl", version = "1.4.1") #dplyr
devtools::install_version("ggsci", version = "2.9") #dplyr
devtools::install_version("showtext", version = "0.9-5") #dplyr
devtools::install_version("extrafont", version = "0.18") #dplyr
devtools::install_version("sysfonts", version = "0.8.8") #dplyr
devtools::install_version("Cairo", version = "3.4.4") #dplyr
devtools::install_version("ggtree", version = "1.6-0") #dplyr
devtools::install_version("stringr", version = "1.5.0") #dplyr
devtools::install_version("amplicon", version = "1.14.2") #dplyr
```

### Linux Package installation

installation time: 0.5\~1 h

```{bash}
#conda install
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
#QIIME2
wget https://data.qiime2.org/distro/core/qiime2-2021.8-py38-linux-conda.yml
conda env create -n qiime2-2021.8 --file qiime2-2021.8-py38-linux-conda.yml
rm qiime2-2021.8-py38-linux-conda.yml
#srs for qiime2
conda activate qiime2-2021.8
pip instFourStages git+https://github.com/vitorheidrich/q2-srs.git
#blast
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/ncbi-blast-2.10.0+-x64-linux.tar.gz 
tar zxvpf ncbi-blast-2.10.0+-x64-linux.tar.gz
vim ~/.bashrc
export PATH=$PATH:$HOME/ncbi-blast-2.10.1+/bin:$PATH
source ~.bashrc
#MAFFT
wget https://mafft.cbrc.jp/alignment/software/mafft_7.487-1_amd64.deb
dpkg -i mafft_7.487-1_amd64.deb
exit
rehash (if necessary)
#IQ-TREE2
conda install -c bioconda iqtree=2.1.4
```

# 
