#### 1. loading required libraries #### 
library(dplyr)
library(ape)
library(iCAMP)
library(NST) # need to be NST >=3.0.3
# version 2020.8.23
rm(list=ls())
t0=Sys.time() # to calculate time cost
#### 2. set folder paths and file names ####
data.set.name <- "_jaccard_FourStages-Filtered1"# prefix in the output filenames
wdImport <- c("E:/Study/SCI/Soil Micro/SCI/Figures/16S/QIIME2-DAD2/FourStages-DAD2-F295+R205-Filtered1/Import")
wdOutput <- c("E:/Study/SCI/Soil Micro/SCI/Figures/16S/QIIME2-DAD2/FourStages-DAD2-F295+R205-Filtered1/Assembly/NST/tNST/jaccard")
getwd()
#### 3. key parameter setting ####
rand.time=1000  # randomization time, 1000 is usually enough. For example test, you may set as 100 or less to save time.
nworker=10 # nworker is thread number for parallel computing, which depends on the CPU core number of your computer.
memory.G=50 # to set the memory size as you need (but should be less than the available space in your hard disk), so that calculation of large tree will not be limited by physical memory. unit is Gb.

#### 4. load data and match IDs ####
setwd(wdImport)
ASVCount <- read.table("feature_table_FourStages_DAD2_F295_R205_Filtered1.txt",header=T,row.names = 1,na.strings = c("NA"))
Group <- read.table("FourStages_SampleData.txt",header=T,row.names = 1,na.strings = c("NA"))
ASVCountT <- data.frame(t(ASVCount))
library("dplyr")
taxonomy_FourStages_Filtered1 <- read.table("Sup_taxonomy_FourStages_DAD2_F295_R205_Filtered1.txt",sep=",",header=T,row.names = 1,fill=TRUE, na.strings = "")
dim(ASVCountT)
dim(Group)
dim(taxonomy_FourStages_Filtered1)

samp.ck=NST::match.name(rn.list=list(comm=ASVCountT,group=Group))
comm=samp.ck$comm
ASVCountT=ASVCountT[,colSums(ASVCountT)>0,drop=FALSE]
group=samp.ck$group

#### 5. Grouping way and metacommunity seting ####
groupSysSpeRhizo=group[,1,drop=FALSE]# if you have multiple ways to group samples, select one way each time.
groupSysSpeRhizo
prefixi=paste0(data.set.name,".SysSpeRhizo") # re-define the prefix in output filenames, specific to the grouping type.
meta.groupSysSpeRhizo=NULL # if treatment and control are from different metacommunities, you may set meta.groupi=groupi

#### 6 taxonomic NST ####
### 6.1 # PF # calculate tNST using jaccard distance ###
dist.method="jaccard" # "jaccard" and "jaccard" are preferred.
t1=Sys.time() # to count time cost
tnst=tNST(comm=comm, group=groupSysSpeRhizo, meta.group=meta.groupSysSpeRhizo, meta.com=NULL,
          dist.method=dist.method, abundance.weighted=TRUE, rand=rand.time,
          output.rand=TRUE, nworker=nworker, LB=FALSE, null.model="PF",
          between.group=FALSE, SES=TRUE, RC=TRUE)
View(tnst)
setwd(wdOutput)
getwd()
save(tnst,file = paste0("tNST",".PF",prefixi,".rda")) # save tNST output in R data format
write.table(tnst$index.grp,file = paste0("tNST.summary",".PF",prefixi,".csv"), quote = FALSE,sep = ",",row.names = FALSE)
write.table(tnst$index.pair.grp,file = paste0("tNST.pairwise",".PF",prefixi,".csv"),quote = FALSE,sep = ",",row.names = FALSE)
write.table(tnst$index.pair.grp,file = paste0("tNST.index.between",".PF",prefixi,".csv"),quote = FALSE,sep = ",",row.names = FALSE)
format(Sys.time()-t1)

## 6.2 Bootstrapping test ##
t1=Sys.time()
groupi=groupSysSpeRhizo # or you may specify a different grouping column from group maxtrix.

tnstbt=nst.boot(nst.result=tnst, group=groupi, rand=rand.time, trace=TRUE,
                two.tail=FALSE, out.detail=TRUE, between.group=FALSE, nworker=nworker)
View(tnstbt)
save(tnstbt,file = paste0("tNST.boot",".PF",prefixi,".rda"))
write.table(tnstbt$summary,file = paste0("tNST.boot.summary",".PF",prefixi,".csv"), quote = FALSE,sep = ",",row.names = FALSE)
write.table(tnstbt$compare,file = paste0("tNST.boot.compare",".PF",prefixi,".csv"), quote = FALSE,sep = ",",row.names = FALSE)
write.table(tnstbt$detail$NST.boot,file = paste0("tNST.boot.detail",".PF",prefixi,".csv"), quote = FALSE,sep = ",",row.names = FALSE)
(t=format(Sys.time()-t1))
getwd()
## 6.3 # PERMANOVA # this is alternative way for signficance test. If bootstrapping works well, PERMANOVA may not be necessary.
t1=Sys.time()
tnstpaov=nst.panova(nst.result=tnst, group = groupi, rand = rand.time, trace = TRUE)
write.table(tnstpaov,file = paste0("tNST.PERMANOVA",".PF",prefixi,".csv"), quote = FALSE,sep = ",",row.names = FALSE)
(t=format(Sys.time()-t1))
