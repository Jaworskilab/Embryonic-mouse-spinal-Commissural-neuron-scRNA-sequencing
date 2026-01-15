library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(patchwork)
library(stringr)
library(SeuratWrappers)
library(monocle3)
library(Rmisc)
setwd("G:/Integrated Data")
#eall<-readRDS("eCombined_Object.RData")
eall<-readRDS("Final_Integrated_Object.RData")
Idents(eall)<-"seurat_clusters"
crossColors<-c("preCrossing"="red", "postCrossing"="blue","None"="lightgrey")
clustercolors<-ggthemes_data[["tableau"]][[
  "color-palettes"]][[
    "regular"]][[
      "Classic 20"]][[
        "value"]]
clustercolors<-clustercolors[c(1:14,17:19)]
domainColors<-c("dI1"="#F8766D", "dI2"="#CD9600",
                "dI4/dILa"="#7CAE00", "dI5/dILb"="#00BE67",
                "dI6"="#00BFC4", "v0"="#00A9FF",
                "v1"="#C77CFF", "v3"="#FF61CC")
patchpath<-"G:/Integrated Data/Figures/Patchwork"

#----MISC----
#r123counts<-readRDS("r123countsDT.RData")
r123counts<-readRDS("r123counts_ofFinalObject.RData")

robo3ROC <- function(obj) {
  r3pos<-WhichCells(obj, expression= Robo3>0)
  r3neg<-WhichCells(obj, expression= Robo3==0)
  
  markers<-FindMarkers(obj, ident.1=r3pos, ident.2=r3neg, test.use="roc", verbose=T)
  mDT<-as.data.table(markers, keep.rownames = T)
  colnames(mDT)[1]<-"gene"
  return(mDT)
}

#----Volcano Plot----
#Set Up Data Table
Idents(e12)<-"cross"
PPX.12<-FindMarkers(e12, ident.1 = "preCrossing" , ident.2 = "postCrossing")
PPX.12<-as.data.table(PPX.12, keep.rownames = T)
colnames(PPX.12)[1]<-"gene"

#Compute Y-Axis
PPX.12[,neglog10P:=-log10(p_val_adj)]

#Graph
ggplot(PPX.12, aes(x=avg_log2FC, y=neglog10P))+geom_point(alpha=0.5)

#----"Boomerang" Plot----
#Find percentile of 0
posmin<-abs(PPX.12[,min(avg_log2FC)])
max<-PPX.12[,max(avg_log2FC)]
pct0<-posmin/(posmin+max)
#Graph
PPX.12<-PPX.12[order(abs(avg_log2FC))] #order for point plotting
ggplot(PPX.12, aes(x=pct.1, y=pct.2))+
  geom_point(aes(color=avg_log2FC))+
  scale_color_gradientn(colors = c("blue","white","red"),
                        values = c(0, pct0-0.05,pct0+0.05, 1))+
  geom_text_repel(data=PPX.12[abs(pct.1-pct.2)>0.5], aes(label=gene),
                  size=3.5, max.overlaps = 15)+
  geom_text_repel(data=PPX.12[abs(avg_log2FC)>2 & abs(pct.1-pct.2)<0.5], aes(label=gene),
                  size=3.5, max.overlaps = 15)+
  labs(title="E12.5 Pre/Post-Crossing DEGs by Expression Frequency", 
       subtitle = "Labeled points have either |pct.1-pct.2| > 0.5 OR avg_log2FC>2",)+
  xlab("Proportion of Pre-Crossing Cells with >1 Transcript")+
  ylab("Proportion of Post-Crossing Cells with >1 Transcript")

p.bmg.12<-ggplot(PPX.12, aes(x=pct.1, y=pct.2))+
  geom_point(aes(color=avg_log2FC))+
  scale_color_gradientn(colors = c("blue","white","red"),
                        values = c(0, pct0-0.05,pct0+0.05, 1))+
  geom_text(data=PPX.12[abs(pct.1-pct.2)>0.5], aes(label=gene))+
  geom_text(data=PPX.12[abs(avg_log2FC)>2], aes(label=gene))

#----Assign Crossing Identity----
#E11.5
cross<-rep("None",length(rownames(e11@meta.data)))
#Filter based on 75th and 90th percentiles of Robo3 and Robo2
preCrossingCells<-WhichCells(e11, expression = Robo3>8 & Robo2<55, slot = "counts")
postCrossingCells<-WhichCells(e11, expression = Robo2>116 & Robo3<3, slot = "counts")
e11<-AddMetaData(e11, metadata = data.frame(cross=cross, row.names = rownames(e11@meta.data)))
e11@meta.data["cross"][which(rownames(e11@meta.data) %in% preCrossingCells),"cross"]<-"preCrossing"
e11@meta.data["cross"][which(rownames(e11@meta.data) %in% postCrossingCells), "cross"]<-"postCrossing"
FeatureScatter(e11, feature1 = "Robo3", feature2 = "Robo2", 
               slot = "counts", jitter = T, group.by = "cross")+
  geom_hline(yintercept = c(116, 55), linetype = "dashed")+
  geom_vline(xintercept = c(3,8), linetype = "dashed")+
  xlim(0,60)+
  scale_color_manual(values=c(None = "gray", preCrossing = "red", postCrossing = "blue"))
FeatureScatter(e11, feature1 = "PC_1", feature2 = "PC_2", group.by = "cross", 
               cells = c(preCrossingCells,postCrossingCells))+
  scale_color_manual(values = c(preCrossing="red", postCrossing="blue"))

Idents(e11)<-"cross"
PPX.11<-FindMarkers(e11, ident.1 = "preCrossing" , ident.2 = "postCrossing")
PPX.11<-as.data.table(PPX.11, keep.rownames = T)
colnames(PPX.11)[1]<-"gene"

posmin<-abs(PPX.11[,min(avg_log2FC)])
max<-PPX.11[,max(avg_log2FC)]
pct0<-posmin/(posmin+max)

PPX.11<-PPX.11[order(abs(avg_log2FC))] #order for point plotting
ggplot(PPX.11, aes(x=pct.1, y=pct.2))+
  geom_point(aes(color=avg_log2FC))+
  scale_color_gradientn(colors = c("blue","white","red"),
                        values = c(0, pct0-0.05,pct0+0.05, 1))+
  geom_text_repel(data=PPX.11[abs(pct.1-pct.2)>0.5], aes(label=gene), size=3.5, max.overlaps = 15)+
  geom_text_repel(data=PPX.11[abs(avg_log2FC)>2 & abs(pct.1-pct.2)<0.5], aes(label=gene), size=3.5, max.overlaps=15)+
  labs(title="E10.5 Pre/Post-Crossing DEGs by Expression Frequency", 
       subtitle = "Labeled points have either |pct.1-pct.2| > 0.5 OR avg_log2FC>2",)+
  xlab("Proportion of Pre-Crossing Cells with >1 Transcript")+
  ylab("Proportion of Post-Crossing Cells with >1 Transcript")

p.bmg.11<-ggplot(PPX.11, aes(x=pct.1, y=pct.2))+
  geom_point(aes(color=avg_log2FC))+
  scale_color_gradientn(colors = c("blue","white","red"),
                        values = c(0, pct0-0.05,pct0+0.05, 1))+
  geom_text_repel(data=PPX.11[abs(pct.1-pct.2)>0.5], aes(label=gene), size=3.5, max.overlaps = 15)+
  geom_text_repel(data=PPX.11[abs(avg_log2FC)>2 & abs(pct.1-pct.2)<0.5], aes(label=gene), size=3.5, max.overlaps=15)+
  labs(title="E11.5 Pre/Post-Crossing DEGs by Expression Frequency", 
       subtitle = "Labeled points have either |pct.1-pct.2| > 0.5 OR avg_log2FC>2",)+
  xlab("Proportion of Pre-Crossing Cells with >1 Transcript")+
  ylab("Proportion of Post-Crossing Cells with >1 Transcript")

#E10.5
cross<-rep("None",length(rownames(e10@meta.data)))
#Filter based on 75th and 90th percentiles of Robo3 and Robo2
preCrossingCells<-WhichCells(e10, expression = Robo3>14 & Robo2<2, slot = "counts")
postCrossingCells<-WhichCells(e10, expression = Robo2>7 & Robo3<4, slot = "counts")
e10<-AddMetaData(e10, metadata = data.frame(cross=cross, row.names = rownames(e10@meta.data)))
e10@meta.data["cross"][which(rownames(e10@meta.data) %in% preCrossingCells),"cross"]<-"preCrossing"
e10@meta.data["cross"][which(rownames(e10@meta.data) %in% postCrossingCells), "cross"]<-"postCrossing"
FeatureScatter(e10, feature1 = "Robo3", feature2 = "Robo2", 
               slot = "counts", jitter = T, group.by = "cross")+
  geom_hline(yintercept = c(7,2), linetype = "dashed")+
  geom_vline(xintercept = c(4,14), linetype = "dashed")+
  #xlim(0,60)+
  scale_color_manual(values=c(None = "gray", preCrossing = "red", postCrossing = "blue"))

FeatureScatter(e10, feature1 = "PC_1", feature2 = "PC_2", group.by = "cross", 
               cells = c(preCrossingCells,postCrossingCells))+
  scale_color_manual(values = c(preCrossing="red", postCrossing="blue"))

Idents(e10)<-"cross"
PPX.10<-FindMarkers(e10, ident.1 = "preCrossing" , ident.2 = "postCrossing")
PPX.10<-as.data.table(PPX.10, keep.rownames = T)
colnames(PPX.10)[1]<-"gene"

posmin<-abs(PPX.10[,min(avg_log2FC)])
max<-PPX.10[,max(avg_log2FC)]
pct0<-posmin/(posmin+max)

PPX.10<-PPX.10[order(abs(avg_log2FC))] #order for point plotting
ggplot(PPX.10, aes(x=pct.1, y=pct.2))+
  geom_point(aes(color=avg_log2FC))+
  scale_color_gradientn(colors = c("blue","white","red"),
                        values = c(0, pct0-0.05,pct0+0.05, 1))+
  geom_text_repel(data=PPX.10[abs(pct.1-pct.2)>0.5], 
                  aes(label=gene), 
                  size=3.5, max.overlaps = 15)+
  geom_text_repel(data=PPX.10[abs(avg_log2FC)>2 & abs(pct.1-pct.2)<0.5], 
                  aes(label=gene), 
                  size=3.5, max.overlaps = 15)+
  labs(title="E10.5 Pre/Post-Crossing DEGs by Expression Frequency", 
       subtitle = "Labeled points have either |pct.1-pct.2| > 0.5 OR avg_log2FC>2",)+
  xlab("Proportion of Pre-Crossing Cells with >1 Transcript")+
  ylab("Proportion of Post-Crossing Cells with >1 Transcript")

p.bmg.10<-ggplot(PPX.10, aes(x=pct.1, y=pct.2))+
  geom_point(aes(color=avg_log2FC))+
  scale_color_gradientn(colors = c("blue","white","red"),
                        values = c(0, pct0-0.05,pct0+0.05, 1))+
  geom_text_repel(data=PPX.10[abs(pct.1-pct.2)>0.5], 
                  aes(label=gene), 
                  size=3, max.overlaps = 15)+
  geom_text_repel(data=PPX.10[abs(avg_log2FC)>2 & abs(pct.1-pct.2)<0.5], 
                  aes(label=gene), 
                  size=3, max.overlaps = 15)+
  labs(title="E10.5 Pre/Post-Crossing DEGs by Expression Frequency", 
       subtitle = "Labeled points have either |pct.1-pct.2| > 0.5 OR avg_log2FC>2",)+
  xlab("Proportion of Pre-Crossing Cells with >1 Transcript")+
  ylab("Proportion of Post-Crossing Cells with >1 Transcript")

#E13.5
cross<-rep("None",length(rownames(e13@meta.data)))
#Filter based on 75th and 90th percentiles of Robo3 and Robo2
preCrossingCells<-WhichCells(e13, expression = Robo3>4 & Robo2<54, slot = "counts")
postCrossingCells<-WhichCells(e13, expression = Robo2>106 & Robo3==0, slot = "counts")
e13<-AddMetaData(e13, metadata = data.frame(cross=cross, row.names = rownames(e13@meta.data)))
e13@meta.data["cross"][which(rownames(e13@meta.data) %in% preCrossingCells),"cross"]<-"preCrossing"
e13@meta.data["cross"][which(rownames(e13@meta.data) %in% postCrossingCells), "cross"]<-"postCrossing"
FeatureScatter(e13, feature1 = "Robo3", feature2 = "Robo2", 
               slot = "counts", jitter = T, group.by = "cross")+
  geom_hline(yintercept = c(106,54), linetype = "dashed")+
  geom_vline(xintercept = c(0,4), linetype = "dashed")+
  #xlim(0,60)+
  scale_color_manual(values=c(None = "gray", preCrossing = "red", postCrossing = "blue"))
FeatureScatter(e13, feature1 = "PC_1", feature2 = "PC_2", group.by = "cross", 
               cells = c(preCrossingCells,postCrossingCells))+
  scale_color_manual(values = c(preCrossing="red", postCrossing="blue"))

Idents(e13)<-"cross"
PPX.13<-FindMarkers(e13, ident.1 = "preCrossing" , ident.2 = "postCrossing")
PPX.13<-as.data.table(PPX.13, keep.rownames = T)
colnames(PPX.13)[1]<-"gene"

posmin<-abs(PPX.13[,min(avg_log2FC)])
max<-PPX.13[,max(avg_log2FC)]
pct0<-posmin/(posmin+max)

PPX.13<-PPX.13[order(abs(avg_log2FC))] #order for point plotting
ggplot(PPX.13, aes(x=pct.1, y=pct.2))+
  geom_point(aes(color=avg_log2FC))+
  scale_color_gradientn(colors = c("blue","white","red"),
                        values = c(0, pct0-0.05,pct0+0.05, 1))+
  geom_text_repel(data=PPX.13[abs(pct.1-pct.2)>0.5], 
                  aes(label=gene), 
                  size=3, max.overlaps = 15)+
  geom_text_repel(data=PPX.13[abs(avg_log2FC)>2 & abs(pct.1-pct.2)<0.5], 
                  aes(label=gene), 
                  size=3, max.overlaps = 15)+
  labs(title="E13.5 Pre/Post-Crossing DEGs by Expression Frequency", 
       subtitle = "Labeled points have either |pct.1-pct.2| > 0.5 OR avg_log2FC>2",)+
  xlab("Proportion of Pre-Crossing Cells with >1 Transcript")+
  ylab("Proportion of Post-Crossing Cells with >1 Transcript")

p.bmg.13<-ggplot(PPX.13, aes(x=pct.1, y=pct.2))+
  geom_point(aes(color=avg_log2FC))+
  scale_color_gradientn(colors = c("blue","white","red"),
                        values = c(0, pct0-0.05,pct0+0.05, 1))+
  geom_text_repel(data=PPX.13[abs(pct.1-pct.2)>0.5], 
                  aes(label=gene), 
                  size=3, max.overlaps = 15)+
  geom_text_repel(data=PPX.13[abs(avg_log2FC)>2 & abs(pct.1-pct.2)<0.5], 
                  aes(label=gene), 
                  size=3, max.overlaps = 15)+
  labs(title="E13.5 Pre/Post-Crossing DEGs by Expression Frequency", 
       subtitle = "Labeled points have either |pct.1-pct.2| > 0.5 OR avg_log2FC>2",)+
  xlab("Proportion of Pre-Crossing Cells with >1 Transcript")+
  ylab("Proportion of Post-Crossing Cells with >1 Transcript")

#----Negative Binomial Thresholding----

#-----Pre/Post Thresholding-----
r123counts<-as.data.table(
  FetchData(eall, 
            vars = c("orig.ident","Robo1","Robo2","Robo3","Cntn2","L1cam","Gapdh"), 
            slot = "counts"),
  keep.rownames = T)

q10.robo3<-quantile(r123counts[Robo3>0 & orig.ident=="E10.5"]$Robo3, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
q10.robo2<-quantile(r123counts[Robo2>0 & orig.ident=="E10.5"]$Robo2, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
pre10<-r123counts[orig.ident == "E10.5" & Robo3>=q10.robo3["90%"] & Robo2<q10.robo2["75%"], rn]
post10<-r123counts[orig.ident == "E10.5" & Robo2>=q10.robo2["90%"] & Robo3<q10.robo3["75%"], rn]

q11.robo3<-quantile(r123counts[Robo3>0 & orig.ident=="E11.5"]$Robo3, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
q11.robo2<-quantile(r123counts[Robo2>0 & orig.ident=="E11.5"]$Robo2, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
pre11<-r123counts[orig.ident == "E11.5" & Robo3>=q11.robo3["90%"] & Robo2<q11.robo2["75%"], rn]
post11<-r123counts[orig.ident == "E11.5" & Robo2>=q11.robo2["90%"] & Robo3<q11.robo3["75%"], rn]

q12.robo3<-quantile(r123counts[Robo3>0 & orig.ident=="E12.5"]$Robo3, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
q12.robo2<-quantile(r123counts[Robo2>0 & orig.ident=="E12.5"]$Robo2, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
pre12<-r123counts[orig.ident == "E12.5" & Robo3>=q12.robo3["90%"] & Robo2<q12.robo2["75%"], rn]
post12<-r123counts[orig.ident == "E12.5" & Robo2>=q12.robo2["90%"] & Robo3<q12.robo3["75%"], rn]

q13.robo3<-quantile(r123counts[Robo3>0 & orig.ident=="E13.5"]$Robo3, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
q13.robo2<-quantile(r123counts[Robo2>0 & orig.ident=="E13.5"]$Robo2, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
pre13<-r123counts[orig.ident == "E13.5" & Robo3>=q13.robo3["90%"] & Robo2<q13.robo2["75%"], rn]
post13<-r123counts[orig.ident == "E13.5" & Robo2>=q13.robo2["90%"] & Robo3<q13.robo3["75%"], rn]

pre.all <- c(pre10, pre11, pre12, pre13)
post.all <- c(post10, post11, post12, post13)

r123counts<-r123counts[,cross:=rep("None",.N)]
r123counts[rn %in% pre.all, cross:="preCrossing"]
r123counts[rn %in% post.all, cross:="postCrossing"]
r123counts[,cross:=as.factor(cross)]

cross <- rep("None", length(rownames(eall@meta.data)))
eall <- AddMetaData(eall, metadata = data.frame(cross=cross, row.names = rownames(eall@meta.data)))
eall@meta.data["cross"][which(rownames(eall@meta.data) %in% pre.all), "cross"]<-"preCrossing"
eall@meta.data["cross"][which(rownames(eall@meta.data) %in% post.all), "cross"]<-"postCrossing"

eall@meta.data$cross<-factor(eall@meta.data$cross, levels = c("preCrossing", "postCrossing", "None"))

FeatureScatter(eall, cells = WhichCells(eall, expression = orig.ident=="E10.5"),
               feature1 = "Robo3", feature2 = "Robo2", slot = "counts", group.by = "cross", jitter = T)+
  scale_color_manual(values=c("red","blue","lightgrey"))+
  geom_hline(yintercept = c(q10.robo2["75%"], q10.robo2["90%"]), linetype="dashed")+
  geom_vline(xintercept = c(q10.robo3["75%"], q10.robo3["90%"]), linetype="dashed")+
  labs(title="E10.5")+NoLegend()

FeatureScatter(eall, cells = WhichCells(eall, expression = orig.ident=="E11.5"),
               feature1 = "Robo3", feature2 = "Robo2", slot = "counts", group.by = "cross", jitter = T)+
  scale_color_manual(values=c("red","blue","lightgrey"))+
  geom_hline(yintercept = c(q11.robo2["75%"], q11.robo2["90%"]), linetype="dashed")+
  geom_vline(xintercept = c(q11.robo3["75%"], q11.robo3["90%"]), linetype="dashed")+
  labs(title="E11.5")+NoLegend()

FeatureScatter(eall, cells = WhichCells(eall, expression = orig.ident=="E12.5"),
               feature1 = "Robo3", feature2 = "Robo2", slot = "counts", group.by = "cross", jitter = T)+
  scale_color_manual(values=c("red","blue","lightgrey"))+
  geom_hline(yintercept = c(q12.robo2["75%"], q12.robo2["90%"]), linetype="dashed")+
  geom_vline(xintercept = c(q12.robo3["75%"], q12.robo3["90%"]), linetype="dashed")+
  labs(title="E12.5")+NoLegend()

FeatureScatter(eall, cells = WhichCells(eall, expression = orig.ident=="E13.5"),
               feature1 = "Robo3", feature2 = "Robo2", slot = "counts", group.by = "cross", jitter = T)+
  scale_color_manual(values=c("red","blue","lightgrey"))+
  geom_hline(yintercept = c(q13.robo2["75%"], q13.robo2["90%"]), linetype="dashed")+
  geom_vline(xintercept = c(q13.robo3["75%"], q13.robo3["90%"]), linetype="dashed")+
  labs(title="E13.5")+NoLegend()

#Multivariate quantile method

multQ <- function(q,age,qr3=q){
  totalN <- r123counts[orig.ident==age,.N]
  totalNonZero <- r123counts[orig.ident==age & !(Robo3==0 & Robo2==0),.N]
  q.Robo3<-quantile(r123counts[orig.ident==age & Robo3>0]$Robo3, probs = qr3)
  q.Robo2<-quantile(r123counts[orig.ident==age & Robo2>0]$Robo2, probs = q)
  n.rightquad <- r123counts[orig.ident==age & Robo3>q.Robo3 & Robo2>q.Robo2,.N]
  resAll<-n.rightquad/totalN
  resNonZero<-n.rightquad/totalNonZero
  notes<-paste(c(paste("Age:",age), 
                 paste("Robo2 nonzero quantile=",q),
                 paste("Robo3 nonzero quantile=",qr3),
                 paste("nquad/ALL:", resAll),
                 paste("nquad/NON-DOUBLEZERO:", resNonZero)),
               collapse = "\n")
  ggplot(r123counts[orig.ident==age], aes(x=Robo3, y=Robo2, color=cross))+
    scale_color_manual(values = c("None"="lightgrey",
                                  "preCrossing"="red","postCrossing"="blue"))+
    geom_jitter(alpha=0.5)+
    geom_hline(yintercept = q.Robo2)+
    geom_vline(xintercept = q.Robo3)+
    theme_bw()+
    theme(legend.position = c(0.9, 0.5), 
          legend.background = element_rect(color="black"))+
    annotate(geom = "text", x=Inf, y=Inf,hjust=1, vjust=1,label=notes)
}
#E10.5___Robo3:0.1, Robo2:0.75
#E11.5___Robo3:0.75, Robo2:0.75
#E12.5___Robo3:0.6, Robo2:0.6
#E13.5___Robo3:0.75, Robo2:0.75
r123counts[,crossNew:=rep("None",.N)]

updateCross <- function(qR2, age, qR3) {
  r123counts<-r123counts[orig.ident==age & cross=="None", crossNew:="None"]
  q.Robo3<-quantile(r123counts[orig.ident==age & Robo3>0]$Robo3, probs=qR3)
  q.Robo2<-quantile(r123counts[orig.ident==age & Robo2>0]$Robo2, probs=qR2)
  newPreCells<-r123counts[cross=="None" & orig.ident==age & 
                            Robo3>=q.Robo3 & Robo2<q.Robo2, 
                          rn]
  newPostCells<-r123counts[cross=="None" & orig.ident==age &
                             Robo2>=q.Robo2 & Robo3<q.Robo3,
                           rn]
  newNoneCells<-r123counts[orig.ident==age & cross!="None" &
                             Robo3>=q.Robo3 & Robo2>=q.Robo2,
                           rn]
  print(length(newNoneCells))
  r123counts[rn %in% newPreCells, crossNew:="newPre"]
  r123counts[rn %in% newPostCells, crossNew:="newPost"]
  r123counts[rn %in% newNoneCells, crossNew:="None"]
}

eall<-AddMetaData(eall, 
                  metadata = data.frame(
                    crossNew=r123counts$crossNew, 
                    row.names = rownames(eall@meta.data)))

colnames(r123counts)[8]<-"cross.75.90"
r123counts[,Cross:=cross.75.90]
r123counts[crossNew=="newPre",Cross:="preCrossing"]
r123counts[crossNew=="newPost",Cross:="postCrossing"]
r123counts[crossNew=="None",Cross:="None"]

s.eall@meta.data$cross<-r123counts$Cross

#START HERE - ROBO1 postX
ggplot(r123counts,aes(x=Robo3, y=Robo1, color=Cross))+
  geom_jitter(alpha=0.4)+facet_wrap(~orig.ident, scales = "free")+
  scale_color_manual(values = c("None"="black",
                                "preCrossing"="red",
                                "postCrossing"="blue"))+
  theme_bw()

R1R3multQ <- function(q,age,qr3=q){
  totalN <- r123counts[orig.ident==age,.N]
  totalNonZero <- r123counts[orig.ident==age & !(Robo3==0 & Robo1==0),.N]
  q.Robo3<-quantile(r123counts[orig.ident==age & Robo3>0]$Robo3, probs = qr3)
  q.Robo1<-quantile(r123counts[orig.ident==age & Robo1>0]$Robo1, probs = q)
  n.rightquad <- r123counts[orig.ident==age & Robo3>q.Robo3 & Robo1>q.Robo1,.N]
  n.preXtotal<-r123counts[orig.ident==age & Cross=="preCrossing",.N]
  n.preXabove<-r123counts[orig.ident==age & Cross=="preCrossing" & Robo1>q.Robo1,.N]
  preXprop<-n.preXabove/n.preXtotal
  resAll<-n.rightquad/totalN
  resNonZero<-n.rightquad/totalNonZero
  notes<-paste(c(paste("Age:",age), 
                 paste("Robo1 nonzero quantile=",q),
                 paste("Robo3 nonzero quantile=",qr3),
                 paste("nquad/ALL:", resAll),
                 paste("nquad/NON-DOUBLEZERO:", resNonZero),
                 paste("preX above Robo1 cutoff:", preXprop)),
               collapse = "\n")
  ggplot(r123counts[orig.ident==age], aes(x=Robo3, y=Robo1, color=Cross))+
    scale_color_manual(values = c("None"="lightgrey",
                                  "preCrossing"="red","postCrossing"="blue"))+
    geom_jitter(alpha=0.5)+
    geom_hline(yintercept = q.Robo1)+
    geom_vline(xintercept = q.Robo3)+
    theme_bw()+
    theme(legend.position = c(0.9, 0.5), 
          legend.background = element_rect(color="black"))+
    annotate(geom = "text", x=Inf, y=Inf,hjust=1, vjust=1,label=notes)
}
#Robo1 nonzero quantiles (>95% of preX have lower Robo1 than cutoff):
#E10.5: 0.95
#E11.5: 0.96
#E12.5: 0.955
#E13.5: 0.86

R1R3updateCross <- function(qR1, age, qR3) {
  r123counts<-r123counts[orig.ident==age & Cross=="None", crossNew:="None"]
  q.Robo3<-quantile(r123counts[orig.ident==age & Robo3>0]$Robo3, probs=qR3)
  q.Robo1<-quantile(r123counts[orig.ident==age & Robo1>0]$Robo1, probs=qR1)
  R1highcells<-r123counts[orig.ident==age & crossNew=="None" &
                            Robo3<q.Robo3 & Robo1>q.Robo1, 
                          rn]
  r123counts[rn %in% R1highcells, crossNew:="R1high"]
}
#r123counts[crossNew=="R1high", Cross:="postCrossing"]

###----Subsetting Integrated Dataset----
#start with eCombined_Object.RData, exclude clusters 4, 17, 18
s.eall<-subset(eall, subset = seurat_clusters %in% c(0,1,2,3,5,6,7,8,9,10,11,12,13,14,15,16))
DefaultAssay(s.eall)<-"integrated"
s.eall<-ScaleData(s.eall)
s.eall<-RunPCA(s.eall, npcs = 30)
s.eall<-FindNeighbors(s.eall, dims = 1:20, reduction = "pca")
s.eall<-FindClusters(s.eall, resolution = 0.5)
s.eall<-RunUMAP(s.eall, dims = 1:20, reduction="pca") #run pre neighb/clust?


newLibID <- paste0(unlist(lapply(strsplit((as.character(s.eall@meta.data$orig.ident)),split = "[.]"),"[",1)),
                   paste0("-", s.eall@meta.data$libID))
newLibID<-as.factor(newLibID)
s.eall<-AddMetaData(s.eall, 
                    metadata = data.frame(embryo=newLibID, 
                                          row.names = rownames(s.eall@meta.data)))
saveRDS(s.eall, file = "Final_Integrated_Object.RData")
saveRDS(r123counts, file = "r123counts_ofFinalObject.RData")

DefaultAssay(eall)<-"RNA"
all.genes<-rownames(eall)
eall<-ScaleData(eall, features = all.genes, verbose = T)

###----Domain Assignment----
#----Combo rule method (Defunct)-----
dGenes<-c("Barhl1","Barhl2","Lhx2","Lhx9","Pou4f1",
          "Lhx1","Lhx5","Foxd3","Foxp2",
          "Tlx3","Isl1","Isl2","Prrxl1",
          "Lbx1","Pax2","Gbx1",
          "Tlx1","Phox2a","Lmx1b",
          "Dmrt3","Wt1","Bhlhe22",
          "Evx1","Evx2","Pitx2",
          "En1",
          "Vsx2","Sox14","Sox21",
          "Tal1","Gata2","Gata3",
          "Sox1",
          "Mnx1","Lhx3",
          "Sim1","Nkx2-2","Olig3")

ddt <- as.data.table(FetchData(eall, vars = c("orig.ident","seurat_clusters",
                                              "nCount_RNA","nFeature_RNA",dGenes)),
                     keep.rownames = T)
ddt[,"domain":=rep("None",.N)]
ddt[,"ndomain":=rep(0,.N)]

ddt[(Barhl1>0 | Barhl2>0 | Lhx2>0 | Lhx9>0),
    `:=` (domain = "dI1",
          ndomain = ndomain+1)]
ddt[Pax2==0 & Bhlhe22==0 (Foxd3>0 | Foxp2>0),
    `:=` (domain = "dI2",
          ndomain = ndomain+1)]
ddt[Isl1>0 & (Tlx3>0 | Prrxl1>0 | Pou4f1>0),
    `:=` (domain = "dI3",
          ndomain = ndomain+1)]
ddt[Lbx1>0 & Pax2>0 & Tlx3==0 & Dmrt3==0 & Wt1==0,
    `:=` (domain = "dI4",
          ndomain = ndomain+1)]
ddt[Lbx1>0 & Pax2==0 & Dmrt3==0 & Wt1==0 & Isl1==0 &
      (Tlx1>0 | Tlx3>0 | Lmx1b>0 | Phox2a>0 | Prrxl1>0 | Pou4f1>0),
    `:=` (domain = "dI5",
          ndomain = ndomain+1)]
ddt[Lbx1>0 & Pax2>0 & Tlx3==0 & (Dmrt3>0 | Wt1>0),
    `:=` (domain = "dI6",
          ndomain = ndomain+1)]
ddt[(Evx1>0 | Evx2>0 | Pitx2>0) & (Barhl1==0 & Barhl2==0 & Lhx2==0 & Lhx9==0),
    `:=` (domain = "v0",
          ndomain = ndomain+1)]
ddt[Lbx1==0 & Bhlhe22>0 & (Foxp2>0 | Foxd3>0 | En1>0),
    `:=` (domain = "v1",
          ndomain = ndomain+1)]
ddt[Vsx2>0 | Sox14>0 | Sox21>0 |
      Tal1>0 | Gata2>0 | Gata3>0 | Sox1>0,
    `:=` (domain = "v2",
          ndomain = ndomain+1)]
ddt[Sim1>0 | `Nkx2-2`>0,
    `:=` (domain = "v3",
          ndomain = ndomain+1)]

ddt[,domain:=as.factor(domain)]
#----Cluster-wise method----
domainDT<-as.data.table(FetchData(eall, 
                                  vars = c("seurat_clusters","Wt1")), 
                        keep.rownames = T)

domainDT[,domain:=rep("None",.N)]
domainDT[seurat_clusters %in% c(2,12), domain:="dI1"]
domainDT[seurat_clusters %in% c(5,14), domain:="dI2"]
domainDT[seurat_clusters %in% c(0,1,3,6,7,9), domain:="dI4/dILa"] #9a
domainDT[seurat_clusters %in% c(4,13,16), domain:="dI5/dILb"]
domainDT[seurat_clusters == 8 | (seurat_clusters==9 & Wt1>0), domain:="dI6"] #9b
domainDT[seurat_clusters == 10, domain:="v0"]
domainDT[seurat_clusters == 11, domain:="v1"]
domainDT[seurat_clusters == 15, domain:="v3"]
domainDT[,domain:=factor(domain, 
                         levels = c("dI1","dI2",
                                    "dI4/dILa","dI5/dILb",
                                    "dI6","v0","v1","v3"))]
eall<-AddMetaData(eall, metadata = data.frame(domain=domainDT$domain,
                                              row.names = rownames(eall@meta.data)))
saveRDS(eall, file = "Final_Integrated_Object.RData")