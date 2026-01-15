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

#|-----------------------------------------------------------------------------|
#|-----------------------------------FIGURES-----------------------------------|
#|-----------------------------------------------------------------------------|

#figpath<-"G:/Integrated Data/Figures"

#Dimplot by clusters
p1 <- DimPlot(eall)+scale_color_manual(values = clustercolors)+
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave("ClusterDimPlot.tiff", plot = p1, path = figpath, 
       width=5, height=5, units = "in", bg="transparent")

#N per replicate, grouped by age, colored by sex
ndt<-as.data.table(FetchData(eall, vars = c("orig.ident","embryo","sex")), keep.rownames = T)
ndt.summ<-ndt[,.N,by=embryo]
ndt.summ[,age:=paste0(substr(embryo,1,3),".5")]
sexVector<-c( #determined by Xist levels
  "M","F","F","M", #e10.5
  "F","M","M", #e11.5
  "M","M","F","F", #e12.5
  "F","M","F") #e13.5
ndt.summ[,sex:=sexVector]
ndt.summ[,age:=factor(age, levels = c("E10.5", "E11.5", "E12.5", "E13.5"))]
p2<-ggplot(ndt.summ, aes(x=age, y=N, fill=sex))+
  geom_dotplot(binaxis="y", stackdir = "center", dotsize = 1, stackratio = 2)+
  theme_pubr()+
  theme(axis.title = element_text(face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid.major = element_line(color = "darkgray", linetype = "dashed"),
        legend.position = "none"
        # legend.position = c(0.1,0.85),
        # legend.box.margin = margin(),
        # legend.background = element_blank(),
        # legend.box.background = element_rect(color = "black", fill="white"),
        #axis.text.x = element_text(angle=90, hjust = 0, vjust=0.5)
  )+
  xlab("Embryonic Age")+ylab("Number of Cells")

ggsave("Nbyage.tiff", plot = p2, path = figpath, width=4, height=5, units = "in", bg="transparent")

#Clusters by age (stacked bar chart by proportion)
clDT<-as.data.table(FetchData(eall, 
                              vars = c("seurat_clusters","orig.ident")),
                    keep.rownames = T)
setkey(clDT, seurat_clusters, orig.ident)
clDT.summ <- clDT[CJ(seurat_clusters, orig.ident, unique = T), .N, by = .EACHI]
p3 <- ggplot(clDT.summ, aes(x=seurat_clusters,
                            y=N,
                            fill=orig.ident))+
  geom_bar(stat = "identity",
           position = "fill", 
           alpha = 0.8,
           color = "black")+
  theme_pubr()+
  scale_y_continuous(expand = c(0,0))+
  xlab("Cluster")+ylab("Proportional Composition")+
  labs(fill = "Age")+
  theme(#legend.justification = c("right","top"),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA),
    axis.title = element_text(face = "bold"))

#Xist vln plots by embryo
p4a<- VlnPlot(eall, features = "Xist", group.by = "embryo", 
              pt.size = 0)+
  geom_vline(xintercept = c(4.5, 7.5, 11.5), linetype="longdash", linewidth=0.8)+
  #NoLegend()+
  theme_pubr()+
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, size=11),
        axis.title = element_text(face = "bold"))+
  xlab("Embryo")+ylab("Xist")+
  labs(title=NULL)+
  annotate(geom = "text", x=2.5, y=4, label="E10.5", fontface=2)+
  annotate(geom = "text", x=6, y=4, label="E11.5", fontface=2)+
  annotate(geom = "text", x=9.5, y=4, label="E12.5", fontface=2)+
  annotate(geom = "text", x=13, y=4, label="E13.5", fontface=2)
#Age-normalized sex per cluster
sexDT<-as.data.table(FetchData(eall, 
                               vars = c("sex","orig.ident","seurat_clusters")),
                     keep.rownames = T)
sexDT[,sex:=factor(sex, levels = c("M","F"))] #orig.ident is already factorized
summDT<-sexDT[,sum:=.N, 
              by=orig.ident][,prop:=.N, 
                             by=.(orig.ident, sex)][,prop:=prop/sum][,sum:=NULL]
setkey(summDT, seurat_clusters, orig.ident, sex)
countDT<-summDT[CJ(seurat_clusters, orig.ident, sex, unique = TRUE), .N, by=.EACHI]
colnames(countDT)[4]<-"clusterN.bysex.byage"
countDT[,clusterN.byage:=sum(clusterN.bysex.byage),by=.(seurat_clusters, orig.ident)]
countDT[,N.byage:=sum(clusterN.bysex.byage),by=orig.ident]
countDT[,N.bysex.byage:=sum(clusterN.bysex.byage), by=.(orig.ident, sex)]
countDT[,prop.bysex.byage:=N.bysex.byage/N.byage]
countDT[,clusterExp.bysex.byage:=clusterN.byage*prop.bysex.byage]
countDT[,clusterExp.bysex:=sum(clusterExp.bysex.byage),by=.(seurat_clusters, sex)]
countDT[,clusterObs.bysex:=sum(clusterN.bysex.byage), by=.(seurat_clusters, sex)]
finalDT<-unique(countDT[,c("seurat_clusters","sex","clusterExp.bysex","clusterObs.bysex")])
finalDT[,obsPerExp:=clusterObs.bysex/clusterExp.bysex]
finalDT[,obsOverobsPlusExp:=clusterObs.bysex/(clusterObs.bysex+clusterExp.bysex)]

totalM<-length(which(eall$sex=="M"))
totalF<-length(which(eall$sex=="F"))

findFisherExact<-function(clusterStr, mtotal=totalM, ftotal=totalF){
  expM<-finalDT[seurat_clusters==clusterStr & sex=="M", clusterExp.bysex]
  expF<-finalDT[seurat_clusters==clusterStr & sex=="F", clusterExp.bysex]
  
  obsM<-finalDT[seurat_clusters==clusterStr & sex=="M", clusterObs.bysex]
  obsF<-finalDT[seurat_clusters==clusterStr & sex=="F", clusterObs.bysex]
  
  nullOR<-((expM/totalM)/(1-(expM/totalM)))/((expF/totalF)/(1-(expF/totalF)))
  return(fisher.test(rbind(c(obsM, totalM-obsM), c(obsF, totalF-obsF)), 
                     or=nullOR))
}


clusters<-c(0:16)
res<-c()
for(i in clusters){
  FE<-findFisherExact(clusterStr = as.character(i))
  res<-c(res, FE)
}
sigDT<-data.table(cluster=c(0:16),p=res)
breaks<-c(0,0.001,0.01,0.05,1)
sigDT[,siglevel:=cut(p,breaks=breaks,labels=F)]
heights<-finalDT[sex=="M",obsPerExp]
sigDT[,heights:=heights]
starlabels <- c("***","**","*","")
sigDT[,star:=starlabels[siglevel]]

p4<-ggplot(finalDT[sex=="M"], aes(x=seurat_clusters, y=obsPerExp))+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0, ymax=1), fill=alpha("lightpink",0.3))+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=1, ymax=2), fill=alpha("lightblue1",0.3))+
  geom_bar(stat="identity", fill="transparent", color="black", linewidth=0.8)+
  geom_hline(yintercept = 1, linetype = "dashed", linewidth=0.8)+
  scale_y_continuous(expand=c(0,0), limits = c(0,2))+
  theme_pubr()+
  theme(panel.border = element_rect(color="black", linetype = "solid", fill=NA),
        axis.title = element_text(face = "bold"))+
  xlab("Cluster")+ylab("Sex Index")

p4.2<- p4a/p4 + plot_layout(heights = c(2,3))

ggsave("sexPlots.tiff", plot = p4.2, path = patchpath, width = 8, height = 6, units = "in")

p4star <- p4 + geom_text(data = sigDT, aes(x=clusters+1, y = heights+0.05, 
                                           label=star), fontface="bold")
ggsave("sexPlotsStar.tiff", plot = p4star, path = patchpath, width = 8, height = 3, units = "in")

# ggplot(finalDT[sex=="M"], aes(x=seurat_clusters, y=obsOverobsPlusExp))+
#   geom_rect(aes(xmin=-Inf, xmax=Inf,ymin=0,ymax=0.5), fill=alpha("lightpink",0.3))+
#   geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0.5, ymax=1), fill=alpha("lightblue1",0.3))+
#   geom_hline(yintercept = 0.5, linetype="longdash")+
#   geom_bar(stat="identity", fill="transparent", color="black", linewidth=1)+
#   ylim(0,1)+
#   theme_pubr()


#Neurotransmitter FeaturePlots and Bar chart
p5<-FeaturePlot(eall, features="Gad1", order = T)+
  NoLegend()+
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
p6<-FeaturePlot(eall, features="Slc6a5", order = T)+
  NoLegend()+
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
p7<-FeaturePlot(eall, features="Slc17a6", order = T)+
  NoLegend()+
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
p8<-FeaturePlot(eall, features="Slc18a3", order = T)+
  NoLegend()+
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave("Gad1FP.tiff", plot = p5, path = figpath, width = 5, height = 5, units = "in")
ggsave("Slc6a5FP.tiff", plot = p6, path = figpath, width = 5, height = 5, units = "in")
ggsave("Slc17a6FP.tiff", plot = p7, path = figpath, width = 5, height = 5, units = "in")
ggsave("Slc18a3FP.tiff", plot = p8, path = figpath, width = 5, height = 5, units = "in")

iClusters <- c(0,1,3,6,7,8,9,11)
eClusters <- c(2,4,5,10,12,13,14,15,16)
total13 <- length(WhichCells(eall, expression = orig.ident=="E13.5"))
iCount13 <- length(WhichCells(eall, expression = orig.ident=="E13.5" &
                                seurat_clusters %in% iClusters))
eCount13 <- length(WhichCells(eall, expression = orig.ident=="E13.5" &
                                seurat_clusters %in% eClusters))
iProp <- iCount13/total13
eProp <- eCount13/total13
ieDF <- data.frame("Type" = c("Inhibitory","Excitatory"),
                   "Proportion" = c(iProp, eProp),
                   "Count" = c(iCount13, eCount13))

p8a <- ggplot(ieDF, aes(x=Type, y=Proportion))+
  geom_bar(stat = "identity",
           color="black",
           fill="lightgrey",
           linewidth=1.5,
           width = 0.4)+
  theme_pubr()+
  scale_y_continuous(limits = c(0,1))+
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = "grey70", linetype = "dashed"),
        axis.title = element_text(face = "bold"))+
  xlab("Neurotransmitter Type") + ylab("Proportional Composition")
#Z/N rep per age
zndt<-as.data.table(FetchData(eall, vars = c("Zfhx4","Neurod6","orig.ident")), keep.rownames = T)
zndt.melt<-melt(zndt, id.vars = c("rn", "orig.ident"),
                measure.vars = c("Zfhx4","Neurod6"),
                variable.name = "gene",
                value.name = "Normalized Expression")
p9<-ggplot(zndt.melt, aes(x=orig.ident, y=`Normalized Expression`, fill=gene))+
  geom_boxplot()+
  scale_fill_manual(values = c(
    "Neurod6"=alpha("orange", 0.5), 
    "Zfhx4"=alpha("purple",0.5)))+
  geom_vline(xintercept = c(1.5,2.5,3.5), linetype="longdash")+
  xlab("Embryonic Age")+
  theme_pubr()+
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(face="bold"),
        panel.grid.major.y = element_line(linetype="dashed", color = "grey70"),
        legend.position = "none"
        # legend.title = element_blank(),
        # legend.position = c(0.85,0.885),
        # legend.background = element_blank(),
        # legend.box.background  = element_rect(color = "black")
  )
ggsave("ZNbyage.tiff", plot = p9, path = figpath, width=5, height=5, units = "in")

#Patchwork Formation
p2.8.9 <- p2 + p8a + p9 + plot_layout(widths = c(1,1,2))

p2.8.9.3 <- p2.8.9 / p3 + plot_layout(heights = c(2,1)) #SAVE MAX DIMS

ggsave("basicInfo.tiff", plot = p2.8.9.3, path = patchpath, width = 20, height = 12, units = "in")

#Cluster Markers
eall.markers <- FindAllMarkers(eall, only.pos = T, logfc.threshold = 1)
eall.markers %>%
  dplyr::group_by(cluster) %>%
  slice_head(n=5) %>%
  ungroup() -> top5
p.clusterheat <- DoHeatmap(eall, features = top5$gene)+NoLegend() #change colorbar to clustercolors!!
ggsave("ClusterHeat.tiff", plot = p.clusterheat, path = patchpath, width = 20, height = 12, units = "in")

#DimPlot by Domain
p10<-DimPlot(eall, group.by = "domain")+
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(title = NULL)
ggsave("DimPlotbyDomain.tiff", plot = p10, path=figpath, width = 5, height = 5, units = "in")

#Stacked Domain Prop by Age
dmDT <- as.data.table(FetchData(eall, vars = c("orig.ident", "domain")),
                      keep.rownames = T)
setkey(dmDT, orig.ident, domain)
dmDT.summ <- dmDT[CJ(orig.ident, domain, unique = T), .N, by = .EACHI]
dmDT.summ[,ageTotal := sum(N), by=orig.ident]
dmDT.summ[,prop := N/ageTotal]
dmDT.summ[,perc := prop*100]

dmDT.summ[,numcode := rep(0,.N)]
dmDT.summ[orig.ident=="E10.5", numcode:=1]
dmDT.summ[orig.ident=="E11.5", numcode:=2]
dmDT.summ[orig.ident=="E12.5", numcode:=3]
dmDT.summ[orig.ident=="E13.5", numcode:=4]

p.dmProp <- ggplot(dmDT.summ, aes(x=numcode, y=prop, fill=domain))+
  geom_area(color="black", alpha=0.6)+
  theme_pubr()+
  theme(legend.position = "right",
        panel.background = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        panel.grid.major.x = element_line(linetype="longdash", color = "grey70"),
        axis.title = element_text(face = "bold"))+
  scale_x_continuous(labels=c("1"="E10.5", "2"="E11.5", "3"="E12.5", "4"="E13.5"))+
  xlab("Embryonic Age") + ylab("Proportional Composition")

p.dmRaw <- ggplot(dmDT.summ, aes(x=numcode, y=N, fill=domain))+
  geom_area(color="black", alpha=0.6)+
  theme_pubr()+
  theme(legend.position = "none", #because of patchwork combo with dm.Prop!!
        panel.background = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        panel.grid.major.x = element_line(linetype="longdash", color = "grey70"),
        axis.title = element_text(face = "bold"))+
  scale_x_continuous(labels=c("1"="E10.5", "2"="E11.5", "3"="E12.5", "4"="E13.5"))+
  xlab("Embryonic Age") + ylab("Number of Cells")

p.dmBoth <- p.dmRaw + p.dmProp

ggsave("DomainComp.tiff", plot = p.dmBoth, path = patchpath, width = 8, height = 5, units = "in")

#Dotplot domain genes
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
p11<-DotPlot(eall, group.by = "domain", features = dGenes, scale = F)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  xlab("Gene")+ylab("Domain")+
  guides(color = guide_colorbar(title = "Mean Expr", order = 2),
         size = guide_legend(title = "Frequency", order = 1))
ggsave("DomainGenesDotPlot.tiff", plot = p11, path = patchpath, width = 11, height = 4, units = "in")

FPgenes <- c("Barhl1", "Barhl2", "Lhx2", "Lhx9", "Pou4f1",
             "Lhx1","Lhx5","Foxd3","Foxp2",
             "Lbx1","Pax2","Gbx1",
             "Tlx3","Phox2a","Lmx1b",
             "Dmrt3","Wt1","Bhlhe22",
             "Evx1","Evx2","Pitx2",
             "En1",
             "Sim1","Nkx2-2")

FPWrapper <- function(feature, plotname=feature){
  FeaturePlot(eall, features = feature, order = T)+
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())+
    labs(title = plotname)
}

p.barhl1 <- FPWrapper("Barhl1")
p.barhl2 <- FPWrapper("Barhl2")
p.lhx2 <- FPWrapper("Lhx2")
p.lhx9 <- FPWrapper("Lhx9")
p.brn3a <- FPWrapper("Pou4f1", plotname = "Brn3a")
p.lhx1 <- FPWrapper("Lhx1")
p.lhx5 <- FPWrapper("Lhx5")
p.foxd3 <- FPWrapper("Foxd3")
p.foxp2 <- FPWrapper("Foxp2")
p.lbx1 <- FPWrapper("Lbx1")
p.pax2 <- FPWrapper("Pax2")
p.gbx1 <- FPWrapper("Gbx1")
p.tlx3 <- FPWrapper("Tlx3")
p.phox2a <- FPWrapper("Phox2a")
p.lmx1b <- FPWrapper("Lmx1b")
p.dmrt3 <- FPWrapper("Dmrt3")
p.wt1 <- FPWrapper("Wt1")
p.bhlhe22 <- FPWrapper("Bhlhe22", plotname = "Bhlhb5")
p.evx1 <- FPWrapper("Evx1")
p.evx2 <- FPWrapper("Evx2")
p.pitx2 <- FPWrapper("Pitx2")
p.en1 <- FPWrapper("En1")
p.sim1 <- FPWrapper("Sim1")
p.nkx22 <- FPWrapper("Nkx2-2", plotname = "Nkx2.2")

p.allFP <- p.barhl1+p.barhl2+p.lhx2+p.lhx9+p.brn3a+p.lhx1+p.lhx5+p.foxd3+
  p.foxp2+p.lbx1+p.pax2+p.gbx1+p.tlx3+p.phox2a+p.lmx1b+p.dmrt3+p.wt1+
  p.bhlhe22+p.evx1+p.evx2+p.pitx2+p.en1+p.sim1+p.nkx22 +
  plot_layout(ncol = 6)

ggsave("alldomainFP.tiff", plot = p.allFP, path = patchpath, width = 15, height = 10, units = "in")



#FeaturePlot domain markers
p12<-FeaturePlot(eall, features = "Dmrt3", order = T)+NoLegend()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
p13<-FeaturePlot(eall, features = "Lmx1b", order = T)+NoLegend()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
p14<-FeaturePlot(eall, features = "Nkx2-2", order = T)+NoLegend()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
p15<-FeaturePlot(eall, features = "Evx1", order = T)+NoLegend()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
ggsave("Dmrt3FP.tiff", plot = p12, path = figpath, width = 3.5, height = 3.5, units = "in")
ggsave("Lmx1bFP.tiff", plot = p13, path = figpath, width = 3.5, height = 3.5, units = "in")
ggsave("Nkx2_2FP.tiff", plot = p14, path = figpath, width = 3.5, height = 3.5, units = "in")
ggsave("Evx1FP.tiff", plot = p15, path = figpath, width = 3.5, height = 3.5, units = "in")

#Domain Heatmap
eall@meta.data$subclusters <- eall@meta.data$seurat_clusters
c9acells<-WhichCells(eall, expression = seurat_clusters == 9 & Wt1>0)
eall@meta.data$subclusters <- as.character(eall@meta.data$subclusters)
#levels(eall@meta.data$subclusters) <- c("2","12","5","14","0","1","3","6","7","9",
#                                        "4","13","16","8","9a","10","11","15")
eall@meta.data$subclusters[row.names(eall@meta.data) %in% c9acells] <- "9a"
eall@meta.data$subclusters <- factor(eall@meta.data$subclusters,
                                     levels = c("2","12",
                                                "5","14",
                                                "0","1","3","6","7","9",
                                                "4","13","16",
                                                "8","9a",
                                                "10",
                                                "11",
                                                "15"))
Idents(eall) <- "domain"
domain.markers <- FindAllMarkers(eall, only.pos = T, logfc.threshold = 1)
domain.markers %>% 
  dplyr::group_by(cluster) %>%
  slice_head(n=10) %>%
  ungroup() -> top10
p.dmheat <- DoHeatmap(eall, features = top10$gene, label = F)
ggsave("DomainHeat.tiff", plot = p.dmheat, path = patchpath, width = 20, height = 12, units = "in")

#Domain breakdown heatmaps
Idents(eall) <- "subclusters"
dI1.markers <- FindMarkers(eall, ident.1 = "2", ident.2 = "12", logfc.threshold = 1)

#Z-N by domain

#di2 vs v1 vlnplot
p16<-VlnPlot(eall, features = "Pou4f1", idents = c(5,11,14), pt.size = 0)+
  xlab("Cluster")+ylab("Expression")
p17<-VlnPlot(eall, features = "Bhlhe22", idents = c(5,11,14), pt.size = 0)+
  xlab("Cluster")+ylab("Expression")
p18<-VlnPlot(eall, features = "Pax2", idents = c(5,11,14), pt.size = 0)+
  xlab("Cluster")+ylab("Expression")
ggsave("Pou4f1Vln.tiff", plot = p16, path = figpath, width = 6, height = 2.5, units = "in")
ggsave("Bhlhe22Vln.tiff", plot = p17, path = figpath, width = 6, height = 2.5, units = "in")
ggsave("Pax2Vln.tiff", plot = p18, path = figpath, width = 6, height = 2.5, units = "in")

#Robo3 vs Robo2 scatterplots
r123dt<-as.data.table(FetchData(eall, 
                                vars = c("Robo1","Robo2","Robo3","cross","crossNew","orig.ident"), 
                                slot = "counts"),
                      keep.rownames = T)

multQFinal <- function(age,qR2,qR3=qR2){
  # totalN <- r123dt[orig.ident==age,.N]
  # totalNonZero <- r123dt[orig.ident==age & !(Robo3==0 & Robo2==0),.N]
  q.Robo3<-quantile(r123dt[orig.ident==age & Robo3>0]$Robo3, probs = qR3)
  q.Robo2<-quantile(r123dt[orig.ident==age & Robo2>0]$Robo2, probs = qR2)
  # n.rightquad <- r123dt[orig.ident==age & Robo3>q.Robo3 & Robo2>q.Robo2,.N]
  # resAll<-n.rightquad/totalN
  # resNonZero<-n.rightquad/totalNonZero
  #notes<-paste(c(paste("Age:",age), 
  #               paste("Robo2 nonzero quantile=",q),
  #               paste("Robo3 nonzero quantile=",qr3),
  #               paste("nquad/ALL:", resAll),
  #               paste("nquad/NON-DOUBLEZERO:", resNonZero)),
  #             collapse = "\n")
  robo2String <- paste0(as.character(qR2*100), "th percentile")
  robo3String <- paste0(as.character(qR3*100), "th percentile")
  maxrobo2 <- r123dt[orig.ident==age, max(Robo2)]
  maxrobo3 <- r123dt[orig.ident==age, max(Robo3)]
  p<-ggplot(r123dt[orig.ident==age], aes(x=Robo3, y=Robo2, color=crossNew))+
    scale_color_manual(values = c("None"="lightgrey",
                                  "preCrossing"="red","postCrossing"="blue"))+
    geom_jitter(alpha=0.5)+
    geom_hline(yintercept = q.Robo2, linetype="longdash", color="grey35", linewidth=0.8)+
    geom_vline(xintercept = q.Robo3, linetype="longdash", color="grey35", linewidth=0.8)+
    annotate(geom = "text", x=maxrobo3, y=q.Robo2, 
             label=robo2String, hjust="inward", vjust=0, color="grey35", fontface=3)+
    annotate(geom = "text", x=q.Robo3, y=maxrobo2,
             label=robo3String, hjust="inward", vjust=1, color="grey35", fontface=3,
             angle=90)+
    theme_pubr()+
    ggtitle(age)+
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          panel.border = element_rect(fill=NA, color="black"),
          panel.grid.major = element_line(linetype = "dotted", colour = "grey70"),
          plot.title = element_text(margin = margin(t=10, b= -20), 
                                    hjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold"))+
    guides(color=guide_legend(title=NULL))
  
  return(p)
}

updateCrossFinal <- function(age, qR2, qR3) {
  r123dt[orig.ident==age, crossNew:="None"]
  q.Robo3<-quantile(r123dt[orig.ident==age & Robo3>0]$Robo3, probs=qR3)
  #print(q.Robo3)
  q.Robo2<-quantile(r123dt[orig.ident==age & Robo2>0]$Robo2, probs=qR2)
  #print(q.Robo2)
  r123dt[orig.ident==age & Robo3>=q.Robo3, crossNew:="preCrossing"]
  r123dt[orig.ident==age & Robo2>=q.Robo2, crossNew:="postCrossing"]
  r123dt[orig.ident==age & Robo3>=q.Robo3 & Robo2>=q.Robo2, crossNew:="None"]
}
#FINAL CROSSING ASSIGNMENTS
updateCrossFinal("E10.5", qR2 = 0.75, qR3 = 0.1)
updateCrossFinal("E11.5", qR2 = 0.75, qR3 = 0.75)
updateCrossFinal("E12.5", qR2 = 0.6, qR3 = 0.6)
updateCrossFinal("E13.5", qR2 = 0.75, qR3 = 0.75)
eall$cross<-r123dt$crossNew

p.105 <- multQFinal("E10.5", 0.75, 0.1)
p.115 <- multQFinal("E11.5", 0.75)
p.125 <- multQFinal("E12.5", 0.6)
p.135 <- multQFinal("E13.5", 0.75)

p.scatter <- p.105 + p.115 + p.125 + p.135 + 
  plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "none")

ggsave(filename = "crossScatter.tiff", plot = p.scatter, path = patchpath, 
       width = 6, height = 6, units = "in")

# p19<-multQPlot(0.75, "E10.5", 0.1)
# p20<-multQPlot(0.75, "E11.5")
# p21<-multQPlot(0.6, "E12.5")
# p22<-multQPlot(0.75, "E13.5")

p.Scatter <- p19 + p20 + p21 + p22 + plot_layout(ncol=2, guides = "collect") 
ggsave("Robo3_2Scatters.tiff", plot = p.Scatter, path = patchpath, width = 10, height = 8, units = "in")


ggsave("e10r3r2.tiff", plot = p19, path = figpath, width = 4, height = 4, units = "in")
ggsave("e11r3r2.tiff", plot = p20, path = figpath, width = 4, height = 4, units = "in")
ggsave("e12r3r2.tiff", plot = p21, path = figpath, width = 4, height = 4, units = "in")
ggsave("e13r3r2.tiff", plot = p22, path = figpath, width = 4, height = 4, units = "in")

#Cross Counts
dt.crosscount <- as.data.table(FetchData(eall, 
                                         vars = c("orig.ident","cross")),
                               keep.rownames = T)
dt.crosscount <- dt.crosscount[cross != "None"]
dt.cc.summ <- dt.crosscount[,.N, by=.(orig.ident, cross)]
dt.cc.summ[,cross:=factor(cross, levels = c("preCrossing","postCrossing"))]

p.cc.age <- ggplot(dt.cc.summ, aes(x=orig.ident, y=N, group=cross, fill=cross))+
  geom_bar(stat="identity", position = "stack", width = 0.6, alpha=0.8, color="black")+
  scale_fill_manual(values = c("preCrossing"="red","postCrossing"="blue"))+
  theme_pubr()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, color="black"),
        panel.grid.major.y = element_line(linetype="longdash", color="grey70"),
        panel.grid.minor.y = element_line(linetype="dashed", color="grey80"))+
  geom_vline(xintercept = c(1.5,2.5,3.5), linetype="longdash", color="grey35")+
  ylab("Cell Count")

dt.cc.total <- dt.crosscount[,.N, by=cross]
dt.cc.total[,cross:=factor(cross,levels=c("preCrossing","postCrossing"))]
p.cc.total <- ggplot(dt.cc.total, aes(x=cross,y=N,fill=cross))+
  geom_bar(stat="identity", position="dodge", width=0.6, alpha=0.8, color="black")+
  scale_fill_manual(values = c("preCrossing"="red", "postCrossing"="blue"))+
  theme_pubr()+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, color="black"),
        panel.grid.major.y = element_line(linetype="longdash", color="grey70"),
        panel.grid.minor.y = element_line(linetype="dashed", color="grey80"))+
  ylab("Cell Count")

p.cc.both <- p.cc.age + p.cc.total + plot_layout(ncol=2, widths = c(3,1))
ggsave(filename = "CrossCounts.tiff", plot = p.cc.both, path = patchpath, 
       width = 11, height = 3, units = "in")

#Tag1/L1cam/Gapdh by cross
Idents(eall)<-"cross"
Idents(eall)<-factor(Idents(eall), levels = c("preCrossing","postCrossing","None",
                                              "newPre","newPost","R1high"))
p23<-VlnPlot(eall, features = c("Cntn2","L1cam","Gapdh"), 
             idents = c("preCrossing","postCrossing"), pt.size = 0, cols = c(alpha("red",0.8),alpha("blue",0.8)))
ggsave("Tag1L1vln.tiff", plot = p23, path = figpath, width = 9, height = 4, units = "in")

Idents(eall) <- "cross"
p.tag1 <- VlnPlot(eall, features = "Cntn2", idents = c("preCrossing", "postCrossing"), 
                  pt.size = 0, cols = crossColors[c(1,2)])+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        legend.position = "none")+
  ylab("Cntn2")+
  labs(title = NULL)

p.l1 <- VlnPlot(eall, features = "L1cam", idents = c("preCrossing", "postCrossing"), 
                pt.size = 0, cols = crossColors[c(1,2)])+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        legend.position = "none")+
  ylab("L1cam")+
  labs(title = NULL)

p.gapdh <- VlnPlot(eall, features = "Gapdh", idents = c("preCrossing", "postCrossing"), 
                   pt.size = 0, cols = crossColors[c(1,2)])+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        legend.position = "none")+
  ylab("Gapdh")+
  labs(title = NULL)

p.pc2 <- VlnPlot(eall, features = "PC_2", idents = c("preCrossing", "postCrossing"), 
                 pt.size = 0, cols = crossColors[c(1,2)])+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        legend.position = "none")+
  ylab("PC 2")+
  labs(title = NULL)

p.crossvlns <- p.tag1 + p.l1 + p.gapdh + p.pc2 + plot_layout(ncol = 4)
ggsave("CrossVlns.tiff", plot = p.crossvlns, path = patchpath, width = 16, height = 5, units = "in")

#Cross Vlns Prettier
dt.cl <- as.data.table(FetchData(eall, vars = c("Cntn2", "L1cam", "Gapdh", "PC_2", "cross")),
                       keep.rownames = T)
dt.cl <- dt.cl[cross != "None"]
dt.cl[,cross:=factor(cross, levels = c("preCrossing", "postCrossing"))]


p.cntn2 <- ggplot(dt.cl, aes(x=cross, y=Cntn2, fill=cross))+
  geom_violin(scale="width", alpha=0.8, linewidth=0.8)+
  scale_fill_manual(values = c("preCrossing"="red", "postCrossing"="blue"))+
  theme_pubr()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        panel.grid.major.y = element_line(linetype = "longdash", color = "grey70"),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, color="black"))+
  geom_vline(xintercept = 1.5, linetype="longdash", color="grey70")

p.l1 <- ggplot(dt.cl, aes(x=cross, y=L1cam, fill=cross))+
  geom_violin(scale="width", alpha=0.8, linewidth=0.8)+
  scale_fill_manual(values = c("preCrossing"="red", "postCrossing"="blue"))+
  theme_pubr()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        panel.grid.major.y = element_line(linetype = "longdash", color = "grey70"),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, color="black"))+
  geom_vline(xintercept = 1.5, linetype="longdash", color="grey70")

p.gapdh <- ggplot(dt.cl, aes(x=cross, y=Gapdh, fill=cross))+
  geom_violin(scale="width", alpha=0.8, linewidth=0.8)+
  scale_fill_manual(values = c("preCrossing"="red", "postCrossing"="blue"))+
  theme_pubr()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        panel.grid.major.y = element_line(linetype = "longdash", color = "grey70"),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, color="black"))+
  geom_vline(xintercept = 1.5, linetype="longdash", color="grey70")

p.pc2 <- ggplot(dt.cl, aes(x=cross, y=PC_2, fill=cross))+
  geom_violin(scale="width", alpha=0.8, linewidth=0.8)+
  scale_fill_manual(values = c("preCrossing"="red", "postCrossing"="blue"))+
  theme_pubr()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        panel.grid.major.y = element_line(linetype = "longdash", color = "grey70"),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, color="black"))+
  geom_vline(xintercept = 1.5, linetype="longdash", color="grey70")+
  ylab("PC 2")

p.combo <- p.cntn2 + p.l1 + p.gapdh + p.pc2 +
  plot_layout(ncol=4, guides = "collect") & theme(legend.position = "none")

ggsave(filename = "crossVlns.tiff", plot = p.combo, path = patchpath,
       width = 15, height = 4, units = "in")

#DimPlot by cross
p24<-DimPlot(eall, cols = crossColors, order = c("preCrossing","postCrossing"))+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
ggsave("DimPlotbyCross.tiff", plot = p24, path = figpath, width = 6, height = 6, units = "in")

#PC2 by cross
p25<-VlnPlot(eall, features = c("PC_2"), idents = c("preCrossing","postCrossing"), cols = crossColors, pt.size = 0)+
  NoLegend()+
  labs(title=NULL)+xlab("")+ylab("PC_2")
ggsave("PC2cross.tiff", plot = p25, path = figpath, width = 4, height = 6, units = "in")

#Crossing Markers USING ROC
Idents(eall) <- "cross"
cross.markers <- FindMarkers(eall, ident.1 = "preCrossing", ident.2 = "postCrossing",
                             test.use = "roc")
cross.markers <- as.data.table(cross.markers, keep.rownames = T)

saveRDS(cross.markers, file = "prePostROCMarkersDT.RData")
cross.markers <- readRDS(file = "prePostROCMarkersDT.RData")

# top20post <- head(cross.markers[avg_log2FC<0],20)$rn
# top20pre <- head(cross.markers[avg_log2FC>0],20)$rn
# dt.cross <- as.data.table(FetchData(eall, 
#                                     vars = c("orig.ident","seurat_clusters",
#                                              "cross","domain")),
#                           keep.rownames = T)
# dt.cross <- dt.cross[cross != "None"]
# dt.cross[,cross:=factor(cross, levels = c("preCrossing","postCrossing"))]
# setorder(dt.cross, cross, orig.ident)
# prePostCells <- dt.cross$rn
# # prePostCells <- WhichCells(eall, expression = crossNew %in% c("preCrossing","newPre",
# #                                                               "postCrossing","newPost","R1high")) #figure out why cross doesn't work
# eall@active.ident <- factor(eall@active.ident, levels = c("preCrossing","postCrossing","None"))
# p.crossHeat <- DoHeatmap(eall, features = c(top20pre,top20post), cells = prePostCells,
#                          label = F, group.colors = crossColors) +
#   guides(color="none")
# 
# ggsave("CrossHeat.tiff", plot = p.crossHeat, path = patchpath, width = 15, height = 9, units = "in")

#Lollipop charts of select Pre/Post genes
# robo1
# robo2
# robo3
# nrp2
# plxna1
# plxna2
# plxna3
# plxnb1
# plxnc1
# epha5
# epha8
# epha10
# ephb3
# ephb6
# dcc
# unc5a
# unc5d
# ptch2
# cntn2
# dscam
# l1cam
# ncam1
# ncam2
# ncam
geneList <- c("Robo1","Robo2","Robo3",
              "Slit1","Slit2",
              "Nell1","Nell2",
              "Nrp1","Nrp2",
              "Plxna1","Plxna2","Plxna3","Plxna4","Plxnb1","Plxnc1",
              "Epha3","Epha4","Epha5","Epha6","Epha7","Epha8","Epha10",
              "Ephb1","Ephb2","Ephb3","Ephb6",
              "Dcc","Unc5a","Unc5d","Neo1", "Rgmb",
              "Ptch2","Ptchd4",
              "Fzd2","Fzd3",
              "Cntn1","Cntn2","Cntn3","Cntn4","Cntn5","Cntnap2","Cntnap4","Cntnap5a","Cntnap5b",
              "L1cam","Ncam1","Ncam2","Ncam",
              "Chl1",
              "Nhlh1","Nhlh2")
#ADD: RGMB, NDFIPs, SLITs, NELLs, PTCHD1(?), SRGAP2, SRGAP3, SPSB4, DPYSL4, CHL1, CDHs, PCDHs
#REMOVE: UNCX, TENMs, SNCs

cc.curated <- cross.markers[rn %in% geneList]
cc.curated[,rank:=rep(0,.N)]
cc.curated[,rank:=match(rn, geneList)]
setorder(cc.curated, rank)
cc.curated[,rn:=factor(rn, levels = rev(rn))]
cc.curated[,group:=ifelse(avg_log2FC>0, "pre", "post")]

myLabels <- c(levels(cc.curated$rn)[1:36], "*Robo3","*Robo2","Robo1")

p.curated <- ggplot(cc.curated, aes(x=rn, y=avg_log2FC, fill=group))+
  geom_rect(xmin=0.5, xmax=2.5, ymin=-Inf, ymax=Inf, fill=alpha("gray90",0.9), color=NA)+
  geom_rect(xmin=5.5, xmax=14.5, ymin=-Inf, ymax=Inf, fill=alpha("gray90",0.9), color=NA)+
  geom_rect(xmin=17.5, xmax=21.5, ymin=-Inf, ymax=Inf, fill=alpha("gray90",0.9), color=NA)+
  geom_rect(xmin=28.5, xmax=30.5, ymin=-Inf, ymax=Inf, fill=alpha("gray90",0.9), color=NA)+
  geom_rect(xmin=32.5, xmax=34.5, ymin=-Inf, ymax=Inf, fill=alpha("gray90",0.9), color=NA)+
  geom_rect(xmin=36.5, xmax=39.5, ymin=-Inf, ymax=Inf, fill=alpha("gray90",0.9), color=NA)+
  geom_hline(yintercept = 0, linewidth=1)+
  geom_vline(xintercept = c(1:39), linetype="dashed")+
  geom_hline(yintercept = c(-3,-2,-1,1,2), linetype="longdash")+
  geom_segment(aes(x=rn, xend=rn, y=0, yend=avg_log2FC), color="black", linewidth=0.8)+
  geom_point(size=3.5,shape=21,stroke=1)+
  coord_flip()+
  theme_pubr()+
  theme(#panel.grid.major.y = element_line(linetype = "dashed", color = "black"),
    #panel.grid.major.x = element_line(linetype = "longdash", color = "grey35"),
    #panel.grid.minor.x = element_line(linetype = "longdash", colour = "grey35"),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1.5),
    axis.title = element_text(face="bold"))+
  scale_fill_manual(values = c("pre"="red","post"="blue"))+
  scale_x_discrete(labels=myLabels)+
  xlab("Gene") + ylab("Average Log2 Fold Change") +
  labs(subtitle = "Guidance Genes of Interest")

ggsave(filename = "CuratedLollipop.tiff", plot = p.curated, path = patchpath,
       width=10, height = 10, units = "in")
# p.prepostcombo <- p.crossHeat / p.curated 
# ggsave(filename = "prepostCombo.tiff", plot=p.prepostcombo, path = patchpath, 
#        width=10, height=20, units = "in")

#Top TFs/CSMs Pre/Post
#(Using cross.markers generated above)
preTFs <- c("Tcf4","Nhlh2", "St18", "Prox1", "Tcf12", "Uncx", "Scrt2", "Pou3f2")
preCRs <-  c("Chd7","Cecr2","Smarca5","Chd4", "Id2", "Zfp536", "Dpf3", "Sobp", "Scmh1")
preCSMs <- c("Robo3","Cntn2","Rgmb","Dcc","Thsd7a","Lrrn1", "Tmeff1", "Igdcc3", "Ddr1")

postTFs <- c("Tshz2", "Rorb")
postCRs <- c("Chd5")
postCSMs <- c("Robo2", "Pcdh7", "Lingo2", "Chl1", "Grin2b", "Tmeff2", "Lsamp",
              "Csmd2", "Adgrb3", "Cntn4", "Dpp10", "Nrxn1", "Cd200", "Jph4", "Tmem130",
              "Lrrc4c", "Negr1", "Neto2", "Tenm3", "Ntm", "Cdh8", "Ncam2", "Csmd3")

TFCSM.genes <- c(preTFs, preCRs, preCSMs, postTFs, postCRs, postCSMs)
dt.tfcsm <- cross.markers[rn %in% TFCSM.genes]
setorder(dt.tfcsm, -avg_log2FC)
dt.tfcsm[,cross:=ifelse(avg_log2FC>0, "pre", "post")]
dt.tfcsm[,rn:=factor(rn, levels = rev(rn))]

p.topTFs <- ggplot(dt.tfcsm[rn %in% c(preTFs, preCRs, postTFs, postCRs)], 
                   aes(x=rn, y=avg_log2FC, fill=cross))+
  geom_hline(yintercept = 0, linewidth=1)+
  geom_hline(yintercept = c(-2,-1,1,2), linetype="longdash")+
  geom_segment(aes(x=rn, xend=rn, y=0, yend=avg_log2FC), color="black", linewidth=0.8)+
  geom_point(size=3.5,shape=21,stroke=1)+
  ylim(-2,2)+
  coord_flip()+
  theme_pubr()+
  theme(panel.grid.major.y = element_line(linetype = "dashed", color = "black"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1.5),
        axis.title = element_text(face="bold"))+
  scale_fill_manual(values = c("pre"="red","post"="blue"))+
  xlab("Gene") + ylab("Average Log2 Fold Change")+
  labs(subtitle = "Transcription Factors and Transcriptional Regulators")

dt.csm <- dt.tfcsm[rn %in% c(preCSMs, postCSMs)]
dt.csm[,rn:=factor(rn, levels = rev(rn))]
myLabels <- c("*Robo2", levels(dt.csm$rn)[2:31], "*Robo3")

p.topCSMs <- ggplot(dt.csm, 
                    aes(x=rn, y=avg_log2FC, fill=cross))+
  geom_hline(yintercept = 0, linewidth=1)+
  geom_hline(yintercept = c(-2,-1,1,2), linetype="longdash")+
  geom_segment(aes(x=rn, xend=rn, y=0, yend=avg_log2FC), color="black", linewidth=0.8)+
  geom_point(size=3.5,shape=21,stroke=1)+
  ylim(-4,3)+
  coord_flip()+
  theme_pubr()+
  theme(panel.grid.major.y = element_line(linetype = "dashed", color = "black"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1.5),
        axis.title = element_text(face="bold"))+
  scale_fill_manual(values = c("pre"="red","post"="blue"))+
  xlab("Gene") + ylab("Average Log2 Fold Change")+
  labs(subtitle = "Cell Surface Molecules") +
  scale_x_discrete(labels=myLabels)

p.TFCSM <- p.topTFs / p.topCSMs + plot_layout(heights = c(2,3))
p.allLollipop <- p.curated | p.TFCSM
ggsave(filename = "Lollipops.tiff", plot = p.allLollipop, path = patchpath,
       width = 15, height = 12, units = "in")

#Other genes of interest to talk about:
#Srgap2, Srgap3, Dpysl4, Spsb4, Igfbpl1, Nav2, Elavl4 (5'UTR isoforms not distinguishable),
#Igf2bp1 (local translation in guidance), Phf21b (promotes neuronal differentiation), Msi2,
#Reep3, Rere (nuclear receptor coregulator-mediated Tx regulation), Ppfia2 (Liprin/LAR),
#Dok6 (axon maintenance, transport, survival), Rbfox1, Arpp21 (post-Tx reg, dendritic branching),
#Rbms3, Tnik, Mab21l2, Khdrbs2, Celf5, Mtss1

#Proportion Pre/Post by Domain
dt.pp <- as.data.table(FetchData(eall, 
                                 vars = c("cross","orig.ident","domain")),
                       keep.rownames = T)
dt.pp <- dt.pp[cross!="None"]
dt.pp[,cross:=factor(cross, levels = c("preCrossing","postCrossing"))]
setkey(dt.pp, orig.ident, domain, cross)
dt.pp.summ <- dt.pp[CJ(orig.ident, domain, cross, unique = T), .N, by=.EACHI]
dt.pp.summ[,total:=sum(N), by=.(orig.ident, domain)]
dt.pp.summ <- dt.pp.summ[cross=="preCrossing"]
dt.pp.summ[,pct.pre:=(N/total)]
dt.pp.summ[,domain:=factor(domain, levels=rev(levels(domain)))]
setorder(dt.pp.summ, domain, orig.ident)

p.dotcross <- ggplot(dt.pp.summ, aes(x=orig.ident, y=domain, fill=100*pct.pre))+
  geom_point(shape=21, size=10, stroke=1) +
  scale_fill_gradient2(low="blue", mid = "white", high = "red", midpoint = 50)+
  theme_pubr()+
  theme(axis.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, color="black", linewidth = 1),
        legend.position = "right",
        panel.grid.major = element_line(linetype = "longdash", color = "grey35"))+
  xlab("Embryonic Age") + ylab("Domain")+
  guides(fill=guide_colorbar(title = "Percent \nPre-Crossing"))

#Z/N by Domain
dt.zn <- as.data.table(FetchData(eall, 
                                 vars = c("domain","Neurod6", "orig.ident")),
                       keep.rownames = T)
#dt.zn[Neurod6>0,n6n:=.N, by=domain]
dt.zn.summ <- dt.zn[Neurod6>0 & orig.ident=="E13.5",.N,by=domain]
setorder(dt.zn.summ, domain)
setorder(dt.zn, domain)
dt.zn.summ$total <- dt.zn[,.N,by=domain]$N
dt.zn.summ[,pct.N:=(N/total)]

p.zn <- ggplot(dt.zn.summ, aes(x=pct.N*100, y=domain, fill=domain)) +
  #geom_point(shape=21,fill="lightgrey",stroke=1,size=8)+
  geom_bar(stat = "identity", position = "dodge",
           linewidth=1, color="black")+
  xlim(0,50)+
  theme_pubr()+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA, color = "black", linewidth = 1),
        panel.grid.major = element_line(linetype = "longdash", color = "grey35"),
        panel.grid.minor.y = element_line(linetype = "dashed", color = "grey70"),
        axis.title = element_text(face = "bold"),
        legend.position = "none")+
  coord_flip()+
  xlab("Percent Neurod6+") + ylab("Domain") +
  scale_fill_manual(values = domainColors)

p.crossdot.zn <- (p.dotcross | p.zn) + plot_layout(widths = c(1,2))
ggsave(filename = "CrossDotZN.tiff", plot = p.crossdot.zn, path = patchpath,
       width = 12, height = 5, units = "in")


#Hox Genes
hoxgenes <- c("Hoxa1","Hoxb1","Hoxd1",
              "Hoxa2","Hoxb2",
              "Hoxa3","Hoxb3","Hoxb3os","Hoxd3","Hoxd3os1",
              "Hoxa4","Hoxb4","Hoxc4","Hoxd4","Hoxd4.1",
              "Hoxa5","Hoxb5","Hoxb5os","Hoxc5",
              "Hoxa6","Hoxb6","Hoxc6",
              "Hoxa7","Hoxb7",
              "Hoxb8","Hoxc8","Hoxd8",
              "Hoxa9","Hoxb9","Hoxc9","Hoxd9",
              "Hoxa10","Hoxc10","Hoxd10",
              "Hoxd11",
              "Hoxd13")

#Domain Subpopulation Dotplots (USING ROC)
#DI1
di1.markers <- as.data.table(FindMarkers(eall, 
                                         ident.1 = "2",
                                         ident.2 = "12",
                                         test.use = "roc"),
                             keep.rownames = T)
di1.markers[,diff.pct:=abs(pct.1-pct.2)]
setorder(di1.markers, -diff.pct)
di1.genes <- c("Prag1","Bmpr1b","Chmp2b","Sema3a","Ctnna3",
               "Evx1", "Cbln1","Kcnk9","Synpr","Vstm2l","Grik1","Asic2","Scrn1",
               "Tpbg","Hotairm1") #removed Pcsk5 bc global marker
di1.markers[,Cluster:=ifelse(avg_log2FC>0, 2, 12)]
di1.markers[,maxpct:=pmax(pct.1, pct.2)]
di1.order <- setorder(di1.markers[rn %in% di1.genes], Cluster, -maxpct)
di1.genes <- rev(di1.order$rn)
#Final list in order:
di1.genes <- c("Kcnk9","Grik1","Evx1","Tpbg","Vstm2l","Asic2","Synpr",
               "Cbln1","Hotairm1","Scrn1","Sema3a","Ctnna3","Chmp2b",
               "Bmpr1b","Prag1") 

p.di1 <- DotPlot(eall, idents = c(2,12), features = di1.genes, scale = F, dot.scale = 8)+
  coord_flip() +
  theme(#aspect.ratio = 5,
    panel.grid.major = element_line(linetype = "dashed", color = "grey70"),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(face = "bold"))+
  ylab("Cluster")+
  xlab("Gene")+
  labs(subtitle = "dI1 subpopulations")+
  guides(color = guide_colorbar("Mean Expr", order = 2),
         size = guide_legend("Frequency", order = 1))+
  scale_color_gradient(low = "grey", high = "blue", limits=c(0,3.3))+
  scale_size_continuous(limits = c(0,100),breaks = c(25,50,75))

#DI2
di2.markers <- as.data.table(FindMarkers(eall,
                                         ident.1 = "5",
                                         ident.2 = "14",
                                         test.use = "roc"),
                             keep.rownames = T)
di2.markers[,diff.pct:=abs(pct.1-pct.2)]
setorder(di2.markers, -diff.pct)
di2.genes <- c("Ebf2","St18","Cachd1","Cbln2","Bmpr1b","Lrrc4c","Tenm1",
               "Mast4","Bcl11b","Tfap2a","Unc5d","Ptprk","Kcnip1","Ripor2",
               "Pcdh17","Tcf4")
di2.markers[,maxpct:=pmax(pct.1, pct.2)]
di2.markers[,Cluster:=ifelse(avg_log2FC>0, 5, 14)]
di2.order <- setorder(di2.markers[rn %in% di2.genes], Cluster, -maxpct)
di2.genes <- rev(di2.order$rn)
#Final list in order
di2.genes <- c("Bmpr1b", "Cbln2", "Ripor2", "St18", "Ebf2", "Bcl11b",
               "Tcf4", "Pcdh17", "Tenm1", "Tfap2a", "Kcnip1", "Ptprk",
               "Unc5d", "Mast4", "Lrrc4c", "Cachd1") 

p.di2 <- DotPlot(eall, idents = c(5,14), features = di2.genes, scale = F, dot.scale = 8)+
  coord_flip() +
  theme(#aspect.ratio = 5,
    panel.grid.major = element_line(linetype = "dashed", color = "grey70"),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(face = "bold"))+
  ylab("Cluster")+
  xlab("Gene")+
  labs(subtitle = "dI2 subpopulations")+
  guides(color=guide_colorbar("Mean Expr", order = 2),
         size = guide_legend("Frequency", order = 1))+
  scale_color_gradient(low = "grey", high = "blue", limits=c(0,3.3))+
  scale_size_continuous(limits = c(0,100),breaks = c(25,50,75))

#DI4/DILa
eall.di4 <- subset(eall, subset = domain == "dI4/dILa")
di4.markers <- as.data.table(FindAllMarkers(eall.di4, 
                                            test.use = "roc", verbose = T),
                             keep.rownames = T)
di4.markers[,diff.pct:=abs(pct.1-pct.2)]
setorder(di4.markers, cluster, -diff.pct)
di4.genes <- c("Kank1","Pdzd2","Cpne5", "Fstl1","Sst",
               "Adarb2","Asic4","St18","Neurod1","Zeb2",
               "Pdzrn3","Thsd4","Rgs16","Dcbld1","Kirrel",
               "Ptchd4","Syt16","Mn1","Nrp2","Gfra1",
               "Grm7","Ntn1","Nalcn","Uncx","Ralyl",
               "Dscam","Zfpm2","Chsy3","March1","Ston2")
di4.markers[,maxpct:=pmax(pct.1, pct.2)]
di4.order <- setorder(di4.markers[rn %in% di4.genes], cluster, -maxpct)
di4.genes <- rev(di4.order$rn)
#Final list in order
di4.genes <- c("Ston2", "Dscam", "March1", "Zfpm2", "Chsy3",
               "Nalcn", "Uncx", "Ntn1", "Ralyl", "Grm7",
               "Mn1", "Ptchd4", "Nrp2", "Syt16", "Gfra1",
               "Rgs16", "Thsd4", "Dcbld1", "Kirrel", "Pdzrn3",
               "St18", "Asic4", "Neurod1", "Adarb2", "Cpne5",
               "Pdzd2", "Sst", "Fstl1", "Kank1")

p.di4 <- DotPlot(eall.di4, idents = c(0,1,3,6,7,9), features = di4.genes, 
                 scale = F, dot.scale = 8)+
  coord_flip() +
  theme(#aspect.ratio = 2,
    panel.grid.major = element_line(linetype = "dashed", color = "grey70"),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(face = "bold"))+
  ylab("Cluster")+
  xlab("Gene")+
  labs(subtitle = "dI4/dILa subpopulations")+
  guides(color=guide_colorbar("Mean Expr", order = 2),
         size = guide_legend("Frequency", order = 1))+
  scale_y_discrete(labels = c("0","1","3","6","7","9a"))+
  scale_color_gradient(low = "grey", high = "blue", limits=c(0,3.3))+
  scale_size_continuous(limits = c(0,100),breaks = c(25,50,75))

#DI5/DILb
eall.di5 <- subset(eall, subset = domain == "dI5/dILb")
di5.markers <- as.data.table(FindAllMarkers(eall.di5, 
                                            test.use = "roc", verbose = T),
                             keep.rownames = T)
di5.markers[,diff.pct:=abs(pct.1-pct.2)]
setorder(di5.markers, cluster, -diff.pct)
di5.genes <- c("Myo16","Gfra1","Abat","Nrp1","Plcb1",
               "Pcdh7","Id4","Rab3ip","Dmd","Pcdh17",
               "Rbp1","Rac3","Tppp3","Syt4")
di5.markers[,max.pct:=pmax(pct.1,pct.2)]
di5.order <- setorder(di5.markers[rn %in% di5.genes], cluster, -pct.1)
di5.genes <- rev(di5.order$rn)
#Final list in order:
di5.genes <- c("Tppp3", "Syt4", "Rac3", "Rab3ip", "Id4",
               "Pcdh17", "Dmd", "Rbp1", "Pcdh7", "Abat",
               "Nrp1", "Myo16", "Gfra1", "Plcb1")

p.di5 <- DotPlot(eall.di5, idents = c(4,13,16), features = di5.genes, 
                 scale = F, dot.scale = 8)+
  coord_flip() +
  theme(#aspect.ratio = 2,
    panel.grid.major = element_line(linetype = "dashed", color = "grey70"),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(face = "bold"))+
  ylab("Cluster")+
  xlab("Gene")+
  labs(subtitle = "dI5/dILb subpopulations")+
  guides(color=guide_colorbar("Mean Expr", order = 2),
         size = guide_legend("Frequency", order = 1))+
  scale_color_gradient(low = "grey", high = "blue", limits=c(0,3.3))+
  scale_size_continuous(limits = c(0,100),breaks = c(25,50,75))

#DI6
eall.di6 <- subset(eall, subset = domain=="dI6")
di6.markers <- as.data.table(FindAllMarkers(eall.di6, 
                                            test.use = "roc", verbose = T),
                             keep.rownames = T)
di6.markers[,diff.pct:=abs(pct.1-pct.2)]
setorder(di6.markers, cluster, -diff.pct)
di6.genes <- c("Brinp1","Irx1","Irx3","Enc1","Irx5",
               "Wt1", "Stk32a", "Sema3a", "Ppp2r2b", "Reln")
di6.markers[,max.pct:=pmax(pct.1,pct.2)]
di6.order <- setorder(di6.markers[rn %in% di6.genes], cluster, -max.pct)
di6.genes <- rev(di6.order$rn)
#Final list in order:
di6.genes <- c("Irx3", "Irx5", "Enc1", "Irx1", "Brinp1",
               "Stk32a", "Reln", "Sema3a", "Ppp2r2b", "Wt1")

p.di6 <- DotPlot(eall.di6, idents = c(8,9), features = di6.genes, 
                 scale = F, dot.scale = 8)+
  coord_flip() +
  theme(#aspect.ratio = 2,
    panel.grid.major = element_line(linetype = "dashed", color = "grey70"),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(face = "bold"))+
  ylab("Cluster")+
  xlab("Gene")+
  labs(subtitle = "dI6 subpopulations")+
  guides(color=guide_colorbar("Mean Expr", order = 2),
         size = guide_legend("Frequency", order = 1))+
  scale_color_gradient(low = "grey", high = "blue", limits=c(0,3.3))+
  scale_size_continuous(limits = c(0,100),breaks = c(25,50,75))+
  scale_y_discrete(labels = c("8", "9b"))

#PATCHWORK OBJECT
layout <- "
11333333444
22333333556
"
p.all <- p.di1 + p.di2 + p.di4 + p.di5 + p.di6 + guide_area() +
  plot_layout(design = layout, guides = "collect")
ggsave(filename = "DomainSubpops.tiff", plot = p.all, path = patchpath, 
       width = 15, height = 8, units = "in")

#Cluster Markers
Idents(eall) <- "seurat_clusters"
p.Satb2 <- VlnPlot(eall, features = c("Satb2"), 
                   pt.size = 0, cols = clustercolors)+
  theme_pubr()+
  theme(legend.position = "none",
        title = element_blank(),
        axis.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, color="black", size=1),
        panel.grid.major.y = element_line(linetype = "longdash", color="grey35"))+
  xlab("Cluster")+ylab("Satb2")
#Satb2 directly negatively regulates DCC in interhemispheric cortical neurons

p.Islr2 <- VlnPlot(eall, features = c("Islr2"),
                   pt.size = 0, cols = clustercolors)+
  theme_pubr()+
  theme(legend.position = "none",
        title = element_blank(),
        axis.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, color="black", size=1),
        panel.grid.major.y = element_line(linetype = "longdash", color="grey35"))+
  xlab("Cluster")+ylab("Islr2")

p.Tfap2b <- VlnPlot(eall, features = c("Tfap2b"),
                    pt.size = 0, cols = clustercolors)+
  theme_pubr()+
  theme(legend.position = "none",
        title = element_blank(),
        axis.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, color="black", size=1),
        panel.grid.major.y = element_line(linetype = "longdash", color="grey35"))+
  xlab("Cluster")+ylab("Tfap2b")

layout <- "
####2222
####2222
11112222
11113333
####3333
"
p.markers <- p.Satb2+p.Islr2+p.Tfap2b+plot_layout(design = layout)
ggsave(filename = "markerVlns.tiff", plot = p.markers, path = patchpath,
       width = 12, height = 5, units = "in")

#Local Reclustering Ventral Domains
#V0
v0.eall<-subset(eall, subset = domain=="v0")
DefaultAssay(v0.eall)<-"integrated"
v0.eall<-ScaleData(v0.eall)
v0.eall<-RunPCA(v0.eall, npcs = 20)
#ElbowPlot(v0.eall)
v0.eall<-FindNeighbors(v0.eall, dims = 1:20, reduction = "pca")
v0.eall<-FindClusters(v0.eall, resolution = 1)
v0.eall<-RunUMAP(v0.eall, dims = 1:20, reduction = "pca")
DefaultAssay(v0.eall)<-"RNA"

v0.de<-as.data.table(FindAllMarkers(v0.eall, test.use = "roc", verbose = T),
                     keep.rownames=T)
v0.genes <- c("Skor1","Lingo2","Onecut3","Asic2","Klhl29",
              "Cblb", "Crybg3", "Igfbpl1", "Klhl35", "Ebf2",
              "Pcdh15","Casz1","Fign","Egfem1","Nin",
              "Lrrc4c","Pou3f1","Ubc",
              "Nr5a2","Ext1","Csmd1","Ldlrad4","Tenm4",
              "Foxb1","Pax6","Pnoc","Rcn1","Brinp2")
v0.de[,max.pct:=pmax(pct.1,pct.2)]
v0.order <- setorder(v0.de[rn %in% v0.genes], cluster, -max.pct)
v0.genes <- rev(v0.order$rn)
#Final list in order:
#add later
p.v0 <- DotPlot(v0.eall, features = v0.genes, 
                scale = F, dot.scale = 8)+
  coord_flip() +
  theme(#aspect.ratio = 2,
    panel.grid.major = element_line(linetype = "dashed", color = "grey70"),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(face = "bold"))+
  ylab("Cluster")+
  xlab("Gene")+
  guides(color=guide_colorbar("Mean Expr", order = 2),
         size = guide_legend("Frequency", order = 1))+
  scale_color_gradient(low = "grey", high = "blue", limits=c(0,3.3))+
  scale_size_continuous(limits = c(0,100),breaks = c(25,50,75))

p.v0dimplot <- DimPlot(v0.eall) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color="black", fill = NA))+
  labs(subtitle = "V0 Local Reclustering")

p.Evx1 <- FeaturePlot(v0.eall, features = "Evx1", order = T)+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
p.Slc17a6 <- FeaturePlot(v0.eall, features = "Slc17a6", order = T)+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
p.Pax2 <- FeaturePlot(v0.eall, features = "Pax2", order = T)+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
p.Pitx2 <- FeaturePlot(v0.eall, features = "Pitx2", order = T)+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
p.Chat <- FeaturePlot(v0.eall, features = "Chat", order = T)+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())


p.v0markers <- p.v0dimplot + p.Evx1 + p.Slc17a6 + p.Pax2 + p.Chat + p.Pitx2 +
  plot_layout(ncol = 2)
p.v0all <- p.v0markers | p.v0
ggsave(filename = "v0local.tiff", plot = p.v0all, path = patchpath, 
       width = 12, height = 8, units = "in")

#V1
v1.eall<-subset(eall, subset = domain=="v1")
DefaultAssay(v1.eall)<-"integrated"
v1.eall<-ScaleData(v1.eall)
v1.eall<-RunPCA(v1.eall, npcs = 20)
ElbowPlot(v1.eall)
v1.eall<-FindNeighbors(v1.eall, dims = 1:20, reduction = "pca")
v1.eall<-FindClusters(v1.eall, resolution = 1)
v1.eall<-RunUMAP(v1.eall, dims = 1:20, reduction = "pca")
DefaultAssay(v1.eall)<-"RNA"

v1.de<-as.data.table(FindAllMarkers(v1.eall, test.use = "roc", verbose = T),
                     keep.rownames=T)

#V3
v3.eall<-subset(eall, subset = domain=="v3")
DefaultAssay(v3.eall)<-"integrated"
v3.eall<-ScaleData(v3.eall)
v3.eall<-RunPCA(v3.eall, npcs = 20)
ElbowPlot(v3.eall)
v3.eall<-FindNeighbors(v3.eall, dims = 1:20, reduction = "pca")
v3.eall<-FindClusters(v3.eall, resolution = 1)
v3.eall<-RunUMAP(v3.eall, dims = 1:20, reduction = "pca")
DefaultAssay(v3.eall)<-"RNA"

v3.de<-as.data.table(FindAllMarkers(v3.eall, test.use = "roc", verbose = T),
                     keep.rownames=T) #REDO THIS

#DI5

#Domain-specific Pre/Post Markers

#Robo1 vs Robo2 vs L1cam
dt.cor <- as.data.table(FetchData(eall, vars = c("orig.ident", "cross", "Robo3",
                                                 "Robo2","Robo1","L1cam","Cntn2","Gapdh"),
                                  slot="counts"),
                        keep.rownames = T)
genes <- c("Robo2","L1cam","Robo1","Cntn2","Gapdh")
cors <- c(-0.323, -0.081, 0.193, 0.565, 0.087)
corDT <- data.table(genes=genes, cor=cors)
corDT[,genes:=factor(genes, levels = c("Robo2","L1cam","Robo1","Cntn2","Gapdh"))]

p.cor <- ggplot(corDT, aes(x=genes, y=cor))+
  geom_bar(stat = "identity", position = "dodge", 
           fill="grey70", color="black", linewidth=1,
           width=0.6)+
  ylim(-1,1)+
  theme_pubr()+
  xlab("Gene") + ylab("Spearman Coefficient versus Robo3")+
  theme(axis.title = element_text(face = "bold"),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid.major.y = element_line(linetype = "longdash", color = "grey35"))+
  geom_hline(yintercept = 0, linewidth=1.5)+
  geom_vline(xintercept = 3.5, linetype="dashed", color="grey35", linewidth=0.7)

ggsave(filename = "Spearman.tiff", plot = p.cor, path = patchpath,
       width = 12, height = 4, units = "in")

#Monocle -- Pseudotime
eall.di1 <- subset(eall, subset = domain == "dI1")
eall.di2 <- subset(eall, subset = domain == "dI2")
eall.di4 <- subset(eall, subset = domain == "dI4/dILa")
eall.di5 <- subset(eall, subset = domain == "dI5/dILb")
eall.di6 <- subset(eall, subset = domain == "dI6")
eall.v0 <- subset(eall, subset = domain == "v0")
eall.v1 <- subset(eall, subset = domain == "v1")
eall.v3 <- subset(eall, subset = domain == "v3")

#Convert to Cell Data Set
cds <- as.cell_data_set(eall.v3)

#Factorize orig.ident, sex, and embryo
colData(cds)$orig.ident <- factor(
  colData(cds)$orig.ident, 
  labels = c("E10.5", "E11.5", "E12.5", "E13.5"))
colData(cds)$sex <- as.factor(colData(cds)$sex)
colData(cds)$embryo <- as.factor(colData(cds)$embryo)

#Pre-process and plot PCs
cds_preprocessed <- preprocess_cds(cds = cds, num_dim = 50, verbose = T)
plot_pc_variance_explained(cds_preprocessed)

#Remove batch effects by timepoint
cds_aligned_age <- align_cds(cds_preprocessed, alignment_group = "orig.ident", verbose = T)
#cds_aligned_embryo_only <- align_cds(cds_preprocessed, alignment_group = "libID", verbose = T)
#cds_aligned_age_embryo <- align_cds(cds_aligned_embryo_only, alignment_group = "orig.ident", verbose = T)

#umapnneighbors 20 or 17 most promising
cds_reduced <- reduce_dimension(cds_aligned_age,verbose = T, reduction_method = "Aligned") #-
plot_cells(cds_reduced)

#cds_reduced_e <- reduce_dimension(cds_aligned_embryo_only, umap.n_neighbors = 10,verbose = T)

rowData(cds_reduced)$gene_name <- rownames(cds_reduced)
rowData(cds_reduced)$gene_short_name <- rowData(cds_reduced)$gene_name

cds_clustered<-cluster_cells(cds_reduced) #-
plot_cells(cds_clustered, color_cells_by = "partition")

cds_learned <- learn_graph(cds_clustered)
plot_cells(cds_learned)

cds_ordered <- order_cells(cds_learned)
plot_cells(cds_ordered, color_cells_by = "pseudotime")

cds.di1 <- cds_ordered
cds.di2 <- cds_ordered
# cds.dI4 <- cds_ordered
cds.di5 <- cds_ordered
cds.di6 <- cds_ordered
cds.v0 <- cds_ordered
cds.v1 <- cds_ordered
cds.v3 <- cds_ordered
traj.coord.di1 <- cds_ordered@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
traj.coord.di2 <- cds_ordered@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
traj.coord.di4 <- cds_ordered@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
traj.coord.di5 <- cds_ordered@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
traj.coord.di6 <- cds_ordered@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
traj.coord.v0 <- cds_ordered@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
traj.coord.v1 <- cds_ordered@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
traj.coord.v3 <- cds_ordered@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]

traj.all <- c(traj.coord.di1, traj.coord.di2, traj.coord.di4,
              traj.coord.di5, traj.coord.di6, traj.coord.v0,
              traj.coord.v1, traj.coord.v3)
saveRDS(traj.all, file = "trajCoordAll.RData")
traj.all <- readRDS(file = "trajCoordAll.RData")

pDT <- as.data.table(FetchData(eall, vars = c("domain","cross","orig.ident",
                                              "Robo3","Robo2","Robo1","Cntn2",
                                              "L1cam","Gapdh","Nhlh1","UMAP_1","UMAP_2"), 
                               slot = "counts"),
                     keep.rownames = T)
pDT[stack(traj.all), on=.(rn=ind), traj.all:=i.values][]
pDT[,maxTraj:=max(traj.all),by=domain]
pDT[,propTraj:=traj.all/maxTraj]
eall <- AddMetaData(eall, metadata = data.frame(propTraj=pDT$propTraj, 
                                                row.names = rownames(eall@meta.data)))
pDT[,trajBin:=cut(propTraj, breaks=20, labels=F)]
pDT[,trajBin:=as.factor(trajBin)]
eall <- AddMetaData(eall, metadata = data.frame(trajBin=pDT$trajBin,
                                                row.names = rownames(eall@meta.data)))
pDTnorm <- as.data.table(FetchData(eall, vars = c("domain","cross","orig.ident",
                                                  "Robo3", "Robo2", "Robo1", "Cntn2",
                                                  "L1cam", "Gapdh", "Nhlh1", "PC_2",
                                                  "Neurod6","Zfhx4","Nfib", "Zfhx3",
                                                  "Neurod2")),
                         keep.rownames = T)
pDTnorm[stack(traj.all), on=.(rn=ind), traj.all:=i.values][]
pDTnorm[,maxTraj:=max(traj.all),by=domain]
pDTnorm[,propTraj:=traj.all/maxTraj]
pDTnorm[,trajBin:=cut(propTraj, breaks=20, labels=F)]
pDTnorm[,trajBin:=as.factor(trajBin)]
pDTnorm[,trajBin:=as.numeric(trajBin)]

makeSummaryDT <- function(gene, inputDT=pDTnorm, group="trajBin"){
  res <- as.data.table(Rmisc::summarySE(data = inputDT,
                                        measurevar = gene,
                                        groupvars = group))
  return(res)
}

summ.Robo3 <- makeSummaryDT("Robo3")
summ.Robo2 <- makeSummaryDT("Robo2")
summ.Robo1 <- makeSummaryDT("Robo1")
summ.L1cam <- makeSummaryDT("L1cam")
summ.Cntn2 <- makeSummaryDT("Cntn2")
summ.PC2 <- makeSummaryDT("PC_2")
summ.Nhlh1 <- makeSummaryDT("Nhlh1")
summ.Zfhx3 <- makeSummaryDT("Zfhx3")
summ.Zfhx4 <- makeSummaryDT("Zfhx4")
summ.Neurod2 <- makeSummaryDT("Neurod2")
summ.Neurod6 <- makeSummaryDT("Neurod6")
summ.Nfib <- makeSummaryDT("Nfib")

plotSumm <- function(summDT, gene, ylabel=gene) {
  p <- ggplot(summDT, aes(x=trajBin, y= .data[[gene]],
                          ymin = .data[[gene]]-.data$sd, ymax= .data[[gene]]+.data$sd))+
    geom_ribbon(alpha=0.5, group=1, fill="lightblue")+
    geom_line(group=1, linewidth=2)+
    theme_pubr()+
    scale_x_continuous(breaks=1:20)+
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill=NA, color="black"),
          panel.grid.major = element_line(linetype = "dashed", color = "grey35"),
          axis.title = element_text(face="bold"))+
    xlab("Pseudotime Bin")+
    ylab(ylabel)
  return(p)
}

p.Robo2 <- plotSumm(summ.Robo2, "Robo2")
p.Robo1 <- plotSumm(summ.Robo1, "Robo1")
p.Robo3 <- plotSumm(summ.Robo3, "Robo3")
p.L1cam <- plotSumm(summ.L1cam, "L1cam")
p.Cntn2 <- plotSumm(summ.Cntn2, "Cntn2")
p.Nhlh1 <- plotSumm(summ.Nhlh1, "Nhlh1")
p.PC2 <- plotSumm(summ.PC2, "PC_2", "PC 2")
p.Zfhx3 <- plotSumm(summ.Zfhx3, "Zfhx3")
p.Zfhx4 <- plotSumm(summ.Zfhx4, "Zfhx4")
p.Neurod2 <- plotSumm(summ.Neurod2, "Neurod2")
p.Neurod6 <- plotSumm(summ.Neurod6, "Neurod6")
p.Nfib <- plotSumm(summ.Nfib, "Nfib")

p.pTimeZN <- p.Zfhx3 + p.Neurod2 + p.Zfhx4 + p.Neurod6 + plot_spacer() + p.Nfib +
  plot_layout(ncol = 2)

p.pTimeCross <- p.Robo3 + p.Robo2 + p.Cntn2 + p.Robo1 + p.Nhlh1 + p.L1cam +
  plot_spacer() + p.PC2 + plot_layout(ncol=2)

ggsave("PTimeCross.tiff", plot = p.pTimeCross, path = patchpath,
       width = 12, height = 12, units = "in")
ggsave("PTimeZN.tiff", plot = p.pTimeZN, path = patchpath, 
       width = 12, height = 9, units = "in")

#Percent Pre/Post by Pseudotime Bin
pTime.dt <- as.data.table(FetchData(eall, vars = c("cross", "trajBin", "domain")),
                          keep.rownames = T)
pTime.dt <- pTime.dt[cross != "None"]
pTime.summ <- pTime.dt[,.N, by=.(cross, trajBin)]
pTime.summ[,total:=sum(N),by=trajBin]
setorder(pTime.summ, trajBin)
pTime.summ[,trajBin:=as.numeric(trajBin)]
pTime.summ[,prop:=N/total]

p.crossProp <- ggplot(pTime.summ, aes(x=trajBin, y=prop, fill=cross))+
  geom_area(color="black", alpha=0.5, linewidth=1)+
  theme_pubr()+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA, color = "black"),
        axis.title = element_text(face="bold"),
        legend.position="none",
        panel.grid.major = element_line(linetype="dashed", color = "grey35"),
        panel.grid.minor.y = element_line(linetype = "dashed", color = "grey70"))+
  scale_fill_manual(values = crossColors)+
  scale_x_continuous(breaks = 1:20)+
  xlab("Pseudotime Bin") + ylab("Proportional Composition")+
  annotate(geom = "text", label="Pre-Crossing", x=1.5, y=0.95, hjust=0)+
  annotate(geom = "text", label="Post-Crossing", x=19.5, y=0.05, hjust=1)

p.crossProgressCombo <- p.dotcross + (p.crossProp/p.zn) + 
  plot_layout(ncol = 2, widths = c(1,2))

ggsave(filename = "CrossProgress.tiff", plot = p.crossProgressCombo, path = patchpath,
       width = 12, height = 5, units = "in")

#DimPlot by Pseudotime
p.dimTime <- FeaturePlot(eall, features = "propTraj", order = T)+
  scale_color_viridis_c(name="Normalized\nPseudotime")+
  theme_pubr()+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, color="black"),
        legend.position = "right")+
  labs(title = NULL)

ggsave("DimPTime.tiff", plot = p.dimTime, path = patchpath,
       width = 5, height=5, units = "in")
# p.pseudo <- p.dimTime + p.crossProgressCombo + plot_layout(ncol = 2, widths = c(1,4))

#Useful for pre/post domain markers:
#VlnPlot(eall, features = "Tcf4", group.by="domain", split.by="cross", pt.size=0,
#idents=c("preCrossing","postCrossing"))+scale_fill_manual(values=crossColors)

#Domain Specific Pre/Post Markers
eall.di1 <- subset(eall, subset = domain == "dI1")
eall.di2 <- subset(eall, subset = domain == "dI2")
eall.di4 <- subset(eall, subset = domain == "dI4/dILa")
eall.di5 <- subset(eall, subset = domain == "dI5/dILb")
eall.di6 <- subset(eall, subset = domain == "dI6")
eall.v0 <- subset(eall, subset = domain == "v0")
eall.v1 <- subset(eall, subset = domain == "v1")
eall.v3 <- subset(eall, subset = domain == "v3")

FindCrossMarkers <- function(inputDT, genes=NULL, logthresh=0.25, minpct=0.1){
  Idents(inputDT) <- "cross"
  res <- as.data.table(FindMarkers(inputDT,
                                   ident.1 = "preCrossing", ident.2 = "postCrossing",
                                   test.use = "roc", features=genes, logfc.threshold=logthresh,
                                   min.pct=minpct),
                       keep.rownames = T)
  return(res)
}

pp.di1 <- FindCrossMarkers(eall.di1)
pp.di2 <- FindCrossMarkers(eall.di2)
pp.di4 <- FindCrossMarkers(eall.di4)
pp.di5 <- FindCrossMarkers(eall.di5)
pp.di6 <- FindCrossMarkers(eall.di6)
pp.v0 <- FindCrossMarkers(eall.v0)
pp.v1 <- FindCrossMarkers(eall.v1)
pp.v3 <- FindCrossMarkers(eall.v3)

countFreq <- function(geneNames){
  counts <- c()
  for (geneName in geneNames){
    res<-0
    if(geneName %in% pp.di1[[1]]){res <- res+1}
    if(geneName %in% pp.di2[[1]]){res <- res+1}
    if(geneName %in% pp.di4[[1]]){res <- res+1}
    if(geneName %in% pp.di5[[1]]){res <- res+1}
    if(geneName %in% pp.di6[[1]]){res <- res+1}
    if(geneName %in% pp.v0[[1]]){res <- res+1}
    if(geneName %in% pp.v1[[1]]){res <- res+1}
    if(geneName %in% pp.v3[[1]]){res <- res+1}
    counts <- c(counts, res)
  }
  return(counts)
}

VlnWrapper <- function(gene){
  p <- VlnPlot(eall, features = gene, group.by = "domain", split.by = "cross",
               idents = c("preCrossing","postCrossing"), pt.size = 0)+
    scale_fill_manual(values = crossColors)+
    theme_pubr()+
    theme(legend.position = "none",
          panel.background = element_blank(),
          panel.border = element_rect(fill=NA, color = "black"),
          axis.title = element_text(face="bold"),
          panel.grid.major.y = element_line(linetype="dashed", color = "grey35"))+
    labs(title = NULL) + xlab("Domain") + ylab(gene)
  return(p)
}

pp.di1[,freq:=countFreq(rn)]
setorder(pp.di1, freq)

pp.di2[,freq:=countFreq(rn)]
setorder(pp.di2, freq)

pp.di4[,freq:=countFreq(rn)]
setorder(pp.di4, freq)

pp.di5[,freq:=countFreq(rn)]
setorder(pp.di5, freq)

pp.di6[,freq:=countFreq(rn)]
setorder(pp.di6, freq)

pp.v0[,freq:=countFreq(rn)]
setorder(pp.v0, freq)

pp.v1[,freq:=countFreq(rn)]
setorder(pp.v1, freq)

pp.v3[,freq:=countFreq(rn)]
setorder(pp.v3, freq)

subtypeMarkers <- c("Dcc", "Hmx3", "Clvs2", "Sertad4", "Dcaf6",
                    "Iglon5", "Antxr2", "Ntrk1", "Kcnh4", "Prkd1")

p.dcc <- VlnWrapper("Dcc")
p.hmx3 <- VlnWrapper("Hmx3")
p.clvs2 <- VlnWrapper("Clvs2")
p.sertad4 <- VlnWrapper("Sertad4")
p.dcaf6 <- VlnWrapper("Dcaf6")
#p.iglon5 <- VlnWrapper("Iglon5")
p.antxr2 <- VlnWrapper("Antxr2")
p.ntrk1 <- VlnWrapper("Ntrk1")
p.kcnh4 <- VlnWrapper("Kcnh4")
p.prkd1 <- VlnWrapper("Prkd1")

p.all <- p.dcc + p.hmx3 + p.clvs2 +
  p.sertad4 + p.dcaf6 + p.antxr2 +
  p.ntrk1 + p.kcnh4 + p.prkd1 +
  plot_layout(ncol=3)

ggsave("DomainCross.tiff", plot = p.all, path = patchpath,
       width = 20, height = 15, units = "in")

#Domain-specific Pre/Post Markers: Lollipops
genes.use <- c("Dcc","Hmx3","Clvs2","Sertad4","Dcaf6","Iglon5","Antxr2","Ntrk1","Kcnh4","Prkd1")

m.di1<-FindCrossMarkers(eall.di1, genes = genes.use, logthresh = 0, minpct = 0)
m.di2<-FindCrossMarkers(eall.di2, genes = genes.use, logthresh = 0, minpct = 0)
m.di4<-FindCrossMarkers(eall.di4, genes = genes.use, logthresh = 0, minpct = 0)
m.di5<-FindCrossMarkers(eall.di5, genes = genes.use, logthresh = 0, minpct = 0)
m.di6<-FindCrossMarkers(eall.di6, genes = genes.use, logthresh = 0, minpct = 0)
m.v0<-FindCrossMarkers(eall.v0, genes = genes.use, logthresh = 0, minpct = 0)
m.v1<-FindCrossMarkers(eall.v1, genes = genes.use, logthresh = 0, minpct = 0)
m.v3<-FindCrossMarkers(eall.v3, genes = genes.use, logthresh = 0, minpct = 0)

mList <- list("dI1"=m.di1,
              "dI2"=m.di2,
              "dI4/dILa"=m.di4,
              "dI5/dILb"=m.di5,
              "dI6"=m.di6,
              "v0"=m.v0,
              "v1"=m.v1,
              "v3"=m.v3
)
m.all <- data.table::rbindlist(mList, idcol = "domain")
m.all[,domain:=factor(domain, levels=c("v3","v1","v0","dI6","dI5/dILb",
                                       "dI4/dILa","dI2","dI1"))]
m.all[,pointcolor:=rep("None",.N)]
m.all[avg_log2FC>=0.25, pointcolor:="preCrossing"]
m.all[avg_log2FC<=-0.25, pointcolor:="postCrossing"]

makeLollipop <- function(geneName){
  p<-ggplot(m.all[rn==geneName], aes(x=domain, y=avg_log2FC))+
    geom_segment(aes(x=domain, xend=domain,
                     y=0, yend=avg_log2FC),
                 color="black")+
    geom_point(aes(color=pointcolor), size=3)+
    scale_color_manual(values = crossColors)+
    coord_flip()+
    geom_hline(yintercept = 0, color="black", linewidth=0.8)+
    ylim(-1, 2.25)+
    theme_pubr()+
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill=NA, color="black"),
          panel.grid.major = element_line(color = "grey35", linetype = "dotted"),
          plot.title = element_text(face = "bold"),
          legend.position = "none")+
    xlab(label = NULL) + ylab(label = "Average Log2 Fold Change")+
    labs(title=geneName)
  return(p)
}

l.prkd1 <- makeLollipop("Prkd1")
l.ntrk1 <- makeLollipop("Ntrk1")
l.antxr2 <- makeLollipop("Antxr2")
l.sertad4 <- makeLollipop("Sertad4")
l.dcc <- makeLollipop("Dcc")
l.kcnh4 <- makeLollipop("Kcnh4")
l.dcaf6 <- makeLollipop("Dcaf6")
l.hmx3 <- makeLollipop("Hmx3")

l.all <- l.prkd1 + l.ntrk1 + l.antxr2 + l.sertad4 + l.dcc + l.kcnh4 + l.dcaf6 + l.hmx3 +
  plot_layout(ncol = 4)

ggsave("DomainLollipops.tiff", plot = l.all, path = patchpath,
       width = 15, height = 8, units = "in")

#Pan Commissural Markers
inAll <- pp.di1[freq==8]
setorder(inAll, -avg_log2FC)

cross.markers <- readRDS(file = "prePostROCMarkersDT.RData")

top20post <- head(cross.markers[avg_log2FC<0],20)$rn
top20pre <- head(cross.markers[avg_log2FC>0],20)$rn

dt <- data.table(gene=c(top20pre,top20post))
dt[,freq:=countFreq(gene)]


dt.cross <- as.data.table(FetchData(eall, 
                                    vars = c("orig.ident","seurat_clusters",
                                             "cross","domain")),
                          keep.rownames = T)
dt.cross <- dt.cross[cross != "None"]
dt.cross[,cross:=factor(cross, levels = c("preCrossing","postCrossing"))]
setorder(dt.cross, cross, orig.ident)
prePostCells <- dt.cross$rn
# prePostCells <- WhichCells(eall, expression = crossNew %in% c("preCrossing","newPre",
#                                                               "postCrossing","newPost","R1high")) #figure out why cross doesn't work
eall@active.ident <- factor(eall@active.ident, levels = c("preCrossing","postCrossing","None"))
p.crossHeat <- DoHeatmap(eall, features = dt[freq==8,gene], cells = prePostCells,
                         label = F, group.colors = crossColors) +
  guides(color="none")

ggsave("CrossHeat.tiff", plot = p.crossHeat, path = patchpath, 
       width = 15, height = 9, units = "in")

#Proportional Pseudotime by Domain
FeaturePlot(eall, cells = WhichCells(eall, expression = domain=="dI1"), features = "propTraj")+
  scale_color_viridis_c()+xlim(-12.5,-8.5)+ylim(-8,-1)

makeTrajPlot <- function(domainName, xmin, xmax, ymin, ymax, legendposition="bottom"){
  p<-FeaturePlot(eall, cells = WhichCells(eall, expression = domain==domainName),
                 features = "propTraj")+
    scale_color_viridis_c(name="Normalized Pseudotime")+xlim(xmin, xmax)+ylim(ymin, ymax)+
    theme_pubr()+
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.position = legendposition)+
    labs(title = domainName)
  return(p)
}

traj.di1 <- makeTrajPlot("dI1", -12.5, -8.5, -8, -1)
traj.di2 <- makeTrajPlot("dI2", -12.5, -8.5, 0, 5.5)
traj.di4 <- makeTrajPlot("dI4/dILa", -1, 11, -10.5, 6)
traj.di5 <- makeTrajPlot("dI5/dILb", -4, 2.5, 6, 11.5)
traj.di6 <- makeTrajPlot("dI6", -6, 2, -7, -2)
traj.v0 <- makeTrajPlot("v0", -6.5, -3.5, -2.5, 4)
traj.v1 <- makeTrajPlot("v1", -6.5, -1, -4, -1.5)
traj.v3 <- makeTrajPlot("v3", -7.5, -4.5, 5.5, 9)
traj.v3.2 <- makeTrajPlot("v3", -7.5, -4.5, 5.5, 9, legendposition = "bottom")
domainPlot <- DimPlot(eall, group.by = "domain", cols = domainColors)+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank())+
  labs(title = NULL)

traj.all <- ((traj.di1 / traj.di2) | (traj.di4) | 
               (traj.di5 / traj.di6) | (traj.v0 / traj.v1 / traj.v3)) + 
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("trajByDomain.tiff", plot = traj.all, path = patchpath,
       width = 12, height = 5, units = "in")

p.nhlh1 <- FPWrapper("Nhlh1")
p.robo3 <- FPWrapper("Robo3")
p.both <- p.nhlh1 / p.robo3
ggsave("nhlh1robo3.tiff", plot = p.both, path = patchpath,
       width = 2.5, height = 5, units = "in")

ggsave("domainDimPlot.tiff", plot = domainPlot, path = patchpath, 
       width = 5, height = 5, units = "in")
