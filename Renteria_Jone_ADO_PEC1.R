
## Jone Renteria (jone.renteria25@gmail.com)

## Microarray analysis: GEO study --> GSE85269

## ----CreateFolders, warning=FALSE, eval=FALSE------------------------------------------------------
setwd(".")
dir.create("data")
dir.create("results")
dir.create("figures")


## ----ReadTargets-----------------------------------------------------------------------------------
targets <- read.csv2("./data/targets.csv", header = TRUE,sep=";") 
knitr::kable(
  targets, booktabs = TRUE,
  caption = 'Content of the targets file used for the current analysis')
str(targets)


## ----installPackages, message=FALSE, warning=FALSE, eval=FALSE-------------------------------------
## Instal all the required packages to run the microarrays analysis in Bioconductor. For this concrete analysis, 
## the following packages would be required: 
## BiocManager::install("mogene20sttranscriptcluster.db")
## BiocManager::install("pd.mogene.2.0.st")


## ----ReadCELfiles, message=FALSE, results='hide', warning=FALSE------------------------------------
library(oligo)
celFiles <- list.celfiles("./data", full.names = TRUE)

#### This must be run the first time to change the CSV file to automatically include the processing date of each array. 
#library(affyio) 
#library(dplyr)
#date<-get.celfile.dates(filename=celFiles) #to know the date of proccesing the files 
#targets<- targets %>%
#write.csv("C:/Users/joner/OneDrive/Escritorio/Master_semestre 3/AD_Omicos/pec1/data/targets.csv")
#targets$date<-as.factor(targets$date)
#str(targets)

library(Biobase)
my.targets <-read.AnnotatedDataFrame(file.path("./data","targets.csv"), 
                                     header = TRUE, row.names = 1, sep=";") 
rawData <- read.celfiles(celFiles, phenoData = my.targets)


## ----ChangeName------------------------------------------------------------------------------------
my.targets@data$ShortName->rownames(pData(rawData))
colnames(rawData) <-rownames(pData(rawData))  

## ----QCRaw, message=FALSE, warning=FALSE, eval=FALSE-----------------------------------------------
library(arrayQualityMetrics)
arrayQualityMetrics(rawData,force = TRUE)


## ----QCRawDataRes, fig.cap="Aspect of the summary table, in the index.html file, produced by the arrayQualityMetrics package on the raw data", echo=FALSE----
###knitr::include_graphics("figures/Figure1.png")


## -----r BoxplotRaw, message=FALSE, fig.cap="Boxplot for arrays intensities (Raw Data)"---------------------------------------------------------------------------------------------
library(ggplot2)
pmexp <- pm(rawData)
sampleNames <- vector();logs = vector()
for (i in 1:12){
  sampleNames = c(sampleNames,rep(my.targets@data[i,4],dim(pmexp)[1]))
  logs = c(logs,log2(pmexp[,i]))}


cols<-c(rep("chartreuse3",3), rep("forestgreen",3), rep("chocolate2",3), rep("firebrick3",3))
logData <- data.frame(logInt=logs,sampleName=sampleNames)
dataBox <- ggplot(logData,aes(sampleName,logInt,fill= sampleName)) +
  geom_boxplot()+
  scale_fill_manual(values=cols)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Boxplot for arrays intensities (Raw Data)")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("")+ylab("Intensities (log)")

ggsave(filename = "Figure2_Boxplot_intensities_raw_data.png",dataBox, path = "./figures" )

dataBox 

##-----------r histoRaw-----------------------------------------------------------------------
dataHist2 <- ggplot(logData, aes(logInt, colour = sampleName)) +
  geom_density() +
  scale_color_manual(values=cols)+
  ggtitle("Signal distribution of raw data")+ 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = "Figure3_Histogram_raw_data.png",dataHist2 ,path = "./figures" )

dataHist2 

##----------------------r scatterRaw, echo=FALSE-----------------------------------------
func_scatt <- function (datos, labels, factor) {
  data <- prcomp(t(datos))
  # plot adjustments
  dataDf <- data.frame(data$x)
  Group <- factor
  perct <- round(data$sdev^2/sum(data$sdev^2)*100,1)
  cols2<-c("chartreuse3","forestgreen","darkolivegreen3","chocolate2","firebrick1","firebrick3")
  
  plot_scatt <- ggplot(dataDf,aes(x=PC1, y=PC2)) +
    geom_point(aes(color = Group), alpha = 0.55, size = 3) +
    geom_text(aes(label=labels),size=2,hjust=0, vjust=0)+
    coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5)) +
    ggtitle("Principal Component Analysis for Raw Data")+ 
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = c(paste("PC1",perct[1],"%")),y=c(paste("PC2",perct[2],"%")))+
    theme_classic() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    scale_color_manual(values=cols2)
  
  plot_scatt
}

fig_Scatt<-func_scatt(datos=exprs(rawData),labels = targets$ShortName, factor = targets$Group)
ggsave(filename = "Figure4_PCA_raw_data.png",fig_Scatt, path = "./figures" )
fig_Scatt

## ----Normalization---------------------------------------------------------------------------------
rma_data <- rma(rawData)

## ----QCNorm, message=FALSE, warning=FALSE, eval=FALSE----------------------------------------------
arrayQualityMetrics(rma_data, outdir = file.path("./results", "QCDir.Norm"), force=TRUE)

## ----QCNormDataRes, fig.cap="Aspect of the summary table, in the index.html file, produced by the arrayQualityMetrics package on normalized data", echo=FALSE----
##knitr::include_graphics("figures/Figure5.png")


##------------BoxplotNorm, message=FALSE, fig.cap="Boxplot for arrays intensities (Normalized Data)"------------
rma_data_matrix<-as.matrix(rma_data)
sampleNames = vector();logs = vector()
for (i in 1:12){
  sampleNames = c(sampleNames,rep(my.targets@data[i,4],dim(rma_data_matrix)[1]))
  logs = c(logs,rma_data_matrix[,i])}

cols<-c(rep("chartreuse3",2), rep("forestgreen",2), rep ("darkolivegreen3",2), rep("chocolate2",2), rep("firebrick1",2), rep("firebrick3",2))
logData <- data.frame(logInt=logs,sampleName=sampleNames)
dataBox_norm <- ggplot(logData,aes(sampleName,logInt,fill= sampleName)) +
  geom_boxplot()+
  scale_fill_manual(values=cols)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Boxplot for arrays intensities (Normalized Data)")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("")+ylab("Intensities (log)")

ggsave(filename = "Figure6_Boxplot_intensities_normalized_data.png",dataBox_norm, path = "./figures" )

dataBox_norm

##---------------------------------r histoNorm----------------------------------------------
dataHist_norm <- ggplot(logData, aes(logInt, colour = sampleName)) +
  geom_density() +
  scale_color_manual(values=cols)+
  ggtitle("Signal distribution of normalized data")+ 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = "Figure7_Histogram_normalized_data.png",dataHist_norm ,path = "./figures" )
dataHist_norm
dev.off()

##--------------------------------r scatterNorm, echo=FALSE----------------------------------------------
fig_Scatt_norm<-func_scatt(datos=rma_data_matrix,labels = targets$ShortName, factor = targets$Group)
ggsave(filename = "Figure8_PCA_normalized_data.png",fig_Scatt, path = "./figures" )
fig_Scatt_norm


## ----BatchDetection, message=FALSE, warning=FALSE--------------------------------------------------
library(pvca)
targets<-pData(rma_data) 
pct_threshold <- 0.6 #select the threshold
batch.factors <- c("Genotype", "Breath","date") #select the factors to analyze
pvcaObj <- pvcaBatchAssess(rma_data, batch.factors,pct_threshold) #run the analysis


## ----plotPVCA, fig.cap="Relative importance of the different factors affecting gene expression----
#plot the results
bp <- barplot(pvcaObj$dat, xlab = "Effects",ylab = "Weighted average proportion variance",
              ylim= c(0,1.1),col = c("chocolate1"), las=2,
              main="PVCA estimation")
axis(1, at = bp, labels = pvcaObj$label, cex.axis = 0.75, las=2)
values = pvcaObj$dat
new_values = round(values,3)
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.7)
dev.off()

##-----------------------r saveplot, echo=FALSE--------------------------------------------------------
png(filename = "./figures/Figure9_barplot.png")
bp <- barplot(pvcaObj$dat, xlab = "Effects",ylab = "Weighted average proportion variance",
              ylim= c(0,1.1),col = c("chocolate1"), las=2,
              main="PVCA estimation")
axis(1, at = bp, labels = pvcaObj$label, cex.axis = 0.75, las=2)
values = pvcaObj$dat
new_values = round(values,3)
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.7)
dev.off()

##---------------SDplot, fig.cap="Values of standard deviations allong all samples for all genes ordered from smallest to biggest"-------------
sds <- apply (exprs(rma_data), 1, sd)
sdsO<- sort(sds)
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes",
     sub="Vertical lines represent 90% and 95% percentiles",
     xlab="Gene index (from least to most variable)", ylab="Standard deviation")
abline(v=length(sds)*c(0.9,0.95))
dev.off()

##---------- saveSDplot, echo=FALSE-------------------------------------------------------
png(filename = "./figures/Figure10_SDplot.png")
sds <- apply (exprs(rma_data), 1, sd)
sdsO<- sort(sds)
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes",
     sub="Vertical lines represent 90% and 95% percentiles",
     xlab="Gene index (from least to most variable)", ylab="Standard deviation")
abline(v=length(sds)*c(0.9,0.95))
dev.off()




## ----saveSDplot, echo=FALSE, results='hide'--------------------------------------------------------
tiff("figures/SDplot.tiff", res = 150, width = 5, height = 5, units = 'in')
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes",
     sub="Vertical lines represent 90% and 95% percentiles",
     xlab="Gene index (from least to most variable)", ylab="Standard deviation")
abline(v=length(sds)*c(0.9,0.95))
dev.off()


## ----Filtering1, results='hide', message=FALSE-----------------------------------------------------
library(genefilter)
library(mogene20sttranscriptcluster.db)
annotation(rma_data) <- "mogene20sttranscriptcluster.db"
filtered <- nsFilter(rma_data, 
                     require.entrez = TRUE, remove.dupEntrez = TRUE,
                     var.filter=TRUE, var.func=IQR, var.cutoff=0.75, 
                     filterByQuantile=TRUE, feature.exclude = "^AFFX")


## ----FilterResults2--------------------------------------------------------------------------------
print(filtered$filter.log)
eset_filtered <-filtered$eset

## ----SaveData1, results='hide', message=FALSE------------------------------------------------------
write.csv(exprs(rma_data), file="./results/normalized.Data.csv")
write.csv(exprs(eset_filtered), file="./results/normalized.Filtered.Data.csv")
save(rma_data, eset_filtered, file="./results/normalized.Data.Rda")


## ----LoadSavedData---------------------------------------------------------------------------------
if (!exists("eset_filtered")) load (file="./results/normalized.Data.Rda")


## ----DesignMatrix, message=FALSE-------------------------------------------------------------------
library(limma)
designMat<- model.matrix(~0+Group, pData(eset_filtered))
colnames(designMat) <- c("WTBAS", "KOBAS", "WTVEN", "KOVEN")
print(designMat)


## ----setContrasts----------------------------------------------------------------------------------
cont.matrix <- makeContrasts (KOvsWT.BAS = KOBAS-WTBAS,
                              KOvsWT.VEN = KOVEN-WTVEN,
                              INT = (KOBAS-WTBAS) - (KOVEN-WTVEN),
                              levels=designMat)
print(cont.matrix)


## ---- linearmodelfit-------------------------------------------------------------------------------
fit<-lmFit(eset_filtered, designMat)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)


## ---- topTabs_1-------------------------------------------------------------------------------------
topTab_KOvsWT.BAS <- topTable (fit.main, number=nrow(fit.main), coef="KOvsWT.BAS", adjust="fdr") 
head(topTab_KOvsWT.BAS)


## ---- topTabs_2-------------------------------------------------------------------------------------
topTab_KOvsWT.VEN <- topTable (fit.main, number=nrow(fit.main), coef="KOvsWT.VEN", adjust="fdr") 
topTab_INT <- topTable (fit.main, number=nrow(fit.main), coef="INT", adjust="fdr") 


## ----GeneAnnotation, message=FALSE, warning=FALSE--------------------------------------------------
annotatedTopTable <- function(topTab, anotPackage)
{
  topTab <- cbind(PROBEID=rownames(topTab), topTab)
  myProbes <- rownames(topTab)
  thePackage <- eval(parse(text = anotPackage))
  geneAnots <- select(thePackage, myProbes, c("SYMBOL", "ENTREZID", "ENSEMBL"))
  annotatedTopTab<- merge(x=geneAnots, y=topTab, by.x="PROBEID", by.y="PROBEID")
  return(annotatedTopTab)
}


## ----annotateTopTables-----------------------------------------------------------------------------
topAnnotated_KOvsWT.BAS <- annotatedTopTable(topTab_KOvsWT.BAS,
                                             anotPackage="mogene20sttranscriptcluster.db")
topAnnotated_KOvsWT.VEN <- annotatedTopTable(topTab_KOvsWT.VEN,
                                             anotPackage="mogene20sttranscriptcluster.db")
topAnnotated_INT <- annotatedTopTable(topTab_INT,
                                      anotPackage="mogene20sttranscriptcluster.db")
write.csv(topAnnotated_KOvsWT.BAS, file="./results/topAnnotated_KOvsWT_BAS.csv")
write.csv(topAnnotated_KOvsWT.VEN, file="./results/topAnnotated_KOvsWT_VEN.csv")
write.csv(topAnnotated_INT, file="./results/topAnnotated_INT.csv")


## ----volcanoPlot, fig.cap="Volcano plots for each of the comparisons with the names of the top 3 genes"}
geneSymbols <- AnnotationDbi::select(mogene20sttranscriptcluster.db, rownames(fit.main), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
par(mfrow=c(1,3))
volcanoplot(fit.main, coef=1, highlight=3, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n"))
abline(v=c(-1,1))
volcanoplot(fit.main, coef=2, highlight=3, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[2], sep="\n"))
abline(v=c(-1,1))
volcanoplot(fit.main, coef=3, highlight=3, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[3], sep="\n"))
abline(v=c(-1,1))
dev.off()


## saving the Volcano plot 
png(filename = "./figures/Figure11_volcanoplot.png")
geneSymbols <- AnnotationDbi::select(mogene20sttranscriptcluster.db, rownames(fit.main), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
par(mfrow=c(1,3))
volcanoplot(fit.main, coef=1, highlight=3, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n"))
abline(v=c(-1,1))
volcanoplot(fit.main, coef=2, highlight=3, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[2], sep="\n"))
abline(v=c(-1,1))
volcanoplot(fit.main, coef=3, highlight=3, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[3], sep="\n"))
abline(v=c(-1,1))
dev.off()


## ----decideTests.1---------------------------------------------------------------------------------
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.1, lfc=1)
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,] 
print(summary(res))

## ---- vennDiagram, fig.cap="Venn diagram showing the genes in common between the three comparisons performed"----
vennDiagram (res.selected[,1:3], cex=0.9)
title("Genes in common between the three comparisons\n Genes selected with FDR < 0.1 and logFC > 1")
dev.off()

## ----vennPlot, echo=FALSE, results='hide'----------------------------------------------------------
png(filename = "./figures/Figure12_vennDiagram.png")
vennDiagram (res.selected[,1:3], cex=0.9)
title("Genes in common between the three comparisons\n Genes selected with FDR < 0.1 and logFC > 1")
dev.off()


## ----data4Heatmap----------------------------------------------------------------------------------
probesInHeatmap <- rownames(res.selected)
HMdata <- exprs(eset_filtered)[rownames(exprs(eset_filtered)) %in% probesInHeatmap,]

geneSymbols <- AnnotationDbi::select(mogene20sttranscriptcluster.db, rownames(HMdata), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
rownames(HMdata) <- SYMBOLS
write.csv(HMdata, file = file.path("./results/data4Heatmap.csv"))

##--------heatmapClustering, fig.cap="Heatmap for expression data grouping genes (rows) and samples (columns) by their similarity"-------------------------
library(gplots)
my_palette <- colorRampPalette(c("blue", "red"))(n = 299)
heatmap.2(HMdata,Rowv = TRUE, Colv = TRUE, dendrogram = "both", main = "Differentially expressed genes \n FDR < 0,1, logFC >=1", scale = "row", col = my_palette, zepcolor = "white", sepwidth = c(0.05,0.05), cexRow = 0.5,           cexCol = 0.9, key = TRUE, keysize = 1.5, density.info = "histogram", ColSideColors = c(rep("red",3),rep("blue",3), rep("green",3), rep("yellow",3)), tracecol = NULL, srtCol = 30)
dev.off()

##----------------saveresults,echo=FALSE-------------------
png(filename = "./figures/Figure13_heatmapGroups.png")
heatmap.2(HMdata,Rowv = TRUE, Colv = TRUE, dendrogram = "both", main = "Differentially expressed genes \n FDR < 0,1, logFC >=1", scale = "row", col = my_palette, zepcolor = "white", sepwidth = c(0.05,0.05), cexRow = 0.5,           cexCol = 0.9, key = TRUE, keysize = 1.5, density.info = "histogram", ColSideColors = c(rep("red",3),rep("blue",3), rep("green",3), rep("yellow",3)), tracecol = NULL, srtCol = 30)
dev.off()

## ----selectGenes-----------------------------------------------------------------------------------
listOfTables <- list(KOvsWT.BAS = topTab_KOvsWT.BAS, 
                     KOvsWT.VEN  = topTab_KOvsWT.VEN, 
                     INT = topTab_INT)
listOfSelected <- list()
for (i in 1:length(listOfTables)){
  # select the toptable
  topTab <- listOfTables[[i]]
  # select the genes to be included in the analysis
  whichGenes<-topTab["adj.P.Val"]<0.15
  selectedIDs <- rownames(topTab)[whichGenes]
  # convert the ID to Entrez
  EntrezIDs<- select(mogene20sttranscriptcluster.db, selectedIDs, c("ENTREZID"))
  EntrezIDs <- EntrezIDs$ENTREZID
  listOfSelected[[i]] <- EntrezIDs
  names(listOfSelected)[i] <- names(listOfTables)[i]
}
sapply(listOfSelected, length)


##------------------------- BiologicalSig-----------------------------------------------------------------
library(clusterProfiler)
listOfData <- listOfSelected[1:3]
comparisonsNames <- names(listOfData)


for (i in 1:length(listOfData)){
  genesIn <- listOfData[[i]]
  comparison <- comparisonsNames[i]
  enrich.result <- enrichGO(gene = genesIn,
                            pvalueCutoff = 0.05,
                            keyType       = 'ENTREZID',
                            readable = T,
                            pAdjustMethod = "BH",
                            OrgDb         = org.Mm.eg.db,
                            ont           = "BP")
  
  cat("##################################")
  cat("\nComparison: ", comparison,"\n")
  print(head(enrich.result))
  
  if (length(rownames(enrich.result@result)) != 0) {
    write.csv(as.data.frame(enrich.result), 
              file =paste0("./results/","CusterProfiler.Results.",comparison,".csv"), 
              row.names = FALSE)
    
    pdf(file=paste0("./results/ClusterProfile_Barplot.",comparison,".pdf"))
    print(barplot(enrich.result, showCategory = 15, font.size = 4, 
                  title = paste0("EnrichGO Analysis for ", comparison,". Barplot")))
    dev.off()
    
    pdf(file = paste0("./results/ClusterProfile_cnetplot.",comparison,".pdf"))
    print(cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, 
                   vertex.label.cex = 0.75))
    dev.off()
    
  }
}


##---------- networkINT, fig.cap="Network obtained from the enrichment analysis on the list obtained from the comparison between KO and WT in their interaction"-------------------------
cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, vertex.label.cex = 0.75)
dev.off()

##------------network, echo=FALSE------------------------------------------------------------------------------------
png(filename = "./figures/Figure14_cnetplot_int.png")
cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, vertex.label.cex = 0.75)
dev.off()


##------------enrichplotINT, fig.cap="Enrichment analysis of the more relevant 20 GO terms in the Interaction"-----------
emapplot(enrich.result,pie="count", pie_scale=1.5, layout="nicely",color = "p.adjust",showCategory = 20) + ggtitle("Enrichplot for GO terms in the Int.")
dev.off()

##-------------enrichplotINTsave,echo=FALSE--------------------------------------------------------------------------
png(filename = "./figures/Figure15_enrichplotINT.png")
emapplot(enrich.result,pie="count", pie_scale=1.5, layout="nicely",color = "p.adjust",showCategory = 20) + ggtitle("Enrichplot INT for GO terms")
dev.off()


##-------------tableSUM_1, echo=FALSE---------------------------------------------------------------------------------
Tab.react1 <- read.csv2(file.path("./results/CusterProfiler.Results.KOvsWT.BAS.csv"), 
                        sep = ",", header = TRUE, row.names = 1)

Tab.react1 <- Tab.react1[1:5, 1:5]
##knitr::kable(Tab.react1, booktabs = TRUE, caption = "First rows and columns for ClusterProfiler results on KOvsWT.BAS comparison")


##-------------tableSUM_2, echo=FALSE---------------------------------------------------------------------------------
Tab.react2 <- read.csv2(file.path("./results/CusterProfiler.Results.KOvsWT.VEN.csv"), 
                        sep = ",", header = TRUE, row.names = 1)

Tab.react2 <- Tab.react2[1:5, 1:5]
##knitr::kable(Tab.react2, booktabs = TRUE, caption = "First rows and columns for ClusterProfiler results on KOvsWT.VEN comparison")


##-------------tableSUM_3, echo=FALSE---------------------------------------------------------------------------------
Tab.react3 <- read.csv2(file.path("./results/CusterProfiler.Results.INT.csv"), 
                        sep = ",", header = TRUE, row.names = 1)

Tab.react3 <- Tab.react3[1:5, 1:5]
##knitr::kable(Tab.react3, booktabs = TRUE, caption = "First rows and columns for ClusterProfiler results between the KO and WT interaction")










