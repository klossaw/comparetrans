#!/usr/bin/env Rscript
### To generate input files of Cytoscape to viaulization OrthoClust results.
### Input files
### - GMX_FPKM.csv and ATH_FPKM.csv: these files are expression matrices.
### - ARATH2GLYMA.RBH.example.txt: this file include pairs of RBH genes from two species.
### - Orthoclust_Results.RData: output file from Section3.4.Step1_OrthoClust.R. 
### Output files
### - for Cytoscape: Cytoscape_Input-edge_[ATH|GMX|RBH].csv
### - for User: ExpressionProfile_Module8.csv, ExpressionProfile_Module8.pdf

HOME_PATH=Sys.getenv("HOME")
PROJ_PATH=paste0(HOME_PATH,"/ATH_GMA/results")
setwd(PROJ_PATH)


### Loading RData saved from the previous step. 
RDataFile="Orthoclust_Results.RData"
lnames = load(file = RDataFile);
GMX_modules=orthoclust_results[[1]]
ATH_modules=orthoclust_results[[2]]

### Reading Gene expression profiles to extract expression values for a module
AthGeFile="ATH_FPKM.csv"
ATH_GE=read.table(AthGeFile, sep = ",", header = T, stringsAsFactors= T, row.names =1, check.names=FALSE)

GmxGeFile="GMX_FPKM.csv"
GMX_GE=read.table(GmxGeFile, sep = ",", header = T, stringsAsFactors= T, row.names =1, check.names=FALSE)


### There would be many modules, but we will explore Module 8 as an example. 
ModulesOfInterest= ortho_module_list_sumOrdered[c(1:10),1:6]
print(ModulesOfInterest)
ModuleName=8


### Extracting datasets for Module 8
	
### Extracting a list of genes in Module 8 for soybean and Arabidopsis respectively
GMX_moduleGenes=names(GMX_modules[which(GMX_modules==ModuleName)])
ATH_moduleGenes=names(ATH_modules[which(ATH_modules==ModuleName)])

### Extracting gene expression values for the genes in Module 8
GMX_moduleGene_GEs=GMX_GE[which(rownames(GMX_GE) %in% GMX_moduleGenes),]
GMX_moduleGene_GEs=GMX_moduleGene_GEs[apply(GMX_moduleGene_GEs, 1, function(row) any(row !=0 )), ]
ATH_moduleGene_GEs=ATH_GE[which(rownames(ATH_GE) %in% ATH_moduleGenes),]
ATH_moduleGene_GEs=ATH_moduleGene_GEs[apply(ATH_moduleGene_GEs, 1, function(row) any(row !=0 )), ]

### Scaling and cetering the extracted gene expression values by genes.
GMX_moduleGene_scaledGEs= t(scale(t(GMX_moduleGene_GEs), center = T, scale = T))
GMX_moduleGene_scaledGEs=na.omit(GMX_moduleGene_scaledGEs)
ATH_moduleGene_scaledGEs= t(scale(t(ATH_moduleGene_GEs), center = T, scale = T))
ATH_moduleGene_scaledGEs=na.omit(ATH_moduleGene_scaledGEs)

### Writing csv file with original and scaled expression leve of genes in Module 8.
GEsOutFile="ExpressionProfile_Module8.csv"
write.table(rbind(colnames(GMX_moduleGene_GEs), GMX_moduleGene_GEs), GEsOutFile, append=F,col.names=F, sep=",")
write.table(rbind(colnames(GMX_moduleGene_scaledGEs), GMX_moduleGene_scaledGEs) , GEsOutFile,append = T,col.names=F, sep=",")
write.table(rbind(colnames(ATH_moduleGene_GEs), ATH_moduleGene_GEs) , GEsOutFile,append = T,col.names=F, sep=",")
write.table(rbind(colnames(ATH_moduleGene_scaledGEs), ATH_moduleGene_scaledGEs) , GEsOutFile,append = T,col.names=F, sep=",")


### Preparing datasets for gene expression plots for Module 8. 

### Gene expression profile for genes from co-expression networks
GMX_moduleGene_GEs_CoExNet= GMX_moduleGene_GEs[which(rownames(GMX_moduleGene_GEs) %in% c(as.character(GMX_edgelist[,1]), as.character(GMX_edgelist[,2]))),]
GMX_moduleGene_scaledGEs_CoExNet= t(scale(t(GMX_moduleGene_GEs_CoExNet), center = T, scale = T))			
ATH_moduleGene_GEs_CoExNet= ATH_moduleGene_GEs[which(rownames(ATH_moduleGene_GEs) %in% c(as.character(ATH_edgelist[,1]), as.character(ATH_edgelist[,2]))),]
ATH_moduleGene_scaledGEs_CoExNet= t(scale(t(ATH_moduleGene_GEs_CoExNet), center = T, scale = T))

### Gene expression profile for genes from the RBH list
GMX_moduleGene_GEs_RBH= GMX_moduleGene_GEs[which(rownames(GMX_moduleGene_GEs) %in% c(as.character(GA_orthologs[,1]))),]
GMX_moduleGene_scaledGEs_RBH= t(scale(t(GMX_moduleGene_GEs_RBH), center = T, scale = T))
ATH_moduleGene_GEs_RBH= ATH_moduleGene_GEs[which(rownames(ATH_moduleGene_GEs) %in% c(as.character(GA_orthologs[,2]))),]
ATH_moduleGene_scaledGEs_RBH= t(scale(t(ATH_moduleGene_GEs_RBH), center = T, scale = T))


### Drawing gene expression plots for Module 8. 

### Setting a Y-axis range for both plots 
Ymin=min(c( hist(c(t(ATH_moduleGene_scaledGEs)),plot=F)$breaks, hist(c(t(GMX_moduleGene_scaledGEs)), plot=F)$breaks))
Ymax=max(c( hist(c(t(ATH_moduleGene_scaledGEs)),plot=F)$breaks, hist(c(t(GMX_moduleGene_scaledGEs)), plot=F)$breaks))
Yrange=seq(Ymin,Ymax, 0.5)

### Naming a plot name and saving the plot
OutPdfFile=gsub(".csv", ".pdf", GEsOutFile)
pdf(OutPdfFile, width=15,height=8)

### Setting a margin of plot
par(mfrow=c(1,2))
par(mar=c(11.1, 5.1, 2.1, 0.0), xpd=TRUE)
op <- par(cex = 1.5)	

### Plotting Soybean gene expression data
par(oma=c(0.1, 0, 0, 0.1))
matplot(t(GMX_moduleGene_scaledGEs), ylim=c(Ymin, Ymax), ylab="Normalized Expression", type=c("l"), lty=1, pch=1, col=8, axes=F, lwd=0.5)
matplot(t(GMX_moduleGene_scaledGEs_RBH), ylim=c(Ymin, Ymax), ylab="", type=c("l"), lty=1, pch=20, col=4, axes=F, add=TRUE, lwd=0.5)
matplot(colMeans(GMX_moduleGene_scaledGEs), ylim=c(Ymin, Ymax), ylab="", type=c("l"), lty=1, pch=1, col=1, axes=F, lwd=5, add=TRUE)
matplot((GMX_moduleGene_scaledGEs[rownames(GMX_moduleGene_scaledGEs)=="Glyma.07G102900", ]), ylim=c(Ymin, Ymax), ylab="", type = c("b"), lty=1, pch=4, col = 3, axes=F,add = TRUE,lwd = 3) 	
matplot((GMX_moduleGene_scaledGEs[rownames(GMX_moduleGene_scaledGEs)=="Glyma.04G245100", ]), ylim=c(Ymin, Ymax), ylab="", type = c("b"), lty=1, pch=4, col = 3, axes=F,add = TRUE,lwd = 3) 	
par(las=2)
axis(side=2,at = Yrange)
axis(side=1,at=1:ncol(GMX_moduleGene_scaledGEs),labels=colnames(GMX_moduleGene_scaledGEs), cex=1)

### Plotting Arabidopsis gene expression data
par(oma=c(0.1, 0, 0, 0.1))
matplot(t(ATH_moduleGene_scaledGEs), ylim=c(Ymin, Ymax), ylab="Normalized Expression", type=c("l"), lty=1, pch=1, col=8, axes=F, add=TRUE, lwd=0.5)
matplot(t(ATH_moduleGene_scaledGEs_RBH), ylim=c(Ymin, Ymax), ylab="", type=c("l"), lty=1, pch=20, col=4, axes=F, add=TRUE, lwd=0.5)
matplot(colMeans(ATH_moduleGene_scaledGEs), ylim=c(Ymin, Ymax), ylab="", type=c("l"), lty=1, pch=1, col=1, axes=F, lwd=5, add=TRUE)
matplot((ATH_moduleGene_scaledGEs[rownames(ATH_moduleGene_scaledGEs)=="AT3G17940", ]), ylim=c(Ymin, Ymax), ylab="", type = c("b"), lty=1, pch=4, col = 3, axes=F,add = TRUE,lwd = 3) 	
matplot((ATH_moduleGene_scaledGEs[rownames(ATH_moduleGene_scaledGEs)=="AT5G52560", ]), ylim=c(Ymin, Ymax), ylab="", type = c("b"), lty=1, pch=4, col = 3, axes=F,add = TRUE,lwd = 3) 	
par(las=2)
axis(side=2,at = Yrange)
axis(side=1,at=1:ncol(ATH_moduleGene_scaledGEs),labels=colnames(ATH_moduleGene_scaledGEs), cex=1)

### Adding a legend and titles of plots
par(xpd=NA)
legend(-11.5,-5.5, legend = c("Cluster Average Expression", "Homologous RFO genes", "Homogous Genes", "Non-homologous Genes"), col=c("black", "green", "blue", "grey"), lwd=5, cex=0.7, horiz = TRUE)
legend(-9,4.2, legend = "", col="white",lwd=NULL, cex=2.2, horiz = TRUE, title="Soybean", bty = "n")
legend(1.4,4.2, legend = "", col="white",lwd=NULL, cex=2.2, horiz = TRUE, title="Arabidopsis", bty = "n")
dev.off()


### Extracting edges for Cytoscape input filess

### Extracting edges for GMX, both genes from coexpression
GMX_edgelist_TMP= GMX_edgelist[GMX_edgelist[,1] %in% GMX_moduleGenes,]
GMX_moduleGenes_EdgeList_Restrict= GMX_edgelist_TMP[GMX_edgelist_TMP[,2] %in% GMX_moduleGenes,]
GMX_moduleGenes_EdgeList_Restrict=unique( as.data.frame(GMX_moduleGenes_EdgeList_Restrict))
row.names(GMX_moduleGenes_EdgeList_Restrict)=NULL; colnames(GMX_moduleGenes_EdgeList_Restrict) = NULL

### Extracting edges for ATH, both genes from coexpression
ATH_edgelist_TMP= ATH_edgelist[ATH_edgelist[,1] %in% ATH_moduleGenes,]
ATH_moduleGenes_EdgeList_Restrict= ATH_edgelist_TMP[ATH_edgelist_TMP[,2] %in% ATH_moduleGenes,]
ATH_moduleGenes_EdgeList_Restrict=unique( as.data.frame(ATH_moduleGenes_EdgeList_Restrict))
row.names(ATH_moduleGenes_EdgeList_Restrict)=NULL; colnames(ATH_moduleGenes_EdgeList_Restrict) = NULL

### Extracting edges for RHB, both genes from RBH list
Orth_edgelist_TMP= GA_orthologs[GA_orthologs[,1] %in% GMX_moduleGenes,]
Orth_moduleGenes_EdgeList_Restrict= Orth_edgelist_TMP[Orth_edgelist_TMP[,2] %in% ATH_moduleGenes,]
Orth_moduleGenes_EdgeList_Restrict=unique(as.data.frame(Orth_moduleGenes_EdgeList_Restrict))
row.names(Orth_moduleGenes_EdgeList_Restrict)=NULL; colnames(Orth_moduleGenes_EdgeList_Restrict) = NULL

### Naming edge list files
EdgeOutFile="Cytoscape_Input"

### Writing edge list files for homologous relations from RBH list and co-expression networks from Soybean and Arabidopsis respectively. 
GMXEdgeFile=paste0(EdgeOutFile, "-edge_GMX.csv"); 
ATHEdgeFile=paste0(EdgeOutFile, "-edge_ATH.csv"); 
RBHEdgeFile=paste0(EdgeOutFile, "-edge_RBH.csv"); 
			
write.table(GMX_moduleGenes_EdgeList_Restrict , GMXEdgeFile, row.names=F, col.names = F, quote=F, sep = ",")
write.table(ATH_moduleGenes_EdgeList_Restrict , ATHEdgeFile, row.names=F, col.names = F, quote=F, sep = ",")
write.table(Orth_moduleGenes_EdgeList_Restrict, RBHEdgeFile, row.names=F, col.names = F, quote=F, sep = ",")
