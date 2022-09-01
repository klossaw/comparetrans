#!/usr/bin/env Rscript

### Input fileis PRJNA197379.csv.Output:
### 1) merged read count file (PRJNA197379.RCmatrix.csv), 
###	2) differential expression results (PRJNA197379.DEresult.csv), 
### 3) normalized read count (PRJNA197379.RCmatrix.norm.csv),
###	4) FPKM matrix (PRJNA197379.RCmatrix.norm.FPKM.csv),
### 5) averaged FPKM (PRJNA197379.RCmatrix.norm.FPKM.mean.csv)

library(DESeq2)


# set working directory
args = commandArgs(TRUE)
setwd(args[1])
print(getwd())


### READ
prjCsv = list.files(getwd(), pattern = ".csv")
print(paste0("===== Project ID: ", prjCsv))
sampleInfo = read.csv(prjCsv, header=TRUE)
rownames(sampleInfo) = sampleInfo[,1]


### To make a column of gene ids that will be merged with readcounts from all SRRs
SRRrcTable = (read.table(paste(sampleInfo[1,1],'.readcounts.txt',sep=''), sep='\t', as.is=T, header=T))
prjRcMatrix = SRRrcTable[1]

### To merge readcounts from all SRRs into the column of gene ids.
print("===== A list of SRRs:")
for (i in 1:nrow(sampleInfo))
{
	print(paste(sampleInfo[i,1]))
	SRRrcTable = (read.table(paste(sampleInfo[i,1]), sep='\t', as.is=T, header=T))
	prjRcMatrix = cbind(prjRcMatrix, SRRrcTable[7])
	colnames(prjRcMatrix)[i+1] = paste(sampleInfo[i,1])
}

### To write merged readcounts into XXX_RCmatrix.csv.
rcMatrixName = paste0(gsub(".csv", "", prjCsv), ".RCmatrix.csv")
write.csv(as.data.frame(prjRcMatrix), file=rcMatrixName, , row.names=F)
print(paste0("===== Merged readcounts file: ", rcMatrixName))


### To perform Differential expression analysis using DESeq2 package
rownames(prjRcMatrix) = prjRcMatrix[,1]
prjRcMatrix = prjRcMatrix[,-1]

dds = DESeqDataSetFromMatrix(countData=prjRcMatrix, colData=sampleInfo, design=~Treatment) 
dds = DESeq(dds)
res = results(dds)
print("===== DE analysis summary:")
print(summary(res))
print(head(res))

### To write Differential expression analysis results into XXX_DEresult.csv.
DEresultName = paste0(gsub(".csv", "", prjCsv), ".DEresult.csv")
write.csv(as.data.frame(res), file=DEresultName, row.names=T)
print(paste0("===== DiffExpres result file: ", DEresultName))


### To normalize read count by size factor   
dds.EstSizFac = estimateSizeFactors(dds)
normRC = counts(dds.EstSizFac, normalized=TRUE)
normRC.round = round(normRC, 3)
 
### To write normalized merged readcounts into XXX_RCmatrix.norm.csv.
normRcName = paste0(gsub(".csv", "", prjCsv), ".RCmatrix.norm.csv")
write.csv(as.data.frame(normRC.round), file=normRcName, row.names=T)
print(paste0("===== Normalized readcounts file: ", normRcName))


### To calculate FPKM using gene length
AnnoData = SRRrcTable[,c("Geneid","Length")]
mcols(dds)$basepairs <- AnnoData[, "Length"]
fpkm(dds)
normFPKM = fpkm(dds)
normFPKM.round = round(normFPKM, 3)

### To write normalized FPKM into XXX_RCmatrix.norm.FPKM.csv.
NormFPKMName = paste0(gsub(".csv", "", prjCsv), ".RCmatrix.norm.FPKM.csv")
write.csv(as.data.frame(normRPKM.round), file=NormFPKMName, row.names=T)
print(paste0("===== Normalized FPKM file: ", NormFPKMName))


### To get average FPKM value for each condition
WholeAveFPKM = data.frame(Geneid=as.matrix(rownames(res)))
#WholeAveFPKM = data.frame(Geneid=res[,1])
print(head(WholeAveFPKM))

print("===== A list of conditions in this proejcts:")
for (cndtn in unique(sampleInfo[,2]))
{ 
	print(cndtn) 
	SRRcndtn = sampleInfo[(which(sampleInfo[,2]==cndtn)),]
	print(SRRcndtn)

	if(dim(SRRcndtn)[1] == 1)
	{
		FPKMcndtn = normRPKM.round[,colnames(normRPKM.round) %in% SRRcndtn[,1]]
		aveFPKMcndtn=as.matrix(FPKMcndtn)
		colnames(aveFPKMcndtn) = cndtn
		WholeAveFPKM = cbind(WholeAveFPKM, round(aveFPKMcndtn, 3))
	}
	else
	{
		FPKMcndtn = normRPKM.round[,colnames(normRPKM.round) %in% SRRcndtn[,1]]
		aveFPKMcndtn = as.matrix(rowMeans(FPKMcndtn))
		colnames(aveFPKMcndtn) = cndtn
		WholeAveFPKM = cbind(WholeAveFPKM, round(aveFPKMcndtn, 3))
	}
}

### To write Averaged FPKM into XXX_RCmatrix.norm.FPKM.mean.csv.
MeanFPKMName = paste0(gsub(".csv", "", prjCsv), ".RCmatrix.norm.FPKM.mean.csv")
write.csv(as.data.frame(WholeAveFPKM), file=MeanFPKMName, row.names=F)
print(paste0("===== Averaged norm FPKM file: ", MeanFPKMName))

