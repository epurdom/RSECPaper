## ----- Setup --------
print(Sys.time())
runClus<-FALSE
removeOthers<-FALSE

library(BiocInstaller) #prints out the bioconductor version, though takes surprisingly long time
## some basic information about the machine
library(benchmarkme) #to get information about RAM and CPUs on machine:
print(get_ram())
print(get_cpu())

library(clusterExperiment)
if(packageVersion("clusterExperiment")<'2.1.5') stop("must have most current version of clusterExperiment to avoid bugs")
library(shape)
print(sessionInfo())

###Directory information
inputDataDir<-"dataInput" 
outputDataDir<-"dataOutput"
if(!dir.exists(inputDataDir)) dir.create(inputDataDir,recursive = TRUE)
if(!dir.exists(outputDataDir)) dir.create(outputDataDir,recursive = TRUE)
	
###Figure settings
figureDir<-"FinalFiguresChen"
if(!dir.exists(figureDir)) dir.create(figureDir,recursive = TRUE)
figHeight<-6
figWidth<-6
plotType<-"png"
makeFigure<-function(name,type=plotType,height=figHeight,width=figWidth){
	if(type=="pdf") pdf(file.path(figureDir,paste(name,".pdf",sep="")), height=height,width=width)
	if(type=="png") png(file.path(figureDir,paste(name,".png",sep="")), height=height*100,width=width*100)
}

urls<-c("https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/chen.rds")

if(runClus){
	sceFile<-file.path(outputDataDir,"chenSCE.rds")
	if(!file.exists(sceFile)){
		library(Seurat)
		if(!file.exists(file.path(inputDataDir,"chen.rds"))){
			print(Sys.time())	
			cat("Starting processing the data:\n")
			st<-proc.time()
			download.file(urls[1], file.path(inputDataDir,"chen.rds"))
		}
		chen <- readRDS("dataInput/chen.rds")
		## Use stringent cutoff to cluster 3,319 cells (min.genes = 2000)
		## seu <- CreateSeuratObject(raw.data = counts(chen), min.genes = 2000)

		## Use all 14,437 cells (min.genes = 0)
		if(!removeOthers) seu <- CreateSeuratObject(raw.data = counts(chen), min.genes = 0)
		else seu <- CreateSeuratObject(raw.data = counts(chen[,chen$cell_type1 != "zothers"]), min.genes = 0)

		## normalization
		seu <- NormalizeData(object = seu)

		## highly variable genes
		seu <- FindVariableGenes(object = seu)
		length(seu@var.genes)

		## scale data
		seu <- ScaleData(object = seu)

		## PCA
		seu <- RunPCA(object = seu, pcs.compute = 50)
		seu <- RunTSNE(object = seu)

		## SCE object with dimensionality reduction
		sce <- chen[,colnames(seu@data)]
		reducedDim(sce, "PCA") <- seu@dr$pca@cell.embeddings
		reducedDim(sce, "TSNE") <- seu@dr$tsne@cell.embeddings

		#Note, don't filter genes or cells, but use the embeddings...

		saveRDS(object=sce,file=sceFile)
	}
	else sce<-readRDS(sceFile)

}

## ----  RSEC ----
#Define parameters for RSEC
rsecFile <- file.path(outputDataDir,'chenRSEC.rds')
paramList<-list(k0s=seq(5,35,by=5),alphas = c(0.1,.3), betas = c(0.9), reduceMethod = "PCA",
			nReducedDims = c(50),clusterFunction = "hierarchical01", minSizes=5,
	consensusProportion = 0.7,consensusMinSize=20,mergeMethod = "adjP", mergeCutoff = 0.1, dendroReduce="PCA",dendroNDims=50)


if(runClus) {
	#### !!! System specific 
	#### will default to run with cores=6 but user will want to adjust
	ncores<-try(as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK')))
	if(inherits(ncores,"try-error") || is.na(ncores)) ncores<-6
	if(ncores> get_cpu()$no_of_cores) ncores<-get_cpu()$no_of_cores/2
	paramList<-c(paramList,list(ncores=ncores-1,random.seed=57892125))
	cat("Number of cores:",ncores,"\n")
	print(Sys.time())	
	cat("Time to cluster data (seconds):\n")
	print(system.time(ceObj <- do.call("RSEC",c(paramList,list(x=sce, isCount=TRUE, subsampleArgs = list(resamp.num=100, clusterFunction="kmeans", clusterArgs=list(nstart=1)),  verbose=TRUE )))))
	print(Sys.time())	
	saveRDS(object=ceObj, file = rsecFile)
}else{
	if(file.exists(rsecFile)) ceObj<-readRDS(file=rsecFile) else stop("file ", rsecFile, " is missing")
	if(!validObject(ceObj)){
		ceObj<-updateObject(ceObj)
	}  #check still valid object in case the package has been updated and changed the class definition...		
}

nClusterMany<-sum(clusterTypes(ceObj)=="clusterMany")
whClusterMany<-which(clusterTypes(ceObj)=="clusterMany")
whConsensus<-which(clusterTypes(ceObj)=="makeConsensus")
whMerge<-which(clusterTypes(ceObj)=="mergeClusters")
cat(paste("Number of total clusters from clusterMany:",nClusterMany,"\n"))
cat(paste("Number of samples:",ncol(ceObj),"\n"))
cat(paste("Number of genes:",nrow(ceObj),"\n"))


#-----------
# make information about clusters from Chen paper as part of object
#-----------
# Notice that
# 1) We do not have "Endo1" or "Endo1" but we do have "Epith1" and "Epith2". Assuming these are the same
# 2) We also do NOT have a cluster "NFO" (newly formed oligodendrocyte) and we DO have two clusters they don't mention: "SCO" and "IMO". Based on their supplemental spreadsheet, which mentions "IMO", "IMO" must be at least a subset of "NFO", since spreadsheet lists "Fyn" as a top marker like "NFO" in main paper, so we treat it as the NFO cluster. It is clear that SCO is a Sox9+ cluster, based on expression of Sox9+ in these cells, but unlikely "NFO" we do not include it in our analysis as part of this marker group because we cannot link it with an existing cluster. 
x<-as.character(colData(ceObj)$cell_type1)
x[x=="zothers"]<- "-1"
xCondensed<-x
xCondensed[grep("GABA",x)]<-"GABA"
xCondensed[grep("Glu",x)]<-"Glu"
xCondensed[grep("Epith",x)]<-"Endothelial" #Cldn5+
xCondensed[x %in% c("POPC","OPC","IMO","MO")]<-"Oligodend" #i.e. Oligo1+; see note on "IMO" above
xCondensed[x %in% c("Astro","Ependy","Tany")]<-"Sox9+" #Sox9+
xCondensed[x %in% c("Micro","Macro")]<-"C1qa+" #C1qa+

xNeuronal<-as.character(xCondensed %in% c("GABA","Glu","Hista"))
xNeuronal[xNeuronal=="TRUE"]<-"Neu"
xNeuronal[xNeuronal=="FALSE"]<-"NonNeu"
xNeuronal[x=="-1"]<-"-1"

ceObj<-addClusterings(ceObj,x,clusterLabel="CellType")
ceObj<-addClusterings(ceObj,xCondensed,clusterLabel="CellTypeCondensed")
ceObj<-addClusterings(ceObj,xNeuronal,clusterLabel="CellTypeNeuronal")

mainMarkers<-c("Snap25","Syt1","Olig1","Cldn5","Sox9","C1qa")
mainTypes<-c("Neuronal","Neuronal","Oligodendrocyte","Endothelial","Astro/Ependy/Tandy","Micro/Macro")
neuMarkers<-c("Snap25","Syt1","Slc17a6","Slc32a1","Gad1","Gad2","Sst")
neuTypes<-c("Neuronal","Neuronal","Glu","GABA","GABA","GABA","GABA")
nneuMarkers<-c("Olig1","Top2a","Pdgfra","Mobp","Fyn",
	"Cldn5","Slc38a5","Myh11",
	"Sox9","Agt","Ccdc153","Rax",
	"C1qa","Cx3cr1","Mrc1")
nneuTypes<-c("Oligodendro","POPC","OPC","MO","IMO",
	"Endothelial","Endo1","Endo2",
	"Astro/Ependy/Tandy","Astro","Ependy","Tany",
	"Micro/Macro","Micro","Macro")
endoMarkers<-c("Cldn5","Slc38a5","Myh11")
endoTypes<-c("Endothelial","Endo1","Endo2")

#-----------
#based on markers (below), define my clusters as neuronal or not (note, code sensitive to numbering changes to clusters in RSEC...)
#-----------
myNeur<-paste("m",sort(c(14,17,20,21,25,28,3,29,10,15,19,22,24,27,30,4,9,31)),sep="")
myNeuronal<-rep("",ncol(ceObj))
myNeuronal[clusterMatrixNamed(ceObj)[,"mergeClusters"] %in% myNeur]<-"Neu"
myNeuronal[!clusterMatrixNamed(ceObj)[,"mergeClusters"] %in%  c("-1",myNeur)]<-"NonNeu"
myNeuronal[clusterMatrixNamed(ceObj)[,"mergeClusters"] %in%  c("-1")]<-"-1"
ceObj<-addClusterings(ceObj,myNeuronal,clusterLabel="myNeuronal")

#-----------
# Align colors
#-----------
#align the colors across all of the clusterings, including clusterMany...
#first align merge and combine
rmColors<-sapply(c("black"),grep,massivePalette)
cellTypeColors<-c("GABA"= "#1F78B4", "Glu"="#FF7F00" , "Endothelial"="#6A3D9A", "Oligodend"="#33A02C","Sox9+"="#E31A1C","C1qa+"="#A6CEE3","Hista"="#bd18ea","SCO"="Cyan")
ceObj<-recolorClusters(ceObj,whichCluster="CellTypeCondensed",value=cellTypeColors)
ceObj<-plotClusters(ceObj,resetColors = TRUE,colPalette=massivePalette[-rmColors], whichClusters=c("CellTypeCondensed","CellType","mergeClusters","makeConsensus"),plot=FALSE,existingColors="firstOnly")
ceObj<-plotClusters(ceObj,resetColors = TRUE,colPalette=massivePalette[-rmColors], existingColors="firstOnly", whichClusters=c(whConsensus,whClusterMany),matchToTop=TRUE,plot=FALSE)


#-----------
# Neuronal markers in clusters
#-----------
##Note, this is just slightly circular, because want to order final picture by the neurons versus not...
##But ran this code for determining which neuronal as well

plotMarkerSet<-function(obj,genes,type,whichCluster,ylim=NULL,...){
	par(las=2)
	bxpOut<-list()
	for(i in 1:length(genes)){
		if(!is.null(ylim)){
			bxpOut[[i]]<-plotFeatureBoxplot(obj,feature=genes[i],whichCluster=whichCluster,main=sprintf("%s (%s marker)",genes[i],type[i]),ylim=ylim,...)
			
		}
		else{
			bxpOut[[i]]<-plotFeatureBoxplot(obj,feature=genes[i],whichCluster=whichCluster,main=sprintf("%s (%s marker)",genes[i],type[i]),...)
			
		}
	
	}
	invisible(bxpOut)
}
neuBxp<-plotMarkerSet(ceObj,genes=neuMarkers[1:2],type=neuTypes[1:2],whichCluster="mergeClusters")
#reorder so neuronal on one size; if relabel the clusters from merge, wouldn't have to do this.
neurInd<-match(myNeur,neuBxp[[1]]$names)
nneurInd<-c(1:length(neuBxp[[1]]$names))[-neurInd]
newOrder<-c(neurInd,nneurInd)
makeOrder<-function(bxpObj,order){
	#note that order may be shorter than ngroups...but should be unique
	order<-unique(order) #just in case
	ngroups<-ncol(bxpObj$stats)
	#reorder the matrices/vectors
	newObj<-lapply(bxpObj,function(x){
		if(!is.null(dim(x))) x[,order]
		else return(x[order]) 
	})
	#deal with outliers/ grouping information
	m<-match(bxpObj$group,order)
	newObj$group<-m
	newObj$out<-bxpObj$out
	#change name of "-1"
	whUn<-grep("-1",newObj$names)
	if(length(whUn)>0) newObj$names[whUn]<-"Unassigned"
	#newObj$names<-paste(newObj$names,paste0("(n=",newObj$n,")"),sep="\n")
	return(newObj)

}
neuBxpOrdered<-lapply(neuBxp,makeOrder,order=newOrder)
makeFigure("plotFeatureBoxplot_neuronalMarkers_all",height=8,width=14)
par(mfrow=c(2,1),cex.main=2,cex.axis=1.7,las=2,mar=c(10.1,5.1,4.1,1.1))
bxp(neuBxpOrdered[[1]],boxfill=neuBxpOrdered[[1]]$colors,main=neuMarkers[1])
bxp(neuBxpOrdered[[2]],boxfill=neuBxpOrdered[[2]]$colors,main=neuMarkers[2])
dev.off()

#-----------
# Difference between neuronal and non-neuronal classification
#-----------
whDiffer<-which(clusterMatrixNamed(ceObj)[,"myNeuronal"]!=clusterMatrixNamed(ceObj)[,"CellTypeNeuronal"])
diffObj<-ceObj[,whDiffer]
diffAssigned<-removeUnassigned(removeUnassigned(diffObj,whichCluster="mergeClusters"),whichCluster="CellType")
tableClusters(diffAssigned,whichClusters=c("myNeuronal","CellTypeNeuronal"))
# > tableClusters(diffAssigned,whichClusters=c("myNeuronal","CellTypeNeuronal"))
#           CellTypeNeuronal
# myNeuronal Neu NonNeu
#     Neu      0     56
#     NonNeu  11      0


makeFigure("plotFeatureScatter_diffNeu_neu",height=6,width=9)
plotFeatureScatter(subsetByCluster(diffAssigned,whichCluster="myNeuronal",value="Neu"),features=mainMarkers,whichCluster="myNeuronal",jitterFactor=3,cex=2)
dev.off()
makeFigure("plotFeatureScatter_diffNeu_nneu",height=6,width=9)
plotFeatureScatter(subsetByCluster(diffAssigned,whichCluster="myNeuronal",value="NonNeu"),features=mainMarkers,whichCluster="myNeuronal",jitterFactor=3,cex=2)
dev.off()



#-----------
##Table how many not assigned
#-----------
cat("Number not assigned by Chen:",sum(clusterMatrixNamed(ceObj)[,"CellType"]=="-1"),"\n")
cat("Number not assigned by RSEC:",sum(clusterMatrixNamed(ceObj)[,"mergeClusters"]=="-1"),"\n")
cat("Number not assigned by both RSEC and Chen:",sum(clusterMatrixNamed(ceObj)[,"mergeClusters"]=="-1" & clusterMatrixNamed(ceObj)[,"CellType"]=="-1"),"\n")

#-----------
# Basic Plot Cluster overlap, all
#-----------
makeFigure("plotClusters_major",height=6,width=12)
par(mfrow=c(1,1),mar=c(1,15,2,.1),cex=2,cex.axis=1.7)
whCl<-rev(c("mergeClusters","CellTypeCondensed"))
lab<-rev(c("RSEC","Major Cell Types\n(Chen et al)"))
plotClusters(ceObj, whichClusters=whCl,axisLine=-1,
	clusterLabels=lab,existingColors="all")
dev.off()
makeFigure("legend_MajorCellTypes",height=4,width=4)
par(cex=2,mar=c(.1,.1,.1,.1))
plotClusterLegend(ceObj,
whichCluster="CellTypeCondensed",title="Major Subtypes\n(Chen et al)",
bty="n")
dev.off()

makeFigure("plotClusters_all",height=6,width=12)
par(mfrow=c(1,1),mar=c(1,15,2,.1),cex=2,cex.axis=1.7)
whCl<-c("mergeClusters","CellType")
lab<-c("RSEC","Clusters of Chen et al")
plotClusters(ceObj, whichClusters=whCl,axisLine=-1,
	clusterLabels=lab,existingColors="all")
dev.off()

col<-colorRampPalette(c("white","black"))(12)
cLab<-1.5
makeFigure("plotTableClusters_CellType_margin2",height=6,width=12)
par(mfrow=c(1,1),mar=c(1,15,2,.1),cex=2,cex.axis=1.7)
whCl<-c("CellType","mergeClusters")
plotClustersTable(ceObj, ignoreUnassigned = TRUE,whichClusters=whCl,margin=2,annLegend = FALSE,colorScale=col,cexRow=cLab,cexCol=cLab)
dev.off()
makeFigure("plotTableClusters_CellTypeCondensed_margin2",height=6,width=12)
par(mfrow=c(1,1),mar=c(1,15,2,.1),cex=2,cex.axis=1.7)
whCl<-c("CellTypeCondensed","mergeClusters")
plotClustersTable(ceObj, ignoreUnassigned = TRUE,whichClusters=whCl,margin=2,annLegend = FALSE,colorScale=col,cexRow=cLab,cexCol=cLab)
dev.off()

makeFigure("plotTableClusters_CellType_bubble_margin2",height=8,width=12)
par(mfrow=c(1,1),mar=c(8,8,8,.1),cex=2,cex.axis=1.7)
whCl<-c("CellType","mergeClusters")
plotClustersTable(ceObj, ignoreUnassigned = TRUE, whichClusters=whCl, margin=2, plotType="bubble",ylab="",xlab="",maxCex=5)
dev.off()
makeFigure("plotTableClusters_CellTypeCondensed_bubble_margin2",height=8,width=12)
par(mfrow=c(1,1),mar=c(8,8,8,.1),cex=2,cex.axis=1.7)
whCl<-c("CellTypeCondensed","mergeClusters")
plotClustersTable(ceObj, ignoreUnassigned = TRUE, whichClusters=whCl, margin=2, plotType="bubble",ylab="",xlab="",maxCex=5)
dev.off()


#-----------
# GABA/Glu split
#-----------
neuAllObj<-ceObj[,xNeuronal=="Neu" | myNeuronal=="Neu"]
makeFigure("plotClusters_neuronalAll",height=6,width=12)
par(mfrow=c(1,1),mar=c(1,15,2,.1),cex=2,cex.axis=1.7)
whCl<-c("mergeClusters","CellTypeCondensed")
lab<-c("RSEC","Glu/GABA\n(Chen et al)")
plotClusters(neuAllObj, whichClusters=whCl,axisLine=-1,
	clusterLabels=lab,existingColors="all")
dev.off()
makeFigure("legend_GluGab",height=4,width=4)
par(cex=2,mar=c(.1,.1,.1,.1))
plotClusterLegend(ceObj[,xNeuronal=="Neu"],
whichCluster="CellTypeCondensed",title="Neuronal Subtypes\n(Chen et al)",
bty="n")
dev.off()

#-----------
# neuronal cluster zoom
#-----------
mNeur<-subsetByCluster(ceObj,whichCluster="mergeClusters",value="m30")
tableClusters(mNeur,whichCluster="CellTypeCondensed")
dim(mNeur)

#plot of neuronal markers for cluster
wh<-which(clusterMatrixNamed(mNeur)[,"CellTypeCondensed"] %in% c("Endothelial","Hista"))
for(i in 1:length(neuMarkers)){
	makeFigure(sprintf("plotFeatureBoxplot_mNeurZoom_%s_%s",neuTypes[i],neuMarkers[i]),"png",width=6,height=6)
	par(las=2,mar=c(10.1,5.1,1.1,1.1),cex.axis=1.7,cex.lab=2)
	plotFeatureBoxplot(mNeur[,-wh],whichCluster="CellTypeCondensed",main="",feature=neuMarkers[i],varwidth=TRUE,ylab="log Gene expression")
	dev.off()	
}


#-----------
##Endothelial
#-----------

whEndoCl<-c("m13","m16","m32","m5","m6","m8")
endoAllObj<-ceObj[,xCondensed=="Endothelial" | clusterMatrixNamed(ceObj)[,"mergeClusters"] %in% whEndoCl]
makeFigure("plotClusters_endoAll",height=6,width=12)
par(mfrow=c(1,1),mar=c(1,15,2,.1),cex=2,cex.axis=1.7)
whCl<-c("mergeClusters","CellType")
lab<-c("RSEC","Endo1/Endo2\n(Chen et al)")
plotClusters(endoAllObj, whichClusters=whCl,axisLine=-1,
	clusterLabels=lab,existingColors="all")
dev.off()
makeFigure("legend_Endo",height=4,width=4)
par(cex=2,mar=c(.1,.1,.1,.1))
plotClusterLegend(ceObj[,xCondensed=="Endothelial"],
whichCluster="CellType",title="Endothelial Subtypes\n(Chen et al)",
bty="n")
dev.off()
makeFigure("legend_EndoRSEC",height=4,width=4)
par(cex=2,mar=c(.1,.1,.1,.1))
plotClusterLegend(ceObj[,clusterMatrixNamed(ceObj)[,"mergeClusters"] %in% whEndoCl],
whichCluster="mergeClusters",title="RSEC",
bty="n")
dev.off()

endoObj<-subsetByCluster(ceObj,whichCluster="CellTypeCondensed",value="Endothelial")
tableClusters(endoObj,whichCluster="mergeClusters")
tableClusters(subsetByCluster(ceObj,whichCluster="mergeClusters",value=whEndoCl),whichCluster="mergeClusters")

prop.table(tableClusters(endoObj,whichCluster=c("mergeClusters","CellType")),margin=1)
ylimit<-c(5,3,6)
for(i in 1:length(endoMarkers)){
	makeFigure(sprintf("plotFeatureBoxplot_endo_%s_%s",endoMarkers[i],"theirs"),"png",width=6,height=6)
	par(las=2,mar=c(10.1,5.1,1.1,1.1),cex.axis=1.7,cex.lab=2)
	plotFeatureBoxplot(endoObj,whichCluster="CellType",main="",feature=endoMarkers[i],ylab="log Gene expression",ylim=c(0,ylimit[i]))
	dev.off()	
	makeFigure(sprintf("plotFeatureBoxplot_endo_%s_%s",endoMarkers[i],"ours"),"png",width=6,height=6)
	par(las=2,mar=c(10.1,5.1,1.1,1.1),cex.axis=1.7,cex.lab=2)
	plotFeatureBoxplot(subsetByCluster(endoObj,whichCluster="mergeClusters",value=whEndoCl),whichCluster="mergeClusters",main="",feature=endoMarkers[i],ylab="log Gene expression",ylim=c(0,ylimit[i]))
	dev.off()	
}
