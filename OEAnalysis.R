## ----- Setup --------
print(Sys.time())
runClus<-FALSE #run RSEC
runClus2<-FALSE #run comparison of algorithms/parameters...
runClus3<-FALSE #run comparison of distances...
parallelSubcores<-FALSE
removeDoublets<-FALSE

library(BiocInstaller) #prints out the bioconductor version, though takes surprisingly long time

library(clusterExperiment)
if(packageVersion("clusterExperiment")<'2.1.5') stop("must have most current version of clusterExperiment to avoid bugs")
library(shape)
library(scran) #bioconductor
library(igraph) #CRAN


###Directory information
inputDataDir<-"dataInput" 
outputDataDir<-"dataOutput"
if(!dir.exists(inputDataDir)) dir.create(inputDataDir,recursive = TRUE)
if(!dir.exists(outputDataDir)) dir.create(outputDataDir,recursive = TRUE)
	
###Figure settings
figureDir<-"FinalFiguresOE"
if(!dir.exists(figureDir)) dir.create(figureDir,recursive = TRUE)
figHeight<-6
figWidth<-6
plotType<-"pdf"
makeFigure<-function(name,type=plotType,height=figHeight,width=figWidth){
	if(type=="pdf") pdf(file.path(figureDir,paste(name,".pdf",sep="")), height=height,width=width)
	if(type=="png") png(file.path(figureDir,paste(name,".png",sep="")), height=height*100,width=width*100)
}

###Function for getting data from git repository
getFileFromGit<-function(filename,unzip=FALSE){
	gitadd<-"https://raw.githubusercontent.com/rufletch/p63-HBC-diff/master/ref/"
	if(!file.exists(file.path(inputDataDir,filename))) {
	  download.file(paste(gitadd,filename,sep=""), file.path(inputDataDir,filename))
	  if(unzip) R.utils::gunzip(file.path(inputDataDir,filename))
	}
}
urls = c("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE95601&format=file&file=GSE95601%5FoeHBCdiff%5FCufflinks%5FeSet%2ERda%2Egz",
         "https://raw.githubusercontent.com/rufletch/p63-HBC-diff/master/ref/oeHBCdiff_clusterLabels.txt")

## some basic information about the machine
library(benchmarkme) #to get information about RAM and CPUs on machine:
print(get_ram())
print(get_cpu())
print(sessionInfo())

## ---- data ----
if(runClus){
	print(Sys.time())	
	cat("Starting processing the data:\n")
	st<-proc.time()
	if(!file.exists(file.path(inputDataDir,"GSE95601_oeHBCdiff_Cufflinks_eSet.Rda"))) {
	  download.file(urls[1], file.path(inputDataDir,"GSE95601_oeHBCdiff_Cufflinks_eSet.Rda.gz"))
	  R.utils::gunzip(file.path(inputDataDir,"GSE95601_oeHBCdiff_Cufflinks_eSet.Rda.gz"))
	}
	getFileFromGit("oeHBCdiff_clusterLabels.txt")
	# if(!file.exists(file.path(inputDataDir,"oeHBCdiff_clusterLabels.txt"))) {
# 	  download.file(urls[2], file.path(inputDataDir,"oeHBCdiff_clusterLabels.txt"))
# 	}
	load(file.path(inputDataDir,"GSE95601_oeHBCdiff_Cufflinks_eSet.Rda"))
	# Count matrix
	E <- assayData(Cufflinks_eSet)$counts_table
	# Remove undetected genes
	E <- na.omit(E)
	E <- E[rowSums(E)>0,]
	dim(E)
	# Remove ERCC and CreER genes
	cre <- E["CreER",]
	ercc <- E[grep("^ERCC-", rownames(E)),]
	E <- E[grep("^ERCC-", rownames(E), invert = TRUE), ]
	E <- E[-which(rownames(E)=="CreER"), ]
	dim(E)

	# Extract QC metrics
	qc <- as.matrix(protocolData(Cufflinks_eSet)@data)[,c(1:5, 10:18)]
	qc <- cbind(qc, CreER = cre, ERCC_reads = colSums(ercc))

	# Extract metadata
	batch <- droplevels(pData(Cufflinks_eSet)$MD_c1_run_id)
	bio <- droplevels(pData(Cufflinks_eSet)$MD_expt_condition)
	clusterLabels <- read.table(file.path(inputDataDir,"oeHBCdiff_clusterLabels.txt"),
	                            sep = "\t", stringsAsFactors = FALSE)
	m <- match(colnames(E), clusterLabels[, 1])

	# Create metadata data.frame
	metadata <- data.frame("Experiment" = bio,
	                       "Batch" = batch,
	                       "publishedClusters" = clusterLabels[m,2],
	                       qc)

	# Symbol for cells not assigned to a lineage in original data
	metadata$publishedClusters[is.na(metadata$publishedClusters)] <- -2

	se <- SummarizedExperiment(assays = list(counts = E),
	                           colData = metadata)

	#---------
	#Identify those that original paper considered doublets between Sus and Neuronal
	#Appears to have been done on raw E...
	#---------
	omp_cyp <- assays(se)$counts["Omp",] > 200 & (assays(se)$counts["Cyp2g1",] > 200 | assays(se)$counts["Cyp1a2",] > 200)
	reg3g<-assays(se)$counts["Reg3g",]>=100
	#90 of ~140 not clustered were excluded because they were considered doublets...
	table(colData(se)$publishedClusters,omp_cyp | reg3g)
	
	colData(se)$Doublets<- omp_cyp | reg3g
	
							   
	#---------
	# QC-metric-based sample-filtering
	#---------
	library(scone)
	data("housekeeping")
	hk = rownames(se)[toupper(rownames(se)) %in% housekeeping$V1]
	mfilt <- metric_sample_filter(assay(se),
	                              nreads = colData(se)$NREADS,
	                              ralign = colData(se)$RALIGN,
	                              pos_controls = rownames(se) %in% hk,
	                              zcut = 3, mixture = FALSE,
	                              plot = TRUE)
	# Simplify to a single logical
	mfilt <- !apply(simplify2array(mfilt[!is.na(mfilt)]), 1, any)

	#---------
	# gene filtering
	#---------
    thresh_fail <-  40  # Tophat counts
    num_fail <- 5  # cells
	mfilt.gf.vec <- rowSums(E[,mfilt] > thresh_fail) > num_fail
	se <- se[mfilt.gf.vec, mfilt]
  
	#---------
	#Normalization of filtered data with 1st PCA of quality metrics 
	#(following original paper in choosing normalization)
	#---------
	fq <- FQ_FN(assay(se)) #fullQ normalization
	qcpca <- prcomp(colData(se)[,colnames(qc)], center=TRUE, scale=TRUE)$x #PCA of quality metrics
	design_mat <- make_design(bio=NULL, batch=NULL, W=qcpca[,1])
	norm <- lm_adjust(log1p(fq), design_mat) #remove 
	assay(se) <- norm



	ed<-proc.time()
	cat("Time to process data (seconds):\n")
	print(ed-st)
}

## ----  RSEC ----
#Define parameters for RSEC
fn <- file.path(outputDataDir,'rsecObj.rda')
paramList<-list(k0s=4:15,alphas = c(0.1,.2), betas = c(0.8,.9), reduceMethod = "PCA",
			nReducedDims = c(20,50,100),clusterFunction = "hierarchical01", minSizes=c(1,5,10),
	consensusProportion = 0.7,consensusMinSize=10,mergeMethod = "adjP", mergeCutoff = 0.1, dendroReduce="var",dendroNDims=1000)

if(runClus) {
	#### !!! System specific 
	#### will default to run with cores=6 but user will want to adjust
	# ncores<-try(as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK')))
	# if(inherits(ncores,"try-error") || is.na(ncores)) ncores<-6
	# if(ncores> get_cpu()$no_of_cores) ncores<-get_cpu()$no_of_cores/2
		ncores<-14
	if(!parallelSubcores){
		paramList<-c(paramList,list(ncores=ncores-1,random.seed=57892125))
		subcores<-1
	}
	else{
		subcores<-ncores-1
		paramList<-c(paramList,list(ncores=1))
	}
	cat("Number of cores:",ncores,"\n")
	print(Sys.time())	
	cat("Time to cluster data (seconds):\n")
	if(removeDoublets){
		print(system.time(ceObj <- do.call("RSEC",c(paramList,list(x=se[,which(!colData(se)$Doublets)], isCount=FALSE, subsampleArgs = list(resamp.num=100, clusterFunction="kmeans", clusterArgs=list(nstart=10),ncores=subcores),  verbose=TRUE )))))
		
	}
	else{
		print(system.time(ceObj <- do.call("RSEC",c(paramList,list(x=se, isCount=FALSE, subsampleArgs = list(resamp.num=100, clusterFunction="kmeans", clusterArgs=list(nstart=10),ncores=subcores),  verbose=TRUE )))))
		
	}
	print(Sys.time())	
	save(ceObj, file = fn)
}else{
	if(file.exists(fn)) load(fn) else stop("file ", fn, " is missing")
	if(!validObject(ceObj)) stop("not valid CE object") #check still valid object in case the package has been updated...
		
}

## ---- summaryInformation
nClusterMany<-sum(clusterTypes(ceObj)=="clusterMany")
whClusterMany<-which(clusterTypes(ceObj)=="clusterMany")
whConsensus<-which(clusterTypes(ceObj)=="makeConsensus")
whMerge<-which(clusterTypes(ceObj)=="mergeClusters")
cat(paste("Number of total clusters from clusterMany:",nClusterMany,"\n"))
cat(paste("Number of samples:",ncol(ceObj),"\n"))
cat(paste("Number of genes:",nrow(ceObj),"\n"))


## ---- assignColors 
#align the colors across all of the clusterings, including clusterMany...
#first align merge and combine
rmColors<-sapply(c("black"),grep,massivePalette)
ceObj<-plotClusters(ceObj,resetColors = TRUE,colPalette=massivePalette[-rmColors], whichClusters=c("mergeClusters","makeConsensus"),plot=FALSE)
ceObj<-plotClusters(ceObj,resetColors = TRUE,colPalette=massivePalette[-rmColors], existingColors="firstOnly", whichClusters=c(whConsensus,whClusterMany),matchToTop=TRUE,plot=FALSE)


## ----  plotPCA
makeFigure("plotPCA",type="png",height=6,width=6)
par(mfrow=c(1,1))
plotReducedDims(ceObj,whichCluster=c("makeConsensus"),whichDims=c(1,2),legend=FALSE)
dev.off()



## ----  plotClusters_makeConsensusBottom
makeFigure("plotClusters_makeConsensusBottom",height=6,width=6)
nCM<-20
nBlank<-50
blankRows<-c(nCM+nBlank,nCM)
par(mar=c(1,15,2,.1),cex.axis=1.7)
plotClustersWorkflow(ceObj,whichClusters=c("makeConsensus"),axisLine=-1,ylim=c(1,nClusterings(ceObj)),clusterManyLabels=FALSE,clusterLabels="Consensus Clustering",existingColors="all",sortBy="highlighted",highlightOnTop=FALSE,nBlank=nBlank,nSize=nCM)
require(shape)
Arrows(x0=mean(par("usr")[1:2]),y0=blankRows[1]-3*nBlank/10,x1=mean(par("usr")[1:2]),y1=blankRows[2]+3*nBlank/10,lwd=3,arr.type="curved")
mtext("Samples",3,line=-.5,cex=2)
title(ylab="Clusterings (different parameters)",line=0,cex.lab=2)
dev.off()

## ---- plotHier_HierContrasts
#pick 50 from each contrast
hierContrasts<-getBestFeatures(ceObj,contrastType = "Dendro",isCount=FALSE, p.value=0.05, number=50)
rmNodeColors<-match(clusterLegend(ceObj)[["mergeClusters"]][,"color"],massivePalette)
nodeColors<-massivePalette[-c(rmNodeColors,rmColors)][c(1,2,7,8,9,10)]
names(nodeColors)<-as.character(unique(hierContrasts$ContrastName))
makeFigure("plotHeatmap_HierContrasts",width=4)
plotContrastHeatmap(ceObj,hierContrasts,breaks=.99,contrastColors=nodeColors,capBreaksLegend=TRUE,cexRow=0,cexCol=0)
dev.off()
makeFigure("plotHeatmap_HierContrasts_centeredscaled",width=4)
plotContrastHeatmap(ceObj,hierContrasts,contrastColors=nodeColors,visualizeData="centeredAndScaled",breaks=.99,capBreaksLegend=TRUE,cexRow=0,cexCol=0)
dev.off()

## ----  plotMergeClusters
makeFigure("plotMergeClusters",width=4,height=6)
par(mar=c(.1,.1,5.1,.1),cex=2.5)
	plotDendrogram(ceObj,whichClusters=c("makeConsensus","mergeClusters"),mergeInfo="mergeMethod",plotType="colorblock",leafType="sample",legend="none",main="",show.node.label=FALSE,
	clusterLabelAngle=90,use.edge.length =TRUE,cex=.7,nodeColors=nodeColors)
dev.off()

## ---- plotHier_FContrasts
#nmf.options(grid.patch=TRUE)
#origClusterColors<-bigPalette[1:nlevels(colData(ceObj)$publishedClusters)]
nunique<-length(unique(hierContrasts$Feature))
FContrasts<-getBestFeatures(ceObj,contrastType = "F",isCount=FALSE,p.value=0.05, number=nunique)
FContrasts<-FContrasts[order(FContrasts$F,decreasing=TRUE),]
makeFigure("plotHeatmap_FContrasts",width=4)
plotHeatmap(ceObj,clusterFeaturesData=FContrasts$IndexInOriginal,breaks=.99,capBreaksLegend=TRUE,cexRow=0,cexCol=0)
dev.off()
makeFigure("plotHeatmap_FContrasts_centeredscaled",width=4)
plotHeatmap(ceObj,clusterFeaturesData=FContrasts$IndexInOriginal,visualizeData="centeredAndScaled",breaks=.99,capBreaksLegend=TRUE,cexRow=0,cexCol=0)
dev.off()
makeFigure("plotHeatmap_FContrasts_noDendro",width=4)
plotHeatmap(ceObj,clusterFeaturesData=FContrasts$IndexInOriginal,treeheight=c(0,50),breaks=.99,capBreaksLegend=TRUE,cexRow=0,cexCol=0)
dev.off()
makeFigure("plotHeatmap_FContrasts_centeredscaled_noDendro",width=4)
plotHeatmap(ceObj,clusterFeaturesData=FContrasts$IndexInOriginal,treeheight=c(0,50),visualizeData="centeredAndScaled",breaks=.99,capBreaksLegend=TRUE,cexRow=0,cexCol=0)
dev.off()


## ----  comparePropNull
propMat<-data.matrix(ceObj@merge_nodeProp[,-c(1:2)])
ordMethods<-c(2,4,1,3,5,6)
colMethods<-palette()[2:7]
wh<-which(!rownames(propMat) %in% names(nodeColors))
allNodeColors<-massivePalette[-c(rmNodeColors,rmColors)][-c(1,2,7,8,9,10)][1:length(wh)]
names(allNodeColors)<-rownames(propMat)[wh]
allNodeColors<-c(nodeColors,allNodeColors)
makeFigure("comparePropNull",type="png",width=10,height=6)
par(cex.axis=1.3,cex.lab=1.7,mar=c(6,6,4,.1))
out<-barplot(t(propMat)[ordMethods,],beside=TRUE,col=colMethods,legend=FALSE,las=2,ylab="Proportion estimated DE")# ,xaxt="n")
# gpCenter<-colMeans(out)
# sapply(1:nrow(propMat),function(ii){
# 	nm<-rownames(propMat)[ii]
# 	axis(1,at=gpCenter[ii],labels=nm,col.axis=allNodeColors[nm],las=2,tick=TRUE)
# })
legend(x=mean(par("usr")[1:2]),y=par("usr")[4]+.1,colnames(propMat)[ordMethods],fill=colMethods,ncol=3,xpd=NA,xjust=.5,title="Method",cex=1.7)
dev.off()
## ----  plotDendrogramNodes
makeFigure("plotDendrogramNodes",type="png",width=12,height=6)
par(mar=c(.1,.1,5.1,.1))
plotDendrogram(ceObj,whichClusters=c("makeConsensus","mergeClusters"),plotType="colorblock",leafType="sample",legend="none",main="",use.edge.length =TRUE,cex=1.2,nodeColors=allNodeColors)
dev.off()
makeFigure("NodeLegend",type="png",width=8,height=6)
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))
legend("center",names(allNodeColors),fill=allNodeColors,ncol=3,title="Node Legend",cex=2)
dev.off()

## ----  plotCoClustering
makeFigure("plotCoClustering")
par(mar=c(.1,.1,5.1,8.1))
plotCoClustering(ceObj,whichClusters=c("mergeClusters","makeConsensus"),annLegend = FALSE)
dev.off()


## ---- pick subsets of parameters
cmParams<-getClusterManyParams(ceObj)

# PC=50,alpha=0.2,beta=0.9,minSize=1, all ks
subset1<-which(cmParams$alpha==0.2 & cmParams$nReducedDims==50 & cmParams$beta==0.9 & cmParams$minSize==1)

# PC=50,alpha=0.1,0.2,beta=0.9,minSize=1, all ks
subset2<-which( cmParams$nReducedDims==50 & cmParams$beta==0.9 & cmParams$minSize==1)


# all PCs,alpha=0.1,0.2,beta=0.8,0.9,minSize=1, all ks
subset3<-which(cmParams$minSize==1)



redoSubset<-function(object,subset,labelTag){
	ceObjSm<-makeConsensus(object,whichClusters= rownames(cmParams)[subset],proportion=paramList$consensusProportion, minSize=paramList$consensusMinSize,clusterLabel=sprintf("makeConsensus, %s",labelTag))
	ceObjSm<-makeDendrogram(ceObjSm,whichCluster=sprintf("makeConsensus, %s",labelTag),reduceMethod=paramList$dendroReduce, nDims=paramList$dendroNDims)
	ceObjSm<-mergeClusters(ceObjSm,clusterLabel=sprintf("mergeClusters, %s",labelTag),cutoff=paramList$mergeCutoff,mergeMethod=paramList$mergeMethod)
	return(ceObjSm)	
}
ceObjSm<-redoSubset(ceObj,subset=subset1,"subset 1")
ceObjSm<-redoSubset(ceObjSm,subset=subset2,"subset 2")
ceObjSm<-redoSubset(ceObjSm,subset=subset3,"subset 3")
clusterLabels(ceObjSm)[clusterLabels(ceObjSm)=="makeConsensus.1"]<-"makeConsensus, full"
clusterLabels(ceObjSm)[clusterLabels(ceObjSm)=="mergeClusters.1"]<-"mergeClusters, full"

makeFigure("plotClusters_makeConsensus_Subsets",width=8)
subsetClusters<-c("makeConsensus, full","makeConsensus, subset 3","makeConsensus, subset 2","makeConsensus, subset 1")
nCl<- c(nClusterMany,length(subset3),length(subset2),length(subset1))
par(mar=c(1.1,12.1,1.1,1.1),cex.axis=1.5)
plotClusters(ceObjSm,whichClusters=subsetClusters,clusterLabels=paste(c("All","Subset:\nminSize=1","Subset:\nminSize=1, beta=0.9,\nPC=50","Subset:\nminSize=1, beta=0.9,\nPC=50, alpha=0.1"),"\n(",nCl ,")",sep=""),axisLine=-5)
dev.off()

##########################
## Note in what follows that we run these in two stages, and define new objects, so as to be able to run the clustering and update the plots separately.
## Normal usage would keep a single object `ceObj` and continue to save to it.
##########################


##########################
## Different parameters/functions
##########################
## ---- standard clustering methods
fn2 <- file.path(outputDataDir,'rsecObj2.rda')

## Create SNN clustering function
library(scran)
library(igraph)
SNN_wrap <- function(x, k, steps = 4, ...) {
  snn <- scran::buildSNNGraph(x, k = k, d = NA, transposed = FALSE) #scran
  res <- igraph::cluster_walktrap(snn, steps = steps) #igraph
  return(res$membership)
}
#check valid function
internalFunctionCheck(SNN_wrap, inputType = "X", algorithmType = "K",
                      outputType="vector")
SNN <- ClusterFunction(SNN_wrap, inputType = "X", algorithmType = "K",
                     outputType="vector")
										 
## run different clustering techniques, parameters, etc.
if(runClus2) {
	#since K has different meaning between SNN and the other methods, will run it separately, using clusterSingle, and not with clusterMany
  print(system.time(ceObj2<-clusterSingle(ceObj, reduceMethod="PCA",nDims=50,
		subsample=FALSE,sequential=FALSE,
		mainClusterArgs=list(clusterFunction=SNN, kRange=c(5,10,15,20,25), findBestK=TRUE),subsampleArgs=NULL, seqArgs=NULL, clusterLabel="SNN") ))

	print(system.time(ceObj2<-clusterMany(ceObj2, reduceMethod="PCA",nReducedDims=50,
		clusterFunction=c("pam"),
		ks=4:15, findBestK=c(FALSE))
	))
	print(system.time(ceObj2<-clusterMany(ceObj2, reduceMethod="PCA",nReducedDims=50,
		clusterFunction=c("pam","kmeans","spectral","hierarchicalK"),
		ks=4:15, findBestK=c(TRUE))
	))
	save(ceObj2, file = fn2)
}else{
	if(file.exists(fn2)) load(fn2) else stop("file ", fn2, " is missing")
}


##########################
## ---- ClusterDistances
##########################
fn3 <- file.path(outputDataDir,'rsecObj3.rda')
if(runClus3) {
	corDist<-function(x){(1-cor(t(x),method="pearson"))/2}
	spearDist<-function(x){(1-cor(t(x),method="spearman"))/2}
	eucDist<-function(x){dist(x)}
	print(system.time(ceObj3<-clusterMany(ceObj, dimReduce="mad",nVarDims=1000,
		clusterFunction=c("pam"),ks=4:15, distFunction=c("corDist","spearDist","eucDist"))
		))
	save(ceObj3, file = fn3)
}else{
	if(file.exists(fn3)) load(fn3) else stop("file ", fn3, " is missing")
}

##########################
## ---- Plots of those comparisons
##########################
clusterParams2<-getClusterManyParams(ceObj2,searchAll=TRUE,whichClusters=grep("clusterMany",clusterTypes(ceObj2)))
# > names(clusterParams2)
#  [1] "clusteringIndex" "nReducedDims"    "k"               "alpha"           "findBestK"       "beta"
#  [7] "minSize"         "sequential"      "subsample"       "clusterFunction"

## ----  plotClusters_diffK
makeFigure("plotClusters_diffK",height=4)
wh<-which(clusterParams2[,"clusterFunction"]=="pam" & !clusterParams2[,"k"]=="NA" & !clusterParams2[,"subsample"])
par(mar=c(3,10,.1,.1))
plotClusters(ceObj2, whichClusters=clusterParams2[wh,"clusteringIndex"],colPalette = c(bigPalette, rainbow(100)),axisLine=-2,matchToTop=FALSE,clusterLabels=paste("K",clusterParams2[wh,"k"],sep="="))
title(xlab="Samples",line=0,cex.lab=2)
title(ylab="Number of clusters (K)",cex.lab=2)
dev.off()

## ----  plotClusters_diffFunctions
ceTemp<-makeConsensus(ceObj2,whichClusters=wh,proportion=4/6)
ceTemp<-makeDendrogram(ceTemp,whichCluster="makeConsensus",reduceMethod="mad")
ceTemp<-mergeClusters(ceTemp,mergeMethod="adjP",cutoff=0.05,plot=FALSE,logFCcutoff=1,clusterLabel="mergeClusters.logFC")
ceTemp<-mergeClusters(ceTemp,mergeMethod="adjP",cutoff=0.05,plot=FALSE,clusterLabel="mergeClusters.standard")

clusterParamsTemp<-getClusterManyParams(ceTemp,searchAll=TRUE,whichClusters=grep("clusterMany",clusterTypes(ceTemp)))
wh<-clusterParamsTemp[which(clusterParamsTemp[,"findBestK"]),"clusteringIndex"]
wh<-c(wh,grep("SNN",clusterLabels(ceTemp)))
funNames<-c(as.character(clusterParamsTemp[which(clusterParamsTemp[,"findBestK"]),"clusterFunction"]),"NN")
ord<-sapply(c("kmeans","NN","pam","hierarchicalK","spectral"),grep,funNames)
bestK<-nClusters(ceTemp,ignoreUnassigned=TRUE)[wh]

#assign colors:
ceTemp<-plotClusters(ceTemp,whichClusters=c("mergeClusters.logFC","mergeClusters.standard","makeConsensus",clusterLabels(ceTemp)[wh[ord]]),resetColors = TRUE)
#manually adjust one of spectral's cluster colors to make orange (better match):
specCL<-clusterLegend(ceTemp)[[wh[ord][5]]]
whAdj<-which(specCL[,"clusterIds"]=="4")
specCL[whAdj,"color"]<-"#FF7F00"
clusterLegend(ceTemp)[[wh[ord][5]]]<-specCL

makeFigure("plotClusters_diffFunctions",height=4)
par(mar=c(3,12,.1,.1))
nBlank<-1
nCM<-1
plotClustersWorkflow(ceTemp,
	whichClusters=c("mergeClusters.logFC","mergeClusters.standard","makeConsensus"),
	clusterLabels=c("mergeClusters, log FC","mergeClusters, standard","makeConsensus"),
	whichClusterMany=clusterLabels(ceTemp)[wh[ord]],
	existingColors ="all",nBlank=nBlank,nSize=nCM,
	clusterManyLabels=paste(funNames[ord],", K=",bestK[ord],sep=""),
	axisLine=-2,highlightOnTop = TRUE)
blankRows<-c(nCM+nBlank,nCM)
title(xlab="Samples",line=0,cex.lab=2)
title(ylab="Clustering Method",cex.lab=2,line=10)
dev.off()

clusterParams3<-getClusterManyParams(ceObj3)

## ----  plotClusters_diffDist
makeFigure("plotClusters_diffDist",height=3)
wh<-which(clusterParams3[,"k"]==12)
wh<-wh[match(c("eucDist","corDist","spearDist"),clusterParams3[wh,"distFunction"])]
par(mar=c(3,10,.1,.1))
plotClusters(ceObj3, whichClusters=clusterParams3[wh,"clusteringIndex"],colPalette = c(bigPalette, rainbow(100)),axisLine=-2,matchToTop=FALSE,
clusterLabels=c("Euclidean","Pearson Corr.","Spearman Rho"))
title(xlab="Samples",line=0,cex.lab=2)
title(ylab="Distance Function",cex.lab=2,line=6)
dev.off()

print(Sys.time())

