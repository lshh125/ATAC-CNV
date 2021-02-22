library(CopyscAT)

output_folder = "results/"
dataset_name = "GSM4138898_scATAC_MPAL1_T1"
  
##### REGULAR WORKFLOW #####
#initialize the environment 
initialiseEnvironment(genomeFile="data/hg19_chrom_sizes.tsv",
                      cytobandFile="data/hg19_1e+06_cytoband_densities_granges.tsv",
                      cpgFile="data/hg19_1e+06_cpg_densities.tsv",
                      binSize=1e6,
                      minFrags=1e4,
                      cellSuffix=c("-1"),
                      lowerTrim=0.5,
                      upperTrim=0.8)

#for this tutorial we will use the sample data included in scData
#to create your own use the process_fragment_file.py script included in the package and run it on a fragments.tsv.gz file of your choosing

setOutputFile(output_folder, dataset_name)
#step 1 normalize the matrix
#USING SAMPLE DATA FROM PACKAGE
#option: if using your own file replace below with the following
scData<-readInputTable("data/GSM4138898_scATAC_MPAL1_T1.tsv")
#scData<-scDataSamp

# Write "GSM4138898_scATAC_MPAL1_T1_signal_distribution.pdf"
scData_k_norm <- normalizeMatrixN(scData,logNorm = FALSE,maxZero=2000,imputeZeros = FALSE,blacklistProp = 0.8,blacklistCutoff=125,dividingFactor=1,upperFilterQuantile = 0.95)

#collapse into chromosome arm level
summaryFunction<-cutAverage
#powval 0.7 - 0.75
scData_collapse<-collapseChrom3N(scData_k_norm,summaryFunction=summaryFunction,binExpand = 1,minimumChromValue = 100,logTrans = FALSE,tssEnrich = 1,logBase=2,minCPG=300,powVal=0.73) 
#apply additional filters
scData_collapse<-filterCells(scData_collapse,minimumSegments = 40,minDensity = 0.1)

#show unscaled chromosome list
graphCNVDistribution(scData_collapse,outputSuffix = "test_violinsn2")
write.csv(scData_collapse, paste0(output_folder, dataset_name, "_intermediate0.csv"), row.names = FALSE)

#compute centers
median_iqr <- computeCenters(scData_collapse,summaryFunction=summaryFunction)
#identify chromosome-level amplifications
candidate_cnvs<-identifyCNVClusters(scData_collapse,median_iqr,useDummyCells = TRUE,propDummy=0.25,minMix=0.01,deltaMean = 0.03,deltaBIC2 = 0.25,bicMinimum = 0.1, subsetSize=600,fakeCellSD = 0.08, uncertaintyCutoff = 0.55,summaryFunction=summaryFunction,maxClust = 4,mergeCutoff = 3,IQRCutoff= 0.2,medianQuantileCutoff = 0.4)
write.csv(candidate_cnvs[[1]], paste0(output_folder, dataset_name, "_intermeiate1.csv"), row.names = FALSE)

#cleanup step
candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.5) #= 1.5)
write.csv(candidate_cnvs_clean[[1]], paste0(output_folder, dataset_name, "_intermeiate2.csv"), row.names = FALSE)

#final results and annotation
final_cnv_list<-annotateCNV4(candidate_cnvs_clean, saveOutput=TRUE,outputSuffix = "clean_cnv",sdCNV = 0.5,filterResults=TRUE,filterRange=0.8)

#identify double minutes / amplifications
#note: this is slow, and may take ~5 minutes
#option to compile this code
library(compiler)
dmRead<-cmpfun(identifyDoubleMinutes)

dm_candidates<-dmRead(scData_k_norm,minCells=100,qualityCutoff2 = 100,minThreshold = 4)
write.table(x=dm_candidates,file=str_c(scCNVCaller$locPrefix,"samp_dm.csv"),quote=FALSE,row.names = FALSE,sep=",")

#assess putative LOH regions
loh_regions<-getLOHRegions(scData_k_norm,diffThreshold = 3,lossCutoff = -0.75,minLength = 2e6,minSeg=2,targetFun=IQR,lossCutoffCells = 200,quantileLimit=0.2,cpgCutoff=100,dummyQuantile=0.6,dummyPercentile=0.4,dummySd=0.1)
write.table(x=loh_regions[[1]],file=str_c(scCNVCaller$locPrefix,"samp_loss.csv"),quote=FALSE,row.names = FALSE,sep=",")