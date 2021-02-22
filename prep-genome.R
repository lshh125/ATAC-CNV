#load the package
library(CopyscAT)

#load your favourite genome
library(BSgenome.Hsapiens.UCSC.hg19)
#use it to save references - this will create a genome chrom.sizes file, a cytobands file and a CpG file
generateReferences(BSgenome.Hsapiens.UCSC.hg19,genomeText = "hg19",tileWidth = 1e6,outputDir = "data/")

