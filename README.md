# ATAC-CNV

## Prepare data

Download [data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4138898) from this Nature Biotech [paper](https://www.nature.com/articles/s41587-019-0332-7) and put into `data/`

## Configure Environment

1. Download [R](https://www.r-project.org/)
2. Install [CopyscAT](https://github.com/spcdot/CopyscAT) as directed in the [tutorial](https://github.com/spcdot/CopyscAT/blob/master/copyscat_tutorial.R).

[Rstudio](https://rstudio.com/) is one of the most widely used IDE. The community version is open source. Jupyter also provides R kernel, but I do not have personal expeirence with it.

## Use CopyscAT to process data
1. Run `prep-genome.R`. You may need to install [BSgenome.Hsapiens.UCSC.hg19](http://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html)

2. Run `process_fragment_file.cmd`. You will need Python.

3. Run `obtain-cnv.R`.
