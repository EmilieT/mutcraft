## ---- message=FALSE, warning=FALSE---------------------------------------
library(mutcraft)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(BSgenome)
ref_genome = "BSgenome.Hsapiens.NCBI.GRCh37"
library(ref_genome, character.only = T)

## ---- message=FALSE, warning=FALSE---------------------------------------
my.files <- list.files(system.file("extdata", package="mutcraft"),
                          pattern = ".vcf", full.names = TRUE)
s.names <- c("sample1", "sample2")
mutnet <- mc.loadVcfs(my.files, s.names, ref.genome=ref_genome)

## ---- message=FALSE, warning=FALSE, fig.width=7,fig.height=5-------------
mutspec <- lapply(mutnet,mc.mutSpectrum,"ref.allele","alt.allele")
mc.plotSpectrum(mutspec,"prop",print.num=F)

## ---- message=FALSE, warning=FALSE---------------------------------------
mut.c <- lapply(mutnet, mc.mutContext, ref_genome)

## ---- message=FALSE, warning=FALSE, fig.width=7,fig.height=5-------------
mc.plotContext(mut.c)

