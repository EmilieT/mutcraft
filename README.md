# Getting started
Mutcraft is a package of funtions to read somatic mutation information from VCF files and visualize spectra, rainplots, strand bias and mutational signatures.

## Installation
```{r , message=FALSE, warning=FALSE}
library(mutcraft)
```



## Reference genome
In order to extract mutation contexts and mutation strands you need to load a reference genome

```{r, message=FALSE, warning=FALSE}
library(BSgenome)
ref_genome = "BSgenome.Hsapiens.NCBI.GRCh37"
library(ref_genome, character.only = T)
```



## Load files
```{r, message=FALSE, warning=FALSE}
my.files <- list.files(system.file("extdata", package="mutcraft"),
                          pattern = ".vcf", full.names = TRUE)
s.names <- c("sample1", "sample2")
mutnet <- mc.loadVcfs(my.files, s.names, ref.genome=ref_genome)
```

## Mutation Spectrum
```{r, message=FALSE, warning=FALSE, fig.width=7,fig.height=5}
mutspec <- lapply(mutnet,mc.mutSpectrum,"ref.allele","alt.allele")
mc.plotSpectrum(mutspec,"prop",print.num=F)
```

## Add mutations context
```{r, message=FALSE, warning=FALSE}
mut.c <- lapply(mutnet, mc.mutContext, ref_genome)
```

## Plot context histogram
```{r, message=FALSE, warning=FALSE, fig.width=7,fig.height=5}
mc.plotContext(mut.c)
```