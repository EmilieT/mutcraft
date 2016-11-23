## ----eval=FALSE----------------------------------------------------------
#  library(mutcraft)

## ----eval=FALSE----------------------------------------------------------
#  library(BSgenome)
#  ref_genome = "BSgenome.Hsapiens.NCBI.GRCh37"
#  library(ref_genome, character.only = T)

## ----eval=FALSE----------------------------------------------------------
#  my.files <- c("myfile1.vcf","myfile2.vcf")
#  s.names <- c("sample1", "sample2")
#  mutnet <- mc.loadVcfs(my.files,s.names, ref.genome=ref_genome)

## ----eval=FALSE----------------------------------------------------------
#  mutspec <- lapply(mutnet,mc.mutSpectrum,"col.ref"="ref.allele","col.alt"="alt.allele")
#  mc.plotSpectrum(mutspec,"prop",print.num=F)

