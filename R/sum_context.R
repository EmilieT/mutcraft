
#' Compute frequencies and proportions of mutations with context
#'
#' @param mutnet.obj A dataframe of mutation info 
#' @param base.ref.col Name of the column with the reference allele in dataframe
#' @param base.alt.col Name of the column with the alternative allele in dataframe
#' @param base.3p.col Name of the column with the base after the mutation in dataframe
#' @param base.5p.col Name of the column with the base before the mutation in dataframe
#' @param trinuc.pc Array of the proportion of each trinucleotide in the genome (for normalisation)
#' @param symetrize Boolean to symetrize reversable mutations or not
#'
#' @return A dataframe with each mutation and context, with frequencies and proportions
#'
#' @examples
mc.sumContext <- function(mutnet.obj, base.ref.col="ref.allele", base.alt.col="alt.allele", base.3p.col = "base.3p", base.5p.col="base.5p", trinuc.pc, symetrize=c(TRUE,FALSE)[1]){
  
  if(!is.data.frame(mutnet.obj))
  {
    for(i in 1:length(mutnet.obj)){
      mutnet.obj[[i]]$sample.name <- names(mutnet.obj)[i]
      mutnet <- do.call("rbind",mutnet)  
    } 
  }
  
  if(is.null(mutnet.obj$sample.name)){
    mutnet.obj$sample.name <- "sample"
  } 
  
  base.ref <- mutnet.obj[,base.ref.col]
  base.alt <- mutnet.obj[,base.alt.col]
  base.5p <- mutnet.obj[,base.5p.col]
  base.3p <- mutnet.obj[,base.3p.col]
  sample.name <- mutnet.obj$sample.name
  
  
  ### symetrization
  
  if(symetrize){
    wS <- which(base.ref %in% c("A","G"))
    base.ref[wS] <- mc.baseComp(base.ref[wS])
    base.alt[wS] <- mc.baseComp(base.alt[wS])
    base.5p[wS] <- mc.baseComp(base.5p[wS])
    base.3p[wS] <- mc.baseComp(base.3p[wS])
  }
  
  mut <- paste(base.ref,base.alt,sep=">")
  context <- paste(base.5p,base.3p,sep="x")
  
  mutab <- as.data.frame(table(as.data.frame(cbind(mut, context, base.5p, base.3p, sample.name))))

  #####
  

  mutab$Freq <- as.numeric(mutab$Freq)
  for(i in 1:ncol(mutab)) if(is.factor(mutab[,i])) mutab[,i] <- as.character(mutab[,i])
  
  mutab$prop <- mutab$Freq/sum(mutab$Freq)
  
  mutab$freqN <- NA
  mutab$trin <- paste(mutab$base.5p,substr(mutab$mut,1,1),mutab$base.3p,sep="")
  
  for(mu in unique(mutab$mut)){
    w <- which(mutab$mut == mu)
    mutab$freqN[w] <- mutab$Freq[w]/(trinucpc[mutab$trin[w]]*2)
  }
  
  mutab$propN <- mutab$freqN/sum(mutab$freqN)
  
  mutab$log10.prop <- log10(mutab$propN*100)
  mutab$log10.prop[mutab$log10.prop %in% -Inf] <- 0
    
  
  return(mutab)
}


