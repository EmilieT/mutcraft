
#' Plot histogram of the 96 possible mutations
#'
#' @param mutnet A list of dataframe with mutations infos for each sample
#' @param base.ref.col Name of the column with the reference allele in dataframe
#' @param base.alt.col Name of the column with the alternative allele in dataframe
#' @param base.3p.col Name of the column with the base after the mutation in dataframe
#' @param base.5p.col Name of the column with the base before the mutation in dataframe
#' @param trinuc.pc Array of the proportion of each trinucleotide in the genome (for normalisation)
#' @param symetrize Boolean to symetrize reversable mutations or not
#'
#' @return A ggplot object to plot
#' 
#' @import ggplot2
#' @import ggthemes
#' 
#' @export
#'
#' @examples
mc.plotContext <- function(mutnet, base.ref.col="ref.allele", base.alt.col="alt.allele", base.3p.col = "base.3p", base.5p.col="base.5p", trinuc.pc, symetrize=c(TRUE,FALSE)[1]){
  
  mutab <- lapply(mutnet,mc.sumContext, base.ref.col, base.alt.col, base.3p.col, base.5p.col, trinuc.pc, symetrize)
  
  if(!is.data.frame(mutab))
  {
    for(i in 1:length(mutab)) mutab[[i]]$sampleName <- names(mutab)[i]
    mutab <- do.call("rbind",mutab)  
  }
  
  p <- ggplot(mutab,aes(x=factor(context),y=propN,fill=mut)) + geom_bar(stat="identity") + facet_grid(sample.name ~ mut) + theme_linedraw() + scale_fill_tableau() + theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1)) + labs(x="Context",y="Proportion (%)", title="Mutation spectra with nucleotide context")
  return(p)
}