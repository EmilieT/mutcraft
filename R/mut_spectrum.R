
#' Make mutations table
#'
#' @param mutnet.obj A dataframe of mutation info
#' @param col.ref Name of the column with the reference alleles
#' @param col.alt Name of the column with the alt alleles
#' @param nuc Array of used nucleotides
#' @param symetrize Boolean to symetrize reversable mutations or not
#'
#' @return A table with the counts and proportions of each mutation type
#' @export
#'
#' @examples
mc.mutSpectrum <- function(mutnet.obj, col.ref="ref.allele", col.alt="alt.allele", nuc=c("A","C","G","T"),  symetrize = T){
  
  mutnet.obj[,col.ref] <- factor(mutnet.obj[,col.ref],levels=nuc)
  mutnet.obj[,col.alt] <- factor(mutnet.obj[,col.alt],levels=nuc)
  
  
  mut.sum <- table(mutnet.obj[,c(col.ref,col.alt)])
  mut.sum <- as.data.frame(mut.sum)
  names(mut.sum) <- c("ref","alt","count")
  mut.sum$mut <- paste(mut.sum$ref,mut.sum$alt,sep=">")
  
  nm <- which(mut.sum$ref == mut.sum$alt)
  if(sum(mut.sum$count[nm])!=0){
    warning(paste("there are identical mutations :",mut.sum$mut[nm][mut.sum$count[nm] != 0]))
  } 
  
  mut.sum <- mut.sum[-nm,]
  
  if(symetrize){
    m0 <- mut.sum[mut.sum$ref %in% c("C","T"),]
    for(i in 1:nrow(m0)) m0$count[i] <- as.numeric(m0$count[i] + mut.sum$count[mc.baseComp(mut.sum$ref) == m0$ref[i] & mc.baseComp(mut.sum$alt) == m0$alt[i]])
  }
  
  m0$prop <- m0$count/sum(m0$count)
  mut.sum <- m0
  
  mut.sum <- mut.sum[,c("mut","ref","alt","count","prop")]
  mut.sum <- mut.sum[order(mut.sum$mut),]
  rownames(mut.sum) <- NULL
  
  return(mut.sum)
}
