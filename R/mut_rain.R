
#' Compute distance between mutations for rainfall plots
#'
#' @param mutnet.obj A dataframe of mutation info
#' @param col.chrom Name of the column with the chromosomes
#' @param col.pos Name of the column with the mutations positions
#' @param col.ref Name of the column with the reference allele
#' @param col.alt Name of the column with the alternative allele
#' @param symetrize Boolean to symetrize reversable mutations or not
#'
#' @return A table with relative distance for each mutation
#' @export
#'
#' @examples
mc.mutRain <- function(mutnet.obj, col.chrom="chrom", col.pos="pos", col.ref="ref.allele", col.alt="alt.allele", symetrize=c(TRUE,FALSE)[1]){
  
  ## genomicOrder
  m0 <- mc.sortByChromosom(mutnet.obj[,c(col.chrom,col.pos,col.ref,col.alt)], col.chrom, col.pos)
  
  ## distance
  m0$d <- c(NA,m0[-1,col.pos] - m0[-nrow(m0),col.pos]) 
  w <- unlist(lapply(unique(m0[,col.chrom]),function(chr) min(which(m0[,col.chrom] == chr))))
  m0$d[w] <- NA
  
  if(!is.null(col.ref) & !is.null(col.alt)){
    
    ## symetrize
    if(symetrize){
      ws <- which(m0[,col.ref] %in% c("A","G"))
      m0[ws,col.ref] <- mc.baseComp(m0[ws,col.ref]) 
      m0[ws,col.alt] <- mc.baseComp(m0[ws,col.alt]) 
    }
    
    m0$mut <- paste(m0[,col.ref],m0[,col.alt],sep=">")
    ##
    m0 <- m0[,c(col.chrom, col.pos, "d", col.ref, col.alt, "mut")]
    names(m0) <- c("chrom", "pos", "d", "ref", "alt", "mut")
  }else{
    m0 <- m0[,c(col.chrom, col.pos, "d")]
    names(m0) <- c("chrom", "pos", "d")
  }
  
  return(m0)
  
}