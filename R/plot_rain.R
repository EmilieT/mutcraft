
#' Plot mutation rainfalls
#'
#' @param mutrain Table of distance between successive mutations in genomic order
#' @param chrom.to.plot List of the chromosomes to include
#' @param color.col Column of the table to use for colors
#'
#' @return A ggplot2 object to plot
#' 
#' @import ggplot2 
#' @import ggthemes
#' 
#' @export
#'
#' @examples
mc.plotRain <- function(mutrain, chrom.to.plot = NULL, color.col="mut"){
  
  if(!is.data.frame(mutrain)){
    for(i in 1:length(mutrain)) {
      mutrain[[i]]$sample.name <- names(mutrain)[i]
    }
    mutrain <- do.call("rbind",mutrain)  
  }
  
  if(!is.null(chrom.to.plot)){
    mutrain <- mutrain[mutrain$chrom %in% chrom.to.plot,]
  }else{
    chrom.to.plot <- unique(mutrain$chrom)
  }
  
  ### couleur
  if(is.null(color.col) | (!color.col %in% names(mutrain))) mutrain$mut <- 1
  if(color.col != "mut") names(mutrain)[names(mutrain) == color.col] <- "mut"
  ###
  
  p <- ggplot(mutrain,aes(x=pos,y=d,col=mut)) + theme_linedraw() + geom_point() + scale_y_log10() + scale_color_tableau(guide = guide_legend(title = "Mutation Type", nrow =  3)) + labs(title="Rainfall Plot", x="Position", y="Distance") + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
  
  if(length(chrom.to.plot > 1)){
    if(length(table(mutrain$sample.name)) <= 1){
      p <- p + facet_grid(. ~ chrom, scales = "free", space="free_x")
    }else{
      p <- p + facet_grid(sample.name ~ chrom , scales = "free", space="free_x")
    }
  }else{
    if(length(table(mutrain$sample.name)) > 1)
    {
      p <- p + facet_grid(sampleName ~ . , scales = "free", space="free_x")
    }
  }
  return(p)
}