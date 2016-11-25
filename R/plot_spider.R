
#' Plot spider spectrum of mutations
#'
#' @param mut.sum Table of the counts and proportions of each mutation type
#'
#' @import ggplot2 
#' @import ggthemes
#' @import scales
#' 
#' @return A ggplot object to plot
#' @export
#'
#' @examples
mc.plotSpider <- function(mut.sum){
  
  if(!is.data.frame(mut.sum)){
    for(i in 1:length(mut.sum)){
      mut.sum[[i]]$sample.name <- names(mut.sum)[i]
    } 
    mut.sum <- do.call("rbind",mut.sum)  
  }
  
  if(is.null(mut.sum$sample.name)){
    mut.sum$sample.name <- "sample"
  } 
  p <- ggplot(mut.sum,
           aes(
             x = mut,
             y = prop,
             colour = sample.name,
             group = sample.name
           )) + geom_line() + coord_polar(theta = "x", direction = -1) + 
    scale_y_continuous(labels = percent_format()) + theme_linedraw() + 
    geom_point(cex = 5, alpha = 0.5) + labs(title = "Spider Spectrum", x = "", y = "") + 
    scale_color_tableau(guide = guide_legend(title = "Samples", nrow =  8))
  return(p)
  
}