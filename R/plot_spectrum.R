


#' Plot mutation spectrum
#'
#' @param mut.sum Table of the counts and proportions of each mutation type
#' @param yaxis Variable to plot: counts or proportions
#' @param print.num Boolean to print values on the plot or not
#'
#' @import ggplot2 
#' @import ggthemes
#' @import scales
#' 
#' @return A ggplot object to plot
#' @export
#'
#' @examples
mc.plotSpectrum <-function(mut.sum, yaxis = c("count", "prop")[1],
           print.num = c(TRUE, FALSE)[1]) {
  
    if (!is.data.frame(mut.sum)) {
      for (i in 1:length(mut.sum))
        mut.sum[[i]]$sample.name <- names(mut.sum)[i]
      mut.sum <- do.call("rbind", mut.sum)
    }
    
    
    if (yaxis == "count") {
      p <-ggplot(mut.sum, aes(
          x = factor(mut),
          y = count,
          fill = mut
        )) + geom_bar(stat = "identity")  + theme_linedraw() + labs(title="Mutation Spectrum (count)", x="Mutation Type", y="Count") + scale_fill_tableau(guide = guide_legend(title = "Mutation Type", nrow=3))
      if (print.num)
        p <- p + geom_text(label = paste(round(mut.sum$count / 1000, 1), "k", sep =""))
    } else if (yaxis == "prop") {
      p <- ggplot(mut.sum, aes(
          x = factor(mut),
          y = prop,
          fill = mut
        )) + geom_bar(stat = "identity") + scale_y_continuous(labels = percent_format()) + theme_linedraw() + labs(title="Mutation Spectrum (%)", x="Mutation Type", y="Percentage") + scale_fill_tableau(guide = guide_legend(title = "Mutation Type", nrow=3))
      if (print.num)
        p <- p + geom_text(label = paste(round(mut.sum$prop * 100, 1), "%", sep = ""))
    }
    
    if (!is.null(mut.sum$sample.name)) {
      p <- p + facet_wrap(~ sample.name)
    }
    
    return(p)
  }