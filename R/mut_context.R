
#' Add context bases to mutations
#'
#' @param mutnet.obj A dataframe of mutations
#' @param ref.genome Name of the loaded reference genome ; if NULL the function will use the "context" column of the dataframe
#'
#' @return The dataframe with context bases added
#' 
#' @export
#'
#' @examples
mc.mutContext <- function(mutnet.obj, ref.genome=NULL){
  
  # extract context of each mutation
  if(!is.null(ref.genome)){
    context <-
      as.character(getSeq(
        get(ref.genome),
        Rle(mutnet.obj$chrom),
        mutnet.obj$pos.1 - 3,
        mutnet.obj$pos.1 + 3
      ))
  }else if(is.null(ref.genome) & "context" %in% names(mutnet.obj)){
    context <- mutnet.obj$context
  }else{
    stop("please provide an avalable ref.genome object")}
  
  # format
  new.context <- base.3p <- base.5p <- rep(NA, length(context))
  
  for (i in 1:length(context)) {
    sp <- strsplit(context[i], "")[[1]]
    sp[4] <- "x"
    new.context[i] <- paste(sp, collapse = '')
    base.5p[i] <- sp[3]
    base.3p[i] <- sp[5]
  }
  
  mutnet.obj$context <- new.context
  mutnet.obj$base.3p <- base.3p
  mutnet.obj$base.5p <- base.5p
  
  return(mutnet.obj)
}
