
mc.baseComp <- function(base = c("A","C","G","T")[1]){
  b0 <- rep(NA,length(base))
  
  b0[base %in% "A"] <- "T"
  b0[base %in% "C"] <- "G"
  b0[base %in% "G"] <- "C"
  b0[base %in% "T"] <- "A"
  
  return(b0)
}