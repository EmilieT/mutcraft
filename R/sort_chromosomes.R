
mc.sortByChromosom <- function(mutnet.obj, col.chrom= "chrom", col.pos="pos"){
  
  mutnet.obj[,col.pos] <- as.numeric(mutnet.obj[,col.pos])
  mutnet.obj[, col.chrom] <- as.character(mutnet.obj[, col.chrom])
  
  chromosomes <- unique(mutnet.obj[,col.chrom])
  
  max.num <- max(as.numeric(chromosomes[grepl("\\d", chromosomes)]))
  
  
  s.1 <- which(mutnet.obj[,col.chrom] == "X")
  s.2 <- which(mutnet.obj[,col.chrom] == "Y")
  mutnet.obj[s.1, col.chrom] <- max.num+1
  mutnet.obj[s.2, col.chrom] <- max.num+2
  
  # other chromosomes to NA
  mutnet.obj[!grepl("\\d", mutnet.obj[, col.chrom]), col.chrom] <- NA
  
  mutnet.obj[, col.chrom] <- as.numeric(mutnet.obj[,col.chrom])
  
  o <- order(mutnet.obj[,col.chrom], mutnet.obj[,col.pos], na.last=TRUE)
  
  mut.sort <- mutnet.obj[o,]
  
  return(mut.sort)
}