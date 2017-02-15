

#' Read VCF files as a list of data frames
#'
#' @param vcf.files A character vector of vcf file names
#' @param sample.names A character vector of sample names
#' @param ref.genome Name of the loaded reference genome
#' @param ... mclapply optional parameters
#'
#' @return A list of dataframes, one for each vcf file
#' @export
#'
#' @import VariantAnnotation
#'  
#' @examples
mc.loadVcfs <- function(vcf.files, sample.names, ref.genome,...) {
  if (length(vcf.files) != length(sample.names)) {
    stop("Number of sample names different than number of files")
  }
  
  if (is.null(ref.genome)) {
    stop("Please provide a reference genome")
  }
  
  print("Loading VCFs files ...")
  
  mut <- mclapply(vcf.files, function(f)
    mc.loadVcf(f, ref.genome),...)
  names(mut) <- sample.names
  for (i in 1:length(mut))
    mut[[i]]$sample.name <- sample.names[i]
  
  return(mut)

}


#' Read one VCF file and create a dataframe
#'
#' @param vcf.file Name of the VCF file
#' @param ref.genome Name of the reference genome to use
#'
#' @return A dataframe with the VCF infos
#' @export
#'
#' @import VariantAnnotation
#' @import SummarizedExperiment
#'
#' @examples
mc.loadVcf <- function(vcf.file, ref.genome) {
  
  if (is.null(vcf.file)) {
    stop("Please provide a valid VCF file")
  }
  
  if (is.null(ref.genome)) {
    stop("Please provide a reference genome")
  }
  
  # load vcf
  vcf <- readVcf(vcf.file, genome = ref.genome)
  
  # restrict to SSM
  vcf.m <- vcf[isSNV(vcf)]
  
  
  ranges <- rowRanges(vcf.m, fixed = T)
  
  # convert to dataframe for usage
  gr <- data.frame(ranges, stringsAsFactors = FALSE)
  
  # need converting for alt allele
  gr$ALT <- as.character(CharacterList(gr$ALT))
  
  # extract genotype in normal
  genotypes <- geno(vcf.m)$GT
  gr$geno.norm <- genotypes[, 2]
  
  # extract allelic counts in normal and tumor
  allelic.counts <- geno(vcf.m)$AD 
  gr$ref.tum.count <- gr$alt.tum.count <- gr$ref.norm.count <- gr$alt.norml.count <- rep(0, nrow(gr))
  
  if(!is.null(allelic.counts)){
    AD.tum <- do.call("rbind", allelic.counts[, 1])
    AD.norm <- do.call("rbind", allelic.counts[, 2])
    gr$ref.tum.count <- AD.tum[, 1]
    gr$alt.tum.count <- AD.tum[, 2]
    gr$ref.norm.count <- AD.norm[, 1]
    gr$alt.norml.count <- AD.norm[, 2]
  }else{
    warning("missing allelic counts info")
  }
  
  mutnet.obj <-
    data.frame(
      chrom = gr$seqnames,
      pos.1 = gr$start,
      pos = gr$end,
      ref.allele = gr$REF,
      alt.allele = gr$ALT,
      ref.count.tumor = gr$ref.tum.count,
      alt.count.tumor = gr$alt.tum.count,
      ref.count.normal = gr$ref.norm.count,
      alt.count.normal = gr$alt.norml.count,
      normal.gt = gr$geno.norm,
      stringsAsFactors = FALSE
    )
  
  mutnet.obj <- mc.formatGenotype(mutnet.obj)
  
  
  return(mutnet.obj)
  
}


mc.formatGenotype <- function(mutnet.obj) {
  new.gt <- rep(NA, length(mutnet.obj$normal.gt))
  
  for (g in 1:length(mutnet.obj$normal.gt)) {
    if (mutnet.obj$normal.gt[g] == "0/0") {
      new.gt[g] <-
        paste(mutnet.obj$ref.allele[g], mutnet.obj$ref.allele[g], sep = "/")
      
    } else if (mutnet.obj$normal.gt[g] == "0/1") {
      new.gt[g] <-
        paste(mutnet.obj$ref.allele[g], mutnet.obj$alt.allele[g], sep = "/")
      
    } else if (mutnet.obj$normal.gt[g] == "1/0") {
      new.gt[g] <-
        paste(mutnet.obj$alt.allele[g], mutnet.obj$ref.allele[g], sep = "/")
      
    } else{
      new.gt[g] <-
        paste(mutnet.obj$alt.allele[g], mutnet.obj$alt.allele[g], sep = "/")
      
    }
  }
  
  mutnet.obj$normal.gt <- new.gt
  
  return(mutnet.obj)
}
