#' A function that help you make get the overlap of SE and permuted CNA excluding chromosome X and Y 
#'
#' This function allows you to make get the overlap of SE and permuted CNA
#' @param data1 Super Enhancer Sample One information about chromosome name, start and end position with columns "REGION_ID","CHROM","START","STOP","NUM_LOCI","CONSTITUENT_SIZE","Intensity","Intensity2","EnhancerRank1","Overlap_Genes","Proximal_Genes","Closest_Gene","EnhancerRank2","IsSuper" 
#' @param data2 Super Enhancer Sample Two information about chromosome name, start and end position with columns "REGION_ID","CHROM","START","STOP","NUM_LOCI","CONSTITUENT_SIZE","Intensity","Intensity2","EnhancerRank1","Overlap_Genes","Proximal_Genes","Closest_Gene","EnhancerRank2","IsSuper" 
#' @param cnadata A list of permutation CNA results
#' @keywords cnVpermtest
#' @export
#' @examples
#' @import GenomicRanges
#' permutation_sum()


permutation_sum <- function(data1,data2,cnadata){
  
  col <- c("REGION_ID","CHROM","START","STOP",
         "NUM_LOCI","CONSTITUENT_SIZE","Intensity","Intensity2",
         "EnhancerRank1","Overlap_Genes","Proximal_Genes",
         "Closest_Gene","EnhancerRank2","IsSuper")
  colnames(data1) <- col
  colnames(data2) <- col
  data1 <- subset(data1,!(data1$CHROM %in% c("chrX","chrY","chrx","chry")))
  data2 <- subset(data2,!(data2$CHROM %in% c("chrX","chrY","chrx","chry")))
  
  chrname.1s <- as.character(data1$CHROM)
  start.1s <- as.numeric(data1$START)
  end.1s <- as.numeric(data1$STOP)
  gr.1s <- GRanges(seqnames=chrname.1s,ranges=IRanges(start=start.1s,end=end.1s))
  
  chrname.2s <- as.character(data2$CHROM)
  start.2s <- as.numeric(data2$START)
  end.2s <- as.numeric(data2$STOP)
  gr.2s <- GRanges(seqnames=chrname.2s,ranges=IRanges(start=start.2s,end=end.2s))
  
  hits <- findOverlaps(gr.1s, gr.2s)
  overlap12 <- intersect(gr.1s, gr.2s)
  overlap <- lapply(cnadata, function(x) lapply(x, function(x) intersect(overlap12, x)))
  num.base <- lapply(overlap, function(x) lapply(x, function(x) sum(width(x))))
  n1 <- sapply(num.base, function(x) sum(unlist(x)))
  n1
  # return(list("n1"=n1, "overlap12"=overlap12))
  
  
}