#' A CNA permutation function for CNA
#'
#' @description This function allows you to do permutation of CNA segments for CNA.
#' @param cnvdatalist A list of CNA data with information about Chromosome, start and end position and segmentmean.
#' @param segmentmean_name The name of CNA segment mean column in the cnvdatalist
#' @param segment_cutoff Cut off value (positive or negative) of CNA segments in CNA data lists
#' @param permtime Permutation times
#' @param core Go parallel number of cores using
#' @export
#' @examples
#' @import doParallel
#' @import foreach
#' @import ggplot2
#' @import GenomicRanges


tcgacnv<-function(cnvdatalist,segmentmean_name,segment_cutoff,permtime,core){


  CNV <- lapply(cnvdatalist,function(x) abs(x[,segmentmean_name])>=segment_cutoff)
  listdata<-mapply(FUN=function(x,CNV){cbind(x,CNV)},cnvdatalist, CNV,SIMPLIFY = FALSE)
  ## Only consider CNV=1
  tcga.cnv <- lapply(listdata, function(x) subset(x,x$CNV==TRUE))
  tcga.cnv <- unname(tcga.cnv[sapply(tcga.cnv, function(x) nrow(x)>0)])
  tcga.cnv <- lapply(tcga.cnv, function(x) subset(x,!(x$Chromosome %in% c("chrX","chrY","chrx","chry"))))

  ## Find available chromosome names for each sample
  chr.name.p<-lapply(tcga.cnv, function(x) as.character(unique(x$Chromosome)))

  ## Calculate the range of permutation: range.p
  range.total <- mapply(FUN=function(x,chr){x[x$Chromosome %in% chr,]},listdata, chr.name.p,SIMPLIFY = FALSE)
  range.p <- lapply(range.total, function(x) by(x,x$Chromosome,function(x) c(min(x$Start),max(x$End))))
  ## Delete the NULL chromosome and size
  range.p <- mapply(FUN = function(x,name){x[name]},range.p, chr.name.p,
                  SIMPLIFY = FALSE)
  ## Calculate the number of observations for each chromosome and delete NA chromosome
  size.p <- lapply(tcga.cnv, function(x) by(x,x$Chromosome,function(x) length(x$Start)))
  size.p <- mapply(FUN = function(x,name){x[name]},size.p, chr.name.p,
                 SIMPLIFY = FALSE)

  len.p <- sapply(size.p, length)
  set.seed(100)
  cl <- makeCluster(core)
  registerDoParallel(cl)

  start.p <- foreach(i=1:length(tcga.cnv)) %:%
    foreach(j=1:len.p[i]) %:%
    foreach(k=1:unlist(size.p[[i]])[j]) %dopar% {
      cnv.p<-tcga.cnv[[i]][tcga.cnv[[i]]$Chromosome==chr.name.p[[i]][j],][k,]
      range<-unlist(range.p[[i]][j])
      start<-sample(range[1]:(range[2]-cnv.p$End+cnv.p$Start),permtime,replace=TRUE)
      end<-start+cnv.p$End-cnv.p$Start
      Chrname<-rep(as.character(chr.name.p[[i]][j]),permtime)
      data.frame(start,end,Chrname)


    }


  df <- foreach(i=1:permtime) %:%
    foreach(l=1:length(tcga.cnv)) %:%
    foreach(j=1:len.p[l],.combine='rbind') %:%
    foreach(k=1:unlist(size.p[l])[j],.combine='rbind') %dopar% {
      data.frame(start.p[[l]][[j]][[k]])[i,]
    }

  df<-lapply(df,function(x) lapply(x, function(x) GRanges(as.character(x$Chrname),ranges=IRanges(x$start,x$end))))

  on.exit(stopCluster(cl))

  return(list("df"=df, "tcga.cnv"=tcga.cnv))
}



