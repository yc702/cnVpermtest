#' A function to make histogram for CNA permutation overlap with Super Enhancers
#'
#' This function allows you to make histogram for CNA permutation.
#' @param permoverlap A dataframe with first column as vector bases of permutated CNA overlaps with Super Enhancers.
#' @param obsoverlap The observed bases of CNA segment overlap with Super Enhancers 
#' @param title Title of the histogram
#' @param position Position of the lable of observed overlap value, right or left of the marked dashed line
#' @param num Number of bases to be away from the marked dashed line
#' @keywords cnVpermtest
#' @export
#' @examples
#' @import ggplot2
#' permhist()

permhist<-function(permoverlap,obsoverlap,title,position="right",num){
  
  if (position =="right"){
    ggplot(data.frame(permoverlap[,1]),aes(permoverlap[,1]))+geom_histogram()+geom_vline(xintercept=obsoverlap,colour="red",linetype="dashed")+
      geom_text(y=20,x=obsoverlap-num,cex=6,label=paste("observed value is",obsoverlap),colour="red")+
      labs(x="Permutation Overlap Bases")+ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5,size=22,face="bold"))+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"))
  }
  else if (position=="left"){
    
    ggplot(data.frame(permoverlap[,1]),aes(permoverlap[,1]))+geom_histogram()+geom_vline(xintercept=obsoverlap,colour="red",linetype="dashed")+
      geom_text(y=20,x=obsoverlap-num,cex=6,label=paste("observed value is",obsoverlap),colour="red")+
      labs(x="Permutation Overlap Bases")+ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5,size=22,face="bold"))+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"))
    
  }
  
}


