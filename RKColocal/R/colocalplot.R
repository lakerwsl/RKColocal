#' @title Colocalization Analysis Plot
#' 
#' @description \code{colocalplot} is used to make colocalization analysis plot, including scatter plot and Li's plot. 
#' Scatter plot is a plot where x axis is pixel intensity of green channel and y axis is pixel intensity of red channel
#' Li's plot provide two plots, the common x axis is the product value and y axis of two plots are pixel intensity of red and green channel, respectively.
#' 
#' @details This function provides the colocalization analysis plot via plot function in \code{ggplot2}. 
#' 
#' @seealso \code{\link{DualImageVisual}}
#' @import grid
#' @import ggplot2
#' @import graphics
#' 
#' @param X A numerical matrix, the intensity matrix of the first channel/red channel. 
#' @param Y A numerical matrix, the intensity matrix of the first channel/green channel.
#' @param mask An logical matrix, where TRUE means pixel is involved into analysis and vice versa. Dimension of mask must be matched with one of image.
#' @param method A string, indicating which plot is made. The plot can only be from \code{c('scatter','li')}.
#' 
#' @return Invisible \code{NULL}
#' 
#' @references Li, Qi, et al. "A syntaxin 1, Gαo, and N-type calcium channel complex at a presynaptic nerve terminal: analysis by quantitative immunocolocalization." Journal of Neuroscience 24.16 (2004): 4070-4081.
#' 
#' @export
#' 
#' @examples 
#' X<-matrix(runif(10*10),nrow=10,ncol=10)
#' Y<-matrix(runif(10*10),nrow=10,ncol=10)
#' colocalplot(X,Y)
#' 
#' @author Shulei Wang
colocalplot <- function(X, Y, mask = NULL, method = 'scatter')
{
  if (!is.matrix(X) ||!is.numeric(X))
    stop("X must be numeric matrix")
  if (!is.matrix(Y) ||!is.numeric(Y))
    stop("Y must be numeric matrix")
  if (nrow(X) != nrow(Y) || ncol(X) != ncol(Y))
    stop("size of X and Y doesn't match")
  if (!is.null(mask))
  {
    if (!is.matrix(mask) || !is.logical(mask))
      stop("mask must be logical matrix")
    if (nrow(X) != nrow(mask) || ncol(X) != ncol(mask))
      stop("size of data and mask doesn't match")
  }
  if (max(X)>1 || min(X)<0)
    stop("the range of X must be [0,1]")
  if (max(Y)>1 || min(Y)<0)
    stop("the range of Y must be [0,1]")
  
  if (is.null(mask))
  {
    AX=as.vector(X)
    AY=as.vector(Y)
  } else {
    AX=X[mask]
    AY=Y[mask]
    AX=as.vector(AX)
    AY=as.vector(AY)
  }
  
  if (method == 'scatter')
  {
    qplot(AX,AY,main="Scatter Plot in ROI",xlab="Intensity X",ylab = "Intensity Y",margins=TRUE) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_point(alpha=0.05)
  } else if (method == 'li') {
    meanX <- mean(AX)
    meanY <- mean(AY)
    prod <- (AX - meanX) * (AY - meanY)
    prodmax <- max(max(prod),-min(prod))
    plot1 <- qplot(x=prod,y=AX,main="Li's Plot in ROI",xlab="(X-x)(Y-y)",ylab = "Intensity X",margins=TRUE) +
      geom_vline(xintercept = 0) +
      xlim(-prodmax, prodmax) +
      geom_point(alpha=0.05)
    plot2 <- qplot(x=prod,y=AY,main="Li's Plot in ROI",xlab="(X-x)(Y-y)",ylab = "Intensity Y",margins=TRUE) +
      geom_vline(xintercept = 0) +
      xlim(-prodmax, prodmax) +
      geom_point(alpha=0.05)
    multiplot(plot1, plot2, cols=2)
  } else {
    stop("method should be one of \"scatter\", \"li\" ")
  }
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}