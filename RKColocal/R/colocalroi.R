#' @title Region based Colocalization Analysis
#'
#' @description \code{colocalroi} can evaluate a statistical quantity for colocalization/correlation in a given region. The colocalization measures \code{colocalroi} support include 
#' Pearson correlation coefficient, Manders' split colocalization coefficients, Spearman’s rank correlation coefficient, Kendall tau correlation coefficient, intensity correlation quotient and maximum truncated Kendall tau correlation coefficient, a robust colocalization measure.
#' 
#' @details The Manders' split colocalization coefficients and maximum truncated Kendall tau correlation coefficient need the input of thresholds. If it is missing, the default thresholds are caculated by Otsu's method.
#' 
#' @seealso \code{\link{colocalroitest}}
#' @importFrom EBImage otsu
#' @importFrom pcaPP cor.fk
#' @importFrom stats cor
#' 
#' @param X A numerical matrix, the intensity matrix of the first channel/red channel. 
#' @param Y A numerical matrix, the intensity matrix of the first channel/green channel.
#' @param mask An logical matrix, where TRUE means pixel is involved into analysis and vice versa. Dimension of mask must be matched with one of image.
#' @param Thresholds, A vector of length 2, the input threshold for each channel when Manders' split colocalization coefficients and maximum truncated Kendall tau correlation coefficient are evaluated.
#' @param method A string, indicating which plot is made. The plot can only be from \code{c('pearson', 'kendall', 'spearman', 'icq', 'mandersM1', 'mandersM2', 'trunkendall')}.
#' 
#' @return The value of colocalization measure.
#' 
#' @references Manders, E. M., et al. "Dynamics of three-dimensional replication patterns during the S-phase, analysed by double labelling of DNA and confocal microscopy." Journal of cell science 103.3 (1992): 857-862.
#' @references Manders, E. M. M., F. J. Verbeek, and J. A. Aten. "Measurement of co‐localization of objects in dual‐colour confocal images." Journal of microscopy 169.3 (1993): 375-382.
#' @references Adler, J., S. N. Pagakis, and I. Parmryd. "Replicate‐based noise corrected correlation for accurate measurements of colocalization." Journal of microscopy 230.1 (2008): 121-133.
#' @references Li, Qi, et al. "A syntaxin 1, Gαo, and N-type calcium channel complex at a presynaptic nerve terminal: analysis by quantitative immunocolocalization." Journal of Neuroscience 24.16 (2004): 4070-4081.
#' @references French, Andrew P., et al. "Colocalization of fluorescent markers in confocal microscope images of plant cells." Nature protocols 3.4 (2008): 619.
#' 
#' @export
#' 
#' @examples
#' X=matrix(runif(32*32),32,32)
#' Y=matrix(runif(32*32),32,32)
#' colocalroi(X,Y)
#' 
#' @author Shulei Wang
colocalroi <- function(X, Y, mask = NULL, Thresholds=NULL, method = 'pearson')
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
  
  if (method == 'mandersM1'||method == 'mandersM2'||method == 'trunkendall')
  {
    if (!is.null(Thresholds))
    {
      if (!is.vector(Thresholds) ||!is.numeric(Thresholds))
        stop("Thresholds must be numeric vector")
      if (length(Thresholds)<2)
        stop("At least two elements in Thresholds")
    } else {
      if (nrow(X)<=1)
      {
        stop("At least two rows for the input image")
      }
      else 
      {
        Thresholds = c(0,0)
        Thresholds[1]<-otsu(X)
        Thresholds[2]<-otsu(Y)
      }
    }
  }
  
  colocalindex = colocalvec(AX, AY, Thresholds, method)
  
  return(colocalindex)
}

colocalvec <- function(x, y, Thresholds, method = 'pearson')
{
  if (!is.vector(x) || !is.numeric(x))
    stop("x must be numeric vector")
  if (!is.vector(y) || !is.numeric(y))
    stop("y must be numeric vector")
  if (length(x) != length(y))
    stop("length of x and y doesn't match")
  if (max(x)>1 || min(x)<0)
    stop("the range of x must be [0,1]")
  if (max(y)>1 || min(y)<0)
    stop("the range of y must be [0,1]")
  if (method == 'pearson')
  {
    colocalindex = cor(x,y,method="pearson")
  }
  else if (method == 'spearman')
  {
    colocalindex = cor(x,y,method="spearman")
  }
  else if (method == 'kendall')
  {
    colocalindex = cor.fk(x,y)
  }
  else if (method == 'icq')
  {
    meanX = mean(x)
    meanY = mean(y)
    prod = (x - meanX) * (y - meanY)
    ratio = sum(prod > 0) / length(x)
    colocalindex = ratio - 0.5
  }
  else if (method == 'mandersM1')
  {
    ThresholdX <- Thresholds[1]
    ThresholdY <- Thresholds[2]
    
    IndexX=x>ThresholdX
    IndexY=y>ThresholdY

    Index=IndexX&IndexY
    colocalindex=sum(Index)/sum(IndexX)
  } else if (method == 'mandersM2')
  {
    ThresholdX <- Thresholds[1]
    ThresholdY <- Thresholds[2]
    
    IndexX=x>ThresholdX
    IndexY=y>ThresholdY
    
    Index=IndexX&IndexY
    colocalindex=sum(Index)/sum(IndexY)
  } else if (method == 'trunkendall') {
    ThresholdX <- Thresholds[1]
    ThresholdY <- Thresholds[2]
    
    colocalindex = colocalrobustvec(x,y,Thresholds)
  }
  else {
    stop("method should be one of \"pearson\", \"kendall\", \"spearman\", \"icq\", \"mandersM1\", \"mandersM2\", \"trunkendall\" ")
  }
  return(colocalindex)
}

colocalrobustvec <- function(x, y, Thresholds)
{
  CutRX<-sum(x<=Thresholds[1])
  CutRY<-sum(y<=Thresholds[2])
  
  RX=rank(x, na.last = TRUE,ties.method = "random")
  RY=rank(y, na.last = TRUE,ties.method = "random")
  n=length(RX)
  Lower=2
  if(n-CutRX>n/Lower)
    CutRX=n-n/Lower
  if(n-CutRY>n/Lower)
    CutRY=n-n/Lower
  IndexX=RX>CutRX
  IndexY=RY>CutRY
  Index=IndexX&IndexY
  RX=RX[Index]
  RY=RY[Index]
  MaxTau=list(0,0)
  names(MaxTau)<-c("Value","Tau")
  MaxTau$Value=0
  alpha=1+1/log(log(n));
  #   iterationX=floor(log(n-CutX)/log(alpha));
  #   iterationY=floor(log(n-CutY)/log(alpha));
  prodX=1;
  while (prodX*alpha+CutRX<n)
  {
    prodX=prodX*alpha
    IndexX=RX>(n-prodX)
    prodY=1;
    while (prodY*alpha+CutRY<n)
    {
      prodY=prodY*alpha
      IndexY=RY>(n-prodY)
      Index=IndexX&IndexY
      l=sum(Index)
      if (l>1)
      {
        RXS=RX[Index]
        RYS=RY[Index]
        # R<-cor(RXS,RYS, method="kendall");
        R<-cor.fk(RXS,RYS)
        Var=sqrt(2*(2*l+5)/9/l/(l-1))
        Re=R/Var
        if(Re>MaxTau$Value)
        {
          MaxTau$Value=Re
          MaxTau$Tau=R
        }
      }
    }
  }
  
  Re = MaxTau$Value
  return (Re)
}

#' @title Statistical Test for Region based Colocalization Analysis
#'
#' @description \code{colocalroitest} can conduct a statistical test for colocalization/correlation in a given region. 
#' \code{colocalroitest} supports both pixel by pixel and block by block permutation test.
#' The colocalization measures \code{colocalroitest} support include 
#' Pearson correlation coefficient, Manders' split colocalization coefficients, Spearman’s rank correlation coefficient, Kendall tau correlation coefficient, intensity correlation quotient and maximum truncated Kendall tau correlation coefficient, a robust colocalization measure.
#' 
#' @details When \code{bsize} is 1 or mask is given, the pixel by pixel permutation test is conducted. 
#' When \code{bsize} is larger than 1 and mask is \code{NULL}, the block by blcok permutation test is conducted. 
#' 
#' @seealso \code{\link{colocalroi}}
#' @seealso \code{\link{plot.colocalroitest}}
#' 
#' @param X A numerical matrix, the intensity matrix of the first channel/red channel. 
#' @param Y A numerical matrix, the intensity matrix of the first channel/green channel.
#' @param mask A logical matrix, where TRUE means pixel is involved into analysis and vice versa. Dimension of mask must be matched with one of image.
#' @param Thresholds, A vector of length 2, the input threshold for each channel when Manders' split colocalization coefficients and maximum truncated Kendall tau correlation coefficient are evaluated.
#' @param method A string, indicating which plot is made. The plot can only be from \code{c('pearson', 'kendall', 'spearman', 'icq', 'mandersM1', 'mandersM2', 'trunkendall')}.
#' @param bsize An interger, the size of permuted block.
#' @param times An interger, the times of permutation for null distribution approximation
#' @param is.parallel A boolean value representing if the permutation is done by parallel computation
#' @param numcore An integer which specifies the number of core in parallel computation if is.parallel is TRUE. 
#' 
#' @import doParallel
#' 
#' @return An object which belongs to the class colocalroitest. 
#' \enumerate{
#' \item \code{Index_Value} The observed value of chosen colocalization measure without any permutation.
#' \item \code{P_value} The P-value calculated from permutation test.
#' \item \code{Null_Distribution} The value of chosen colocalization measure for each permutation.
#' }
#' 
#' @references Costes, Sylvain V., et al. "Automatic and quantitative measurement of protein-protein colocalization in live cells." Biophysical journal 86.6 (2004): 3993-4003.
#' 
#' @export
#' 
#' @examples
#' X=matrix(runif(64*64),64,64)
#' Y=matrix(runif(64*64),64,64)
#' colocalroitest(X,Y)
#' 
#' @author Shulei Wang
colocalroitest <- function(X, Y, mask = NULL, Thresholds=NULL, method = 'pearson', bsize = 1, times=100, is.parallel=FALSE, numcore=2)
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
  if (!is.integer(times))
    times = as.integer(times)
  if (times <= 0)
    stop("times must be positive")
  if (!is.logical(is.parallel))
    stop("is.logical must be logical")
  if (!is.integer(numcore))
    numcore = as.integer(numcore)
  if (numcore <= 0)
    stop("numcore must be positive")
  if (!is.integer(bsize))
    bsize = as.integer(bsize)
  if (bsize <= 0)
    stop("bsize must be positive")
  R = colocalroi(X, Y, mask, method)
  if (!is.null(mask))
  {
    AX=X[mask]
    AY=Y[mask]
    AX=as.vector(AX)
    AY=as.vector(AY)
    if (is.parallel)
    {
      cl <- makeCluster(numcore)
      registerDoParallel(cl)
      Re<-foreach (i=1:times,.combine='c') %dopar%
      {
        TempR<-colocalvec(sample(AX,length(AX)), AY, Thresholds, method)
        TempR
      }
      stopCluster(cl)
    }
    else {
      Re<-rep(0,times)
      for (i in 1: times)
      {
        TempR<-colocalvec(sample(AX,length(AX)), AY, Thresholds, method)
        Re[i]=TempR
      }
    }
  } else {
    if (bsize<=1)
    {
      AX=as.vector(X)
      AY=as.vector(Y)
      if (is.parallel)
      {
        cl <- makeCluster(numcore)
        registerDoParallel(cl)
        Re<-foreach (i=1:times,.combine='c') %dopar%
        {
          TempR<-colocalvec(sample(AX,length(AX)), AY, Thresholds, method)
          TempR
        }
        stopCluster(cl)
      }
      else {
        Re<-rep(0,times)
        for (i in 1: times)
        {
          TempR<-colocalvec(sample(AX,length(AX)), AY, Thresholds, method)
          Re[i]=TempR
        }
      }
    } else {
      if (is.parallel)
      {
        cl <- makeCluster(numcore)
        registerDoParallel(cl)
        Re<-foreach (i=1:times,.combine='c') %dopar%
        {
          TempR<-colocalroi(BlockPermu(X,bsize), Y, Thresholds, method)
          TempR
        }
        stopCluster(cl)
      }
      else
      {
        Re<-rep(0,times)
        for (i in 1: times)
        {
          TempR<-colocalroi(BlockPermu(X,bsize), Y, Thresholds, method)
          Re[i]=TempR
        }
      }
    }
  }
  Pvalue=(sum(Re>R)+1)/(times+1)
  ReList=list(0,0,c(0,0))
  names(ReList)<-c("Index_Value","P_value","Null_Distribution")
  ReList$Index_Value <- R
  ReList$P_value <- Pvalue
  ReList$Null_Distribution <- Re
  
  class(ReList) <- "colocalroitest"
  return (ReList)
}


#' @title Plot Permutation Test Result for an \code{colocalroitest} Object
#' 
#' @description \code{plot.colocalroitest} makes a plot for permutation test result, including a histgram of null distribution, a red dashed line representing oberserved value and P-value at the top.
#' 
#' @details The hitogram of null distribution relies on \code{ggplot2}.
#' 
#' @importFrom reshape2 melt
#' @import ggplot2
#' 
#' @param x An object of class \code{\link{colocalroitest}} 
#' @param ... Arguments to be passed to methods.
#' 
#' @return Invisible \code{NULL}
#' 
#' @examples  
#' X=matrix(runif(64*64),64,64)
#' Y=matrix(runif(64*64),64,64)
#' co <- colocalroitest(X,Y)
#' plot(co)  
#' 
#' @seealso \code{\link{colocalroitest}} 
#' 
#' @export
#' 
#' @author Shulei Wang
plot.colocalroitest <- function(x, ...)
{
  if (class(x)!="colocalroitest")
    stop("colocalroitestPlot is designed for class colocalroitest")
  DF=melt(x$Null_Distribution)
  n=nrow(DF)
  minDF=min(DF)
  maxDF=max(DF)
  binwidth=(maxDF-minDF)/n*5
  ggplot(DF, aes(x=value)) +
    geom_histogram(binwidth=binwidth, colour="black", fill="white") + 
    geom_vline(aes(xintercept=x$Index_Value),color="red", linetype="dashed", size=1) +
    labs(title=paste('Pvalue:',x$P_value))
}

#' @title Block permutation of matrix 
#'
#' @description \code{BlockPermu} permutes single channel of image/matrix block-wisely.
#' 
#' @details This function can randomly permute the matrix block-wisely. 
#' 
#' @seealso \code{\link{colocalroitest}}
#' 
#' @param X A matrix representing the channel of image/matrix.
#' @param BSize An integer that is the size of permuted block.
#' 
#' @return The matrix after block-wise permutation.
#' 
#' @export
#' 
#' @examples
#' X=matrix(runif(64*64),64,64)
#' BlockPermu(X,8)
#' 
#' @author Shulei Wang
BlockPermu <- function(X,BSize)
{
  Y=X
  nr<-nrow(X)
  nc<-ncol(X)
  nbr<-floor(nr/BSize)
  nbc<-floor(nc/BSize)
  for (i in 1:nbr)
  {
    SR<-(i-1)*BSize+1
    ER<-i*BSize
    Temp<-Y[SR:ER,]
    for (j in 1:nbc)
    {
      Perm<-sample(1:nbc, nbc)
      SC<-(j-1)*BSize+1
      EC<-j*BSize
      SNC<-(Perm[j]-1)*BSize+1
      ENC<-Perm[j]*BSize
      Temp[,SC:EC]=Y[SR:ER,SNC:ENC]
    }
    Y[SR:ER,]=Temp
  }
  for (i in 1:nbc)
  {
    SC<-(i-1)*BSize+1
    EE<-i*BSize
    Temp<-Y[,SC:EC]
    for (j in 1:nbr)
    {
      Perm<-sample(1:nbr, nbr)
      SR<-(j-1)*BSize+1
      ER<-j*BSize
      SNR<-(Perm[j]-1)*BSize+1
      ENR<-Perm[j]*BSize
      Temp[SR:ER,]=Y[SNR:ENR,SC:EC]
    }
    Y[,SC:EC]=Temp
  }
  return (Y)
}