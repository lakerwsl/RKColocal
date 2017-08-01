#' @title Dual-Channel Image Visualization
#'
#' @description \code{DualImageVisual} is used to display dual channel image.
#' Each channel can be displayed individually or two channels can be displayed after merging together.
#' This function can also increase contrast ratio of plot by scaling.
#' If only region of interest is displayed, mask can be used.
#'
#' @details This function depends on \code{\link[EBImage]{display}} mainly.
#' First, image is preprocessed by scaling or masking.
#' Then, two plots are made when each channel is diplayed individually, and only one is displayed otherwise.
#' In addition, after scaling, the largest intensity in each channel is always 1.
#'
#' @seealso \code{\link{DualImageSave}}
#' @importFrom EBImage display
#' @importFrom abind abind
#' @importClassesFrom EBImage Image
#'
#' @param X A numerical matrix, the intensity matrix of the first channel/red channel.
#' @param Y A numerical matrix, the intensity matrix of the first channel/green channel.
#' @param isSplit An logical value. If it's \code{TRUE}, then two channels dispaly individually, and if \code{FALSE}, a image with merging two channels dispaly.
#' @param isScale An logical value to determine whether displayed image is scaled.
#' @param mask An logical matrix, where TRUE means pixel is involved into analysis and vice versa. Dimension of mask must be matched with \code{X} and \code{Y}.
#'
#' @return Invisible \code{NULL}
#'
#' @export
#'
#' @examples
#' X<-matrix(runif(10*10),nrow=10,ncol=10)
#' Y<-matrix(runif(10*10),nrow=10,ncol=10)
#' DualImageVisual(X,Y)
#'
#' @author Shulei Wang
DualImageVisual <- function (X, Y, isSplit=FALSE, isScale=TRUE, mask=NULL)
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
  ZeroMatrix=matrix(0,nrow = nrow(X),ncol = ncol(X))

  if (isSplit)
  {
    if (isScale)
    {
      red = X / max(X)
    }
    else
    {
      red = X
    }
    if (!is.null(mask))
    {
      red[!mask] = 0
    }
    red <- Image(abind(red,ZeroMatrix,ZeroMatrix,along = 3),colormode="Color")
    if (isScale)
    {
      green=Y/max(Y)
    }
    else
    {
      green = Y
    }
    if (!is.null(mask))
    {
      green[!mask] = 0
    }
    green <- Image(abind(ZeroMatrix,green,ZeroMatrix,along = 3),colormode="Color")
    par(mfrow=c(1,2))
    display(green,title = "green channel")
    display(red,title = "red channel")
  }
  else
  {
    if (isScale)
    {
      red = X / max(X)
      green=Y / max(Y)
    } else {
      red = X
      green = Y
    }
    if (!is.null(mask))
    {
      red[!mask] = 0
      green[!mask] = 0
    }
    imaged <- Image(abind(red,green,ZeroMatrix,along = 3),colormode="Color")
    display(imaged,title = "merged image")
  }
}

#' @title Dual-Channel Image Writing
#'
#' @description \code{DualImageSave} is used to save dual channel image.
#' Each channel can be saved individually or two channels can be saved as single overlayed image.
#' The same with \code{\link{DualImageVisual}}, this function can also increase contrast ratio of plot by scaling.
#' If only region of interest is saved, mask can be used.
#'
#' @details This function depends on \code{\link[EBImage]{writeImage}} mainly.
#' \code{X} and \code{Y} are preprocessed by scaling or masking firstly.
#' Then, tao images name of which end with 'red' and 'green' are saved if two channels are saved individually, or only one image is saved.
#' In addition, after scaling, the largest intensity in each channel is always 1.
#'
#' @seealso \code{\link{DualImageVisual}}
#' @importFrom EBImage writeImage
#' @importFrom abind abind
#' @importClassesFrom EBImage Image
#'
#' @param X A numerical matrix, the intensity matrix of the first channel/red channel.
#' @param Y A numerical matrix, the intensity matrix of the first channel/green channel.
#' @param isSplit An logical value. If it's \code{TRUE}, then two channels dispaly individually, and if \code{FALSE}, a image with merging two channels dispaly.
#' @param isScale An logical value to determine whether displayed image is scaled.
#' @param mask An logical matrix, where TRUE means pixel is involved into analysis and vice versa. Dimension of mask must be matched with one of image.
#' @param filepath A string indicating where the file shall be saved.
#' @param name A string for the name of saved image.
#'
#' @return Invisible \code{NULL}
#'
#' @export
#'
#' @examples
#' X<-matrix(runif(10*10),nrow=10,ncol=10)
#' Y<-matrix(runif(10*10),nrow=10,ncol=10)
#' DualImageSave(X,Y, name = 'sample.jpeg')
#'
#' @author Shulei Wang
DualImageSave <- function (X, Y, filepath=NULL, name, isSplit=FALSE, isScale=TRUE, mask=NULL)
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
  ZeroMatrix=matrix(0,nrow = nrow(X),ncol = ncol(X))

  if (isSplit)
  {
    if (isScale)
    {
      red = X / max(X)
    }
    else
    {
      red = X
    }
    if (!is.null(mask))
    {
      red[!mask] = 0
    }
    red <- Image(abind(red,ZeroMatrix,ZeroMatrix,along = 3),colormode="Color")
    if (isScale)
    {
      green=Y/max(Y)
    }
    else
    {
      green = Y
    }
    if (!is.null(mask))
    {
      green[!mask] = 0
    }
    green <- Image(abind(ZeroMatrix,green,ZeroMatrix,along = 3),colormode="Color")
    if (is.null(filepath))
    {
      greenname = paste("green_",name,sep = "")
      redname = paste("red_",name,sep = "")
    } else {
      greenname = paste(filepath,"green_",name,sep = "")
      redname = paste(filepath,"red_",name,sep = "")
    }
    writeImage(green, greenname, quality = 85)
    writeImage(red, redname, quality = 85)
  }
  else
  {
    if (isScale)
    {
      red = X / max(X)
      green=Y / max(Y)
    } else {
      red = X
      green = Y
    }
    if (!is.null(mask))
    {
      red[!mask] = 0
      green[!mask] = 0
    }
    imaged <- Image(abind(red,green,ZeroMatrix,along = 3),colormode="Color")
    if (!is.null(filepath))
    {
      name = paste(filepath,name,sep = "")
    }
    writeImage(imaged, name, quality = 85)
  }
}
