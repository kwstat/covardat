# covar.r
# Time-stamp: c:/x/notes/covar.r
# Copyright 2017, Kevin Wright

# todo: toeplitz

rvar <- function(obj, ...) UseMethod("rvar")

rvar.asreml <- function(obj){
  s2 <- obj$sigma2
  s2 * rcor(obj)
}

# ----------------------------------------------------------------------------

# Extract the R-side correlation matrix
#'
#' Extract the R-side variance components from a fitted model (asreml or lme)
#' and format them as matrices.
#' 
#' @param obj A fitted model object.
#' 
#' @return A list of matrices.
#' @export
#' @docType methods
#' @rdname rcor-methods
#' @examples \dontrun{
#' require('agridat')
#' dat <- hanks.sprinkler
#' dat <- transform(dat, subf=factor(subplot),
#'                  irrf=factor(irr))
#' dat <- dat[order(dat$block, dat$gen, dat$subplot),]
#' m1 <- asreml(yield ~ gen + dir + irrf + gen:dir + gen:irrf + dir:irrf,
#'              data=dat,
#'              random= ~ block + block:dir + block:irrf,
#'              rcov= ~ block:gen:corb(subf, k=3))
#' require("lucid")
#' lucid(rcor(m1)$subf[1:8,1:8],dig=2)
#' }
rcor <- function(obj) UseMethod("rcor")

# ----------------------------------------------------------------------------

#' @method rcor gls
#' @importFrom nlme corMatrix
#' @export
#' @rdname rcor-methods
rcor.gls <- function(obj){
  # Extract the R-side correlation matrix from lme object
  rc <- obj$modelStruct$corStruct
  # corMatrix returns a list, one item per group
  rc <- corMatrix(rc)[[1]]
  return(rc)
}

# ----------------------------------------------------------------------------

#' @method rcor asreml
#' @export
#' @rdname rcor-methods
rcor.asreml <- function(obj){
  # asreml stores variance parameters in lists.
  # The top-level list is for each at() item.
  # Within each item of the top-level list, each item in the second-level list
  # is for one of the terms in the Kronecker product
  if(length(obj$R.param) > 1) {
    stop("Not available when 'at()' is part of rcov.  Try something like: \n kw:::make.ar1(obj$R.param[[1]][[2]]")
  }

  rparam <- obj$R.param$R
  nmat <- length(rparam)
  out <- vector(length=nmat, mode="list") # Initialize a list
  for(ii in 1:nmat){
    # For each item in list, turn the parameters into a matrix
    rpi <- rparam[[ii]]

    if(names(rpi)[1]=="s2")
      oi <- make.var(rpi)
    else if(rpi$model=="ante")
      oi <- make.ante(rpi)
    else if(rpi$model=="ar1")
      oi <- make.ar1(rpi)
    else if(rpi$model=="id")
      oi <- make.id(rpi)
    else if(rpi$model=="cor")
      oi <- make.cor(rpi)
    else if(rpi$model=="corb") # aka toeplitz
       oi <- make.corb(rpi)
    else if(rpi$model=="corg")
       oi <- make.corg(rpi)
    else if(rpi$model=="facv")
       oi <- make.facv(rpi)
    else if(rpi$model=="us")
       oi <- make.us(rpi)
    out[[ii]] <- oi
  }
  names(out) <- names(rparam)
  class(out) <- c("asremlrmat", class(out))

  attr(out, "formula") <- obj$call$rcov

  return(out)
}

#' @method print asremlrmat
#' @param x Object to print
#' @param n Size of matrix corner to print
#' @param dig Number of digits for printing
#' @param ... Other arguments (not used)
#' @import lucid
#' @importFrom stats cov2cor
#' @export
#' @rdname rcor-methods
print.asremlrmat <- function(x, n=6, dig=3, ...){
  # Use n=0 to print the full matrix
  print(attr(x, "formula"))
  for(ii in 1:length(x)){
    xi <- x[[ii]]

    if(n>0 & n<ncol(xi))
      tmp <- paste(" [1:",n,"] ",sep="")
    else tmp <- ""
    cat(names(x)[[ii]], tmp, ":\n", sep="")

    if(n==0L)
      print(lucid(xi, dig=dig), quote=FALSE)
    else
      print(lucid(xi[1:min(nrow(xi),n) ,
               1:min(ncol(xi),n)], dig=dig), quote=FALSE)

  }
  invisible(x)
}

make.ante <- function(mi){
  mat <- diag(length(mi$levels))
  mat[upper.tri(mat,diag=TRUE)] <- mi$initial
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  colnames(mat) <- rownames(mat) <- mi$levels
  cov2cor(mat)
}

make.ar1 <- function(mi){
  mat <- diag(length(mi$levels))
  mat <- mi$initial ^ abs(row(mat) - col(mat))
  colnames(mat) <- rownames(mat) <- mi$levels
  mat
}

make.facv <- function(mi){
  # Right now only works with facv( , 1)
  fa1 <- mi$initial[grep("\\.fa1", names(mi$initial))]
  sv <- mi$initial[grep("\\.var", names(mi$initial))]
  mat <- fa1 %*% t(fa1) + diag(sv)
  colnames(mat) <- rownames(mat) <- mi$levels
  mat
}

make.var <- function(mi){
  mat <- matrix(mi$s2)
}

make.id <- function(mi){
  mat <- diag(length(mi$levels))
  colnames(mat) <- rownames(mat) <- mi$levels
  mat
}

make.cor <- function(mi){
  # uniform correlation, aka compound symmetry
  mat <- diag(length(mi$levels))
  mat[lower.tri(mat)] <- mi$initial
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  colnames(mat) <- rownames(mat) <- mi$levels
  return(mat)
}

make.corb <- function(mi){
  # banded correlation
  mat <- diag(length(mi$levels))
  mat[lower.tri(mat)] <-
    mi$initial[abs(row(mat)-col(mat))[lower.tri(mat)]]
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  colnames(mat) <- rownames(mat) <- mi$levels
  return(mat)
}

make.corg <- function(mi){
  # general correlation, i.e. unstructured off-diagonal
  mat <- diag(length(mi$levels))
  # Fill the upper tri, column-wise
  mat[upper.tri(mat)] <- mi$initial
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  colnames(mat) <- rownames(mat) <- mi$levels
  return(mat)
}

make.toep <- function(mi){
  # toeplitz
  #return(mat)
}

make.us <- function(mi){
  #browser()
  # unstructured
  mat <- diag(length(mi$levels))
  # Fill the upper tri, column-wise
  #mat[upper.tri(mat)] <- mi$initial
  mat[upper.tri(mat, diag=TRUE)] <- mi$initial
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  colnames(mat) <- rownames(mat) <- mi$levels
  return(mat)
}
