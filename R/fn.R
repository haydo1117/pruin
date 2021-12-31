fn_gs0 <- function(Fn,x,n=7,...){
    ##
    stopifnot(typeof(Fn)=="closure")
    stopifnot(x>0)
    stopifnot(n>=1)
    ##
    temp_f <- function(.k123456,...){
        akn(.k123456,n)*Fn(.k123456*log(2)/x,...)/x
    }
    temp <- seq(from=1,to=2*n) |> purrr::map_dbl(temp_f,...) |> purrr::reduce(sum)
    ##
    temp*log(2)
}

#' Compute Inverse Laplace transform using Gaver-Stehfest algorithm
#'
#' Note Gaver-Stehfest algorithm is theoretically convergent but practically high precision arithmetic is required.
#' This is not achieved in R as R4.0+ uses "long-double" in C, which is 64bit floating point. As a result, using high
#' \code{n} would yield results divergent to infinity. Therefore, small \code{n} should be used. Using default arguments,
#' values for different \code{n} would be printed along with a plot which assists in determining an appropriate \code{n}.
#'
#' @param Fn Function. Function to be inversed.
#' @param x Numeric (positive). Location of the function to be approximated.
#' @param n Numeric. Number of terms used in the Gaver-Stehfest algorithm.
#' @param try_n Logical. Whether all n values up to \code{n} are evaluated. Useful to determine \code{n}.
#' @param plot Logical. In case \code{try_n} is \code{TRUE}, a plot is also generated.
#' @param ... Other arguments to be passed to \code{Fn}. (".k123456" is reserved. An event of argument crash should be unlikely...)
#'
#' @return Numeric. Vector of length \code{n} if \code{try_n} is \code{TRUE} (Default). Otherwise, length is 1.
#' @export
#'
#' @examples
#'
#' library(pruin)
#'
#' ## constant function
#' Fconst <- function(x,k)k/x
#'
#' fn_gs(Fconst,2,n=15,k=2)
#'
fn_gs <- function(Fn,x,n=10,try_n=TRUE,plot=TRUE,...){
    stopifnot(n>=1)
    if(!try_n){
        return(fn_gs0(Fn,x,n,...))
    }
    temp <- seq(1,n) |> purrr::map_dbl(fn_gs0,Fn=Fn,x=x,...)
    names(temp) <- sprintf("n=%s",1:n)
    ##
    if(plot){
        plot(1:n,temp,ylab="function value",xlab="values of `n`")
    }
    ##
    temp
}
