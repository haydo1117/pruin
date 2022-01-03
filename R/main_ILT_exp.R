get_R <- function(delta, c, lambda, beta){
    ##
    stopifnot(delta>=0)
    stopifnot(c>0)
    stopifnot(lambda>0)
    stopifnot(beta>0)
    ##
    A <- c
    B <- c*beta-lambda-delta
    C <- -beta*delta
    Det <- B^2-4*A*C
    stopifnot(Det>=0)
    ##
    if(delta==0){
        return(B/A)
    }
    ##
    -1*(-B-sqrt(Det))/(2*A)
}

Lprob_ut <- function(delta,c,lambda,beta,u){
    ##
    stopifnot(delta>0)
    stopifnot(u>=0)
    r <- get_R(delta,c,lambda,beta)
    ##
    (1-r/beta)*exp(-r*u)/delta
}

#' Calculate the finite ruin probability of Poisson-Exponential process via Inverse Laplace transform using
#' Gaver-Stehfest algorithm
#'
#' Note Gaver-Stehfest algorithm is theoretically convergent but practically high precision arithmetic is required.
#' This is not achieved in R as R4.0+ uses "long-double" in C, which is 64bit floating point. As a result, using high
#' \code{n} would yield results divergent to infinity. Therefore, small \code{n} should be used. Using default arguments,
#' values for different \code{n} would be printed along with a plot which assists in determining an appropriate \code{n}.
#'
#' All arguments must be of length 1. For vectorized arguments, either use a loop or \code{apply} family, or \code{purrr}
#' package.
#'
#' Net profit condition is assumed, i.e. \code{c-lambda/beta > 0}.
#'
#' @param u Numerical. Initial surplus.
#' @param t Numerical. Terminal time.
#' @param c Numerical. Premium rate.
#' @param lambda Numerical. Claim arrival rate ("lambda" parameter of the Poisson process).
#' @param beta Numerical. Rate of the exponential distribution for the severity (mean of the distribution is \code{1/beta}).
#' @param ... Optional argument to be passed to \code{\link{fn_gs}}.
#'
#' @return Numeric. Vector of length \code{n} if \code{try_n} is \code{TRUE} (Default). Otherwise, length is 1.
#' @export
#'
#' @examples
#'
#' library(pruin)
#'
#' c <- 1
#' lambda <- 1
#' beta <- 1.2
#'
#' ruin_prob_exp_gs(1,1,c,lambda,beta)
#'
ruin_prob_exp_gs <- function(u,t,c,lambda,beta,...){
    ##
    stopifnot(is.numeric(u))
    stopifnot(is.numeric(t))
    stopifnot(is.numeric(c))
    stopifnot(is.numeric(lambda))
    stopifnot(is.numeric(beta))
    ##
    stopifnot(length(u)==1)
    stopifnot(length(t)==1)
    stopifnot(length(c)==1)
    stopifnot(length(lambda)==1)
    stopifnot(length(beta)==1)
    ##
    stopifnot(u>=0)
    stopifnot(t>=0)
    if(t==0){return(1)}
    stopifnot((c-lambda/beta) > 0)
    ##
    fn_gs(Lprob_ut,x=t,c=c,lambda=lambda,beta=beta,u=u,...)
}








