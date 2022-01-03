
#' Create a function to generate generalised inverse gaussian family
#'
#'
#' A function to create "family" objects is returned. The "family" object is to be used in the \code{family}
#' argument in \code{\link{ruin_prob_ls}}. It contains all essential information to calculate the
#' ruin probability, e.g. calculation of \eqn{\Theta_{f,k}} for various functions \eqn{f}.
#'
#' Appendix A3 in article "Finite-time ruin probabilities using bivariate Laguerre series" is used.
#' The sum to infinity terms in (A18), (A19) are truncated at \code{M2} and \code{M3} (\code{M2>M3}).
#'
#' Ideally, large \code{M2} and \code{M3} should be used. However, due to limit in numerical precision in R
#' (long double in C, i.e. 64bit) calculation of \eqn{\Theta_{f,k}} for large \eqn{k} is unstable. As a result,
#' \code{M2} and \code{M3} has to be moderate.
#'
#' Computer programs with high precision arithmetic built-in (e.g. Mathematica or Matlab) should be used for
#' reliable results.
#'
#' @param M2 Numerical of length 1.
#' @param M3 Numerical of length 1.
#'
#' @return A function to create a family object with parameters:
#' \itemize{
#'   \item \code{a}: Numerical of length 1.
#'   \item \code{b}: Numerical of length 1.
#'   \item \code{alpha}: Numerical of length 1.
#' }
#' See \url{https://en.wikipedia.org/wiki/Generalized_inverse_Gaussian_distribution}.
#' \code{alpha} is \eqn{p} in the page.
#'
#' @export
#'
#' @examples
#'
#' library(pruin)
#'
#' a <- 5
#' b <- 0.02
#' alpha <- 2.5
#'
#' family <- gig()(a=a,b=b,alpha=alpha) # a list object
#'
#' family[c("name","par","mean")] # extract some information
#'
gig <- function(M2=45,M3=32) {
    ##
    validate_M(M2,1)
    validate_M(M3,1)
    stopifnot(M2>M3)
    ##
    create_family(
        name = "gig",
        par_list = rlang::pairlist2(a=,b=,alpha=),
        f_validate = validate_gig,
        f_E = E_gig,
        f_theta = Theta_gig,
        f_scale = scale_gig,
        M2=M2,M3=M3
    )
}

validate_gig <- function(a,b,alpha){
    stopifnot(length(a)==1)
    stopifnot(length(b)==1)
    stopifnot(length(alpha)==1)
    stopifnot(a>0)
    stopifnot(b>0)
    stopifnot(is.numeric(alpha))
    invisible(0)
}

E_gig <- function(a,b,alpha){
    sqrt(b)*besselK(x=sqrt(a*b), nu=alpha+1)/
        sqrt(a)/besselK(x=sqrt(a*b), nu=alpha)
}

scale_gig <- function(par, s) {
    validate_s(s)
    par$a <- par$a/s
    par$b <- par$b*s
    ##
    par
}

Theta_gig <- function(M,a,b,alpha){
    ##
    validate_gig(a,b,alpha)
    validate_M(M)
    ## j
    fj <- function(j){
        stopifnot(j>=0)
        t1a <- (a/b)^(alpha/2)
        t1b <- ((a+1)/b)^((alpha+j)/2)
        t2a <- besselK(x=sqrt((a+1)*b),nu=alpha+j)
        t2b <- besselK(x=sqrt(a*b),nu=alpha)
        ##
        t1a/t1b*t2a/t2b/factorial(j)
    }
    ## k
    fk <- function(k){
        stopifnot(k>=0)
        jv <- purrr::map_dbl(0:k, fj)
        tv <- (-1)^(0:k)*choose(k,0:k)
        ##
        sum(jv*tv)
    }
    ##
    purrr::map_dbl(0:M,fk)
}


