
#' Create a function to generate exponential family
#'
#' A special case of \code{\link{comb_exponential}}.
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
#'   \item \code{beta}: Numerical of length 1. Rate parameter of the exponential distribution,
#'   i.e. mean is \code{1/beta}.
#' }
#' @export
#'
#' @examples
#'
#' library(pruin)
#'
#' beta <- 1
#'
#' family <- exponential()(beta=beta) # a list object
#'
#' family[c("name","par","mean")] # extract some information
#'
exponential <- function(M2=45,M3=32) {
    ##
    validate_M(M2,1)
    validate_M(M3,1)
    stopifnot(M2>M3)
    ##
    create_family(
        name = "exponential",
        par_list = rlang::pairlist2(beta=),
        f_validate = validate_exponential,
        f_E = E_exp,
        f_theta = Theta_exp,
        f_scale = scale_exp,
        M2=M2,M3=M3
    )
}

validate_exponential <- function(beta){
    stopifnot(length(beta)==1)
    stopifnot(beta>0)
    invisible(0)
}

E_exp <- function(beta){
    1/beta
}

scale_exp <- function(par,s){
    validate_s(s)
    par$beta <- par$beta/s
    ##
    par
}

Theta_exp <- function(M,beta){
    ##
    validate_M(M)
    validate_exponential(beta)
    ##
    fk <- function(k) ((1-1/(beta+0.5))^k)*beta/(beta+0.5)
    purrr::map_dbl(0:M,fk)
}







