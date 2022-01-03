
#' Create a function to generate combination of exponential family
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
#'   \item \code{w}: Numerical vector sum to 1 describing the weights for different exponentials. Negative weights can be used.
#'   \item \code{beta}: Numerical vector. Rate of the individual exponentials.
#' }
#' @export
#'
#' @examples
#'
#' library(pruin)
#'
#' w <- c(2,-1)
#' beta <- c(1,2)
#'
#' family <- comb_exponential()(w=w,beta=beta) # a list object
#'
#' family[c("name","par","mean")] # extract some information
#'
comb_exponential <- function(M2=45,M3=32) {
    ##
    validate_M(M2,1)
    validate_M(M3,1)
    stopifnot(M2>M3)
    ##
    create_family(
        name = "comb_exponential",
        par_list = rlang::pairlist2(w=,beta=),
        f_validate = validate_comb_exponential,
        f_E = E_comb_exp,
        f_theta = Theta_comb_exp,
        f_scale = scale_comb_exp,
        M2=M2,M3=M3
    )
}

validate_comb_exponential <- function(w, beta){
    stopifnot(length(w)>0)
    stopifnot(length(w)==length(beta))
    stopifnot(all(beta>0))
    stopifnot(sum(w)==1)
    invisible(0)
}

E_comb_exp <- function(w, beta){
    sum(w*1/beta)
}

scale_comb_exp <- function(par,s){
    stopifnot(s>0)
    par$beta <- par$beta/s
    ##
    par
}

Theta_comb_exp <- function(M,w,beta){
    ##
    validate_M(M)
    validate_comb_exponential(w,beta)
    ##
    fk <- function(k,w) sum(w*((1-1/(beta+0.5))^k)*beta/(beta+0.5))
    purrr::map_dbl(0:M,fk,w=w)
}







