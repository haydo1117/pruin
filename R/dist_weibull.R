
#' Create a function to generate weibull family
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
#'   \item \code{alpha}: Numerical of length 1.
#'   \item \code{beta}: Numerical of length 1.
#' }
#' See \url{https://en.wikipedia.org/wiki/Generalized_inverse_Gaussian_distribution}.
#' \code{alpha} is \eqn{k} in the page.
#' \code{beta} is \eqn{\lambda} in the page.
#'
#' @export
#'
#' @examples
#'
#' library(pruin)
#'
#' alpha <- 2
#' beta <- 1
#'
#' family <- weibull()(alpha=alpha, beta=beta) # a list object
#'
#' family[c("name","par","mean")] # extract some information
#'
weibull <- function(M2=45,M3=32) {
    ##
    validate_M(M2,1)
    validate_M(M3,1)
    stopifnot(M2>M3)
    ##
    create_family(
        name = "weibull",
        par_list = rlang::pairlist2(alpha=,beta=),
        f_validate = validate_weibull,
        f_E = E_weibull,
        f_theta = Theta_weibull,
        f_scale = scale_weibull,
        M2=M2,M3=M3
    )
}


validate_weibull <- function(alpha, beta){
    stopifnot(length(alpha)==1)
    stopifnot(length(beta)==1)
    stopifnot(alpha>0)
    stopifnot(beta>0)
    invisible(0)
}

E_weibull <- function(alpha, beta) beta*gamma(1+1/alpha)

scale_weibull <- function(par,s){
    validate_s(s)
    par$beta <- par$beta*s
    ##
    par
}

Theta_weibull <- function(M,alpha,beta,nterms=40){
    ##
    # via moments
    ##
    validate_M(M)
    validate_weibull(alpha,beta)
    validate_M(nterms,1)
    ##
    # part1: evaluate $E[Y^j/j! exp(-Y/2)]$
    exp_ye_j <- function(n,j){
        ## n from 0 to infinity
        ## individual term
        (beta^j)*((-0.5*beta)^n)*gamma(1+(n+j)/alpha)/factorial(n)/factorial(j)
    }
    ##
    exp_ye <- function(j){
        purrr::map_dbl(0:nterms, exp_ye_j, j=j) |> sum()
    }
    # part2: evaluate $E[L_k(Y) exp(-Y/2)]
    el <- function(k){
        sum(
            purrr::map_dbl(0:k, exp_ye) * (-1)^(0:k) * choose(k, 0:k)
        )
    }
    ##
    purrr::map_dbl(0:M, el)
}




