
#' Create a function to generate truncated normal family
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
#'   \item \code{mu}: Numerical of length 1.
#'   \item \code{sigma}: Numerical of length 1.
#' }
#' @export
#'
#' @examples
#'
#' library(pruin)
#'
#' mu <- 1
#' sigma <- 1
#'
#' family <- truncated_normal()(mu=mu,sigma=sigma) # a list object
#'
#' family[c("name","par","mean")] # extract some information
#'
truncated_normal <- function(M2=45,M3=32) {
    ##
    validate_M(M2,1)
    validate_M(M3,1)
    stopifnot(M2>M3)
    ##
    create_family(
        name = "truncated_normal",
        par_list = rlang::pairlist2(mu=,sigma=),
        f_validate = validate_tnorm,
        f_E = E_tnorm,
        f_theta = Theta_tnorm,
        f_scale = scale_tnorm,
        M2=M2,M3=M3
    )
}

validate_tnorm <- function(mu, sigma){
    stopifnot(length(mu)==1)
    stopifnot(length(sigma)==1)
    stopifnot(sigma>0)
    invisible(0)
}

E_tnorm <- function(mu, sigma){
    mu+sigma*stats::dnorm(-mu/sigma)/(1-stats::pnorm(-mu/sigma))
}

scale_tnorm <- function(par,s){
    validate_s(s)
    par$mu <- par$mu*s
    par$sigma <- par$sigma*s
    ##
    par
}


eval_moments_tnorm <- function(M,mu,sigma){
    ##
    validate_M(M)
    if(M==0){return(1)}
    validate_tnorm(mu,sigma)
    ##
    temp <- rep(0,M+2)
    temp[[2]] <- 1
    for (k in 1:M){
        if(k==1){
            temp[[k+2]] <- (k-1)*(sigma^2)*temp[[k]] + mu*temp[[k+1]]-sigma*((-stats::dnorm(-mu/sigma))/(1-stats::pnorm(-mu/sigma)))
        }else{
            temp[[k+2]] <- (k-1)*(sigma^2)*temp[[k]] + mu*temp[[k+1]]
        }
    }
    ##
    temp[-1]
}


Theta_tnorm <- function(M,mu,sigma,nterms=40){
    ##
    # via moments
    ##
    validate_M(M)
    validate_tnorm(mu,sigma)
    validate_M(nterms,1)
    ##
    # part1: evaluate $E[Y^j/j! exp(-Y/2)]$
    mm <- eval_moments_tnorm(nterms+M,mu,sigma)
    exp_ye_j <- function(n,j){
        ## n from 0 to infinity
        ## individual term
        ((-0.5)^n)*mm[[n+j+1]]/factorial(n)/factorial(j)
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




