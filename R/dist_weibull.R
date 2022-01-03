
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




