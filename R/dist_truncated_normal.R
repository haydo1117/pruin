
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
    mu+sigma*dnorm(-mu/sigma)/(1-pnorm(-mu/sigma))
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
            temp[[k+2]] <- (k-1)*(sigma^2)*temp[[k]] + mu*temp[[k+1]]-sigma*((-dnorm(-mu/sigma))/(1-pnorm(-mu/sigma)))
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




