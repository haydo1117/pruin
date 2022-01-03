

create_family <- function(name,
                          par_list,
                          f_validate,
                          f_E,
                          f_theta,
                          f_scale,
                          M2,M3){
    ##
    par_name <- names(par_list)
    args0 <- c(par_list,rlang::pairlist2(M2=M2,M3=M3))
    ##
    rlang::new_function(
        args0,
        rlang::expr({
            ##
            (!!f_validate)(!!!rlang::syms(par_name))
            validate_M(M2,min_M=1)
            ##
            E <- (!!f_E)(!!!rlang::syms(par_name))
            p0 <- purrr::map(!!(names(args0)), \(x) rlang::env_parent()[[x]])
            names(p0) <- !!(names(args0))
            ##
            par <- purrr::map(!!(names(par_list)), \(x) rlang::env_parent()[[x]])
            names(par) <- !!(names(par_list))
            f_theta2 <- function(M,par){
                par$M <- M
                ## update the following line
                do.call(!!f_theta,par)
            }
            ##
            list(
                name = !!name,
                par = p0,
                mean = E,
                Theta = \(M) f_theta2(M,par),
                ThetaDD = \(M) ThetaDD(M,f_theta2,par),
                Theta_pe = \(M) Theta_pe(M,M2,f_theta2,E,par),
                Theta_pebar = \(M) Theta_pebar(M,M2,M3,f_theta2,E,par),
                Theta_psi = \(M,c,lambda) Theta_psi(M,c,lambda,M2,M3,f_theta2,E,par),
                #optional
                scale = !!f_scale
            )

            ##
        })
    )
}

ThetaDD <- function(M,f_theta,par){
    ##
    validate_M(M,2)
    ##
    theta <- f_theta(M,par)
    thetaD <- theta[-1]-theta[-length(theta)]
    ##
    thetaD[-1]-thetaD[-length(thetaD)]
}

Theta_pe <- function(M,M2,f_theta,E,par){
    validate_M(M)
    validate_M(M2,min_M = M+1)
    validate_M(E)
    ##
    Theta <- f_theta(M2,par)
    ##
    res <- Theta[0:M+1]
    for (k in 0:M){
        temp_v <- Theta[(k+2):length(Theta)]
        alternate1 <- (-1)^(1:length(temp_v))
        res[[k+1]] <- res[[k+1]] + 2*sum(alternate1*temp_v)
    }
    ##
    res*2/E
}

Theta_pebar <- function(M,M2,M3,f_theta,E,par){
    ##
    validate_M(M3,min_M=M+1)
    validate_M(M2,min_M=M3+1)
    ##
    Theta_pe <- Theta_pe(M=M3,M2,f_theta,E,par)
    ##
    res <- Theta_pe[0:M+1]
    for (k in 0:M){
        temp_v <- Theta_pe[(k+2):length(Theta_pe)]
        alternate1 <- (-1)^(1:length(temp_v))
        res[[k+1]] <- res[[k+1]] + 2*sum(alternate1*temp_v)
    }
    ##
    2*res
}

Theta_psi <- function(M,c,lambda,M2,M3,f_theta,E,par){
    ##
    validate_M(M)
    validate_clambda(c,lambda)
    # compute up to M, return vector of size M+1
    ##
    stopifnot((c-lambda*E) > 0)
    ##
    theta <- c/(lambda*E)-1
    ##
    res <- rep(0,M+1)
    Theta_pe <- Theta_pe(M,M2,f_theta,E,par)
    Theta_pebar <- Theta_pebar(M,M2,M3,f_theta,E,par)
    divisor <- 1+theta-Theta_pe[[1]]
    res[[1]] <- Theta_pebar[[1]]/divisor
    if(M==0){
        return(res[[1]])
    }
    ##
    for (m in 1:M){
        v <- Theta_pe[(m+1):1] # m:0
        d <- v[-length(v)] - v[-1]
        r <- res[1:m]
        res[[m+1]] <- (Theta_pebar[[m+1]] + sum(r*d))/divisor
    }
    ##
    res
}


phi <- function(x,k){
    ##
    stopifnot(k>=0)
    stopifnot(x>=0)
    ##
    L <- as.function(mpoly::laguerre(k), silent=TRUE)
    ##
    L(x)*exp(-x/2)
}

phi_coef <- function(k){
    ##
    stopifnot(k>=0)
    ##
    f <- mpoly::laguerre(k)
    ##
    purrr::map_dbl(1:length(f),\(i)f[[i]][["coef"]])
}

validate_M <- function(M,min_M=0){
    stopifnot(M>=min_M)
    stopifnot(length(min_M)==1)
    stopifnot(length(M)==1)
    invisible(0)
}

validate_clambda <- function(c,lambda){
    stopifnot(length(c)==1)
    stopifnot(length(lambda)==1)
    stopifnot(c>0)
    stopifnot(lambda>0)
    invisible(0)
}

validate_s <- function(s){
    stopifnot(s>0)
    stopifnot(length(s)==1)
    invisible(0)
}
