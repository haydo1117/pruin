

pad_zero <- function(v,n){
    ##
    stopifnot(n>=0)
    stopifnot(length(v)>n)
    ##
    c(rep(0,n),v[1:(length(v)-n)])
}

vec2lowtM <- function(v){
    ##
    stopifnot(length(v)>1)
    ##
    purrr::map(1:length(v)-1,pad_zero,v=v) |> purrr::reduce(c) |> matrix(nrow = length(v))
}


G1 <- function(c,lambda,family){
    ##
    stopifnot(c>0)
    stopifnot(lambda>0)
    ##
    (c-1)/2-lambda+lambda*family$Theta(0)
}

G2 <- function(c,lambda,family){
    ##
    stopifnot(c>0)
    stopifnot(lambda>0)
    ##
    th <- family$Theta(1)
    ##
    (c+1)/2+lambda+lambda*(th[[2]]-2*th[[1]])
}

G3 <- function(c,lambda,family){
    ##
    stopifnot(c>0)
    stopifnot(lambda>0)
    ##
    (c+1)/2-lambda+lambda*family$Theta(0)
}

G4 <- function(c,lambda,family){
    ##
    stopifnot(c>0)
    stopifnot(lambda>0)
    ##
    th <- family$Theta(1)
    ##
    (c-1)/2+lambda+lambda*(th[[2]]-2*th[[1]])
}



Q1 <- function(c,lambda,family,M){
    ##
    stopifnot(M>2)
    ##
    v <- c(
        G1(c,lambda,family),
        G2(c,lambda,family),
        lambda*family$ThetaDD(M)
    )
    ##
    mat <- vec2lowtM(v)
    mat[1,] <- mat[1,]-c
    ##
    mat
}

Q2 <- function(c,lambda,family,M){
    ##
    stopifnot(M>2)
    ##
    v <- c(
        G3(c,lambda,family),
        G4(c,lambda,family),
        lambda*family$ThetaDD(M)
    )
    ##
    mat <- vec2lowtM(v)
    mat[1,] <- mat[1,]-c
    ##
    mat
}


Q3 <- function(c,lambda,family,M,check=FALSE){
    res <- solve(
        Q1(c,lambda,family,M),
        Q2(c,lambda,family,M)
    )
    if(check){
        ev <- eigen(res)$values
        cat(sprintf("Maximum magnitude of eigen value of Q3 is %s \n", max(abs(ev))))
    }
    ##
    res
}

Q4 <- function(c,lambda,family,M,t){
    ##
    stopifnot(t>0)
    ##
    q1 <- Q1(c,lambda,family,M)
    q3 <- Q3(c,lambda,family,M)
    ##
    i <- diag(rep(1,M+1))
    m0 <- 0.5*i+solve(i-q3,q3)
    e <- expm::expm(m0*(-t))
    ##
    e
}


Det <- function(c,lambda,family,M,check=FALSE,include_error=FALSE){
    ##
    stopifnot(M>=2)
    ##
    theta <- c/(lambda*family$mean)-1
    ##
    i <- diag(rep(1,M+1))
    q1 <- Q1(c,lambda,family,M)
    q3 <- Q3(c,lambda,family,M,check)
    ##
    v_M <- family$Theta_psi(M,c,lambda)
    v <- c(
        family$Theta_psi(0,c,lambda),
        v_M[-1]-v_M[-length(v_M)]
    )
    ##
    res <- as.numeric(solve(i-q3,solve(q1,v)))
    err <- 1/(1+theta)+sum(res)
    if(check){
        cat(sprintf("Approximation error for phi(0,t=inf) is %s \n",err))
    }
    if(include_error){
        return(list(res=res,err=err))
    }
    ##
    res
}


ruin_prob_ls <- function(u,t,c,lambda,family,M=20,psi_u=NULL,u_scale=1,t_scale=1,check=FALSE,include_error=FALSE){
    ##
    stopifnot(is.numeric(u))
    stopifnot(is.numeric(t))
    stopifnot(is.numeric(c))
    stopifnot(is.numeric(lambda))
    stopifnot(is.numeric(M))
    stopifnot(is.numeric(u_scale))
    stopifnot(is.numeric(t_scale))
    if(!is.null(psi_u)) stopifnot(is.numeric(psi_u))
    ##
    stopifnot(length(u)==1)
    stopifnot(length(t)==1)
    stopifnot(length(c)==1)
    stopifnot(length(lambda)==1)
    stopifnot(length(M)==1)
    stopifnot(length(u_scale)==1)
    stopifnot(length(t_scale)==1)
    if(!is.null(psi_u)) stopifnot(length(psi_u)==1)
    ##
    if(!is.null(family$scale)){
        u <- u*u_scale
        c <- c*u_scale
        family <- do.call(
            # family$name,
            do.call(family$name, list()),
            family$scale(family$par, u_scale))
        # browser()
    }else{
        if(u_scale!=1){
            warning(sprintf("No scale transform available for family: `%s`", family$name))
        }
    }
    ##
    t <- t*t_scale
    lambda <- lambda/t_scale
    c <- c/t_scale
    ##
    e <- Q4(c,lambda,family,M,t)
    dd <- Det(c,lambda,family,M,check,include_error)
    d <- if(!include_error) dd else dd$res
    ##
    phi_vec <- purrr::map_dbl(0:M,phi,x=u)
    ##
    h <- -1*sum(
        as.numeric(e %*% d) * phi_vec
    )
    ##
    if(is.null(psi_u)){
        psi_u <- sum(
            family$Theta_psi(M,c,lambda) * phi_vec
        )
    }
    ##
    res <- psi_u-h
    if(include_error){
        return(list(res=res,err=dd$err))
    }
    res
}

