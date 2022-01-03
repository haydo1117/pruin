
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


