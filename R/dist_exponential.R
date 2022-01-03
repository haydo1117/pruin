
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







