
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







