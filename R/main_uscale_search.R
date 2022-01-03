uscale_search <- function(n=NULL,.cand=NULL,...){
    ##
    # create .cand if not supplied
    ##
    if(is.null(.cand)){
        validate_M(n,2)
        n_half <- floor(n/2)
        n2 <- n_half
        n1 <- n-n_half
        ##
        .cand <- 10^c(seq(-1,0,length.out=n1), seq(0,0.5,length.out=n2+1)[-1])
    }
    validate_M(length(.cand),2)
    ##
    f <- function(x){
        ruin_prob_ls(u_scale = x, include_error = TRUE,...)
    }
    res <- lapply(.cand, f)
    l_res <- purrr::map_dbl(res, \(x)x$res)
    l_err <- purrr::map_dbl(res, \(x)x$err)
    ##
    data.frame(
        u_scale = .cand,
        res = l_res,
        err = l_err
    ) |>
    dplyr::arrange(abs(err))
}
