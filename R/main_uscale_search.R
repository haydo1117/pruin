#' Searching for the best values of u_scale for ruin_prob_ls
#'
#' It turns out with limited precision, scale change in \code{\link{ruin_prob_ls}} has a huge impact on the calculation.
#' With good \code{u_scale}, reliable \code{ruin_prob_ls} is computed.
#'
#' This function tries different values of \code{u_scale} and attach the corresponding \code{ruin_prob_ls} value and the
#' error of \code{ruin_prob_ls(u=0,t=Inf)}. A table is then return to assist to either chosing \code{u_scale} for
#' a separate run of \code{ruin_prob_ls}, or decision for a value for the ruin probability
#' (if the top values are close to each other).
#'
#' @param n Numerical of length 1, valued at least 2. Default to be 100. Ignore if \code{.cand} is used.
#' @param .cand Numerical vector of length at least 2. Candidates for \code{u_scale} in \code{ruin_prob_ls}.
#' Automatically generated using \code{n} if not provided.
#' @param ... Other arguments to be passed to \code{\link{ruin_prob_ls}}, i.e. \code{u}, \code{t}, \code{c},
#' \code{lambda} and \code{family}.
#'
#' @return A data.frame object with columns "u_scale", "res" and "err", sorted in the absolute value of "err".
#' @export
#'
#' @examples
#'
#' library(pruin)
#'
#' ## reliable example (small "err", top "res" are close and reasonable)
#' e <- exponential()(beta=1)
#' uscale_search(n=20,u=10,t=2,c=1.1,lambda=1,family=e)
#'
#' ## unreliable example (large "err", "res" unreasonable)
#' g <- gig()(a=0.32497,b=0.61543,alpha=-0.75)
#' uscale_search(n=20,u=4,t=4,c=1.1,lambda=1,family=g)
#'
#' ## default n is 100. It takes some time to run
#' \dontrun{
#'   uscale_search(u=4,t=4,c=1.1,lambda=1,family=g) # no need to specify n
#' }
#'
uscale_search <- function(n=100,.cand=NULL,...){
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
