

term_j <- function(j,k,n){
    temp <- j^(n+1)*choose(n,j)/factorial(n)*choose(j,k-j)*choose(2*j,j)
    sign_temp <- (-1)^(n+k)
    ##
    temp*sign_temp
}

akn <- function(k,n){
    ##
    stopifnot(n>=1)
    stopifnot(k>=1)
    stopifnot(k<=2*n)
    ##
    seq_temp <- seq(from=floor((k+1)/2), to=min(k,n))
    ##
    seq_temp |> purrr::map_dbl(term_j,k=k,n=n) |> purrr::reduce(sum)
}

akn_collection <- function(n){
    k_seq <- seq(1,2*n)
    k_seq |> purrr::map_dbl(akn,n=n)
}
