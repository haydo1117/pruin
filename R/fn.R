fn_gs0 <- function(Fn,x,n=7){
    ##
    stopifnot(typeof(Fn)=="closure")
    stopifnot(x>0)
    stopifnot(n>=1)
    ##
    temp_f <- function(k){
        akn(k,n)*Fn(k*log(2)/x)/x
    }
    temp <- seq(from=1,to=2*n) |> purrr::map_dbl(temp_f) |> purrr::reduce(sum)
    ##
    temp*log(2)
}

fn_gs <- function(Fn,x,n=20,try_n=TRUE,plot=TRUE){
    stopifnot(n>=1)
    if(!try_n){
        return(fn_gs0(Fn,x,n))
    }
    temp <- seq(1,n) |> purrr::map_dbl(fn_gs0,Fn=Fn,x=x)
    names(temp) <- sprintf("n=%s",1:n)
    ##
    if(plot){
        plot(1:n,temp)
    }
    ##
    temp
}
