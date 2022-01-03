#
#
#
# ## should have scale method, validate method
#
# validate_exponential <- function(beta){
#     stopifnot(beta>0)
#     invisible(0)
# }
#
# exponential <- function(beta){
#     ##
#     validate_exponential(beta)
#     ##
#     list(
#         name = "exponential",
#         par = list(beta=beta),
#         mean = 1/beta,
#         Theta = \(M) Theta_exp(M,beta),
#         ThetaDD = \(M) Theta_expDD(M,beta),
#         Theta_pe = \(M) Theta_pe_exp(M,beta),
#         Theta_pebar = \(M) Theta_pebar_exp(M,beta),
#         Theta_psi = \(M,c,lambda) Theta_psi_exp(M,c,lambda,beta),
#         #optional
#         scale = \(par, s) {par$beta <- par$beta/s;par},
#         psi_u = \(u,c,lambda) {
#             theta <- c/(lambda/beta)-1
#             r <- theta*beta/(1+theta)
#             ##
#             exp(-r*u)/(1+theta)
#         }
#     )
# }
#
#
#
#
#
# Theta_exp <- function(M,beta){
#     ##
#     # basically $$\int_0^\infty \beta exp(-\beta x) phi_k(x) dx$$
#     ##
#     stopifnot(beta>0)
#     stopifnot(M>=0)
#     ##
#     fk <- function(k) ((1-1/(beta+0.5))^k)*beta/(beta+0.5)
#     purrr::map_dbl(0:M,fk)
# }
#
# Theta_expDD <- function(M,beta){
#     ##
#     stopifnot(beta>0)
#     stopifnot(M>=2)
#     ##
#     theta <- Theta_exp(M,beta)
#     thetaD <- theta[-1]-theta[-length(theta)]
#     ##
#     thetaD[-1]-thetaD[-length(thetaD)]
# }
#
# Theta_pe_exp <- function(M,beta){
#     Theta_exp(M,beta)
# }
#
# Theta_pebar_exp <- function(M,beta){
#     Theta_pe_exp(M,beta)/beta
# }
#
# Theta_psi_exp <- function(M,c,lambda,beta){
#     ##
#     stopifnot(M>=0)
#     # compute up to M, return vector of size M+1
#     ##
#     stopifnot((c-lambda/beta) > 0)
#     stopifnot(c>0)
#     stopifnot(lambda>0)
#     stopifnot(beta>0)
#     ##
#     theta <- c/lambda*beta-1
#     ##
#     res <- rep(0,M+1)
#     Theta_pe <- Theta_pe_exp(M,beta)
#     Theta_pebar <- Theta_pebar_exp(M,beta)
#     divisor <- 1+theta-Theta_pe[[1]]
#     res[[1]] <- Theta_pebar[[1]]/divisor
#     if(M==0){
#         return(res[[1]])
#     }
#     ##
#     for (m in 1:M){
#         v <- Theta_pe[(m+1):1] # m:0
#         d <- v[-length(v)] - v[-1]
#         r <- res[1:m]
#         res[[m+1]] <- (Theta_pebar[[m+1]] + sum(r*d))/divisor
#     }
#     ##
#     res
# }
#
#
#
#
#
#
#
#
#
