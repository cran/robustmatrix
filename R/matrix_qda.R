# matrix_qda <- function(X, MU, SIGMA, OMEGA, priors){
#   p <- dim(X)[2]
#   q <- dim(X)[3]
#
#   k <- dim(MU)[1]
#
#   dets <- numeric()
#   for(l in 1:k){
#     dets[l] <- (determinant(SIGMA[l,,])$modulus*q+determinant(OMEGA[l,,])$modulus*p)
#   }
#
#   y_classes <- character()
#   dfs <- matrix(nrow = dim(X)[1], ncol = k)
#   for(i in 1:dim(X)[1]){
#     for(l in 1:k){
#       dfs[i,l] <- -0.5*dets[l] - 0.5*mmd(X = X[i,,], Mu = MU[l,,], Sigma = SIGMA[l,,], Omega = OMEGA[l,,]) + log(priors[l])
#     }
#
#     y_classes[i] <- which.max(dfs[i,])
#   }
#   return(list("y" = y_classes, "dfs" = dfs))
# }
