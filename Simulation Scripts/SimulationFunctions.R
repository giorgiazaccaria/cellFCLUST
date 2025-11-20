######## FUNCTIONS #######
LabelSwitchCFC <- function(mu, mu.th, model){
  G = nrow(mu.th)
  permut = pracma::perms(1:G)
  
  rmse <- c()
  
  for (i in 1:nrow(permut)){
    rmse[i] <- sqrt(mean((mu.th - mu[permut[i, ], ])^2))
  }
  
  idx = permut[which.min(rmse), ]
  best_model = model
  for (g in 1:G){
    best_model$mu[g, ] = model$mu[idx[g], ]
    best_model$sigma[[g]] = model$sigma[[idx[g]]]
    best_model$pp[g] = model$pp[idx[g]]
    best_model$U[, g] = model$U[, idx[g]]
    best_model$X.imp[,,g] = model$X.imp[,,idx[g]]
  }
  best_model$label = apply(best_model$U, 1, which.max)
  return(best_model)
}

LabelSwitchFKM <- function(mu, mu.th, model){
  G = nrow(mu.th)
  permut = pracma::perms(1:G)
  
  rmse <- c()
  
  for (i in 1:nrow(permut)){
    rmse[i] <- sqrt(mean((mu.th - mu[permut[i, ], ])^2))
  }
  
  idx = permut[which.min(rmse), ]
  best_model = model
  for (g in 1:G){
    best_model$H[g, ] = model$H[idx[g], ]
    best_model$U[, g] = model$U[, idx[g]]
    best_model$F[,,g] = model$F[,,idx[g]]
  }
  best_model$clus[, 1] = apply(best_model$U, 1, which.max)
  best_model$clus[, 2] = apply(best_model$U, 1, max)
  
  return(best_model)
}

LabelSwitchFTC <- function(mu, mu.th, model){
  G = nrow(mu.th)
  permut = pracma::perms(1:G)
  
  rmse <- c()
  
  for (i in 1:nrow(permut)){
    rmse[i] <- sqrt(mean((mu.th - mu[permut[i, ], ])^2))
  }
  
  idx = permut[which.min(rmse), ]
  best_model = model
  for (g in 1:G){
    best_model$centers[, g] = model$centers[, idx[g]]
    best_model$z[, g] = model$z[, idx[g]]
    best_model$cov[,,g] = model$cov[,,idx[g]]
    best_model$weights[g] = model$weights[idx[g]]
  }
  best_model$cluster[rowSums(best_model$z)!=0] = apply(best_model$z, 1, which.max)[rowSums(best_model$z)!=0]
  
  return(best_model)
}


split.along.dim <- function(a, n){
  setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                  array, dim = dim(a)[-n], dimnames(a)[-n]),
           dimnames(a)[[n]])
}

LabelSwitchUFTCP <- function(mu, mu.th, model){
  G = nrow(mu.th)
  permut = pracma::perms(1:G)
  
  rmse <- c()
  
  for (i in 1:nrow(permut)){
    rmse[i] <- sqrt(mean((mu.th - mu[permut[i, ], ])^2))
  }
  
  idx = permut[which.min(rmse), ]
  best_model = model
  for (g in 1:G){
    best_model$centers[g, ] = model$centers[idx[g], ]
    best_model$covariance[[g]] = model$covariance[[idx[g]]]
    best_model$pp[g] = model$pp[idx[g]]
    best_model$U[, g] = model$U[, idx[g]]
  }
  best_model$label = apply(best_model$U, 1, which.max)
  return(best_model)
}
