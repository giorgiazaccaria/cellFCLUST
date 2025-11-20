# Simulations cellFCLUST - Baseline 2
# Paper "Robust fuzzy clustering with cellwise outliers" by Zaccaria, G., 
# Benzakour, L., García-Escudero, L.A., Greselin, F., Mayo-Íscar, A.

rm(list = ls())

# Load the necessary packages.
required_packages <- c("mclust", "fclust", "doParallel", "parallel", "xtable", "ggplot2", "tidyr", "dplyr", "patchwork")
for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    message(paste("Package", package, "is not installed. Attempting to install now."))
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

# Load cellFCLUST
source("cellFCLUST.R", echo = TRUE)

set.seed(12345)

# Generation of the parameters depending on the scenario.
nsample <- 100
n.obs <- 500
G <- 4
p <- 10
alpha.out <- c(0)
alpha.mis <- 0
out.val <- alpha.out*n.obs*p
NA.val <- alpha.mis*n.obs*p

# Tuning parameter setting
nrep <- 40
nstart <- 10
niter <- 10
zero.tol <- 1e-16
rndstart <- 3
tuning_param_init <- matrix(0, nrow = length(alpha.out), ncol = 8)
for (set.out in 1:length(alpha.out)){
  alpha_tclust <- alpha.out[set.out]*2
  alpha_1 = alpha_2 <- alpha.out[set.out]
  alpha.A1 <- alpha.out[set.out]
  alpha.A2 <- alpha.out[set.out]*2
  tuning_param_init[set.out, ] <- c(alpha_tclust, alpha_1, alpha_2, alpha.A1, alpha.A2, 
                                    nrep, nstart, niter)
}

# Parameter setting
# Weights
pp <- c(0.2, 0.2, 0.3, 0.3)
# Labels
label.th <- NULL
for (g in 1:G) {
  label.th <- c(label.th, rep(g, n.obs*pp[g]))
}
# Means 
muval <- 1*pracma::rand(G-1, p-2)
mu <- matrix(0, G, p)
mu[2, 1:2] <- c(-0.5, 0.5)
mu[2, 3:p] <- sign((-1)^(1:(p-2)))*0.5
mu[3, ] <- 1:p%%3-1
mu[4, ] <- (2:(p+1)%%4-1)/2
# Var-cov matrices
Sigma <- vector (mode = "list", length = G)
Sigma[[1]] <- matrix(0, p, p)
for (j in 1:(p-1)) {
  for (jj in (j+1):p) {
    Sigma[[1]][j ,jj] <- (0.6)^(abs(j-jj))
  }
}
Sigma[[2]] <- matrix(0, p, p)
for (j in 1:(p-1)) {
  for (jj in (j+1):p) {
    Sigma[[2]][j ,jj] <- (-0.6)^(abs(j-jj))
  }
}
Sigma[[3]] <- matrix(0, p, p)
for (j in 1:(p-1)) {
  for (jj in (j+1):p) {
    Sigma[[3]][j ,jj] <- (0.7)^(abs(j-jj))
  }
}
Sigma[[4]] <- matrix(0, p, p)
for (j in 1:(p-1)) {
  for (jj in (j+1):p) {
    Sigma[[4]][j ,jj] <- (-0.7)^(abs(j-jj))
  }
}
Sigma[[1]] <- Sigma[[1]] + t(Sigma[[1]]) + diag(p)
Sigma[[2]] <- Sigma[[2]] + t(Sigma[[2]]) + diag(p)
Sigma[[3]] <- Sigma[[3]] + t(Sigma[[3]]) + diag(p)
Sigma[[4]] <- Sigma[[4]] + t(Sigma[[4]]) + diag(p)
Sigma.inv <- lapply(Sigma, solve)
maxfact.th <- max(unlist(lapply(Sigma, function(x){eigen(x)$values})))/min(unlist(lapply(Sigma, function(x){eigen(x)$values})))
maxfact <- 24

# Scaling the sigma
s = 16
Sigma[[1]] <- Sigma[[1]]/s
Sigma[[2]] <- Sigma[[2]]/s
Sigma[[3]] <- Sigma[[3]]/s
Sigma[[4]] <- Sigma[[4]]/s

# Data generation 
X <- vector(mode = "list", length = nsample)
Xout <- vector(mode = "list", length = length(alpha.out))
Wth <- vector(mode = "list", length = length(alpha.out))
for (set.out in 1:length(alpha.out)) {
  Xout[[set.out]] <- vector(mode = "list", length = nsample)
  Wth[[set.out]] <- vector(mode = "list", length = nsample)
  for (samp in 1:nsample) {
    # Data generation
    for (g in 1:G) {
      if (g == 1) {
        X[[samp]] <- MASS::mvrnorm(n.obs*pp[g], mu[g, ], Sigma[[g]], tol = .Machine$double.xmin)
      } else {
        X[[samp]] <- rbind(X[[samp]], MASS::mvrnorm(n.obs*pp[g], mu[g, ], Sigma[[g]], tol = .Machine$double.xmin))
      }
    }
    ## Contamination
    # Outlier generation
    Xout[[set.out]][[samp]] <- X[[samp]]
    Wth[[set.out]][[samp]] <- matrix(1, nrow = n.obs, ncol = p)
    for (j in 1:p){
      repl <- sample(n.obs, out.val[set.out]/p)
      Xout[[set.out]][[samp]][repl, j] <- runif(out.val[set.out]/p, min = -30, max = 30)
      Wth[[set.out]][[samp]][repl, j] <- 0
    }
    # Check outlyingness
    if (alpha.out[set.out] > 0) {
      out.rows <- which(rowSums(Wth[[set.out]][[samp]]) < p)
      for (i in 1:length(out.rows)) {
        out <- which(Wth[[set.out]][[samp]][out.rows[i], ] == 0)
        D <- sweep(mu, 2, Xout[[set.out]][[samp]][out.rows[i], ], FUN = "-")
        while (any(!as.matrix(lapply(1:G, function(ii) {D[ii, , drop = FALSE] %*% Sigma.inv[[ii]] %*% t(D[ii, , drop = FALSE])})) > qchisq(p = 0.99, df = p))) {
          Xout[[set.out]][[samp]][out.rows[i], out] <- runif(length(out), min = -30, max = 30)
          D <- sweep(mu, 2, Xout[[set.out]][[samp]][out.rows[i], ], FUN = "-")
        }
      }
    }
    # Missing generation
    miss <- 1:(n.obs*p)
    miss <- miss[!(miss %in% repl)]
    replNA <- sample(miss, NA.val)
    Xout[[set.out]][[samp]][replNA] <- NA
    Wth[[set.out]][[samp]][replNA] <- 0
  }
}

# Compute the memberships
m = 2
Uth <- vector(mode = "list", length = length(alpha.out))
for (set.out in 1:length(alpha.out)) {
  Uth[[set.out]] <- vector(mode = "list", length = nsample)
  for (samp in 1:nsample) {
    dif.pattern <- unique(Wth[[set.out]][[samp]])
    pat.unit <- match(data.frame(t(Wth[[set.out]][[samp]])), data.frame(t(dif.pattern))) 
    Uth[[set.out]][[samp]] <- update_U_new(Xout[[set.out]][[samp]], m, dif.pattern, pat.unit, pp, mu, Sigma)
  }
}

# Run the models

# cellFCLUST
result.cellFCLUST <- vector(mode = "list", length = length(alpha.out))
for (set.out in 1:length(alpha.out)){
  
  tuning_param_init_set.out = data.frame(
    alpha_tclust = tuning_param_init[set.out, 1], 
    alpha_1 = tuning_param_init[set.out, 2],
    alpha_2 = tuning_param_init[set.out, 3], 
    alpha.A1 = tuning_param_init[set.out, 4], 
    alpha.A2 = tuning_param_init[set.out, 5], 
    nrep = tuning_param_init[set.out, 6], 
    nstart = tuning_param_init[set.out, 7], 
    niter = tuning_param_init[set.out, 8])
  
  result.cellFCLUST[[set.out]] <- vector(mode = "list", length = nsample)
  
  ncores <- detectCores()
  cl <- makeCluster(ncores - 2)
  registerDoParallel(cl)
  
  result.cellFCLUST[[set.out]] <- foreach (i = 1:100, .packages = c("MASS", "tclust", "purrr"), .errorhandling = "pass") %dopar% {
    cellFCLUST(Xout[[set.out]][[i]], G = G, m = m, tuning_param_init = tuning_param_init_set.out, 
               alpha = 1-alpha.out[set.out], maxfact = maxfact, rndstart = 3)         
  }
  
  stopCluster(cl)
}
print("Finished.")

# F-TCLUST
set.seed(12345)
result.FTclust <- vector(mode = "list", length = length(alpha.out))

for (set.out in 1:length(alpha.out)){
  
  result.FTclust[[set.out]] <- vector(mode = "list", length = nsample)
  
  ncores <- detectCores()
  cl <- makeCluster(ncores - 2)
  registerDoParallel(cl)
  
  result.FTclust[[set.out]] <- foreach (i = 1:100, .packages = c("MASS", "tclust", "Rcpp"), .errorhandling = "pass") %dopar% {
    tclust:::.tclust.int(Xout[[set.out]][[i]], k = G, alpha = 0, restr.fact = maxfact, fuzzy = TRUE, m = 2)         
  }
  
  stopCluster(cl)
}
print("Finished.")

# FKM
set.seed(12345)
result.FKM <- vector(mode = "list", length = length(alpha.out))

for (set.out in 1:length(alpha.out)){
  
  result.FKM[[set.out]] <- vector(mode = "list", length = nsample)
  
  ncores <- detectCores()
  cl <- makeCluster(ncores - 2)
  registerDoParallel(cl)
  
  result.FKM[[set.out]] <- foreach (i = 1:100, .packages = c("MASS", "tclust", "purrr", "fclust"), .errorhandling = "pass") %dopar% {
    FKM(Xout[[set.out]][[i]], k = G, m = 1.35, RS = 50, conv = 1e-6, maxit = 500)         
  }
  
  stopCluster(cl)
}
print("Finished.")

# FKM.noise
set.seed(12345)
result.FKM.noise <- vector(mode = "list", length = length(alpha.out))

for (set.out in 1:length(alpha.out)){
  
  result.FKM.noise[[set.out]] <- vector(mode = "list", length = nsample)
  
  ncores <- detectCores()
  cl <- makeCluster(ncores - 2)
  registerDoParallel(cl)
  
  result.FKM.noise[[set.out]] <- foreach (i = 1:100, .packages = c("MASS", "tclust", "purrr", "fclust"), .errorhandling = "pass") %dopar% {
    FKM.noise(Xout[[set.out]][[i]], k = G, m = 1.35, RS = 50, conv = 1e-6, maxit = 500, delta = 10)         
  }
  
  stopCluster(cl)
}
print("Finished.")

# FKM.gkb
set.seed(12345)
result.FKM.gkb <- vector(mode = "list", length = length(alpha.out))

for (set.out in 1:length(alpha.out)){
  
  result.FKM.gkb[[set.out]] <- vector(mode = "list", length = nsample)
  
  ncores <- detectCores()
  cl <- makeCluster(ncores - 2)
  registerDoParallel(cl)
  
  result.FKM.gkb[[set.out]] <- foreach (i = 1:100, .packages = c("MASS", "tclust", "purrr", "fclust"), .errorhandling = "pass") %dopar% {
    FKM.gkb(Xout[[set.out]][[i]], k = G, m = 1.5, RS = 50, conv = 1e-6, maxit = 500, mcn = maxfact, vp = c(24, 24, 30, 30), gam = 0.1)         
  }
  
  stopCluster(cl)
}
print("Finished.")

# FKM.gkb.noise
set.seed(12345)
result.FKM.gkb.noise <- vector(mode = "list", length = length(alpha.out))

for (set.out in 1:length(alpha.out)){
  
  result.FKM.gkb.noise[[set.out]] <- vector(mode = "list", length = nsample)
  
  ncores <- detectCores()
  cl <- makeCluster(ncores - 2)
  registerDoParallel(cl)
  
  result.FKM.gkb.noise[[set.out]] <- foreach (i = 1:100, .packages = c("MASS", "tclust", "purrr", "fclust"), .errorhandling = "pass") %dopar% {
    FKM.gkb.noise(Xout[[set.out]][[i]], k = G, m = 1.5, RS = 50, conv = 1e-6, maxit = 500, mcn = maxfact, vp = c(24, 24, 30, 30), gam = 0.1, delta = 10)         
  }
  
  stopCluster(cl)
}
print("Finished.")

# UFTCP
source("UFTCP.R")
set.seed(12345)
result.UFTCP <- vector(mode = "list", length = length(alpha.out))

for (set.out in 1:length(alpha.out)){
  
  result.UFTCP[[set.out]] <- vector(mode = "list", length = nsample)
  
  ncores <- detectCores()
  cl <- makeCluster(ncores - 2)
  registerDoParallel(cl)
  
  result.UFTCP[[set.out]] <- foreach (i = 1:100, .packages = c("MASS", "tclust", "purrr", "fclust"), .errorhandling = "pass") %dopar% {
    UFtkmeans(0, Xout[[set.out]][[i]], G, 1.35, 0, maxiter = 500, tol = 1e-6, rndstart = 50)         
  }
  
  stopCluster(cl)
}
print("Finished.")
