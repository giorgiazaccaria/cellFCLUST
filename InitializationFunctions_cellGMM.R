## Functions for the initialization of cellGMM #################################

################################################################################
w_initial <- function (y, 
                       alpha_tclust = 0.2, 
                       alpha_1 = 0.1, 
                       alpha_2 = 0.1, 
                       K, 
                       met = 1) {
  n <- nrow(y)
  p <- ncol(y)
  # UNIVARIATE tclust
  tc1 <- list()
  for (j1 in 1:p) {
    tc1[[j1]] <- tclust::tclust(y[!is.na(y[, c(j1)]), c(j1)], k = K, alpha = alpha_tclust)
  }
  # BIVARIATE tclust
  tc2 <- list()
  for (j1 in (1:(p-1))) {
    tc2[[j1]] <- list()
    for (j2 in ((j1+1):p)) {
      tc2[[j1]][[j2]] <- tclust::tclust(y[apply(is.na(y[, c(j1,j2)]), 1, sum) == 0 , c(j1,j2)], k = K, alpha = alpha_tclust)
    }
  }
  # UNIVARIATE identification
  f1 <- array(NA, c(n, p))
  for (j1 in (1:p)) {
    dm1 <- array(NA, c(n, K))
    for (k in 1:tc1[[j1]]$k){
      y_ <- NA
      y_ <- y[, c(j1)] - tc1[[j1]]$centers[, k]
      dm1[, k] <- apply(y_%*% solve(tc1[[j1]]$cov[,,k], tol = .Machine$double.xmin)*y_, 1, sum)
    }
    f1[, j1] <- apply(dm1, 1, min, na.rm = TRUE)
  }
  qq1 <- quantile(f1, 1-alpha_1, na.rm = TRUE)
  ww1 <- (f1<=qq1) + 0
  # BIVARIATE identification
  f2 <- array(NA,c(n, p, p))
  for (j1 in (1:(p-1))) {
    for (j2 in ((j1+1):p)) {
      dm2 <- array(NA,c(n,K))
      for (k in 1:tc2[[j1]][[j2]]$k){
        y_ <- y[, c(j1, j2)] - matrix(1, nrow = n, ncol = 1)%*%rbind(tc2[[j1]][[j2]]$centers[, k])
        dm2[, k] <- apply(y_%*% solve(tc2[[j1]][[j2]]$cov[,,k])*y_, 1, sum)
      }
      f2[, j1, j2] <- apply(dm2, 1, min, na.rm = TRUE)
    }
  }
  for (j in 1:p) {f2[ww1[, j] == 0, j, ] <- NA}
  for (j in 1:p) {f2[ww1[, j] == 0, , j] <- NA}
  
  ff2 <- array(NA, c(n, p))
  if (met==1) for (j in 1:p) {ff2[, j] <- apply(f2[, j, ], c(1), sum, na.rm = TRUE) + apply(f2[, , j], c(1), sum, na.rm = TRUE)}
  if (met==2) for (j in 1:p) {ff2[, j] <- max(apply(f2[, j, ],c(1), max, na.rm = TRUE), apply(f2[, , j], c(1), max, na.rm = TRUE), na.rm = TRUE)}
  ff2[ww1 == 0] = NA
  
  qq2 <- quantile(ff2, 1-alpha_1/(1-alpha_2), na.rm = TRUE)
  ww2 <- (ff2<=qq2) + 0
  ww2[is.na(ww2)] <- 1
  ww <- ww1*ww2
  ww[is.na(y)] <- 0
  www = list()
  www$ww1 <- ww1
  www$ww2 <- ww2
  www$ww <- ww
  return(www)
}

################################################################################
prepars_initial <- function (x = x, 
                             K = K, 
                             alpha = alpha, 
                             q = floor(p/2) + 1, 
                             nrep = 40, 
                             nstart = 10, 
                             niter = 10) {
  
  p <- ncol(x)
  s <- list()
  rws <- list()
  for (i in 1:nrep) {
    s[[i]] <- sample(1:p, q, replace = FALSE)
    ind_s <- (apply(!is.na(x[, s[[i]]]), 1, sum) == q)
    y <- x[ind_s, c(s[[i]])]
    tcl <- try(tclust::tclust(x = y, k = K, alpha = alpha, drop.empty.clust = FALSE), silent = TRUE) # Added "drop.empty.clust = FALSE" !!!
    if (!is.numeric(tcl[[1]])) {
      tcl <- tclust::tkmeans(x = y, k = K)
      ss <- array(NA, c(length(s[[i]]), ncol(tcl$center)))
      for (kk in 1:ncol(tcl$center)) {ss[, kk] = apply(rbind(y[tcl$cluster==kk, ]), 2, sd)^2}
      s2 <- mean(ss, na.rm = TRUE)
      tcl$cov <- array(NA, c(length(s[[i]]), length(s[[i]]), ncol(tcl$center)))
      for (kk in 1:ncol(tcl$center)) {tcl$cov[, , kk] <- s2*diag(rep(1, length(s[[i]])))}
    }
    
    tc <- t(tcl$centers)
    ctr <- matrix(NA, nrow = nrow(tc), ncol = p)
    ctr[, s[[i]]] <- tc
    cts <- array(NA, c(p, p, nrow(tc)))
    cts[s[[i]], s[[i]], ] <- tcl$cov
    for (k in 1:nrow(tc)) {if (k == 1) {cts__ <- cts[, , 1]} else {cts__ <- rbind(cts__, cts[, , k])}}
    
    if (i == 1) {ctr_ <- ctr; rws[[i]] <- 1:nrow(ctr)} else {
      rws[[i]] <- (nrow(ctr_) + 1):(nrow(ctr_) + nrow(ctr))
      ctr_ <- rbind(ctr_,ctr)
    }
    
    if (i == 1) cts_ <- cts__ else cts_ <- rbind(cts_,cts__)
  }
  pars <- list()
  pars$ctr_ <- ctr_ 
  pars$cts_ <- cts_ 
  pars$rws <- rws 
  return(pars)
}

################################################################################
pars_initial <- function (x = x,
                          K = K,
                          alpha = alpha,
                          nstart = 10,
                          niter = 10,
                          ctr_,
                          cts_,
                          rws){
  p <- ncol(x)
  xx <- ctr_
  nrep <- length(rws)
  best_obj <- Inf
  for (start in 1:nstart) {
    ssaa <- sample(nrep, 1)
    nu <- xx[rws[[ssaa]], , drop = FALSE]
    ssaa2 <- sample(nrep, 1)
    if (nrow(nu) < K)  {nu_ <- xx[rws[[ssaa2]], ][1:(K-nrow(nu)), ] + 1e-10; nu <- rbind(nu, nu_)}
    for (iter in 1:niter) {
      dd <- array(NA, c(nrow(xx), K))
      for (k in 1:K) dd[,k] <- apply((xx - matrix(1, nrow = nrow(xx), ncol = 1)%*%nu[k, ])^2, 1, sum, na.rm = TRUE)
      ee <- array(NA, c(nrow(xx), K))
      for (k in 1:K) ee[, k] <- apply(!is.na((xx-matrix(1, nrow = nrow(xx), ncol = 1)%*%nu[k, ] ) ^2 ), 1, sum)
      dd[ee==0] <- Inf
      dist <- apply(dd, 1, min)
      group <- apply(dd, 1, which.min)
      qq <- quantile(dist, 1-alpha)
      group[dist>qq] <- 0
      for (k in 1:K) nu[k,] <- apply(rbind(xx[group==k, ]), 2, mean, na.rm = TRUE)
    }
    for (kk in 1:K) for (jj in 1:p)  if (sum(is.nan(nu[kk, jj]))>0)  nu[kk, jj] <- mean(nu[, jj], na.rm = TRUE)
    obj = sum(dd)
    if (obj < best_obj )  {best_obj = obj; best_init_mean = nu; best_nu_assign = group}
  }
  cts_a <- array(NA, c(p, p, nrow(ctr_)))
  cts_K <- array(NA,c(p, p, K))
  for (ii in 1:nrow(ctr_)) cts_a[,,ii] <- cts_[(p*(ii-1)+1):(ii*p), ]
  for (kk in 1:K) cts_K[, , kk] = apply(cts_a[, , best_nu_assign == kk] , c(1,2), mean, na.rm = TRUE)
  for (kk in 1:K) for (jj in 1:p)  if (sum(is.nan(cts_K[jj, jj, kk]))>0)  cts_K[jj, jj, kk] = mean(cts_K[1:p, 1:p, kk], na.rm = TRUE)
  for (kk in 1:K) for (jj1 in 1:p) for (jj2 in 1:p) if (sum(is.nan(cts_K[jj1, jj2, kk]))>0)  cts_K[jj1, jj2, kk] = runif(1)*(cts_K[jj1, jj1, kk]*cts_K[jj2, jj2, kk])^0.5
  
  best <- list()
  best$init_mean <- best_init_mean
  best$nu_assign <- best_nu_assign
  
  best$init_cov <- cts_K
  
  init_pi_ <- runif(K)
  best$init_pi_ <- init_pi_/sum(init_pi_)
  group <- informal_eval(x, K, best_init_mean)
  best$init_pi_ <- prop.table(table(group))
  return(best)
}

################################################################################
informal_eval <- function(X,
                          K,
                          nu) {
  dd <- array(NA, c(nrow(X), K))
  for (k in 1:K) dd[, k] <- apply((X - matrix(1, nrow = nrow(X), ncol = 1)%*%nu[k, ])^2,1, sum, na.rm = TRUE)
  group <- apply(dd, 1, which.min)
  return(group)
}

################################################################################
## Alternative initialization

initfunc_stepwise <- function(X, 
                              G,
                              m,
                              alpha, 
                              equal.weights) {
  n.obs <- nrow(X)
  p <- ncol(X)
  init <- vector(mode = "list", length = 0)
  invisible(capture.output(Ximp <- cellWise::DI(X)$Ximp))
  tclust.param <- tclust::tclust(Ximp, k = G, alpha = 1 - alpha, restr.fact = 1000)
  init$mu <- matrix(t(tclust.param$centers), nrow = G, ncol = p)
  init$Sigma <- vector(mode = "list", length = G)
  Wg <- array(data = 0, dim = c(n.obs, p, G))
  for (g in 1:G) {
    init$Sigma[[g]] <- restr_eigen(tclust.param$cov[, , g], ceigen = 1e-6)
    out <- cellWise:: cellHandler(X, mu = init$mu[g, ], Sigma = init$Sigma[[g]])
    Wout <- matrix(1, n.obs, p)
    Wout[c(out$indcells, out$indNAs)] <- 0
    Wg[, , g] <- Wout
  }
  nout <- apply(Wg, 3, rowSums)
  init$label <- apply(nout, 1, function(x){
    if (length(unique(x)) != G) {
      return(sample(1:G, 1))
    } else{
      return(which.max(x))
    }
  })
  init$W <- matrix(0, nrow = n.obs, ncol = p)
  for (i in 1:n.obs) {
    init$W[i, ] <- Wg[i, , init$label[i]]
  }
  init$U <- diag(G)
  init$U <- init$U[init$label, , drop = FALSE]
  init$pp <- update_weights(init$U, m, equal.weights) 
  return(init)
}



