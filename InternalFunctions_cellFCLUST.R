## Internal functions for cellGMM ###########

################################################################################
tr <-
  function(x)
    sum(diag(as.matrix(x)))

################################################################################
norm <-
  function(x,
           type = c("standard", "center", "range", "robust"),
           alpha = NULL) {
    type <- match.arg(type,
                      choices = eval(formals(norm)$type),
                      several.ok = FALSE)
    x <- as.matrix(x)
    switch(
      type,
      "standard" = scale(x, center = apply(x, 2, mean, na.rm = TRUE), scale = apply(x, 2, sd, na.rm = TRUE)),
      "center"   = scale(x, center = apply(x, 2, mean, na.rm = TRUE), scale = FALSE),
      "range"    = apply(x, 2, function(x)
        (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))),
      "robust" = {init_MCD <- cellMCD(x, alpha = alpha, checkPars = list(silent = TRUE))
      x <- sweep(x, 2, init_MCD$mu, FUN = "-")
      x <- sweep(x, 2, sqrt(diag(init_MCD$S)), FUN = "/")
      (x - init_MCD$mu) / sqrt(diag(init_MCD$S))}
    )
    return(x)
  }

################################################################################
rand.member <-
  function(n.obs,
           G) {
    if (G > n.obs)
      stop('The number of groups is larger than the number of observations')
    if (G == n.obs) {
      U <- diag(n.obs)
    } else if (G == 1) {
      U <- c(rep(1, n.obs))
    }
    else {
      U <- matrix(0, n.obs, G)
      U[1:G,] = diag(G)
      U[(G + 1):n.obs, 1] <- 1
      for (i in (G + 1):n.obs) {
        U[i,] <- U[i, sample(G)]
      }
      U <- U[sample(n.obs),]
    }
    return(as.matrix(U))
  }

################################################################################
initfunc_new <- function(X, 
                         G,
                         m,
                         maxfact,
                         zero.tol,
                         alpha_tclust, 
                         alpha_1,
                         alpha_2,
                         alpha.A1, 
                         alpha.A2,
                         nrep,
                         nstart,
                         niter) {
  p <- ncol(X)
  W <- w_initial(y = X, alpha_tclust = alpha_tclust, alpha_1 = alpha_1, alpha_2 = alpha_2, K = G, 1)$ww
  x <- X
  x[W==0] <- NA
  prpar <- prepars_initial(x = x, K = G, alpha = alpha.A1, q = floor(p/2) + 1, nrep = nrep, nstart = nstart, niter = niter)
  par <- pars_initial(x = x, K = G, alpha = alpha.A2, nstart = nstart, niter = niter, prpar$ctr_, prpar$cts_, prpar$rws)
  group <- informal_eval(X = x, K = G, par$init_mean)
  Sigma <- purrr::array_tree(par$init_cov, 3)
  if (any(unlist(lapply(lapply(Sigma, function(x){eigen(x)$values}), is.complex)))) {
    refg <- which(unlist(lapply(lapply(Sigma, function(x) {eigen(x)$values}), is.complex)) == TRUE)
    for (g in refg) {
      decSigma <- eigen(Sigma[[g]])
      Sigma[[g]] <- Re(decSigma$vectors)%*%diag(Re(decSigma$values))%*%t(Re(decSigma$vectors))
    }
  }
  par$init_cov_constr <- vector(mode = "list", length = G)
  par$init_cov_constr <- restr_diffax(Sigma, maxfact = maxfact, zero.tol = zero.tol, csize = tabulate(group))
  return(list(mu = par$init_mean,
              Sigma = par$init_cov_constr,
              pp = as.vector(par$init_pi_),
              W = W,
              label = group))
}

################################################################################
impute_mis <-
  function(X,
           difpat,
           pat,
           pp,
           mu,
           Sigma) {
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
    G <- dim(mu)[1]
    Ximp <- array(rep(X, G), dim = c(n.obs, p, G))
    sigma.imp <- array(0, dim = c(p, p, n.obs, G))
    for (t in 1:nrow(difpat)) {
      var.mis <- difpat[t, , drop = FALSE] == 0
      for (g in 1:G) {
        if (any(var.mis)) {
          if (all(var.mis)) {
            Ximp[pat == t, , g] <- mu[g, , drop = FALSE]
            un.t <- which(pat==t)
            for (i in 1:length(un.t)) {
              sigma.imp[, , un.t[i], g] <- Sigma[[g]]
            }
          } else {
            var.obs <- !var.mis
            mu_m <- mu[g, var.mis, drop = FALSE]
            mu_o <- mu[g, var.obs, drop = FALSE]
            sigma_mo <- Sigma[[g]][var.mis, var.obs, drop = FALSE]
            sigma_om <- Sigma[[g]][var.obs, var.mis, drop = FALSE]
            sigma_mm <- Sigma[[g]][var.mis, var.mis, drop = FALSE]
            sigma_oo_inv <- solve(Sigma[[g]][var.obs, var.obs], tol = .Machine$double.xmin)
            Ximp[pat == t, var.mis, g] <- 
              sweep(t(sigma_mo %*% sigma_oo_inv %*% t(sweep(X[pat == t, var.obs, drop = FALSE], 2, mu_o, "-", check.margin = FALSE))), 2, mu_m, "+")
            sigma.mis <- sigma_mm - sigma_mo %*% sigma_oo_inv %*% sigma_om
            sigma.imp[var.mis, var.mis, pat==t, g] <- sigma.mis
          }
        }
      }
    }
    return(list(Ximp = Ximp,
                sigma.imp = sigma.imp))
  }

################################################################################
impute_mis_best <-
  function(X,
           Ximp,
           W,
           label) {
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
    Ximp.final <- X
    for (i in which(rowSums(W) != p)) {
      var.mis <- W[i, , drop = FALSE] == 0
      Ximp.final[i, var.mis] <- Ximp[i, var.mis, label[i], drop = FALSE]
    }
    return(Ximp.final)
  }

################################################################################
update_W <- 
  function(X,
           G,
           W,
           U,
           m,
           pp,
           mu,
           Sigma,
           alpha) {
    n.obs <- dim(W)[1]
    p <- dim(W)[2]
    ordering <- 1:p # Alternatively, for instance, sample(1:p, p)
    for (j in 1:p) {
      unit.obs <- which(!is.na(X[, ordering[j]]))
      Delta <- rep(0, length(unit.obs))
      for (i in 1:length(unit.obs)) {
        var.obs <- which(W[unit.obs[i], ] == 1)
        var.obs <- var.obs[var.obs != ordering[j]]
        xC.hat <- matrix(as.double(NA), G, 2)
        if (length(var.obs) == 0) {
          for (g in 1:G) {
            xC.hat[g, 1] <- mu[g, ordering[j]]
            xC.hat[g, 2] <- Sigma[[g]][ordering[j], ordering[j]]
            if (xC.hat[g, 2] < .Machine$double.xmin) {
              xC.hat[g, 2] <- (.Machine$double.xmin)
            }
          }
        } else {
          for (g in 1:G) {
            smo_sooinv = Sigma[[g]][ordering[j], var.obs, drop = FALSE] %*% 
              .Internal(La_chol2inv(.Internal(La_chol(Sigma[[g]][var.obs, var.obs,  drop = FALSE], pivot = FALSE, tol = .Machine$double.xmin)), length(var.obs)))
            
            xC.hat[g, 1] <- mu[g, ordering[j]] + smo_sooinv%*%(X[unit.obs[i], var.obs] - mu[g, var.obs])
            xC.hat[g, 2] <- Sigma[[g]][ordering[j], ordering[j]] - smo_sooinv%*%Sigma[[g]][var.obs, ordering[j], drop = FALSE]
            if (xC.hat[g, 2] < .Machine$double.xmin) {
              xC.hat[g, 2] <- (.Machine$double.xmin)
            }
          }
        }
        Delta[i] <- (-0.5 * sum((U[unit.obs[i], ]^m)*(log(2*pi) + log(xC.hat[, 2]) + ((X[unit.obs[i], ordering[j]] - xC.hat[, 1])^2)/xC.hat[, 2])))
      }
      # cutoff <- sort(Delta, decreasing = TRUE)[ceiling(alpha * n.obs)]
      cutoff <- sort(Delta, decreasing = TRUE)[ceiling(alpha * length(unit.obs))]
      clean <- unit.obs[Delta >= cutoff]
      out <- unit.obs[Delta < cutoff]
      # clean <- which(Delta >= cutoff) # If missing are considered in the alpha %
      # out <- which(Delta < cutoff) # If missing are considered in the alpha %
      W[clean, ordering[j]] <- 1
      W[out, ordering[j]] <- 0
    }
    return(W)
  }

################################################################################
get_alpha <- function(Delta){
  sorted_delta <- sort(Delta, decreasing = FALSE)
  n.obs <- length(sorted_delta)
  slope = (max(sorted_delta) - min(sorted_delta))/(1 - 1/n.obs)
  intercept = max(sorted_delta) - slope
  dp = abs(sorted_delta - (slope * (1:n.obs/n.obs) + intercept))/sqrt(1 + slope^2)
  return(which(dp == max(dp)))
}

################################################################################
update_U_new <-
  function(X,
           m,
           difpat,
           pat,
           pp,
           mu,
           Sigma) {
    G <- dim(mu)[1]
    ll <- matrix (NA, nrow(X), G) ## to hold the density for each point in the sample and w.r.to each component, n rows, K columns
    for (tt in 1:nrow(difpat)) {  
      var.obs <- difpat[tt, , drop = FALSE] == 1
      if (any(var.obs)) {
        for (g in 1:G) {	
          mu_o <- mu[g, var.obs, drop = FALSE]
          sigma_oo <- Sigma[[g]][var.obs, var.obs, drop = FALSE]
          ll[pat == tt, g] <- pp[g]*exp(dmnorm2(X[pat == tt, var.obs, drop = FALSE], drop(mu_o), sigma_oo))  
        }
      } else if (!any(var.obs)) {
        ll[pat == tt, ] <-  t(replicate(length(which(pat == tt)), pp))
      }
    }
    if (m > 1) {
      ll <- ll*(ll>=0) # Retain only non negative values of the densities, evaluated at each point i (row) with reference to group g (col)
      assig <- apply(ll, 1, which.max)
      U <- matrix (0, nrow = nrow(X), ncol = G)		
      max.ll <- apply(ll, 1, max)
      max.ll.greater1 <- max.ll >= 1
      iid <- cbind(which(max.ll.greater1), assig[max.ll.greater1])
      U[iid] <- 1
      if (any(!max.ll.greater1)) {
        log.ll <- ll * 0
        log.ll[ll > 0] <- -log(ll[ll > 0]) * (ll[ll > 0] < 1)
        ll2 <- log.ll %*% matrix(rep(diag(G), G), nrow = G)
        yy <- matrix (0, nrow = G, ncol = G * G)
        for (ii in 1:G) {
          yy[ii, ((ii-1)*G+1):(ii*G)] = rep(1, G)
        }
        ll3 <- log.ll %*% yy
        U2_ij__ <- (((ll3 * (ll2 > 0) / (ll2 + (ll2 == 0)))^(1/(m-1))) %*% t(yy))
        U2 <- 1*(U2_ij__  > 0) / (U2_ij__ + (U2_ij__ == 0))
        U[!max.ll.greater1, ] <- U2[!max.ll.greater1, ] 
      }
    } else if (m == 1) {
      assig <- apply(ll, 1, which.max)
      U <- diag(G)
      U <- U[assig, ]
    }
    U[which(U < sqrt(.Machine$double.eps), arr.ind = TRUE)] <- sqrt(.Machine$double.eps)
    return (U)
  }

################################################################################
update_weights <-
  function(U, m, equal.weights) {
    if (equal.weights) {
      pp <- rep(1/ncol(U), ncol(U))
    } else {
      pp <- colSums(U^m) / sum(U^m)
    }
    return(pp)
  }

################################################################################
update_param_fuzzy <-
  function(Ximp,
           U,
           m,
           sigma.imp,
           maxfact,
           zero.tol,
           label) {
    n.obs <- dim(Ximp)[1]
    p <- dim(Ximp)[2]
    G <- dim(U)[2]
    mu.new <- matrix(0, G, p)
    Sigma.new <- list()
    for (g in 1:G){
      mu.new[g, ] <- sweep(t(U[, g]^m) %*% Ximp[, , g], 2, 1 / sum(U[, g]^m), "*")
      Xo <- array(0, dim = c(p, p, n.obs))
      for (i in 1:n.obs) {
        Xo[, , i] <- (U[i, g]^m)*(tcrossprod((Ximp[i, , g] - mu.new[g, ])) + sigma.imp[, , i, g])
      }
      Xo <- apply(Xo, 1:2, sum)
      Sigma.temp <-  Xo / sum(U[, g]^m)
      Sigma.new[[g]] <- Sigma.temp
    }
    constr <- restr_diffax(Sigma.new, maxfact = maxfact, zero.tol = zero.tol, csize = colSums(U^m))
    Sigma.new <- constr
    return(list(mu = mu.new,
                Sigma = Sigma.new))
  }

################################################################################
restr_diffax <- 
  function (Sigma, 
            maxfact, 
            zero.tol,
            csize) {  
    G <- length(Sigma)
    p <- dim(Sigma[[1]])[1]
    u <- array (NA, dim = c(p, p, G))
    d <- array (NA, dim = c(p, G))
    for (g in 1:G) {
      ev <- eigen(Sigma[[g]])
      u[, , g] <- ev$vectors
      d[, g] <- ev$values
    }
    d[d < 0] <- 0
    d <- restr2_eigenv(d, csize, maxfact, zero.tol)
    if (!(max(d) > zero.tol)) {
      return(Sigma)
    }
    for (g in 1:G) {
      Sigma[[g]] <- u[, , g] %*% diag(d[, g], nrow = p) %*% t(u[, , g])
    }
    return(Sigma)
  }

################################################################################
restr2_eigenv <- 
  function(autovalues,
           csize,
           maxfact,
           zero.tol) {
    c <- maxfact
    d <- t(autovalues)
    p <- dim(autovalues)[1]
    G <- dim(autovalues)[2]
    n.obs <- sum(csize)
    nis <- matrix(data = csize, nrow = G, ncol = p)
    d_ <- sort(c(d, d/c))
    dim <- length(d_)
    d_1 <- d_
    d_1[dim + 1] <- d_[dim]*2
    d_2 <- c(0, d_)
    ed <- (d_1 + d_2)/2
    dim <- dim + 1;
    if ((max(d[nis>0]) <= zero.tol))
      return(matrix(0, nrow = p, ncol = G))
    if (abs(max(d[nis>0])/min(d[nis>0]))<=c) {
      d[nis == 0] <- mean(d[nis>0])
      return(t(d))
    }
    t <- s <- r <- array(0, dim = c(G, dim))
    sol <- sal <- array(0, dim = c(dim))
    for (mp_ in 1:dim){
      for (i in 1:G) {
        r[i, mp_] <- sum((d[i, ]<ed[mp_])) + sum((d[i, ]>ed[mp_]*c))
        s[i, mp_] <- sum(d[i, ]*(d[i, ]<ed[mp_]))
        t[i, mp_] <- sum(d[i, ]*(d[i, ]>ed[mp_]*c))
      }
      sol[mp_] <- sum(csize/n.obs*(s[, mp_] + t[, mp_]/c))/(sum(csize/n.obs*(r[, mp_])))
      e <- sol[mp_]*(d<sol[mp_]) + d*(d>=sol[mp_])*(d<=c*sol[mp_]) + (c*sol[mp_])*(d>c*sol[mp_])
      o <- -1/2*nis/n.obs*(log(e) + d/e)
      sal[mp_] <- sum(o)
    }
    eo <- which.max(c(sal))
    m <- sol[eo]
    t(m*(d<m) + d*(d>=m)*(d<=c*m) + (c*m)*(d>c*m))
  }

################################################################################
restr_eigen <- 
  function(Sigma,
           ceigen) {
    if (min(eigen(Sigma, symmetric = TRUE)$values) < ceigen) {
      Sigma <- Sigma + ceigen * diag(dim(Sigma)[1])
    }
    return(Sigma)
  }

################################################################################
dmnorm2 <- 
  function(X, 
           mu,
           sigma){
    return(log(((2 * pi)^(- length(mu) / 2))) + log(det(sigma)^(-1/2)) -0.5 * mahalanobis(X, mu, sigma))
  }

################################################################################
obj_cellFCLUST <-
  function(X,
           difpat,
           pat,
           U,
           m,
           pp,
           mu,
           Sigma) {
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
    G <- length(pp)
    lf <- matrix(as.double(NA), n.obs, G)
    for (t in 1:nrow(difpat)){
      var.obs <- difpat[t, , drop = FALSE] == 1
      if (any(var.obs)) {
        for (g in 1:G) {
          mu_o <- mu[g, var.obs, drop = FALSE]
          sigma_oo <- Sigma[[g]][var.obs, var.obs, drop = FALSE]
          lf[pat == t, g] <- (U[pat == t, g]^m)*(dmnorm2(X[pat == t, var.obs, drop = FALSE], drop(mu_o), sigma_oo) + log(pp[g]))
        }
      } else {
        lf[pat == t, ] <- (U[pat == t, ]^m)*log(pp)
      }
    }
    loglik <- sum(lf, na.rm = TRUE) 
    return(loglik)
  }

