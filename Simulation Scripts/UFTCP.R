# Name: Unsupervised Fuzzy Trimmed C Prototypes (FTPC)
UFtkmeans <- function(d_alpha,
                      X, 
                      G, 
                      m, 
                      max_alpha = 0.50,
                      normalization = NULL,
                      maxiter = 500,
                      tol = 1e-6,
                      rndstart = 20,
                      showprogress = TRUE){
  source("FTCP.R")
  
  if (max_alpha==0){
    return(Ftkmeans(X, G, m, 0, normalization, maxiter, tol, rndstart))
  }
  VDpa <- vector("double", length(seq(0, max_alpha, d_alpha)))
  model <- vector("list", length(seq(0, max_alpha, d_alpha)))
  best.model <- list()
  k = 1
  for(alpha in seq(0, max_alpha, d_alpha)){
    model[[k]] <- Ftkmeans(X, G, m, alpha, normalization, maxiter, tol, rndstart)
    centers = model[[k]]$centers
    U = model[[k]]$U
    Dpa = 0
    for (i in 1:G){
      Fj = cov.wt(X, U[, i]^(m), center = centers[i, ], method = "ML")$cov
      Dpa = Dpa + sum(U[, i])/det(Fj)^(0.5)
    }
    
    VDpa[k] <- Dpa
    k = k+1
  }
  a_star = which.max(diff(VDpa)) + 1
  best.model <- model[[a_star]]
  best.model$alpha <- seq(0, max_alpha, d_alpha)[a_star]
  return(best.model)
}