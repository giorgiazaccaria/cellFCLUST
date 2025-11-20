# Name: Fuzzy Trimmed C Prototypes (FTPC) 
# Competitor for cellFCLUST
#
## Description
#
## Arguments:
## X:                  (n.obs x p) numeric data matrix or data frame 
## G:                  Number of the clusters (numeric)
## m:                  Fuzzifier parameter (numeric)
## alpha:              Trimming ratio (default: 0.50, numeric)
## normalization:      Type of normalization for the data: NULL, "standard", "center", "range" "SVD" (default: NULL; string)
## maxiter:            Maximum number of iterations (default: 500; numeric)
## tol:                Tolerance value for convergence (default: 1e-6; numeric)
## rndstart:           Number of random starts (default: 20; numeric)
## showprogress:       Progress of the code (default: TRUE; logical)

## Values:
## call:              Matched call.
## X:                 Input data matrix
## X.imputed:         Imputed data matrix for each component
## X.imputed.best:    Imputed data matrix considering the unit-component assignment
## W:                 (n.obs x p) matrix with zeros corresponding to missing and/or outlying values
## label:             (n.obs x 1) vector of cluster membership, based on the Maximum A Posteriori
## U                  (n.obs x G) membership matrix             
## pp:                (1 x G) vector of weights
## mu:                (G x p) matrix of the component mean vectors
## sigma:             List of G (p x p) component covariance matrices
## obj:               Vector of the objective computed for each iteration of the best loop
## obj.final:         Final value of the objective function
## init.param:        A list of initial parameters (W, pp, mu, Sigma, U, label) or NULL
## loop:              Loop corresponding to the best model
## iter:              Actual number of iterations needed to reach convergence
##################################################################################################################################################
Ftkmeans <- function(X, 
                     G, 
                     m, 
                     alpha = 0.50,
                     normalization = NULL,
                     maxiter = 500,
                     tol = 1e-6,
                     rndstart = 20,
                     showprogress = TRUE) {
  call <- mget(names(formals())[-c(1, 9)],sys.frame(sys.nframe()))
  # PRE-PROCESSING 
  if (!is.null(normalization)) {
    X <- norm(X, normalization, alpha)
  } else {
    X <- as.matrix(X)
  }
  n.obs <- dim(X)[1]
  p <- dim(X)[2]
  if (G < 1 || G > n.obs || G%%1 != 0) {
    stop("G is not properly fixed: G must be an integer chosen into [1, dim(X)[1]]", call. = TRUE)
  }
  # STARTING
  obj = Inf
  for (loop in 1:rndstart) {
    vobj = c()
    count <- 0
    conv <- 1
    # INITIALIZATION (as well as tclust)
    X_init <- X[sample(1:n.obs, G*(p + 1)), ]
    label_init <- rep(1:G, each = p + 1)
    centers <- as.matrix(aggregate(X_init, by = list(label = label_init), FUN = mean)[-1], G, p)
    covariance <- lapply(1:G, function(g) {
      Xg <- as.matrix(X[label_init == g, ])
      crossprod(scale(Xg, scale = FALSE)) / nrow(Xg)})
    U <- matrix(runif(n.obs*G, 0, 1), nrow = n.obs, ncol = G)  
    U <- U/apply(U, 1, sum)
    
    model = list()
    model$init = list()
    model$init$centers = centers
    model$init$covariance = covariance
    model$init$U = U
    model$init$X_init
    model$init$label_init

    # Iterations
    while (count < maxiter & conv > tol){
      
      # Calculate distances
      d2 <- matrix(NA, n.obs, G)
      D2 <- matrix(NA, n.obs, G)
      a <- matrix(NA, 1, G)
      for (i in 1:G){
        d2[, i] <- rowSums(sweep(X, 2, centers[i, ], "-")^2)
        #D2[, i] <- d2[, i] - rowSums((sweep(X, 2, centers[i, ], "-") %*% eigen(covariance[[i]])$vector)^2)
        #D2[, i] <- d2[, i] - (sweep(X, 2, centers[i, ], "-") %*% eigen(covariance[[i]])$vector[, 1])^2
        #a[1, i] = 1 - min(eigen(covariance[[i]])$values)/max(eigen(covariance[[i]])$values)
      }
      #Z2 = sweep(D2, 2, a, "*") + sweep(d2, 2, (1-a), "*")
      Z2 = d2
      # Compute the loss
      h2 = rowSums(Z2^(1/(1-m)))^(1-m)
      # Trim
      cutoff = sort(h2, decreasing = T)[max(ceiling(alpha*n.obs), 1)]
      t.idx = h2 > cutoff
      rel.h2 = h2[!t.idx]
      # Update U
      U[!t.idx, ] <- Z2[!t.idx, ]^(1/(1-m))/rowSums(Z2[!t.idx, ]^(1/(1-m)))
      U[t.idx, ] <- 0
      
      # Compute centers
      old_centers <- centers
      centers <- crossprod(U^m, X)/colSums(U^m)
      for (i in 1:G){
        covariance[[i]] <- cov.wt(X, U[, i]^(m), center = centers[i, ], method = "ML")$cov
      }
      
      conv = sum(abs(centers-old_centers))
      
      count = count + 1
      
      vobj = c(vobj, sum(rel.h2))
      #vobj = c(vobj, conv)
      if (conv < tol | count >= maxiter){
        model$centers <- centers
        model$U <- U
        model$label <- rep(0, n.obs)
        model$label[!t.idx] <- apply(model$U[!t.idx, ,drop=F], 1, which.max)
        model$covariance = covariance
        model$call = call
        model$X = X
        model$W = matrix(1, n.obs, p)
        model$W[t.idx, ] <- 0
        model$pp <- colSums(U^m)/sum(U^m)
        model$loop = loop
        model$iter = count
        model$obj = vobj
        model$obj.final = sum(rel.h2)
      }
    }
    if (sum(rel.h2) < obj){
      obj <- sum(rel.h2)
      best.model <- model
    }
  }
  return(best.model)
}


