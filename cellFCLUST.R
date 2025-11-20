source('InternalFunctions_cellFCLUST.R')
source('InitializationFunctions_cellFCLUST.R')

###################################################################################################################################################
# Name: cellFCLUST
# Fuzzy clustering with cellwise outlier detection with constraints and allowing missing data
#
## Description
#
## Arguments:
## X:                  (n.obs x p) numeric data matrix or data frame allowing missing data
## G:                  Number of the clusters (numeric)
## m:                  Fuzzifier parameter (numeric)
## tuning_param_init:: A data frame of tuning parameters for the initialization: alpha_tclust, alpha_1, alpha_2, alpha.A1, alpha.A2, nrep, nstart, niter (numeric)
## alpha:              Fraction of cells per variable that must remain unflagged (default: 0.75, which represents a lower bound; numeric)
## quant:              Tuning constant to flag cells (default: 0.99; numeric)
## equal.weights:      Value indicating if the weights must be equal (TRUE) or not (FALSE) (default: FALSE; logical)
## maxfact:            Constant value greater than 1. Larger values imply larger differences of component covariance matrices; a value of 1 specifies the strongest restriction (default: 30; numeric)
## zero.tol:           Tolerance for the eigenvalues (default: 1e-16; numeric)
## normalization:      Type of normalization for the data: NULL, "standard", "center", "range" "SVD" (default: NULL; string)
## maxiter:            Maximum number of iterations (default: 500; numeric)
## tol:                Tolerance value for convergence (default: 1e-6; numeric)
## rndstart:           Number of random starts (default: 20; numeric)
## init:               Type of initialization: "cellGMM", "stepwise" (default: cellgmm; character)
## manual_initparam:   A list of initial parameters (W, pp, mu, Sigma, U, label) (default: NULL, i.e., it is ignored)
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
cellFCLUST <- function(X,
                       G,
                       m,
                       tuning_param_init,
                       alpha = 0.75,
                       quant = 0.99,
                       equal.weights = FALSE,
                       maxfact = 30,
                       zero.tol = 1e-16, 
                       normalization = NULL,
                       maxiter = 500,
                       tol = 1e-6,
                       rndstart = 20,
                       init_flag = "cellGMM",
                       manual_initparam = NULL,
                       showprogress = TRUE) {
  call <- mget(names(formals())[-c(1, 16)],sys.frame(sys.nframe()))
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
  for (loop in 1:rndstart) {
    count <- 0
    conv <- 1
    wp <- 0
    # INITIALIZATION 
    if (is.null(manual_initparam)) {
      if (init_flag == "cellGMM") {
        init <- initfunc_new(X = X, G = G, m = m, maxfact = maxfact, zero.tol = zero.tol, 
                             alpha_tclust = tuning_param_init$alpha_tclust, 
                             alpha_1 = tuning_param_init$alpha_1,
                             alpha_2 = tuning_param_init$alpha_2,
                             alpha.A1 = tuning_param_init$alpha.A1, 
                             alpha.A2 = tuning_param_init$alpha.A2,
                             nrep = tuning_param_init$nrep,
                             nstart = tuning_param_init$nstart,
                             niter = tuning_param_init$niter)
      } else if (init_flag == "stepwise") {
        init <- initfunc_stepwise(X = X, G = G, m = m, alpha = alpha, equal.weights = equal.weights)
      } 
    } else {
      init <- manual_initparam
    }
    model <- init[-which(names(init) == "label")]
    label <- init$label
    dif.pattern <- unique(model$W)
    pat.unit <- match(data.frame(t(model$W)), data.frame(t(dif.pattern))) 
    model$U <- update_U_new(X, m, dif.pattern, pat.unit, model$pp, model$mu, model$Sigma)
    if (showprogress) {
      pb <- txtProgressBar(min = 0, max = rndstart, style = 3, 
                           width = 75, char = "=")
    }
    # ITERATIONS  
    while (conv > tol && count <= maxiter) {
      if (count == 0) {
        if (sum(unique(unlist(label))) != sum(1:G)) {
          #which(! sort(unique(unlist(label))) %in% 1:G)
          wp <- 1
          obj <- -.Machine$double.xmax
          break
        }
      }
      count <- count + 1 
      ## E-STEP 
      # UPDATE W
      if (alpha < 1) {
        model$W <- update_W(X = X, G = G, W = model$W, U = model$U, m = m, pp = model$pp, mu= model$mu, Sigma = model$Sigma, alpha = alpha)
        dif.pattern <- unique(model$W)
        pat.unit <- match(data.frame(t(model$W)), data.frame(t(dif.pattern)))  
      }
      # UPDATE U AND X[W^{c}]   
      model$U <- update_U_new(X, m, dif.pattern, pat.unit, model$pp, model$mu, model$Sigma)
      label <- apply(model$U, 1, which.max)
      if (sum(unique(unlist(label))) != sum(1:G)) {
        wp <- 1
        obj.c <- -.Machine$double.xmax
        break
      }
      val.impute <- impute_mis(X, dif.pattern, pat.unit, model$pp, model$mu, model$Sigma)
      model$X.imputed <- val.impute$Ximp
      # M-STEP   
      model$pp <- update_weights(model$U, m, equal.weights)
      comp.param <- update_param_fuzzy(Ximp = val.impute$Ximp, U = model$U, m = m, sigma.imp = val.impute$sigma.imp, maxfact = maxfact, zero.tol = zero.tol, label = label)
      model$mu <- comp.param$mu
      model$Sigma <- comp.param$Sigma
      # OBJECTIVE FUNCTION 
      obj.c <- obj_cellFCLUST(X, dif.pattern, pat.unit, model$U, m, model$pp, model$mu, model$Sigma) 
      if (count > 1) {
        vobj <- c(vobj, obj.c)
        conv <- obj.c - obj
      } else if (count == 1) {
        vobj <- obj.c
      }
      obj <- obj.c
    }
    if (loop == 1) {
      W.best <- model$W
      label.best <- label
      model.best <-  model[-which(names(model) == "W")]
      obj.best <- vobj
      obj.final.best <- obj.c
      init.best <- init
      loop.best <- loop
      iter.best <- count
      if (wp == 0) {
        X.imputed.best <- impute_mis_best(X = X, Ximp = model.best$X.imputed, W = W.best, label = label.best)
      } else if (wp == 1) {
        X.imputed.best <- NULL
      }
    } else {
      if (obj.c > obj.final.best) {
        W.best <- model$W
        label.best <- label
        model.best <- model[-which(names(model) == "W")]
        obj.best <- vobj
        obj.final.best <- obj.c
        init.best <- init
        loop.best <- loop
        iter.best <- count 
        if (wp == 0) {
          X.imputed.best <- impute_mis_best(X = X, Ximp = model.best$X.imputed, W = W.best, label = label.best)
        } else if (wp == 1) {
          X.imputed.best <- NULL
        }
      }
    }
    if (G == 1) {
      break
    }
    if (showprogress) {
      setTxtProgressBar(pb, loop)
      cat("Loop", loop, "/", rndstart)
    }
  } # END LOOP
  if (showprogress) {
    close(pb)
  }
  if (obj.final.best == -.Machine$double.xmax) {
    print("The solution has a number of cluster < G and the objective is not computed.")
  }
  return(
    list(
      call = call,
      X = X,
      X.imputed = model.best$X.imputed,
      X.imputed.best = X.imputed.best,
      W = W.best,
      label = label.best,
      U = model.best$U,
      pp = model.best$pp,
      mu = model.best$mu,
      sigma = model.best$Sigma,
      obj = obj.best,
      obj.final = obj.final.best,
      init.param = init.best,
      loop = loop.best,
      iter = iter.best
    )
  )
} # END FUNCTION