solve_core_given_factors <- function(X, U) {
  d <- length(U)
  U_pinv <- lapply(1:d, function(k) {
    UtU    <- crossprod(U[[k]])
    lambda <- max(eigen(UtU, symmetric = TRUE, 
                        only.values = TRUE)$values) * 1e-10
    solve(UtU + diag(lambda, nrow(UtU))) %*% t(U[[k]])
  })
  ttl(X, U_pinv, ms = 1:d)
}

compute_objective_with_normalization <- function(X, G, U) {
  result <- list()
  result$core <- G
  result$factors <- U
  result_normalized <- normalize_factors(result)
  X_recon <- reconstruct_unimodal(result_normalized)
  X_recon@data <- X_recon@data / sum(X_recon@data)
  abs_error <- fnorm(X - X_recon) 
  # X_recon <- ttl(G, U, ms = 1:length(U))
  # abs_error <- fnorm(X - X_recon)
  return(abs_error)
}

compute_objective <- function(X, G, U) {
  X_recon <- ttl(G, U, ms = 1:length(U))
  return(fnorm(X - X_recon))
}

# This uses solve_core_given_factors where we update each U_k using only the gradient direction (no solve) and
# after all factors are updated, recompute G exactly via pseudoinverses
update_mode_with_descent <- function(X, G, U, mode_idx, current_obj) {
  
  k <- mode_idx
  d <- length(U)
  
  X_unfold <- k_unfold(X, k)@data
  G_unfold <- k_unfold(G, k)@data
  
  other_modes <- setdiff(d:1, k)
  if (length(other_modes) > 0) {
    kron_prod <- U[[other_modes[1]]]
    if (length(other_modes) > 1)
      for (m in other_modes[-1])
        kron_prod <- kronecker(kron_prod, U[[m]])
  } else {
    kron_prod <- matrix(1, 1, 1)
  }
  
  V       <- tcrossprod(G_unfold, kron_prod)
  VVt     <- tcrossprod(V)
  VVt_reg <- VVt + diag(1e-8 * mean(diag(VVt)), nrow(VVt))
  
  # Proper unconstrained least squares solution at current point
  # This is a genuine descent direction for the unconstrained problem
  U_unconstrained <- tcrossprod(X_unfold, V) %*% solve(VVt_reg)
  
  # Project each column onto unimodal constraint
  U_projected <- U_unconstrained
  for (j in 1:ncol(U_projected)) {
    U_projected[, j] <- best_unimodal(U_unconstrained[, j])
    # if (sum(abs(U_projected[, j])) < 1e-10){
    #   cat("reinitializing ... \n")
    #   U_projected[, j] <- best_unimodal(rnorm(nrow(U_projected)))
    # }
  }
  
  U_new      <- U
  U_new[[k]] <- U_projected
  
  # Core always computed exactly via pseudoinverses — this is where
  # the ill-conditioning fix lives, not in the factor update
  G_new   <- solve_core_given_factors(X, U_new)
  obj_new <- compute_objective(X, G_new, U_new)
  
  if (obj_new <= current_obj)
    return(list(U = U_new, G = G_new, obj = obj_new, accepted = TRUE))
  
  return(list(U = U, G = G, obj = current_obj, accepted = FALSE))
}

# This is the main algorithm implementing unimodal Tucker decomposition using ALS
tucker_unimodal_als <- function(X, ranks, 
         max_iter = 500, 
         tol = 1e-4, 
         random_init = FALSE,
         verbose = TRUE) 
{
  # Convert to tensor if needed
  if (inherits(X, "table")) {
    X <- as.tensor(unclass(X))
  } else if (!inherits(X, "Tensor")) {
    X <- as.tensor(X)
  }
  
  d <- X@num_modes
  dims <- X@modes
  
  if (length(ranks) != d) {
    stop("Length of ranks must equal number of modes")
  }
  
  if (verbose) {
    cat("Tucker ALS with Unimodal Constraints - GUARANTEED DESCENT\n")
    cat("Tensor dimensions:", dims, "\n")
    cat("Ranks:", ranks, "\n\n")
  }
  
  # HOSVD initialization
  U <- list()
  for (k in 1:d) {
    X_unfold <- k_unfold(X, k)@data
    max_rank <- min(ranks[k], nrow(X_unfold), ncol(X_unfold))
    if (random_init){
      if (verbose) cat("Initializing randomly ...\n")
      U_init <- matrix(rnorm(dims[k] * ranks[k]), nrow = dims[k], ncol = ranks[k]) 
    }else{
      if (verbose) cat("Initializing with HOSVD...\n")
      svd_result <- svd(X_unfold, nu = max_rank, nv = 0)
      U_init <- svd_result$u[, 1:ranks[k], drop = FALSE]
    }
    # Project to unimodal
    for (j in 1:ranks[k]) {
      U_init[, j] <- best_unimodal(U_init[, j])
      if (sum(abs(U_init[, j])) < 1e-10) {
        U_init[, j] <- best_unimodal(rnorm(dims[k]))
      }
    }
    U[[k]] <- U_init
  }
  
  # Initialize core
  # the original version
  # G <- ttl(X, lapply(1:d, function(k) t(U[[k]])), ms = 1:d)
  # the new version for pseudoinverse
  G <- solve_core_given_factors(X, U)
  
  
  # Compute initial error
  current_obj <- compute_objective(X, G, U)
  
  if (verbose) {
    cat(sprintf("Initial error: %.6f (%.2f%%)\n", 
                current_obj, 100 * current_obj / fnorm(X)))
    cat("Starting ALS iterations...\n")
  }
  
  # Track history
  history <- list(
    reconstruction_error = numeric(max_iter),
    mode_accepted = matrix(FALSE, max_iter, d)
  )
  
  # Initialise diagnostics storage
  diagnostics <- list(
    cond_numbers   = vector("list", d),  # one vector per mode
    max_cosine_sim = vector("list", d),
    min_entropy    = vector("list", d)
  )
  for (k in 1:d) {
    diagnostics$cond_numbers[[k]]   <- numeric(0)
    diagnostics$max_cosine_sim[[k]] <- numeric(0)
    diagnostics$min_entropy[[k]]    <- numeric(0)
  }
  
  # Main ALS loop
  for (iter in 1:max_iter) {
    obj_start_iter <- current_obj
    
    # Update each mode
    for (k in 1:d) {
      result <- update_mode_with_descent(X, G, U, k, current_obj)
      # result <- update_mode_with_descent_diagnostic(X, G, U, k, current_obj, diagnostics)
      if (result$accepted) {
        U <- result$U
        G <- result$G
        current_obj <- result$obj
        # diagnostics <- result$diagnostics
        history$mode_accepted[iter, k] <- TRUE
      } else {
        history$mode_accepted[iter, k] <- FALSE
      }
    }
    
    # Store error
    history$reconstruction_error[iter] <- current_obj
    relative_error <- current_obj / fnorm(X)
    
    # Compute change
    obj_change <- abs(current_obj - obj_start_iter) / (obj_start_iter + 1e-10)
    
    if (verbose && (iter <= 10 || iter %% 10 == 0)) {
      n_accepted <- sum(history$mode_accepted[iter, ])
      cat(sprintf("Iter %4d: Error = %.6f (%.2f%%), Change = %.2e, Accepted = %d/%d\n",
                  iter, current_obj, 100 * relative_error, obj_change, n_accepted, d))
    }
    
    # Check convergence
    converged <- FALSE
    
    # Criterion 1: Objective change is tiny
    if (obj_change < tol) {
      converged <- TRUE
      convergence_reason <- "objective stagnation"
    }
    
    # Criterion 2: No modes accepted (stuck)
    if (iter > 5 && all(history$mode_accepted[iter, ] == FALSE)) {
      converged <- TRUE
      convergence_reason <- "no modes can improve"
    }
    
    # Criterion 3: Very small error
    if (relative_error < tol) {
      converged <- TRUE
      convergence_reason <- "target error achieved"
    }
    
    if (converged) {
      if (verbose) {
        cat(sprintf("\nConverged at iteration %d (%s)\n", iter, convergence_reason))
        cat(sprintf("Final error: %.6f (%.2f%%)\n", 
                    current_obj, 100 * relative_error))
      }
      history$reconstruction_error <- history$reconstruction_error[1:iter]
      history$mode_accepted <- history$mode_accepted[1:iter, , drop = FALSE]
      break
    }
  }
  
  if (iter == max_iter && verbose) {
    cat(sprintf("\nReached max iterations (%d)\n", max_iter))
    cat(sprintf("Final error: %.6f (%.2f%%)\n", current_obj, 100 * current_obj / fnorm(X)))
  }
  
  # Trim history
  if (iter < max_iter) {
    history$reconstruction_error <- history$reconstruction_error[1:iter]
    history$mode_accepted <- history$mode_accepted[1:iter, , drop = FALSE]
  }
  
  result <- list(
    core = G,
    factors = U,
    history = history,
    ranks = ranks,
    diagnostics = diagnostics,
    converged = converged
  )
  
  class(result) <- "tucker_unimodal"
  return(result)
}

# Unimodal projection function uses Bro's (1998) algorithm for unimodal projection with improved convergence
project_to_unimodal_fast <- function(v, sign="+") {
  n <- length(v)
  if (n <= 2) {
    if (sign == "+"){ 
      return(pmax(v, 0))
    }else{
      return(pmin(v, 0))
    }
  }else{
    cn <- rep(0,n)
    best_error <- sum((v - cn)^2)
    best_v <- cn
    for (p in 1:n){
      if (p==1){ 
        if (sign == "+"){ 
          p2 <- pmax(rev(isoreg(rev(v[(p+1):n]))$yf), 0)
          cn <- c(max(v[p],0), p2) 
        }else{
          p2 <- pmin(isoreg(v[(p+1):n])$yf, 0)
          cn <- c(min(v[p],0), p2) 
        }
      }else{
        if (p==n){
          if (sign == "+"){ 
            p1 <- pmax(isoreg(v[1:(p-1)])$yf, 0)
            cn <- c(p1, max(v[p],0))
          }else{
            p1 <- pmin(rev(isoreg(rev(v[1:(p-1)]))$yf), 0)
            cn <- c(p1, min(v[p],0))
          }
        }else{
          if (sign == "+"){ 
            p1 <- pmax(isoreg(v[1:(p-1)])$yf, 0)
            p2 <- pmax(rev(isoreg(rev(v[(p+1):n]))$yf), 0)
            cn <- c(p1, max(v[p],0), p2)
          }else{
            p1 <- pmin(rev(isoreg(rev(v[1:(p-1)]))$yf), 0)
            p2 <- pmin(isoreg(v[(p+1):n])$yf, 0)
            cn <- c(p1, min(v[p],0), p2)
          }
        }
      }
      if (sign == "+"){
        if (all(v[p]>=cn)){
          error <- sum((v - cn)^2)
          if (error < best_error) {
            best_error <- error
            best_v <- cn
          }
        }
      }else{
        if (all(v[p]<=cn)){
          error <- sum((v - cn)^2)
          if (error < best_error) {
            best_error <- error
            best_v <- cn
          }
        }
      }
    }  
    return(best_v)
  }
}

best_unimodal <- function(v){
  vp <- project_to_unimodal_fast(v, sign = "+")
  error.v <- sum((v - vp)^2)
  
  up <- project_to_unimodal_fast(v, sign = "-")
  error.u <- sum((v - up)^2)
  
  if (error.v <= error.u){
    best_v <- vp
  }else{
    best_v <- up
  }
  return(best_v)
}

normalize_factors <- function(result){
  # Handle both naming conventions
  if (!is.null(result$U)) {
    factors <- result$U
    core <- result$Z
  } else if (!is.null(result$factors)) {
    factors <- result$factors
    core <- result$core
  } else {
    stop("Result must have either $U/$Z or $factors/$core")
  }
  
  d <- length(factors)
  ranks <- sapply(factors, ncol)
  
  # Create normalized factor matrices
  factors_normalized <- factors
  core_normalized <- core
  
  for (k in 1:d) {
    min_vals <- apply(factors[[k]], 2, min)
    max_vals <- apply(factors[[k]], 2, max)
    norm.const <- apply(factors[[k]], 2, function(x) max(abs(x)))
    min.max <- rbind(min_vals, max_vals)
    which.norm.const <- apply(min.max, 2, function(x) which.max(abs(x)))
    sign.val <- sign(min.max[cbind(which.norm.const, 1:ncol(min.max))])
    
    scaling_factors <- sign.val * ifelse(norm.const == 0, 0, 1/norm.const)
    factors_normalized[[k]] <- sweep(factors[[k]], MARGIN=2, scaling_factors, `*`)
    
    # Adjust core tensor to compensate for scaling
    scaling_matrix <- diag(sign.val * ifelse(norm.const == 0, 0, norm.const), nrow = ranks[k])
    
    # Apply scaling in mode k
    core_normalized <- ttm(core_normalized, scaling_matrix, k)
  }
  
  # Create new result object preserving original naming
  result_normalized <- result
  if (!is.null(result$U)) {
    result_normalized$U <- factors_normalized
    result_normalized$Z <- core_normalized
  } else {
    result_normalized$factors <- factors_normalized
    result_normalized$core <- core_normalized
  }
  result_normalized$ranks <- ranks
  
  return(result_normalized)
}

# Reconstruct tensor
reconstruct_unimodal <- function(result) {
  # Handle both naming conventions
  if (!is.null(result$U)) {
    # rTensor convention (tucker() output)
    core <- result$Z
    factors <- result$U
  } else if (!is.null(result$factors)) {
    # Custom convention
    core <- result$core
    factors <- result$factors
  } else {
    stop("Result must have either $U/$Z or $factors/$core")
  }
  ttl(core, factors, ms = 1:length(factors))
}

# plot convergence and factor matrices
plot_unimodal_results <- function(result, ALS=FALSE) {
  # Handle both naming conventions
  if (!is.null(result$U)) {
    factors <- result$U
  } else if (!is.null(result$factors)) {
    factors <- result$factors
  } else {
    stop("Result must have either $U or $factors")
  }
  
  # Convergence plots
  if (ALS){
    plot_als_convergence(result)
  }
  
  # Determine which factors to plot (skip trivial constant factors)
  d <- length(factors)
  factors_to_plot <- c()
  
  for (k in 1:d) {
    # Check if factor is trivial (1×1 matrix with constant value)
    if (nrow(factors[[k]]) == 1 && ncol(factors[[k]]) == 1) {
      # Skip this trivial factor
      next
    }
    factors_to_plot <- c(factors_to_plot, k)
  }
  
  # Plot only non-trivial factors
  n_plots <- length(factors_to_plot)
  
  if (n_plots > 0) {
    n_rows <- ceiling(n_plots / 2)
    n_cols <- min(n_plots, 2)
    par(mfrow = c(n_rows, n_cols))
    
    for (k in factors_to_plot) {
      # matplot(factors[[k]], type = 'l', lwd = 3, lty = 1,
      #         # main = paste("Mode", k, "Factor Matrix"),
      #         main = "Basic Functions",
      #         xlab = "Ordinal Value", ylab = "Factor Loading",
      #         ylim = c(min(factors[[k]]),max(factors[[k]])),
      #         col = col[1:ncol(factors[[k]])])
      
      # par(mar = c(bottom, left, top, right)), 
      # the default value for mar is c(5.1, 4.1, 4.1, 2.1).
      par(mar=c(5.1,5.1,2.1,2.1))
      
      matplot(factors[[k]], type = 'l', lwd = 3, lty = 1,
              # main = paste("Mode", k, "Factor Matrix"),
              main = "",
              xlab = "Original values", ylab = "Factor values", cex.axis = 1.5, cex.lab = 1.5,
              ylim = c(min(factors[[k]]),max(factors[[k]])),
              col = col[1:ncol(factors[[k]])])
      grid()
    }
    par(mfrow = c(1, 1))
  } else {
    message("No non-trivial factors to plot")
  }
}

#' Compare marginal distributions of original and reconstructed tensors
#' @param X_original Original tensor
#' @param result Result from tucker_unimodal_admm
#' @param normalize If TRUE, normalize reconstruction to sum to 1
#' @return List with marginal L2 errors for each mode
plot_marginals_comparison <- function(X_original, result, normalize = TRUE) {
  # Handle both naming conventions
  if (!is.null(result$U)) {
    factors <- result$U
  } else if (!is.null(result$factors)) {
    factors <- result$factors
  } else {
    stop("Result must have either $U or $factors")
  }
  
  # Reconstruct tensor
  X_recon <- reconstruct_unimodal(result)
  
  # Normalize if requested
  if (normalize) {
    X_recon@data <- X_recon@data / sum(X_recon@data)
  }
  
  d <- length(factors)
  tensor_dims <- dim(X_original@data)
  
  # Determine which modes to plot (skip trivial factors and singleton dimensions)
  modes_to_plot <- c()
  for (k in 1:d) {
    # Skip if factor is trivial (1×1) or tensor dimension is 1
    if ((nrow(factors[[k]]) == 1 && ncol(factors[[k]]) == 1) || tensor_dims[k] == 1) {
      next
    }
    modes_to_plot <- c(modes_to_plot, k)
  }
  
  n_plots <- length(modes_to_plot)
  
  if (n_plots > 0) {
    n_rows <- ceiling(n_plots / 2)
    n_cols <- min(n_plots, 2)
    par(mfrow = c(n_rows, n_cols))
    
    for (k in modes_to_plot) {
      # Compute marginals by summing over all other dimensions
      marginal_orig <- apply(X_original@data, k, sum)
      marginal_recon <- apply(X_recon@data, k, sum)
      
      # Plot comparison
      x_vals <- 1:length(marginal_orig)
      y_max <- max(c(marginal_orig, marginal_recon)) * 1.1
      
      # plot(x_vals, marginal_orig, type = 'b', col = 'blue', lwd = 3, pch = 16,
      #      # main = paste("Mode", k, "Marginal Distribution"),
      #      main = "Marginal Distribution",
      #      xlab = "Ordinal Value", ylab = "Probability",
      #      ylim = c(0, y_max))
      
      # par(mar = c(bottom, left, top, right)), 
      # the default value for mar is c(5.1, 4.1, 4.1, 2.1).
      par(mar=c(5.1,5.1,2.1,2.1))
      
      plot(x_vals, marginal_orig, type = 'b', col = 'blue', lwd = 3, pch = 16,
           # main = paste("Mode", k, "Marginal Distribution"),
           main = "",
           xlab = "Original Values", ylab = "Probability values", cex.axis = 1.5, cex.lab = 1.5,
           ylim = c(0, y_max))
      lines(x_vals, marginal_recon, type = 'b', col = 'red', lwd = 3, pch = 17)
      grid()
    }
    
    par(mfrow = c(1, 1))
  } else {
    message("No non-trivial modes to plot")
  }
  
  # Return marginal errors for all modes (including trivial ones for completeness)
  marginal_errors_abs <- sapply(1:d, function(k) {
    marginal_orig <- apply(X_original@data, k, sum)
    marginal_recon <- apply(X_recon@data, k, sum)
    sqrt(sum((marginal_orig - marginal_recon)^2))
  })
  
  marginal_errors_rel <- sapply(1:d, function(k) {
    marginal_orig <- apply(X_original@data, k, sum)
    marginal_recon <- apply(X_recon@data, k, sum)
    sqrt(sum((marginal_orig - marginal_recon)^2)) / sqrt(sum(marginal_orig^2))
  })
  
  names(marginal_errors_abs) <- paste("Mode", 1:d)
  names(marginal_errors_rel) <- paste("Mode", 1:d)
  
  result_list <- list(
    absolute_L2 = marginal_errors_abs,
    relative_L2 = marginal_errors_rel
  )
  
  invisible(result_list)
}

