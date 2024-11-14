# Package
library(fda)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(reshape2)
library(splines)
library(snowfall)
library(fields)
library(limSolve)
library(igraph)
library(Rfast)
library(ggnetwork)
library(network)
library(mvtnorm)
library(Matrix)
library(extraDistr)
library(xtable)
library(mgcv)
library(pracma)
library(clusterGeneration)

################################################################################################################
# Main function
#' Graphical Functional Principal Component Analysis (GFPCA)
#' 
#' @description
#' A function to implement graphical FPCAs.
#'
#' @details
#' This is a generic function to implement graphical versions of DFPCA 
#' or SFPCA for multivariate functional time series (MFTS). The function requires 
#' that the MFTS data are densely and regularly observed. 
#' @param MFTS A Z*p*J array containing the regularly observed MFTS data, where
#' Z is the number of observed time points, p is the dimension of MFTS, and
#' J is the time length of MFTS.
#' @param Lt A vector containing the observed time points of MFTS. We require
#' the time points to be contained in the interval [0,1]. 
#' @param Dynamic A logical evaluating to TRUE or FALSE indicating whether 
#' the graphical DFPCA or SFPCA should be conducted.
#' @param Max.comp The maximum number of components.
#' @param FVE A numeric in [0,1] indicating the fraction of variance explained 
#' for determining the number of components. If FVE = 0, the number of component 
#' is selected by the ratio of variance explained.
#' @param Mean.zero A logical evaluating to TRUE or FALSE indicating whether 
#' the mean functions of MFTS are zero.
#' @return A list with components
#' \item{comp_num}{Number of components.}
#' \item{xi_dyn_IN or xi_sta_IN}{A list of FPC scores computed by integration.}
#' \item{xi_dyn_CE or xi_sta_CE}{A list of FPC scores computed by conditional expectation.}
#' \item{mean_function}{A list of mean functions of MFTS, where Lt contains the time grid of the estimated values.}
#' \item{Functional filters or eigenfunctions}{A list of the estimated functional filters (or eigenfunctions), where Lt contains the time grid of the estimated values.}
#' \item{eigen_matrix}{A p*p*(J/2)*comp_num array containing the estimated eigenmatrices for different frequences.}
#' \item{Phi}{A p*p*(J/2)*comp_num array containing the estimated inverse of eigenmatrices by incorporating graph constraints.}
#' \item{mea_error}{A vector containing the estimated variances of the measurement errors.}
#' \item{dmean_smo_fda}{The pre-smoothed coefficients after removing mean trends.}
GFPCA <- function(MFTS, Lt, Dynamic, Max.comp, FVE, Mean.zero){
  
  subject_length <- dim(MFTS)[2] # Dimension
  time_length <- dim(MFTS)[3] # Temporal Length
  
  X_basis <- create.bspline.basis(c(0, 1), breaks = Lt, norder = 4) # Basis functions for pre-smoothing
  X_mat <- eval.basis(X_basis, seq(0, 1, 0.0001))
  tra_mat <- t(sqrtm(crossprod(X_mat) * 0.0001)$B) # Transformation matrix
  true_basis <- eval.basis(X_basis, seq(0, 1, 0.01)) %*% solve(tra_mat) # Orthogonal basis functions
  
  smo_fda <- smooth_fda(MFTS, Lt, subject_length, time_length, X_basis, tra_mat)
  sigma_est <- smo_fda[[1]] # Estimated variances of measurement errors
  smo_fda <- smo_fda[[2]]
  
  # Mean function
  if(Mean.zero == T){
    smo_fda_mean <- matrix(0, dim(smo_fda)[1], dim(smo_fda)[2])
  }else{
    smo_fda_mean <- apply(smo_fda, c(1,2), mean)
  }
  
  mean_function <- list(Lt = seq(0, 1, 0.01),
                        val = true_basis %*% smo_fda_mean)
  
  dmean_smo_fda <- sapply(1:time_length, function(j){
    sapply(1:subject_length, function(i){
      smo_fda[,i,j] - smo_fda_mean[,i]
    }, simplify = "array")
  }, simplify = "array")
  
  # Estimating functional filters
  grid_freq <- seq(-pi, pi, length.out = 1001)[-1001]
  q <- time_length ^ 0.4
  r <- floor(q) 
  
  cov_matrix <- eigen_pre(subject_length, time_length, dmean_smo_fda, q, r)
  sp_matrix <- sp_tran(r, grid_freq, cov_matrix)
  result <- eigen_func_dyn(sp_matrix, grid_freq, FVE, Max.comp, time_length)
  
  comp_num <- result$comp_num
  eigen_value_freq <- result$eigen_value_freq
  eigen_vector_freq <- result$eigen_vector_freq
  eigen_vector_time <- result$eigen_vector_time
  
  functional_filter <- list(Lt = seq(0, 1, 0.01),
                            val = lapply(1:comp_num, function(k){
                              true_basis %*% eigen_vector_time[[k]]
                            }))
  
  time_length_eig <- time_length / 2
  Freq_eig <- 2 * pi / time_length * (1:time_length_eig)
  
  if(Dynamic == T){
    # Estimating eigenmatrices
    eigen_matrix <- eigen_matrix_con(time_length, comp_num, subject_length,
                                     dmean_smo_fda, q, r, cov_matrix, Freq_eig, time_length_eig)
    
    ## Incorporating graph constraints
    lambda <- Sel_lam(eigen_matrix)
    
    eigen_matrix_t <- lapply(1:dim(eigen_matrix)[4], function(k){
      gra_selection(eigen_matrix[,,,k], lambda[k])
    })
    Phi <- array(0, dim(eigen_matrix))
    for(k in 1:dim(Phi)[4]){
      Phi[,,,k] <- eigen_matrix_t[[k]][[1]]
    }
    
    # Score extraction
    basis <- eval.basis(X_basis, Lt)
    basis <- basis %*% solve(tra_mat)
    
    ## Preparation
    xi <- lapply(1:comp_num, function(k){
      matrix(0, subject_length, ncol(eigen_vector_time[[k]]) - 1 + time_length)
    })
    
    fre_tra <- lapply(1:comp_num, function(k){
      time_length_k <- ncol(eigen_vector_time[[k]]) - 1 + time_length
      return(sapply(Freq_eig, function(theta){
        exp(complex(1, 0, 1) * theta * (1:time_length_k)) / sqrt(2 * pi * time_length_k)
      })) 
    })
    
    squ_fre_tra <- lapply(1:comp_num, function(k){
      time_length_k <- ncol(eigen_vector_time[[k]]) - 1 + time_length
      return(sapply(1:time_length_eig, function(theta){
        a <- fre_tra[[k]][,theta]
        return(tcrossprod(a, Conj(a)))
      }, simplify = "array"))
    })
    
    des_mat <- lapply(1:comp_num, function(k){
      sapply(1:time_length, function(j){
        pre_mat <- matrix(0, nrow(eigen_vector_time[[k]]), dim(xi[[k]])[2])
        pre_mat[,j:(j-1+dim(eigen_vector_time[[k]])[2])] <- eigen_vector_time[[k]]
        pre_mat <- basis %*% pre_mat 
        return(pre_mat)
      }, simplify = "array")
    })
    
    dmean_fda <- sapply(1:time_length, function(j){
      sapply(1:subject_length, function(i){
        c(MFTS[,i,j] - basis %*% smo_fda_mean[,i])
      }, simplify = "array")
    }, simplify = "array")
    
    squ_des_mat <- lapply(1:comp_num, function(k){
      mat <- sapply(1:time_length, function(j){
        crossprod(des_mat[[k]][,,j])
      }, simplify = "array")
      mat <- apply(mat, c(1, 2), sum)
      return(mat)
    })
    
    ## Score extraction by integration
    xi_dyn_IN <- lapply(1:comp_num, function(k){
      lag_t <- (ncol(eigen_vector_time[[k]]) - 1) / 2
      return(sapply(1:(ncol(eigen_vector_time[[k]]) - 1 + time_length), function(j){
        sapply(1:subject_length, function(i){
          if(j <= lag_t){
            val <- 0
          }else if(j > time_length + lag_t){
            val <- 0
          }else{
            mat <- cbind(matrix(0, dim(dmean_smo_fda)[1], lag_t),
                         dmean_smo_fda[,i,],
                         matrix(0, dim(dmean_smo_fda)[1], lag_t))
            val <- sum(diag(t(mat[,(j+lag_t):(j-lag_t)]) %*% 
                              eigen_vector_time[[k]]))
          }
          return(val)
        })
      }))
    }) 
    
    ## Score extraction by conditional expectation
    xi_dyn_CE <- score_est(xi_dyn_IN, fre_tra, des_mat, dmean_fda, Phi, squ_des_mat,
                           squ_fre_tra, sigma_est, subject_length, time_length, comp_num)[[1]]
    
    return(
      list(
        comp_num = comp_num,
        xi_dyn_IN = xi_dyn_IN,
        xi_dyn_CE = xi_dyn_CE,
        mean_function = mean_function,
        functional_filter = functional_filter,
        Phi = eigen_matrix,
        mea_error = sigma_est,
        dmean_smo_fda = dmean_smo_fda
      )
    )
  }else{
    # Estimating eigenfunctions and eigenmatrices
    cov_matrix_0 <- eigen_pre(subject_length, time_length, dmean_smo_fda, 0, 0)
    result <- eigen_func_sta(cov_matrix_0, FVE, Max.comp)
    comp_num_sta <- result$comp_num_sta
    eigen_vector_sta <- result$eigen_vector_sta
    
    eigenfunction <- list(Lt = seq(0, 1, 0.01),
                          val = lapply(1:comp_num, function(k){
                            true_basis %*% eigen_vector_sta[k,]
                          }))
    
    eigen_matrix <- eigen_matrix_con(time_length, comp_num_sta, subject_length,
                                     dmean_smo_fda, q, r, cov_matrix, Freq_eig, time_length_eig)
    
    # Estimating eigenmatrices
    eigen_matrix <- eigen_matrix_con(time_length, comp_num, subject_length,
                                     dmean_smo_fda, q, r, cov_matrix, Freq_eig, time_length_eig)
    
    ## Incorporating graph constraints
    lambda <- Sel_lam(eigen_matrix)
    
    eigen_matrix_t <- lapply(1:dim(eigen_matrix)[4], function(k){
      gra_selection(eigen_matrix[,,,k], lambda[k])
    })
    Phi <- array(0, dim(eigen_matrix))
    for(k in 1:dim(Phi)[4]){
      Phi[,,,k] <- eigen_matrix_t[[k]][[1]]
    }
    
    # Score extraction
    basis <- eval.basis(X_basis, Lt)
    basis <- basis %*% solve(tra_mat)
    
    ## Preparation
    xi <- lapply(1:comp_num_sta, function(k){
      matrix(0, subject_length, time_length)
    })
    
    fre_tra <- lapply(1:comp_num_sta, function(k){
      time_length_k <- time_length
      return(sapply(Freq_eig, function(theta){
        exp(complex(1, 0, 1) * theta * (1:time_length_k)) / sqrt(2 * pi * time_length_k)
      })) 
    })
    
    squ_fre_tra <- lapply(1:comp_num_sta, function(k){
      return(sapply(1:time_length_eig, function(theta){
        a <- fre_tra[[k]][,theta]
        a <- tcrossprod(a, Conj(a))
        a <- (a + t(Conj(a))) / 2
        return(a)
      }, simplify = "array"))
    })
    
    des_mat <- lapply(1:comp_num_sta, function(k){
      sapply(1:time_length, function(j){
        pre_mat <- matrix(0, nrow(eigen_vector_sta), dim(xi[[k]])[2])
        pre_mat[,j] <- eigen_vector_sta[k,]
        pre_mat <- basis %*% pre_mat 
        return(pre_mat)
      }, simplify = "array")
    })
    
    dmean_fda <- sapply(1:time_length, function(j){
      sapply(1:subject_length, function(i){
        c(MFTS[,i,j] - basis %*% smo_fda_mean[,i])
      }, simplify = "array")
    }, simplify = "array")
    
    squ_des_mat <- lapply(1:comp_num_sta, function(k){
      mat <- sapply(1:time_length, function(j){
        crossprod(des_mat[[k]][,,j])
      }, simplify = "array")
      mat <- apply(mat, c(1, 2), sum)
      return(mat)
    })
    
    ## Score extraction by integration 
    xi_sta_IN <- lapply(1:comp_num_sta, function(k){
      sapply(1:time_length, function(j){
        t(dmean_smo_fda[,,j]) %*% eigen_vector_sta[k,]
      })
    }) 
    
    ##  Score extraction by conditional expectation
    xi_sta_CE <- score_est(xi_sta_IN, fre_tra, des_mat, dmean_fda, Phi, squ_des_mat,
                           squ_fre_tra, sigma_est, subject_length, time_length, comp_num_sta)[[1]]
    
    return(
      list(
        comp_num = comp_num_sta,
        xi_sta_IN = xi_sta_IN,
        xi_sta_CE = xi_sta_CE,
        mean_function = mean_function,
        eigenfunction = eigenfunction,
        eigen_matrix = eigen_matrix,
        Phi = Phi,
        mea_error = sigma_est,
        dmean_smo_fda = dmean_smo_fda
      )
    )
  }
}

################################################################################################################
# Basic function
## Summation operation for an array
array_last_sum <- function(array){
  dim_array <- dim(array)
  array <- matrix(array, prod(dim_array[1:2]))
  array <- rowsums(array)
  return(matrix(array, dim_array[1]))
}

################################################################################################################
# Functions for data generation
## Generating scores
score_gen <- function(subject_length, time_length, comp_length, rho, covariance_matrix, lag){
  score <- array(0, c(comp_length, subject_length, (time_length + 2)))
  init <- rmvnorm(1, rep(0, comp_length * subject_length), covariance_matrix)
  score[,,1] <- sapply(1:subject_length, function(i){
    init[((i - 1) * comp_length + 1):(i * comp_length)]
  })
  for(j in 2:(time_length + 2)){
    noise <- rmvnorm(1, rep(0, comp_length * subject_length), covariance_matrix)
    score[,,j] <- sapply(1:subject_length, function(i){
      rho * score[,i,j-1] + noise[((i - 1) * comp_length + 1):(i * comp_length)]
    })
  }
  return(score)
}

## Generating mean functions
ture_mean <- function(t, mean_sample){
  fourier(t, nbasis = 3) %*% mean_sample
}

## Generating eigenfunctions
basis_f <- function(x, lag, k, comp_length){
  use_num <- (k - 1) * 2 * (k <= lag) + lag * 2 * (k > lag) 
  
  real_lag <- 1 * (k <= lag) + 0 * (k > lag)
  weight <- sapply((-1):1, function(g){
    ifelse(abs(g) <= real_lag, exp(- 2 * abs(g)), 0)
  })
  weight <- sqrt(weight / sum(weight))
  main <- fourier(x, nbasis = comp_length + lag * 2, period = 1)[,1:(comp_length + lag * 2)]
  main <- sapply((-1):1, function(g){
    if(abs(g) > real_lag){
      rep(0, length(x))
    }else if(g == 0){
      main[,k]
    }else if(g  < 0){
      main[,g + real_lag + comp_length + use_num + 1]
    }else{
      main[,g + real_lag + comp_length + use_num]
    }
  })
  main <- t(main) * weight
  return(t(main))
}

##  Generating functional data
if(exists("sep")){
  if(sep == T){
    fda_gen <- function(x, i, j, mean_sample, comp_length, score, lag, subject_length){
      as.numeric(ture_mean(x, mean_sample) +
                   rowsums(sapply(1:comp_length, function(k){
                     basis_f(x, lag, k, comp_length) %*% score[k,i,j:(j + 2)]
                   })))
    }
  }else{
    fda_gen <- function(x, i, j, mean_sample, comp_length, score, lag, subject_length){
      as.numeric(ture_mean(x, mean_sample) + 
                   rowsums(sapply(1:comp_length, function(k){
                     (basis_f(x, lag, k, comp_length) * (5 * sin(i * x / subject_length) + 1)) %*% score[k,i,j:(j + 2)]
                   })))
    }
  }
}

################################################################################################################
# Functions for GDFPCA & GSFPCA
#' Pre-smoothing of Multivariate Functional Time Series (MFTS)
#'  
#' @description
#' A function to pre-smooth MFTS data.
#' 
#' @param fda An array containing the observed MFTS data.
#' @param x_fda The observed time points of MFTS.
#' @param subject_length The dimension of MFTS.
#' @param time_length The time length of MFTS.
#' @param X_basis A basis object for pre-smoothing.
#' @param tra_mat A transformation matrix for the basis functions.
#' @return A list with the pre-smoothed coefficients and the estimated variances of measurement errors.
smooth_fda <- function(fda, x_fda, subject_length, time_length, X_basis, tra_mat){
  
  fit <- lapply(1:subject_length, function(i){
    return(lapply(1:time_length, function(j){
      fit <- smooth.spline(x = x_fda, y = fda[,i,j], cv = F)
      return(fit)
    }))
  })
  
  error <- sapply(1:subject_length, function(i){
    sum(sapply(1:time_length, function(j){
      fit[[i]][[j]]$pen.crit / (dim(fda)[1] - fit[[i]][[j]]$df)
    })) / time_length
  })
  
  coef <- sapply(1:time_length, function(j){
    sapply(1:subject_length, function(i){
      return(as.numeric(tra_mat %*% fit[[i]][[j]]$fit$coef))
    }, simplify = "array")
  }, simplify = "array")
  
  return(list(error, coef))
}

#' Estimating Cross-covariance
#'
#' @description
#' A function to estimate cross-covariance for different time lags.
#' 
#' @param subject_length The dimension of MFTS.
#' @param time_length The time length of MFTS.
#' @param dmean_smo_fda The pre-smoothed coefficients after removing mean trends.
#' @param p,r The time lags for the estimated covariance functions.
#' @return An array of the estimated coefficients for the cross-covariance functions.
eigen_pre <- function(subject_length, time_length, dmean_smo_fda, q, r) {
  
  # Create spectral density
  Result <- lapply(0:r, function(g){
    matrix(rowsums(sapply(1:subject_length, function(i){
      tcrossprod(dmean_smo_fda[,i,(1+g):time_length], dmean_smo_fda[,i,1:(time_length-g)]) / time_length
    })), dim(dmean_smo_fda)[1])
  })
  
  cov_matrix <- array(0, c(2 * r + 1, dim(dmean_smo_fda)[1], dim(dmean_smo_fda)[1]))
  if(r == 0){
    cov_matrix[1,,] <- Result[[1]] / 2 / pi
  }else{
    for(i in r:1){
      cov_matrix[r-i+1,,] <- t(Result[[i+1]]) * (1 - i / q) / 2 / pi
    }
    for(i in 0:r){
      cov_matrix[r+i+1,,] <- Result[[i+1]] * (1 - i / q) / 2 / pi
    }
  }
  
  return(cov_matrix)
}

#' Constructing Spectral Density Kernels
#'
#' @description
#' A function to construct spectral densities given cross-covariance.
#' 
#' @param r The time lags for estimating the spectral density kernels.
#' @param grid_freq The frequencies of the estimated spectral density kernels.
#' @param cov_matrix An array of the estimated coefficients for the cross-covariance functions.
#' @return An array of the estimated coefficients for spectral density kernels.
sp_tran <- function(r, grid_freq, cov_matrix) {
  exp_freq <- exp(complex(1, 0, 1) * tcrossprod((-r):r, grid_freq))
  return(sapply(1:dim(cov_matrix)[3], function(t){
    crossprod(exp_freq, cov_matrix[,,t])
  }, simplify = "array"))
}

#' Estimating Eigenfunctions
#'
#' @description
#' A function to estimate eigenfunctions for static FPCA.
#' 
#' @param cov_matrix_0 The estimated coefficients for lag-zero covariance functions.
#' @param sel_eig The fraction of variance explained.
#' @param max_comp The maximum number of components.
#' @return Estimations of the number of components and the coefficients of eigenfunctions.
eigen_func_sta <- function(cov_matrix_0, sel_eig, max_comp){
  eigen_decom_sta <- eigen(cov_matrix_0[1,,], symmetric = T)
  eigen_value_sta <- eigen_decom_sta$values 
  eigen_vector_sta <- t(eigen_decom_sta$vectors)
  
  tol_var_comp_s <- eigen_value_sta[eigen_value_sta >= 0]
  if(sel_eig == 0){
    comp_num_sta <- which.max((tol_var_comp_s[-length(tol_var_comp_s)] / tol_var_comp_s[-1])[1:max_comp])
  }else{
    comp_num_sta <- sum(cumsum(tol_var_comp_s / sum(tol_var_comp_s)) < sel_eig) + 1
  }
  
  return(list(comp_num_sta = comp_num_sta, eigen_vector_sta = eigen_vector_sta))
}

#' Estimating Functional Filters
#'
#' @description
#' A function to estimate functional filters for dynamic FPCA.
#' 
#' @param sp_matrix An array of the estimated coefficients for spectral density kernels.
#' @param grid_freq The frequencies of the estimated spectral density kernels.
#' @param sel_eig The fraction of variance explained.
#' @param max_comp The maximum number of components.
#' @param time_length The time length of MFTS.
#' @return Estimations of the number of components and the coefficients of functional filters.
eigen_func_dyn <- function(sp_matrix, grid_freq, sel_eig, max_comp, time_length){
  
  ### Eigen-decomposition of spectral density matrix for different frequencies
  eigen_decom <- lapply(1:length(grid_freq), function(k){
    eigen_decom_freq <- eigen(sp_matrix[k,,], symmetric = T)
    eigen_value_freq <- eigen_decom_freq$values 
    eigen_vector_freq <- t(Conj(eigen_decom_freq$vectors))
    return(list(eigen_value_freq, eigen_vector_freq))
  })
  
  eigen_value_freq <- sapply(1:length(grid_freq), function(k){
    eigen_decom[[k]][[1]]
  })
  
  eigen_vector_freq <- sapply(1:length(grid_freq), function(k){
    eigen_decom[[k]][[2]]
  }, simplify = "array")
  
  ### Determine the number of component
  tol_var_comp <- rowMeans(eigen_value_freq) * 2 * pi
  tol_var_comp <- tol_var_comp[tol_var_comp >= 0]
  if(sel_eig == 0){
    comp_num <- which.max((tol_var_comp[-length(tol_var_comp)] / tol_var_comp[-1])[1:max_comp])
  }else{
    comp_num <- sum(cumsum(tol_var_comp / sum(tol_var_comp)) < sel_eig) + 1
  }
  
  # Determine the eigen-function on time domain
  eigen_vector_time <- lapply(1:comp_num, function(k){
    l <- 0
    phi <- Re(eigen_vector_freq[k,,] %*% exp(- complex(1, 0, 1) * l * grid_freq)) / length(grid_freq)
    while(l < time_length) {
      l <- l + 1
      phi_l_1 <- Re(eigen_vector_freq[k,,] %*% exp(- complex(1, 0, 1) * l * grid_freq)) / length(grid_freq)
      phi_l_2 <- Re(eigen_vector_freq[k,,] %*% exp(complex(1, 0, 1) * l * grid_freq)) / length(grid_freq)
      phi <- cbind(phi_l_2, phi, phi_l_1)
    }
    
    mark <- which.max(colsums(phi ^ 2)) - time_length - 1
    l <- 0
    phi <- Re(eigen_vector_freq[k,,] %*% exp(- complex(1, 0, 1) * (l + mark) * grid_freq)) / length(grid_freq)
    while(((sum(phi ^ 2)) < 0.99) & (l < time_length / 5)){
      l <- l + 1
      phi_l_1 <- Re(eigen_vector_freq[k,,] %*% exp(- complex(1, 0, 1) * (mark + l) * grid_freq)) / length(grid_freq)
      phi_l_2 <- Re(eigen_vector_freq[k,,] %*% exp(- complex(1, 0, 1) * (mark - l) * grid_freq)) / length(grid_freq)
      phi <- cbind(phi_l_2, phi, phi_l_1)
    }
    return(phi)
  })
  
  return(list(comp_num = comp_num, 
              eigen_value_freq = eigen_value_freq,
              eigen_vector_time = eigen_vector_time
  ))
}

#' Estimating Eigen-matrices
#'
#' @description
#' A function to estimate eigen-matrices.
#' 
#' @param time_length The time length of MFTS.
#' @param comp_num The number of components.
#' @param subject_length The dimension of MFTS.
#' @param dmean_smo_fda The pre-smoothed coefficients after removing mean trends.
#' @param p,r The time lags for the estimated covariance functions.
#' @param cov_matrix An array of the estimated coefficients for the cross-covariance functions.
#' @param Freq_eig The frequencies for the estimated eigenmatrices.
#' @param time_length_eig The number of the evaluated frequencies.
#' @return An array of the estimated eigenmatrices.
eigen_matrix_con <- function(time_length, comp_num, subject_length,
                             dmean_smo_fda, q, r, cov_matrix, 
                             Freq_eig, time_length_eig){
  
  sp_matrix <- sp_tran(r, Freq_eig, cov_matrix)
  
  # Eigen-decomposition of spectral density matrix for different frequencies
  eigen_decom <- lapply(1:length(Freq_eig), function(k){
    eigen_decom_freq <- eigen(sp_matrix[k,,], symmetric = T)
    eigen_value_freq <- eigen_decom_freq$values 
    eigen_vector_freq <- t(Conj(eigen_decom_freq$vectors))
    return(list(eigen_value_freq, eigen_vector_freq))
  })
  
  eigen_vector_freq <- sapply(1:length(Freq_eig), function(k){
    eigen_decom[[k]][[2]]
  }, simplify = "array")
  
  cov_matrix <- array(0, c(2 * r + 1, prod(dim(dmean_smo_fda)[1:2]), prod(dim(dmean_smo_fda)[1:2])))
  dmean_smo_fda_1 <- sapply(1:dim(dmean_smo_fda)[3], function(k){
    c(dmean_smo_fda[,,k])
  })
  Result <- lapply(0:r, function(g){
    tcrossprod(dmean_smo_fda_1[,(1+g):time_length], dmean_smo_fda_1[,1:(time_length-g)]) / time_length
  })
  
  if(r == 0){
    cov_matrix[1,,] <- Result[[1]] / 2 / pi
  }else{
    for(i in r:1){
      cov_matrix[r-i+1,,] <- t(Result[[i+1]]) * (1 - i / q) / 2 / pi
    }
    for(i in 0:r){
      cov_matrix[r+i+1,,] <- Result[[i+1]] * (1 - i / q) / 2 / pi
    }
  }
  
  sp_matrix <- sp_tran(r, Freq_eig, cov_matrix)
  
  eigen_matrix <- sapply(1:comp_num, function(k){
    sapply(1:time_length_eig, function(theta){
      sapply(1:subject_length, function(j){
        sapply(1:subject_length, function(i){
          t(eigen_vector_freq[k,,theta]) %*% sp_matrix[theta,((i-1)*dim(eigen_vector_freq)[1]+1):(i*dim(eigen_vector_freq)[1]),((j-1)*dim(eigen_vector_freq)[1]+1):(j*dim(eigen_vector_freq)[1])] %*% Conj(eigen_vector_freq[k,,theta])
        })
      })
    }, simplify = "array")
  },  simplify = "array")
  
  return(eigen_matrix)
}

#' Testing Dynamic Weak Separability
#'
#' @description
#' A function to illustrate the validity of dynamic weak separability.
#' 
#' @param time_length The time length of MFTS.
#' @param comp_num The number of components.
#' @param subject_length The dimension of MFTS.
#' @param dmean_smo_fda The pre-smoothed coefficients after removing mean trends.
#' @param p,r The time lags for the estimated covariance functions.
#' @param Freq_eig The frequencies for the estimated eigenmatrices.
#' @param time_length_eig The number of the evaluated frequencies.
#' @return An object to illustrate the validity of dynamic weak separability.
test_weaksep <- function(time_length, comp_num, subject_length,
                         dmean_smo_fda, q, r,  
                         Freq_eig, time_length_eig){
  
  cov_matrix <- eigen_pre(subject_length, time_length, dmean_smo_fda, q, r)
  sp_matrix <- sp_tran(r, Freq_eig, cov_matrix)
  
  eigen_decom <- lapply(1:length(Freq_eig), function(k){
    eigen_decom_freq <- eigen(sp_matrix[k,,], symmetric = T)
    eigen_value_freq <- eigen_decom_freq$values 
    eigen_vector_freq <- t(Conj(eigen_decom_freq$vectors))
    return(list(eigen_value_freq, eigen_vector_freq))
  })
  
  eigen_vector_freq <- sapply(1:length(Freq_eig), function(k){
    eigen_decom[[k]][[2]]
  }, simplify = "array")
  
  cov_matrix <- array(0, c(2 * r + 1, prod(dim(dmean_smo_fda)[1:2]), prod(dim(dmean_smo_fda)[1:2])))
  dmean_smo_fda_1 <- sapply(1:dim(dmean_smo_fda)[3], function(k){
    c(dmean_smo_fda[,,k])
  })
  Result <- lapply(0:r, function(g){
    tcrossprod(dmean_smo_fda_1[,(1+g):time_length], dmean_smo_fda_1[,1:(time_length-g)]) / time_length
  })
  
  if(r == 0){
    cov_matrix[1,,] <- Result[[1]] / 2 / pi
  }else{
    for(i in r:1){
      cov_matrix[r-i+1,,] <- t(Result[[i+1]]) * (1 - i / q) / 2 / pi
    }
    for(i in 0:r){
      cov_matrix[r+i+1,,] <- Result[[i+1]] * (1 - i / q) / 2 / pi
    }
  }
  
  sp_matrix <- sp_tran(r, Freq_eig, cov_matrix)
  
  cross_eigen_matrix <- lapply(1:comp_num, function(k){
    sapply(k:comp_num, function(k_1){
      sapply(1:time_length_eig, function(theta){
        sapply(1:subject_length, function(j){
          sapply(1:subject_length, function(i){
            t(eigen_vector_freq[k,,theta]) %*% sp_matrix[theta,((i-1)*dim(eigen_vector_freq)[1]+1):(i*dim(eigen_vector_freq)[1]),((j-1)*dim(eigen_vector_freq)[1]+1):(j*dim(eigen_vector_freq)[1])] %*% Conj(eigen_vector_freq[k_1,,theta])
          })
        })
      }, simplify = "array")
    },  simplify = "array")
  })
  
  cov_mat <- lapply(1:comp_num, function(k){
    sapply(1:dim(cross_eigen_matrix[[k]])[4], function(k_1){
      sum(abs(cross_eigen_matrix[[k]][,,,k_1]) ^ 2)
    })
  })
  
  return(cov_mat)
}

#' Regularized Estimator for Precision Matrices
#'
#' @description
#' Joint regularized estimator for precision matrices by ADMM.
#' 
#' @param S An array of eigenmatrices.
#' @param lambda A tuning parameter to the joint regularized estimator.
#' @return A list of two collections of precision matrices.
gra_selection <- function(S, lambda){
  rho <- 1
  time_length_eig <- dim(S)[3] 
  
  X <- sapply(1:dim(S)[3], function(k){
    diag(1, dim(S)[2])
  }, simplify = "array")
  X_t <- array(0, dim(S))
  Z <- X
  U <- array(0, dim(S))
  init <- 1
  while(((sum((abs(X-X_t))) / sum((abs(X_t))) > 10 ^ (-5)) | (init <= 30) | (sum(abs(X - Z)) > 10 ^ (-10))) &
        (init < 100)
  ){
    init <- init + 1
    X_t <- X
    A <- S + rho * (U - Z)
    X <- sapply(1:dim(A)[3], function(k){
      step_1 <- eigen(A[,,k], symmetric = T)
      step_1$values <- Re((- step_1$values + sqrt(step_1$values ^ 2 + 4 * rho)) / 2 / rho)
      return(
        tcrossprod(crossprod(t(step_1$vectors), diag(step_1$values)), Conj(step_1$vectors))
      )
    }, simplify = "array")
    
    Y <- X + U
    mat <- Re(sapply(1:dim(Y)[1], function(h){
      sapply(1:dim(Y)[1], function(g){
        sqrt(sum(abs(Y[g,h,]) ^ 2))
      })}))
    mat <- (1 - lambda / rho / mat) * ((1 - lambda / rho / mat) > 0)
    diag(mat) <- 1
    
    Z <- sapply(1:dim(Z)[3], function(k){
      mat * Y[,,k]
    }, simplify = "array")
    U <- U + (X - Z)
    rm(Y)
  }
  
  return(list(X, Z))
}

# Incorporating a known graph constraint on the estimated precision matrices
graph_pre <- function(S, adjacency_matrix, subject_length){
  
  adjacency_matrix_t <- adjacency_matrix
  diag(adjacency_matrix_t) <- F
  
  eigen_matrix_test_t <- S
  eigen_matrix_test_tt <- eigen_matrix_test_t + 1000
  init <- 1
  while ((sum(abs(eigen_matrix_test_t - eigen_matrix_test_tt)) / sum(abs(eigen_matrix_test_tt)) > 10 ^ (-3)) & (init < 200))  {
    init <- init + 1
    eigen_matrix_test_tt <-eigen_matrix_test_t
    for(i in 1:subject_length){
      eigen_matrix_test_t[-i, i,] <- sapply(1:(dim(eigen_matrix_test_t)[3]), function(theta){
        beta <- rep(0, dim(S)[1] - 1)
        if(sum(adjacency_matrix[i,-i]) != 0){
          beta[adjacency_matrix_t[i,-i]] <- solve(eigen_matrix_test_t[adjacency_matrix_t[i,],adjacency_matrix_t[i,],theta], S[adjacency_matrix_t[i,],i,theta])
        }
        return(eigen_matrix_test_t[-i,-i,theta] %*% beta)
      })
    }
    error <- sum(sapply(1:(dim(eigen_matrix_test_t)[3]), function(theta){
      sum(abs(solve(eigen_matrix_test_t[,,theta])[adjacency_matrix == F])) 
    }))
  }
  
  inv_eigen_matrix_test <- sapply(1:(dim(eigen_matrix_test_t)[3]), function(theta){
    a <- array(0, dim(eigen_matrix_test_t)[1:2])
    a[adjacency_matrix] <- solve(eigen_matrix_test_t[,,theta])[adjacency_matrix]
    return((a + t(Conj(a))) / 2)
  }, simplify = "array")
  
  return(inv_eigen_matrix_test)
}

# Select lambda by AIC
## Estimating the degree of smoothness for AIC
smooth_est <- function(eigen_matrix, X_basis, Freq_eig){
  time_length_eig <- dim(eigen_matrix)[3]
  subject_length <- dim(eigen_matrix)[2]
  
  fit_dat_re <- matrix(0, subject_length * subject_length, time_length_eig)
  fit_dat_im <- matrix(0, subject_length * subject_length, time_length_eig)
  for(i in 1:subject_length){
    fit_dat_re[((i-1)  * subject_length + 1):(i * subject_length),] <- Re(eigen_matrix[i,,])
    fit_dat_im[((i-1)  * subject_length + 1):(i * subject_length),] <- Im(eigen_matrix[i,,])
  }
  
  fit_dat <- rbind(fit_dat_re, fit_dat_im)
  fit_TOT <- sum((fit_dat - rowMeans(fit_dat)) ^ 2)
  lambda <- 0
  R2 <- 1
  init <- 0
  while((R2 > 0.95) & (init < 1000)){
    lambda <- lambda + 0.01
    init <- init + 1
    fit <- smooth.basisPar(Freq_eig, t(fit_dat), X_basis, Lfdobj = int2Lfd(2), lambda = lambda)
    R2 <- 1 - fit$SSE / fit_TOT
  }
  
  return(fit$df)
} 

AIC <- function(eigen_matrix, time_length_eig, lambda, k, X_basis, Freq_eig){
  
  S <- eigen_matrix[,,,k]
  S_lam <- gra_selection(S, lambda)
  
  df <- smooth_est(S_lam[[1]], X_basis, Freq_eig)
  
  aic <- sum(Re(sapply(1:time_length_eig, function(theta){
    Re(sum(t(eigen_matrix[,,theta,k]) * (S_lam[[1]])[,,theta]) -
         sum(log((eigen((S_lam[[1]])[,,theta], symmetric = T))[[1]])))
  }))) + df * sum(abs(apply(S_lam[[2]], c(1, 2), function(u) sqrt(sum(abs(u) ^ 2)))) != 0)
  
  aic <- is.finite(aic) * aic + (1 - is.finite(aic)) * 10 ^ 10
  
  return(aic)
}

#' Selection of Tuning Parameters for Joint Graphical Lasso
#'
#' @description
#' A function to select tuning parameters for joint graphical Lasso via AIC.
#' 
#' @param eigen_matrix An array of eigenmatrices.
#' @return A selected tuning parameter.
Sel_lam <- function(eigen_matrix){
  Freq_eig <- 2 * pi / (2 * dim(eigen_matrix)[3]) * (1:dim(eigen_matrix)[3])
  X_basis <- create.bspline.basis(c(0, pi), breaks = Freq_eig, norder = 4) 
  
  lambda <- sapply(1:dim(eigen_matrix)[4], function(k){
    upper <- max(apply(eigen_matrix[,,,k], c(1,2), function(h){sqrt(sum(abs(h) ^ 2))})[lower.tri(diag(1:dim(eigen_matrix)[2]))])
    lamb <- optim(par = upper / 2, fn = AIC, 
                  eigen_matrix = eigen_matrix,
                  time_length_eig = dim(eigen_matrix)[3], 
                  k = k,
                  X_basis = X_basis, 
                  Freq_eig = Freq_eig,
                  method = "Brent",
                  lower = 0,
                  upper = upper,
                  control = list(reltol = 1))$par
    return(lamb)
  })
  
  return(lambda)
}

## Conditional likelihood for the score extraction
Lik_xi <- function(xi, fre_tra, des_mat, dmean_fda, eigen_matrix, sigma_est, comp_num, time_length){
  graph_prior <- Re(sum(sapply(1:comp_num, function(k){
    mat <- crossprod(t(xi[[k]]), fre_tra[[k]])
    sum(sapply(1:dim(fre_tra[[k]])[2], function(theta){
      Conj(mat[,theta]) %*% eigen_matrix[,,theta,k] %*% mat[,theta] 
    }))
  })))
  
  obser_lik <- sum(sapply(1:time_length, function(j){
    sum(t(dmean_fda[,,j] - array_last_sum(sapply(1:comp_num, function(k){
      tcrossprod(des_mat[[k]][,,j], xi[[k]])
    }, simplify = "array"))) ^ 2  / 2 / sigma_est)
  })) 
  
  return(graph_prior + obser_lik)
}

#' Score Extraction
#'
#' @description
#' A function to extract FPCA scores for Graphical FPCA.
#' 
#' @param xi An initial value of the scores.
#' @param fre_tra,des_mat,dmean_fda,squ_des_mat,squ_fre_tra Some auxiliary objects for the score extraction.
#' @param eigen_matrix An array of eigenmatrices.
#' @param sigma_est A vector of the estimated variances for measurement errors.
#' @param subject_length The dimension of MFTS.
#' @param time_length The time length of MFTS.
#' @param comp_num The number of components.
#' @return A list containing the estimated scores and other objects relating to the convergence of the algorithm.
score_est <- function(xi, fre_tra, des_mat, dmean_fda, eigen_matrix, squ_des_mat,
                      squ_fre_tra, sigma_est, subject_length, time_length, comp_num){
  Lik <- Lik_xi(xi, fre_tra, des_mat, dmean_fda, eigen_matrix, sigma_est, comp_num, time_length)
  
  gradient <- lapply(1:comp_num, function(k){
    matrix(1, dim(xi[[k]])[1], dim(xi[[k]])[2])
  })
  
  xi_doubel <- xi
  xi_doubel_t <- xi
  xi_t <- lapply(1:comp_num, function(k){
    xi[[k]] + 100
  })
  abs_lik <- 1
  Lik_t <- Lik
  ite <- 1
  half <- 1
  init_step <- 1
  while((sum(abs(unlist(xi) - unlist(xi_t))) / sum(abs(unlist(xi_t))) > 0.005) & ((sum((unlist(gradient)) ^ 2) > 10 ^ (-2)))){
    xi_t <- xi
    xi_doubel_t <- xi_doubel
    gradient <- lapply(1:comp_num, function(k){
      (tcrossprod(xi[[k]], squ_des_mat[[k]]) - array_last_sum(sapply(1:time_length, function(j){
        t(crossprod(des_mat[[k]][,,j], dmean_fda[,,j] - array_last_sum(sapply((1:comp_num)[-k], function(g){
          tcrossprod(des_mat[[g]][,,j], xi[[g]])
        }, simplify = "array"))
        )) 
      }, simplify = "array"))) / sigma_est +
        array_last_sum(sapply(1:dim(eigen_matrix)[3], function(theta) {
          Re(tcrossprod(tcrossprod(eigen_matrix[,,theta,k], t(xi[[k]])), t(squ_fre_tra[[k]][, , theta])))
        }, simplify = "array")) * 2
    })
    
    if(half < 100){
      if(half > 2){
        init_step <- init_step / 1.5
      }
      step_size <- init_step
      half <- 1
      while((Lik_t > Lik - 0.5 * step_size * sum((unlist(gradient)) ^ 2)) &
            (half < 100)){
        step_size <- step_size * 0.8
        half <- half + 1
        xi_doubel <- lapply(1:comp_num, function(k){
          xi_t[[k]] - step_size * gradient[[k]]
        })
        xi <- lapply(1:comp_num, function(k){
          xi_doubel[[k]] + ite / (ite + 3) * (xi_doubel[[k]] - xi_doubel_t[[k]])
        })
        Lik_t <- Lik_xi(xi, fre_tra, des_mat, dmean_fda, eigen_matrix, sigma_est, comp_num, time_length)
      }
    }else{
      step_size <- init_step
      while((Lik_t > Lik - 0.5 * step_size * sum((unlist(gradient)) ^ 2))){
        step_size <- step_size * 0.8
        xi <- lapply(1:comp_num, function(k){
          xi_t[[k]] - step_size * gradient[[k]]
        })
        Lik_t <- Lik_xi(xi, fre_tra, des_mat, dmean_fda, eigen_matrix, sigma_est, comp_num, time_length)
      }
    }
    abs_lik <- abs(Lik - Lik_t)
    Lik <- Lik_t
    ite <- ite + 1
  }
  return(list(xi, sum((unlist(gradient)) ^ 2), abs_lik))
}

################################################################################################################
# Conducting GDFPCA & GSFPCA
## First-step estimation
first_est <- function(x_fda, fda, mean_zero, max_comp, sel_eig){
  
  # Basis 
  subject_length <- dim(fda)[2]
  time_length <- dim(fda)[3]
  
  X_basis <- create.bspline.basis(c(0, 1), breaks = x_fda, norder = 4) # Basis functions for pre-smoothing
  X_mat <- eval.basis(X_basis, seq(0, 1, 0.0001))
  tra_mat <- t(sqrtm(crossprod(X_mat) * 0.0001)$B) # Transformation matrix
  true_basis <- eval.basis(X_basis, seq(0, 1, 0.01)) %*% solve(tra_mat) # Orthogonal basis functions
  
  smo_fda <- smooth_fda(fda, x_fda, subject_length, time_length, X_basis, tra_mat)
  sigma_est <- smo_fda[[1]]
  smo_fda <- smo_fda[[2]]
  
  # Preparation
  ## Mean
  if(mean_zero == T){
    smo_fda_mean <- matrix(0, dim(smo_fda)[1], dim(smo_fda)[2])
  }else{
    smo_fda_mean <- apply(smo_fda, c(1,2), mean)
  }
  
  dmean_smo_fda <- sapply(1:time_length, function(j){
    sapply(1:subject_length, function(i){
      smo_fda[,i,j] - smo_fda_mean[,i]
    }, simplify = "array")
  }, simplify = "array")
  
  ## Static FPCA
  ### Separated FPCA
  Sep_sta <- lapply(1:subject_length, function(i){
    cov_matrix_0 <- eigen_pre(1, time_length, array(dmean_smo_fda[,i,], c(dim(dmean_smo_fda)[1], 1, dim(dmean_smo_fda)[3])), 0, 0)
    return(eigen_func_sta(cov_matrix_0, sel_eig, max_comp))
  })
  
  ### Weakly-separable FPCA
  cov_matrix_0 <- eigen_pre(subject_length, time_length, dmean_smo_fda, 0, 0)
  result <- eigen_func_sta(cov_matrix_0, sel_eig, max_comp)
  comp_num_sta <- result$comp_num_sta
  eigen_vector_sta <- result$eigen_vector_sta
  
  ## Dynamic FPCA
  grid_freq <- seq(-pi, pi, length.out = 1001)[-1001]
  q <- time_length ^ 0.4
  r <- floor(q) 
  
  ### Separated FPCA
  Sep_dyn <- lapply(1:subject_length, function(i){
    cov_matrix <- eigen_pre(1, time_length, array(dmean_smo_fda[,i,], c(dim(dmean_smo_fda)[1], 1, dim(dmean_smo_fda)[3])), q, r)
    sp_matrix <- sp_tran(r, grid_freq, cov_matrix)
    return(eigen_func_dyn(sp_matrix, grid_freq, sel_eig, max_comp, time_length))
  })
  
  ### Weakly-separable FPCA
  cov_matrix <- eigen_pre(subject_length, time_length, dmean_smo_fda, q, r)
  sp_matrix <- sp_tran(r, grid_freq, cov_matrix)
  result <- eigen_func_dyn(sp_matrix, grid_freq, sel_eig, max_comp, time_length)
  
  comp_num <- result$comp_num
  eigen_value_freq <- result$eigen_value_freq
  eigen_vector_freq <- result$eigen_vector_freq
  eigen_vector_time <- result$eigen_vector_time
  
  time_length_eig <- time_length / 2
  Freq_eig <- 2 * pi / time_length * (1:time_length_eig)
  
  eigen_matrix <- eigen_matrix_con(time_length, max(comp_num, comp_num_sta), subject_length,
                                   dmean_smo_fda, q, r, cov_matrix, Freq_eig, time_length_eig)
  
  return(list(
    comp_num = comp_num,
    comp_num_sta = comp_num_sta,
    dmean_smo_fda = dmean_smo_fda,
    smo_fda_mean = smo_fda_mean,
    X_basis = X_basis,
    tra_mat = tra_mat,
    eigen_vector_sta = eigen_vector_sta,
    eigen_vector_time = eigen_vector_time,
    sigma_est = sigma_est,
    true_basis = true_basis,
    eigen_matrix = eigen_matrix,
    dmean_smo_fda = dmean_smo_fda,
    Sep_sta = Sep_sta,
    Sep_dyn = Sep_dyn
  ))
}

## Second-step estimation
sco_est <- function(adjacency_matrix, x_fda, fda, first_step, lambda){
  
  comp_num <- first_step[[1]]
  comp_num_sta <- first_step[[2]]
  dmean_smo_fda <- first_step[[3]]
  smo_fda_mean <- first_step[[4]]
  X_basis <- first_step[[5]]
  tra_mat <- first_step[[6]]
  eigen_vector_sta <- first_step[[7]]
  eigen_vector_time <- first_step[[8]]
  sigma_est <- first_step[[9]]
  true_basis <- first_step[[10]]
  eigen_matrix <- first_step[[11]]
  dmean_smo_fda <- first_step[[12]]
  Sep_sta <- first_step[[13]]
  Sep_dyn <- first_step[[14]]
  rm(first_step)
  
  # Basis 
  subject_length <- dim(fda)[2]
  time_length <- dim(fda)[3]
  q <- time_length ^ 0.4
  r <- floor(q)
  
  ## Eigen-matrix
  time_length_eig <- time_length / 2
  Freq_eig <- 2 * pi / time_length * (1:time_length_eig)
  
  ### Unknown graph constrain
  eigen_matrix_t <- lapply(1:dim(eigen_matrix)[4], function(k){
    gra_selection(eigen_matrix[,,,k], lambda[k])
  })
  eigen_matrix_1 <- array(0, dim(eigen_matrix))
  for(k in 1:dim(eigen_matrix)[4]){
    eigen_matrix_1[,,,k] <- eigen_matrix_t[[k]][[1]]
  }
  
  ### Known graph constrain
  S <- array(0, c(dim(eigen_matrix)[1:2], prod(dim(eigen_matrix)[3:4])))
  for(k in 1:dim(eigen_matrix)[4]){
    S[,,((k - 1) * time_length_eig + 1):(k * time_length_eig)] <- eigen_matrix[,,,k]
  }
  
  eigen_matrix_t <- graph_pre(S, adjacency_matrix, subject_length)
  eigen_matrix_2 <- array(0, dim(eigen_matrix))
  for(i in 1:dim(eigen_matrix)[4]){
    eigen_matrix_2[,,,i] <- eigen_matrix_t[,,((i - 1) * time_length_eig + 1):(i * time_length_eig)] 
  }
  rm(eigen_matrix_t, S)
  
  ## FDA reconstruction
  basis <- eval.basis(X_basis, x_fda)
  basis <- basis %*% solve(tra_mat)
  
  ### Static FPCA
  #### Preparation
  xi <- lapply(1:comp_num_sta, function(k){
    matrix(0, subject_length, time_length)
  })
  
  fre_tra <- lapply(1:comp_num_sta, function(k){
    time_length_k <- time_length
    return(sapply(Freq_eig, function(theta){
      exp(complex(1, 0, 1) * theta * (1:time_length_k)) / sqrt(2 * pi * time_length_k)
    })) 
  })
  
  squ_fre_tra <- lapply(1:comp_num_sta, function(k){
    return(sapply(1:time_length_eig, function(theta){
      a <- fre_tra[[k]][,theta]
      a <- tcrossprod(a, Conj(a))
      a <- (a + t(Conj(a))) / 2
      return(a)
    }, simplify = "array"))
  })
  
  des_mat <- lapply(1:comp_num_sta, function(k){
    sapply(1:time_length, function(j){
      pre_mat <- matrix(0, nrow(eigen_vector_sta), dim(xi[[k]])[2])
      pre_mat[,j] <- eigen_vector_sta[k,]
      pre_mat <- basis %*% pre_mat 
      return(pre_mat)
    }, simplify = "array")
  })
  
  dmean_fda <- sapply(1:time_length, function(j){
    sapply(1:subject_length, function(i){
      c(fda[,i,j] - basis %*% smo_fda_mean[,i])
    }, simplify = "array")
  }, simplify = "array")
  
  squ_des_mat <- lapply(1:comp_num_sta, function(k){
    mat <- sapply(1:time_length, function(j){
      crossprod(des_mat[[k]][,,j])
    }, simplify = "array")
    mat <- apply(mat, c(1, 2), sum)
    return(mat)
  })
  
  #### Score extraction
  xi_sta_sep <- lapply(1:subject_length, function(i){
    lapply(1:Sep_sta[[i]]$comp_num_sta, function(k){
      sapply(1:time_length, function(j){
        dmean_smo_fda[,i,j] %*% Sep_sta[[i]]$eigen_vector_sta[k,]
      })
    })
  }) # Separated FPCA
  
  xi_wg_sta <- lapply(1:comp_num_sta, function(k){
    sapply(1:time_length, function(j){
      t(dmean_smo_fda[,,j]) %*% eigen_vector_sta[k,]
    })
  }) # Weakly-separable FPCA
  
  xi_sta_1 <- score_est(xi_wg_sta, fre_tra, des_mat, dmean_fda, eigen_matrix_1, squ_des_mat,
                        squ_fre_tra, sigma_est, subject_length, time_length, comp_num_sta)
  # GSFPCA
  
  xi_sta_2 <- score_est(xi_wg_sta, fre_tra, des_mat, dmean_fda, eigen_matrix_2, squ_des_mat,
                        squ_fre_tra, sigma_est, subject_length, time_length, comp_num_sta)
  # (G)SFPCA
  
  ### Dynamic FPCA
  #### Preparation
  xi <- lapply(1:comp_num, function(k){
    matrix(0, subject_length, ncol(eigen_vector_time[[k]]) - 1 + time_length)
  })
  
  fre_tra <- lapply(1:comp_num, function(k){
    time_length_k <- ncol(eigen_vector_time[[k]]) - 1 + time_length
    return(sapply(Freq_eig, function(theta){
      exp(complex(1, 0, 1) * theta * (1:time_length_k)) / sqrt(2 * pi * time_length_k)
    })) 
  })
  
  squ_fre_tra <- lapply(1:comp_num, function(k){
    time_length_k <- ncol(eigen_vector_time[[k]]) - 1 + time_length
    return(sapply(1:time_length_eig, function(theta){
      a <- fre_tra[[k]][,theta]
      return(tcrossprod(a, Conj(a)))
    }, simplify = "array"))
  })
  
  des_mat <- lapply(1:comp_num, function(k){
    sapply(1:time_length, function(j){
      pre_mat <- matrix(0, nrow(eigen_vector_time[[k]]), dim(xi[[k]])[2])
      pre_mat[,j:(j-1+dim(eigen_vector_time[[k]])[2])] <- eigen_vector_time[[k]]
      pre_mat <- basis %*% pre_mat 
      return(pre_mat)
    }, simplify = "array")
  })
  
  dmean_fda <- sapply(1:time_length, function(j){
    sapply(1:subject_length, function(i){
      c(fda[,i,j] - basis %*% smo_fda_mean[,i])
    }, simplify = "array")
  }, simplify = "array")
  
  squ_des_mat <- lapply(1:comp_num, function(k){
    mat <- sapply(1:time_length, function(j){
      crossprod(des_mat[[k]][,,j])
    }, simplify = "array")
    mat <- apply(mat, c(1, 2), sum)
    return(mat)
  })
  
  #### Score extraction
  xi_dyn_sep <- lapply(1:subject_length, function(i){
    lapply(1:Sep_dyn[[i]]$comp_num, function(k){
      lag_t <- (ncol(Sep_dyn[[i]]$eigen_vector_time[[k]]) - 1) / 2
      return(sapply(1:(ncol(Sep_dyn[[i]]$eigen_vector_time[[k]]) - 1 + time_length), function(j){
        if(j <= lag_t){
          val <- 0
        }else if(j > time_length + lag_t){
          val <- 0
        }else{
          mat <- cbind(matrix(0, dim(dmean_smo_fda)[1], lag_t),
                       dmean_smo_fda[,i,],
                       matrix(0, dim(dmean_smo_fda)[1], lag_t))
          val <- sum(diag(t(mat[,(j+lag_t):(j-lag_t)]) %*% 
                            Sep_dyn[[i]]$eigen_vector_time[[k]]))
        }
        return(val)
      }))
    })
  }) # Separated FPCA
  
  xi_wg_dyn <- lapply(1:comp_num, function(k){
    lag_t <- (ncol(eigen_vector_time[[k]]) - 1) / 2
    return(sapply(1:(ncol(eigen_vector_time[[k]]) - 1 + time_length), function(j){
      sapply(1:subject_length, function(i){
        if(j <= lag_t){
          val <- 0
        }else if(j > time_length + lag_t){
          val <- 0
        }else{
          mat <- cbind(matrix(0, dim(dmean_smo_fda)[1], lag_t),
                       dmean_smo_fda[,i,],
                       matrix(0, dim(dmean_smo_fda)[1], lag_t))
          val <- sum(diag(t(mat[,(j+lag_t):(j-lag_t)]) %*% 
                            eigen_vector_time[[k]]))
        }
        return(val)
      })
    }))
  }) # Weakly-separable FPCA
  
  xi_dyn_1 <- score_est(xi_wg_dyn, fre_tra, des_mat, dmean_fda, eigen_matrix_1, squ_des_mat,
                        squ_fre_tra, sigma_est, subject_length, time_length, comp_num)
  # GDFPCA
  
  xi_dyn_2 <- score_est(xi_wg_dyn, fre_tra, des_mat, dmean_fda, eigen_matrix_2, squ_des_mat,
                        squ_fre_tra, sigma_est, subject_length, time_length, comp_num)
  # (G)DFPCA
  
  Result <- list(
    score_sta_sep = xi_sta_sep,
    score_sta_wg = xi_wg_sta,
    score_sta = xi_sta_1[[1]],
    score_sta_kg = xi_sta_2[[1]],
    score_dyn_sep = xi_dyn_sep,
    score_dyn_wg = xi_wg_dyn,
    score_dyn = xi_dyn_1[[1]],
    score_dyn_kg = xi_dyn_2[[1]],
    eigen_vector_sta = eigen_vector_sta,
    eigen_vector_dyn = eigen_vector_time,
    exp_tol_mean = smo_fda_mean,
    comp_num = comp_num,
    comp_num_sta = comp_num_sta,
    mea_error = sigma_est,
    basis = true_basis,
    obs_basis = basis,
    dat = fda,
    eigen_matrix = eigen_matrix,
    eigen_matrix_1 = eigen_matrix_1,
    lambda = lambda,
    dmean_smo_fda = dmean_smo_fda,
    Sep_sta = Sep_sta,
    Sep_dyn = Sep_dyn
  )
  
  return(Result)
}

################################################################################################################
# Function to conduct the simulation study
sim_fun <- function(subject_length, time_length, comp_length, covariance_matrix,
                    seed, adjacency_matrix, rho,
                    mean_sample, lag, mea_error_vec, 
                    point_num){
  set.seed(seed * 1000)
  Result <- NULL
  
  while(is.null(Result) == T){
    try(
      {
        x_fda <- seq(0, 1, length.out = point_num + 1)[-c(point_num + 1)]
        score <- score_gen(subject_length, time_length, comp_length, rho, covariance_matrix, lag)
        
        ## Functional data Generation
        fda <- array(NA, c(length(x_fda), subject_length, time_length))
        for(i in 1:subject_length){
          for(j in 1:time_length){
            fda[,i,j] <- fda_gen(x_fda, i, j, mean_sample[[i]], comp_length, score, lag, subject_length) + rnorm(length(x_fda), 0, mea_error_vec[i])
          }
        }
        
        max_comp <- comp_length + 4
        first_step <- first_est(x_fda, fda, mean_zero = T, max_comp, sel_eig = 0)
        lambda <- Sel_lam(first_step[[11]])
        
        Result <- sco_est(adjacency_matrix, x_fda, fda, first_step, lambda)
        
      }, silent = T
    )
  }
  
  Result$mean_sample <- mean_sample
  Result$score <- score
  Result$ture_adjacency_matrix <- adjacency_matrix
  
  return(Result)
}

# Function to calculate the NMSE
res_function <- function(comp_length, time_length, subject_length, Result, lag){
  basis <- Result[[1]]$basis
  grid_time <- seq(0, 1, 0.01)
  
  ture_fda <- sapply(1:time_length, function(j){
    sapply(1:subject_length, function(i){
      sapply(1:length(Result), function(h){
        fda_gen(grid_time, i, j, Result[[h]]$mean_sample[[i]], comp_length, Result[[h]]$score, lag, subject_length)
      }, simplify = "array")
    }, simplify = "array")
  }, simplify = "array")
  
  sum_tol_var <- sapply(1:length(Result), function(h){
    sum(sapply(1:time_length, function(j){
      sapply(1:subject_length, function(i){
        sum(ture_fda[,h,i,j] ^ 2) * (grid_time[2] - grid_time[1])
      })
    }))
  })
  
  res <- sapply(1:time_length, function(j){
    sapply(1:subject_length, function(i){
      sapply(1:comp_length, function(comp_k){
        sapply(1:length(Result), function(h){
          fit_sta_sep <- basis %*% ((rowSums(sapply(1:min(comp_k, Result[[h]]$Sep_sta[[i]]$comp_num_sta), function(k){
            Result[[h]]$Sep_sta[[i]]$eigen_vector_sta[k,] * (Result[[h]]$score_sta_sep[[i]][[k]][j])
          }))) + Result[[h]]$exp_tol_mean[,i])
          fit_wg_sta <- basis %*% ((rowSums(sapply(1:min(comp_k, Result[[h]]$comp_num_sta), function(k){
            Result[[h]]$eigen_vector_sta[k,] * (Result[[h]]$score_sta_wg[[k]][i,j])
          }))) + Result[[h]]$exp_tol_mean[,i])
          fit_sta <- basis %*% ((rowSums(sapply(1:min(comp_k, Result[[h]]$comp_num_sta), function(k){
            Result[[h]]$eigen_vector_sta[k,] * (Result[[h]]$score_sta[[k]][i,j])
          }))) + Result[[h]]$exp_tol_mean[,i])
          fit_kg_sta <- basis %*% ((rowSums(sapply(1:min(comp_k, Result[[h]]$comp_num_sta), function(k){
            Result[[h]]$eigen_vector_sta[k,] * (Result[[h]]$score_sta_kg[[k]][i,j])
          }))) + Result[[h]]$exp_tol_mean[,i])
          
          fit_dyn_sep <- basis %*% ((rowSums(sapply(1:min(comp_k, Result[[h]]$Sep_dyn[[i]]$comp_num), function(k){
            Result[[h]]$Sep_dyn[[i]]$eigen_vector_time[[k]] %*% Result[[h]]$score_dyn_sep[[i]][[k]][j:(j+ncol(Result[[h]]$Sep_dyn[[i]]$eigen_vector_time[[k]])-1)]
          }))) + Result[[h]]$exp_tol_mean[,i])
          fit_wg_dyn <- basis %*% ((rowSums(sapply(1:min(comp_k, Result[[h]]$comp_num), function(k){
            Result[[h]]$eigen_vector_dyn[[k]] %*% (Result[[h]]$score_dyn_wg[[k]][i,j:(j+ncol(Result[[h]]$eigen_vector_dyn[[k]])-1)]) 
          }))) + Result[[h]]$exp_tol_mean[,i])
          fit_dyn <- basis %*% ((rowSums(sapply(1:min(comp_k, Result[[h]]$comp_num), function(k){
            Result[[h]]$eigen_vector_dyn[[k]] %*% (Result[[h]]$score_dyn[[k]][i,j:(j+ncol(Result[[h]]$eigen_vector_dyn[[k]])-1)]) 
          }))) + Result[[h]]$exp_tol_mean[,i])
          fit_kg_dyn <- basis %*% ((rowSums(sapply(1:min(comp_k, Result[[h]]$comp_num), function(k){
            Result[[h]]$eigen_vector_dyn[[k]] %*% (Result[[h]]$score_dyn_kg[[k]][i,j:(j+ncol(Result[[h]]$eigen_vector_dyn[[k]])-1)]) 
          }))) + Result[[h]]$exp_tol_mean[,i])
          return(cbind(fit_sta_sep, fit_wg_sta, fit_sta, fit_kg_sta, 
                       fit_dyn_sep, fit_wg_dyn, fit_dyn, fit_kg_dyn))
        }, simplify = "array")
      }, simplify = "array")
    }, simplify = "array")
  }, simplify = "array")
  
  mse <- sapply(1:comp_length, function(comp_k){
    sapply(1:8, function(k){
      sapply(1:time_length, function(j){
        sapply(1:subject_length, function(i){
          colSums((res[,k,,comp_k,i,j] - ture_fda[,,i,j]) ^ 2) * 
            (grid_time[2] - grid_time[1]) / sum_tol_var
        }, simplify = "array")
      }, simplify = "array")
    }, simplify = "array")
  }, simplify = "array")
  
  mse <- apply(mse, c(1,4,5), sum)
  mse <- apply(mse, c(2,3), mean)
  return(mse)
}
