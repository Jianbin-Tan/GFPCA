# Setting address
rm(list=ls())
# set working directory to be where the current script is located
mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

# Importing data and function
source("Functions.R")
load("Data/fda_dat.rda")

MFTS <- fda_dat$fda # Multivariate Functional Time Series Data
Lt <- seq(0, 1, length.out = dim(MFTS)[1] + 1)[-(dim(MFTS)[1] + 1)] # Regular time girds on [0,1]

# Running GFPCA
Result_dyn <- GFPCA(MFTS, Lt, Dynamic = T, Max.comp = 8, FVE = 0.8, Mean.zero = F)
Result_sta <- GFPCA(MFTS, Lt, Dynamic = F, Max.comp = 8, FVE = 0.8, Mean.zero = F)

# Remarks for the function "GFPCA" 
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

## Curves' reconstruction
subject_length <- dim(MFTS)[2] # Dimension
time_length <- dim(MFTS)[3] # Temporal Length

fit_fda <- sapply(1:subject_length, function(i){
  sapply(54:60, function(j){

    fit_WSFPCA <- rowSums(sapply(1:Result_sta$comp_num, function(k){
      Result_sta$eigenfunction$val[[k]] * Result_sta$xi_sta_IN[[k]][i,j]
    })) + Result_sta$mean_function$val[,i]

    fit_GSFPCA <- rowSums(sapply(1:Result_sta$comp_num, function(k){
      Result_sta$eigenfunction$val[[k]] * Result_sta$xi_sta_CE[[k]][i,j]
    })) + Result_sta$mean_function$val[,i]

    fit_WDFPCA <- rowSums(sapply(1:Result_dyn$comp_num, function(k){
      Result_dyn$functional_filter$val[[k]] %*% (Result_dyn$xi_dyn_IN[[k]])[i,j:(j+ncol(Result_dyn$functional_filter$val[[k]])-1)]
    })) + Result_dyn$mean_function$val[,i]

    fit_GDFPCA <- rowSums(sapply(1:Result_dyn$comp_num, function(k){
      Result_dyn$functional_filter$val[[k]] %*% (Result_dyn$xi_dyn_CE[[k]])[i,j:(j+ncol(Result_dyn$functional_filter$val[[k]])-1)]
    })) + Result_dyn$mean_function$val[,i]

    return(cbind(fit_WSFPCA, fit_GSFPCA, fit_WDFPCA, fit_GDFPCA))
  }, simplify = "array")
}, simplify = "array")

dimnames(fit_fda)[[1]] <- Result_dyn$mean_function$Lt
dimnames(fit_fda)[[2]] <- c("WSFPCA", "GSFPCA", "WDFPCA", "GDFPCA")
dimnames(fit_fda)[[3]] <- 54:60
dimnames(fit_fda)[[4]] <- fda_dat$loc_inf$name
fit_fda <- fit_fda[,,,c(3, 17, 24)]
fit_fda <- melt(fit_fda)
fit_fda <- cbind(fit_fda, rep(NA, nrow(fit_fda)), rep(NA, nrow(fit_fda)))
colnames(fit_fda) <- c("Hour", "Type", "Day", "Station", "Value", "Observed", "City")

for(k in 1:nrow(fit_fda)){
  z <- which(fit_fda[k,]$Hour == round(Lt, 2))
  if(length(z) > 0){
    fit_fda[k,]$Observed <- MFTS[z,which(fit_fda[k,]$Station == fda_dat$loc_inf$name),fit_fda[k,]$Day]
  }
  fit_fda[k,]$City <- as.character(fda_dat$loc_inf$city)[which(fit_fda[k,]$Station == fda_dat$loc_inf$name)]
}

fit_fda$Type <- factor(fit_fda$Type, levels = c("WSFPCA", "GSFPCA", "WDFPCA", "GDFPCA"))

ggplot(fit_fda) +
  geom_line(aes(x = Hour, y = Value, color = Type, linetype = Type), size = 0.8) +
  geom_point(aes(x = Hour, y = Observed), size = 0.8) +
  facet_wrap(~ City + Day, nrow = 3) +
  labs(x = "Hours", y = "PM2.5 concentrations",
       title = "",
       colour = "", fill = "", linetype = "") +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 24, 6)) +
  scale_color_manual(values = c("blue", "blue", "red", "red")) +
  scale_linetype_manual(values = c(2, 1, 2, 1)) +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        text = element_text(size = 15),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0))

## Illustration of dynamic weak separability
q <- time_length ^ 0.4
r <- floor(q) 
dmean_smo_fda <- Result_dyn$dmean_smo_fda
Freq_eig <- seq(-pi, pi, length.out = 101)[-101]
time_length_eig <- length(Freq_eig)
comp_num_pre <- Result_dyn$comp_num + 2

cov_mat <- test_weaksep(time_length, comp_num_pre, subject_length,
                        dmean_smo_fda, q, r,  
                        Freq_eig, time_length_eig)

cor_mat <- matrix(0, comp_num_pre, comp_num_pre)
for(i in 1:comp_num_pre){
  cor_mat[i,i:comp_num_pre] <- cov_mat[[i]]
}
cor_mat <- (cor_mat + t(cor_mat)) / 2

cor_mat <- diag(diag(cor_mat) ^ (-1/2)) %*% cor_mat %*%diag(diag(cor_mat) ^ (-1/2))
cor_mat <- round(cor_mat, 2)
colnames(cor_mat) <- rownames(cor_mat) <- paste0("Component ", 1:4)
cor_dat <- melt(cor_mat)
cor_dat <- data.frame(cor_dat, values = cor_dat$value)
for(i in 1: nrow(cor_dat)){
  if(cor_dat[i,]$values == 1){
    cor_dat[i,]$values <- NA
  }
}

ggplot(cor_dat) + 
  geom_tile(aes(x = Var1, y = Var2, fill = value)) +
  geom_text(aes(x = Var1, y = Var2, label = values), size = 2.5, colour = "black", alpha = 0.7) +
  scale_fill_gradientn(colors = c('white','#bdc9e1','#74a9cf','#2b8cbe','#045a8d')) +
  labs(x = "", y = "",
       title = "",
       colour = "", fill = "", linetype = "") +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        # text = element_text(size = 15),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) 

