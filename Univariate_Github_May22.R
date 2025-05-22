
Rcpp::sourceCpp("/Users/michael/Downloads/Funclustering_cpp/src/RfunclustMain.cpp")
Rcpp::sourceCpp("/Users/michael/Downloads/Funclustering_cpp/src/mfpca.cpp")
source("/Users/michael/Downloads/Funclustering_cpp/R/funclust.r")
source("/Users/michael/Downloads/Funclustering_cpp/R/cpp_data.R")
source("/Users/michael/Downloads/Funclustering_cpp/R/input.R")
source("/Users/michael/Downloads/Funclustering_cpp/R/output.R")



####
base_avg <- read.csv("/Users/michael/OneDrive - University of Leeds/urban_transition-v.0.2/output/multivariate functional data/baseline/aggregated data/baseline_average_income.csv",
                     header = TRUE)
###
### collapse training
ct_avg <- read.csv("/Users/michael/OneDrive - University of Leeds/urban_transition-v.0.2/output/multivariate functional data/collapse training/ct_average_income.csv",
                   header = TRUE)
###
### total collapse 
tc_avg <- read.csv("/Users/michael/OneDrive - University of Leeds/urban_transition-v.0.2/output/multivariate functional data/total collapse/tc_average_income.csv",
                   header = TRUE)
###
### transition
ts_avg <- read.csv("/Users/michael/OneDrive - University of Leeds/urban_transition-v.0.2/output/multivariate functional data/transition/ts_average_income.csv",
                   header = TRUE)
####merge all the observations into one dataset
all_avg <- rbind(base_avg,ct_avg,tc_avg, ts_avg)
mat_avg <- as.matrix(all_avg)
t <- seq(1,677)
set.seed(123)
nbasis_hddc <- 35  
###
#basisobj <- create.fourier.basis(rangeval = range(t), nbasis = nbasis_hddc)
basisobj <- create.bspline.basis(rangeval = range(t), nbasis = nbasis_hddc)
### Smoothing the average income
all_avg_smooth <- Data2fd(t, t(all_avg), basisobj =basisobj )
col_scenario <- c(rep(1,19),rep(2,19), rep(3,19), rep(4,19))
plot.fd(all_avg_smooth, col=col_scenario, xlab = "time",ylab = "average income")
###
### smoothing the proportion of skilled population
psk_all_smooth <- Data2fd(t, t(psk_all), basisobj =basisobj )
plot.fd(psk_all_smooth, col=col_scenario, xlab = "time",ylab = "p.o.s")
###
### smoothing the unemployment rate
une_all_smooth <- Data2fd(t, t(une_all), basisobj =basisobj )
plot.fd(une_all_smooth, col=col_scenario, xlab = "time",ylab = "unemployment rate")
### Descriptive analysis
### Plot the smoothed data (Plot above)
### The smoothed curves for average income, proportion of skilled population, 
### and unemployment rate have been plotted using the corresponding 
### functional data objects (all_avg_smooth, psk_all_smooth, and une_all_smooth) in the above code.
###
###
### Visualize the mean function per indicator per cluster
all_avg_smooth_matrix <- eval.fd(t, all_avg_smooth)
matplot(all_avg_smooth_matrix,type="l", xlab="t", ylab="average income",,col=col_scen,
        ylim=c(min(all_avg_smooth_matrix), max(all_avg_smooth_matrix)))
mean_avg_baseline <- colMeans(all_avg_smooth_matrix[,1:19])
mean_avg_ct <- colMeans(all_avg_smooth_matrix[,20:38])
mean_avg_tc <- colMeans(all_avg_smooth_matrix[,39:57])
mean_avg_ts <- colMeans(all_avg_smooth_matrix[,58:76])
# Define column index groups
col_groups <- list(
  baseline = 1:19,
  ct = 20:38,
  tc = 39:57,
  ts = 58:76
)

# Compute column means 
mean_avg_list <- lapply(col_groups, function(idx) rowMeans(all_avg_smooth_matrix[, idx]))

# Extract values 
mean_avg_baseline <- mean_avg_list$baseline
mean_avg_ct <- mean_avg_list$ct
mean_avg_tc <- mean_avg_list$tc
mean_avg_ts <- mean_avg_list$ts
mean_avg_scen <- t(cbind(mean_avg_baseline,mean_avg_ct,mean_avg_tc,mean_avg_ts))

### visualize the average income mean per scenario
matplot(t(mean_avg_scen),lty=1,lwd=2, type="l",col=c(1,2,3,4), xlab="t", ylab="average-income")

###
### compute the mean of proportion of skilled population
psk_all_smooth_matrix <- eval.fd(t, psk_all_smooth)
matplot(psk_all_smooth_matrix,type="l", xlab="t", ylab="p.o.s",,col=col_scen,
        ylim=c(min(psk_all_smooth_matrix), max(psk_all_smooth_matrix)))
# Compute column means 
mean_psk_list <- lapply(col_groups, function(idx) rowMeans(psk_all_smooth_matrix[, idx]))

# Extract values 
mean_psk_baseline <- mean_psk_list$baseline
mean_psk_ct <- mean_psk_list$ct
mean_psk_tc <- mean_psk_list$tc
mean_psk_ts <- mean_psk_list$ts
mean_psk_scen <- t(cbind(mean_psk_baseline,mean_psk_ct,mean_psk_tc,mean_psk_ts))

### visualize the proportion of skilled population mean
matplot(t(mean_psk_scen), lty=1,lwd=2,type="l",col=c(1,2,3,4), xlab="t", ylab="p.o.s")

###
###
### compute the mean of unemployment rate 
une_all_smooth_matrix <- eval.fd(t, une_all_smooth)
matplot(une_all_smooth_matrix,type="l", xlab="t", ylab="unemployment rate",,col=col_scen,
        ylim=c(min(une_all_smooth_matrix), max(une_all_smooth_matrix)))

# Compute column means 
mean_une_list <- lapply(col_groups, function(idx) rowMeans(une_all_smooth_matrix[, idx]))

# Extract values 
mean_une_baseline <- mean_une_list$baseline
mean_une_ct <- mean_une_list$ct
mean_une_tc <- mean_une_list$tc
mean_une_ts <- mean_une_list$ts
mean_une_scen <- t(cbind(mean_une_baseline,mean_une_ct,mean_une_tc,mean_une_ts))

### visualize the mean of unemployment rate
matplot(t(mean_une_scen), lty=1,lwd=2, type="l",col=c(1,2,3,4), xlab="t", ylab="unemployment rate")

### Outlier detection
### The "MBD" approach is used to calculate the band depth
total_curves <- ncol(all_avg_smooth_matrix)
group_size <- 19
num_groups <- total_curves %/% group_size
scenario_names <- c("baseline", "collase training", "total collapse", "transition")
### plot the functional box plot for average income per scenario
fbox_avg_list <- list()
for (i in 1:num_groups) {
  start_index <- (i - 1) * group_size + 1
  end_index <- i * group_size
  subset_matrix <- all_avg_smooth_matrix[, start_index:end_index]
  
  fbox <- fbplot(subset_matrix, x = t, method = "MBD",
                 outliercol = "green", barcol = "blue",
                 xlab = "Time", ylab = "Average Income",fullout=TRUE,
                 main=scenario_names[i])
  if(length(fbox$outpoint) > 0) {
    matlines(t, subset_matrix[, fbox$outpoint],
             col = "green", lwd = 2, lty = 1)
  }
}
####
###
###plot the proportion of skilled population
fbox_psk_mbd_outlier <- list()
for (i in 1:num_groups) {
  start_index <- (i - 1) * group_size + 1
  end_index <- i * group_size
  subset_matrix <- psk_all_smooth_matrix[, start_index:end_index]
  
  fbox <- fbplot(subset_matrix, x = t, method = "MBD",
                 outliercol = "green", barcol = "blue",
                 xlab = "Time", ylab = "p.o.s",fullout=TRUE,
                 main=scenario_names[i])
  if(length(fbox$outpoint) > 0) {
    matlines(t, subset_matrix[, fbox$outpoint],
             col = "green", lwd = 2, lty = 1)
  }
}

### plot the functional box plot for unemployment rate
for (i in 1:num_groups) {
  start_index <- (i - 1) * group_size + 1
  end_index <- i * group_size
  subset_matrix <- une_all_smooth_matrix[, start_index:end_index]
  
  fbox <- fbplot(subset_matrix, x = t, method = "MBD",
                 outliercol = "red", outlierlwd = 4.5, barcol = "blue",
                 xlab = "Time", ylab = "unemployment rate",fullout=TRUE,
                 main=scenario_names[i])
  if(length(fbox$outpoint) > 0) {
    matlines(t, subset_matrix[, fbox$outpoint],
             col = "green", lwd = 2, lty = 1)}
  
}
### Model fitting
###
###load the "Funclustering" package for funclust
Rcpp::sourceCpp("/Users/michael/Downloads/Funclustering_cpp/src/Rfunclust_export.cpp")
Rcpp::sourceCpp("/Users/michael/Downloads/Funclustering_cpp/src/mfpca.cpp")
source("/Users/michael/Downloads/Funclustering_cpp/R/funclust.r")
source("/Users/michael/Downloads/Funclustering_cpp/R/cpp_data.R")
source("/Users/michael/Downloads/Funclustering_cpp/R/input.R")
source("/Users/michael/Downloads/Funclustering_cpp/R/output.R")

### For funclust, K=2 case:
res_avg_k2 <- funclust(all_avg_smooth, K=2)
### plot the loglikelihood
plot(y=res_avg_k2$loglikTotal, x=c(1:length(res_avg_k2$loglikTotal)), 
     type="l", xlab="iteration", ylab="log-likelihood")
###check the clustering result
res_avg_k2$cls
### check and visualize the posterior probability
res_avg_k2$tik
boxplot(res_avg_k2$tik)
### plot the clustering result
plot.fd(all_avg_smooth, col=res_avg_k2$cls, xlab="Time", ylab="average income" )

#####
##### proportion of skilled population
res_psk_k2 <- funclust(psk_all_smooth, K=2)
### plot the loglikelihood
plot(y=res_psk_k2$loglikTotal, x=c(1:length(res_psk_k2$loglikTotal)), 
     type="l", xlab="iteration", ylab="log-likelihood")
###check the clustering result
res_psk_k2$cls
### check and visualize the posterior probability
res_psk_k2$tik
boxplot(res_psk_k2$tik)
### plot the clustering result
plot.fd(psk_all_smooth, col=res_psk_k2$cls, xlab="Time", ylab="p.o.s" )

#####
##### unemployment rate
res_une_k2 <- funclust(une_all_smooth, K=2)
### plot the loglikelihood
plot(y=res_une_k2$loglikTotal, x=c(1:length(res_une_k2$loglikTotal)), 
     type="l", xlab="iteration", ylab="log-likelihood")
###check the clustering result
res_une_k2$cls
### check and visualize the posterior probability
res_une_k2$tik
boxplot(res_une_k2$tik)
###
### plot the clustering result
plot.fd(une_all_smooth, col=res_une_k2$cls, xlab="Time", ylab="unemplyoment rate" )

####
####
#### calculate the estimated mean and 95% standard deviation bands
#### of each variables: average income, proportion of skilled population
#### and unemployment rate.
###
### helper functions
compute_cluster_summaries <- function(X, tau) {
  # X: n x T matrix of functional data
  # tau: n x K matrix of posterior probabilities
  
  n <- nrow(X)
  T_len <- ncol(X)
  K <- ncol(tau)
  
  cluster_means <- vector("list", K)
  cluster_vars  <- vector("list", K)
  cluster_weights <- numeric(K)
  
  for (k in 1:K) {
    tau_k <- tau[, k]
    W_k <- sum(tau_k)
    cluster_weights[k] <- W_k
    
    # Weighted mean
    mu_k <- colSums(X * tau_k) / W_k
    
    # Centered data
    centered <- sweep(X, 2, mu_k)
    
    # Weighted variance
    var_k <- colSums(centered^2 * tau_k) / W_k
    
    cluster_means[[k]] <- mu_k
    cluster_vars[[k]]  <- var_k
  }
  cluster_weights <- cluster_weights/sum(cluster_weights)
  return(list(
    means = cluster_means,
    variances = cluster_vars,
    weights = cluster_weights
  ))
}
###
plot_multiple_curves_vs_clusters  <- function(
    X,                    # n × T matrix of curves
    cluster_means,        # list of K mean‐vectors (length T)
    cluster_vars,         # list of K var‐vectors  (length T)
    time = NULL,          # length‐T time grid
    curve_indices = NULL, # integer vector: which rows of X to plot
    cluster_colors = NULL,# length‐K vector of colors for clusters
    curve_colors = NULL,  # length‐|curve_indices| for the curves
    function_value = "Functional value"
) {
  K     <- length(cluster_means)
  T_len <- if (is.null(time)) ncol(X) else length(time)
  if (is.null(time)) time <- seq_len(T_len)
  if (is.null(cluster_colors))
    cluster_colors <- rainbow(K)
  if (!is.null(curve_indices) && is.null(curve_colors))
    curve_colors <- rep("darkorange", length(curve_indices))
  
  # 1) Precompute upper/lower envelopes by index
  upper <- vector("list", K)
  lower <- vector("list", K)
  for (k in seq_len(K)) {
    sd_k       <- sqrt(cluster_vars[[k]])
    upper[[k]] <- cluster_means[[k]] + 1.96 * sd_k
    lower[[k]] <- cluster_means[[k]] - 1.96 * sd_k
  }
  
  # 2) Determine overall y‐range
  if (!is.null(curve_indices)) {
    y_min <- min(c(unlist(lower), X[curve_indices, , drop = FALSE]))
    y_max <- max(c(unlist(upper), X[curve_indices, , drop = FALSE]))
  } else {
    y_min <- min(unlist(lower))
    y_max <- max(unlist(upper))
  }
  
  # 3) Empty base plot
  plot(time, cluster_means[[1]], type = "n",
       ylim = c(y_min, y_max),
       xlab = "Time", ylab = function_value,
       main = if (!is.null(curve_indices)) {
         paste("Curves", paste(curve_indices, collapse = ", "),
               "vs. Cluster Means")
       } else {
         "Cluster Means with Confidence Bands"
       })
  
  # 4) Add each cluster’s band + mean
  for (k in seq_len(K)) {
    polygon(c(time, rev(time)),
            c(upper[[k]], rev(lower[[k]])),
            col   = adjustcolor(cluster_colors[k], alpha.f = 0.2),
            border = NA)
    lines(time, cluster_means[[k]], col = cluster_colors[k], lwd = 4)
  }
  
  # 5) Overlay selected curves if any
  if (!is.null(curve_indices)) {
    for (i in seq_along(curve_indices)) {
      idx <- curve_indices[i]
      lines(time, X[idx, ],
            col = curve_colors[i],
            lwd = 2,
            lty = 1)
    }
    
    # 6) Legend for curves
    legend("bottomleft",
           legend = as.character(curve_indices),
           col    = curve_colors,
           lty    = 1,
           lwd    = 2,
           bg     = "white",
           cex = 0.75,
           y.intersp = 0.5,
           seg.len = 1)
  }
}
####
####
#### Average income
### calculate and plot the estimated mean and standard deviation

### average income
avg_smoothed_mat <- eval.fd(t, all_avg_smooth)
cluster_summ_avg_k2 <- compute_cluster_summaries(avg_smoothed_mat, res_avg_k2$tik)
plot_multiple_curves_vs_clusters(avg_smoothed_mat, cluster_summ_avg_k2$means,
                                 cluster_summ_avg_k2$variances,
                                 cluster_colors=c(1,2),
                                 function_value = "average income")

### proportion of skilled population
psk_smoothed_mat <- eval.fd(t, psk_all_smooth)
cluster_summ_psk_k2 <- compute_cluster_summaries(psk_smoothed_mat, res_psk_k2$tik)
plot_multiple_curves_vs_clusters(psk_smoothed_mat,cluster_summ_psk_k2$means,
                                 cluster_summ_psk_k2$variances,
                                 cluster_colors=c(1,2),
                                 function_value = "p.o.s")
###
### unemployment rate
une_smoothed_mat <- eval.fd(t, une_all_smooth)
cluster_summ_une_k2 <- compute_cluster_summaries(une_smoothed_mat, res_une_k2$tik)
plot_multiple_curves_vs_clusters(une_smoothed_mat,cluster_summ_une_k2$means,
                                 cluster_summ_une_k2$variances,
                                 cluster_colors=c(1,2),
                                 function_value = "unemployment rate")
###
###
### Apply to funHDDC
### average income
funhddc_res_avg <- funHDDC(all_avg_smooth,K=2)
plot.fd(all_avg_smooth, col= funhddc_res_avg$class,ylab = "average income")
cluster_summ_avg_k2_hddc <- compute_cluster_summaries(avg_smoothed_mat, funhddc_res_avg$posterior)
plot_multiple_curves_vs_clusters(avg_smoothed_mat, cluster_summ_avg_k2_hddc$means,
                                 cluster_summ_avg_k2_hddc$variances,
                                 cluster_colors=c(1,2),
                                 function_value = "average income")

###
#### p.o.s
funhddc_res_psk <- funHDDC(psk_all_smooth,K=2)
plot.fd(psk_all_smooth, col =funhddc_res_psk$class, ylab = "p.o.s", xlab="time" )
which(funhddc_res_psk$class ==1)
which(funhddc_res_psk$class ==2)
cluster_summ_psk_k2_hddc <- compute_cluster_summaries(psk_smoothed_mat, funhddc_res_psk$posterior)
plot_multiple_curves_vs_clusters(psk_smoothed_mat, cluster_summ_psk_k2_hddc$means,
                                 cluster_summ_psk_k2_hddc$variances,
                                 cluster_colors=c(1,2),
                                 curve_indices = c(40,50),
                                 function_value = "p.o.s")

### unemployment rate
funhddc_res_une <- funHDDC(une_all_smooth,K=2)
funhddc_res_une$class
plot.fd(une_all_smooth, col=funhddc_res_une$class)
which(funhddc_res_une$class ==2)
which(funhddc_res_une$class ==1)
cluster_sum_une_k2_hddc <- compute_cluster_summaries(t(une_smoothed_mat),funhddc_res_une$posterior)
plot_multiple_curves_vs_clusters(une_smoothed_mat, cluster_sum_une_k2_hddc$means,
                                 cluster_sum_une_k2_hddc$variances,
                                 cluster_colors=c(1,2),
                                 curve_indices = c(38),
                                 function_value = "unemployment rate")

#####
########K=3
set.seed(112233)
###
### funclust
### average income
res_avg_k3 <- funclust(all_avg_smooth, K=3)
### p.o.s
res_psk_k3 <- funclust(psk_all_smooth, K=3)
### unemployment rate
res_une_k3 <- funclust(une_all_smooth, K=3)

### funHDDC
### average income
funhddc_res_avg_k3 <- funHDDC(all_avg_smooth,K=3)
save(funhddc_res_avg_k3, file = "/Users/michael/OneDrive - University of Leeds/Multivariate functional data/Model-based clustering/Rcode for mbclustering/multivariate_simulation_new/smartcity/plot/k=3/funHDDC/funhddc_res_avg_k3.RData")
funhddc_res_avg_k3$class
which(funhddc_res_avg_k3$class == 1)
which(funhddc_res_avg_k3$class == 2)
which(funhddc_res_avg_k3$class == 3)
plot.fd(all_avg_smooth, col=funhddc_res_avg_k3$class, xlab="time",ylab="average income")
### proportion of skilled population
funhddc_res_psk_k3 <- funHDDC(psk_all_smooth,K=3,init = "random")
save(funhddc_res_psk_k3, file = "/Users/michael/OneDrive - University of Leeds/Multivariate functional data/Model-based clustering/Rcode for mbclustering/multivariate_simulation_new/smartcity/plot/k=3/funHDDC/p_o_s/funhddc_res_psk_k3.RData")
funhddc_res_psk_k3$posterior
which(funhddc_res_psk_k3$class == 1)
which(funhddc_res_psk_k3$class == 2)
which(funhddc_res_psk_k3$class == 3)
plot(y=funhddc_res_psk_k3$loglik_all, x= c(1:length(funhddc_res_psk_k3$loglik_all)),
     type="l", ylab="log-likelihood", xlab="iteration")
#####
levels(funhddc_res_psk_k3$class) <- c("3","1","2")
# (Now internally 1→3, 2→1, 3→2)
# If you then need numeric:
funhddc_res_psk_k3$class <- as.integer(as.character(funhddc_res_psk_k3$class))
which(funhddc_res_psk_k3$class == 1)
which(funhddc_res_psk_k3$class == 2)
which(funhddc_res_psk_k3$class == 3)
plot.fd(psk_all_smooth, col= funhddc_res_psk_k3$class, xlab="time",
        ylab="p.o.s")
###
post_psk_k3 <- cbind(funhddc_res_psk_k3$posterior[,2],funhddc_res_psk_k3$posterior[,3], funhddc_res_psk_k3$posterior[,1])
cluster_sum_psk_k3 <- compute_cluster_summaries(t(psk_smoothed_mat),post_psk_k3)
#city 58: originate from transition scenario, but being clustered to
### the cluster consists of collapse training scenario
plot_multiple_curves_vs_clusters(psk_smoothed_mat, cluster_sum_psk_k3$means,
                                 cluster_sum_psk_k3$variances,
                                 cluster_colors=c(1,2,3),
                                 curve_indices = c(58),
                                 function_value = "p.o.s")

# then select those entries strictly between 19 and 39:
### these cities are clustered to the cluster that consist of baseline
### and transition scenario.

mis_psk_index <- psk_k3_class2[
  psk_k3_class2 > 19 & psk_k3_class2 < 39
]
mis_psk_index
plot_multiple_curves_vs_clusters(psk_smoothed_mat, cluster_sum_psk_k3$means,
                                 cluster_sum_psk_k3$variances,
                                 cluster_colors=c(1,2,3),
                                 curve_indices = mis_psk_index[1:5],
                                 curve_colors = 3,
                                 function_value = "p.o.s")

plot_multiple_curves_vs_clusters(psk_smoothed_mat, cluster_sum_psk_k3$means,
                                 cluster_sum_psk_k3$variances,
                                 cluster_colors=c(1,2,3),
                                 curve_indices = mis_psk_index[6:9],
                                 curve_colors = 3,
                                 function_value = "p.o.s")

###unemployment rate
funhddc_res_une_k3 <- funHDDC(une_all_smooth,K=3)
save(funhddc_res_une_k3, file ="/Users/michael/OneDrive - University of Leeds/Multivariate functional data/Model-based clustering/Rcode for mbclustering/multivariate_simulation_new/smartcity/plot/k=3/funHDDC/unemployment rate/funhddc_res_une_k3.RData")
funhddc_res_une_k3$class
which(funhddc_res_une_k3$class == 1)
which(funhddc_res_une_k3$class == 2)
which(funhddc_res_une_k3$class == 3)
### switch the label
funhddc_res_une_k3$class <- with(funhddc_res_une_k3,
                                 ifelse(class == 2, 3,
                                        ifelse(class == 3, 2,
                                               class  # leave 1’s as 1
                                        )
                                 )
)
plot.fd(une_all_smooth, col=funhddc_res_une_k3$class,xlab="time",
        ylab="unemployment rate")

### 58,67
post_une_k3 <- cbind(funhddc_res_une_k3$posterior[,1],
                     funhddc_res_une_k3$posterior[,3],
                     funhddc_res_une_k3$posterior[,2])
cluster_sum_une_k3 <- compute_cluster_summaries(t(une_smoothed_mat),post_une_k3)
###
plot_multiple_curves_vs_clusters(une_smoothed_mat, cluster_sum_une_k3$means,
                                 cluster_sum_une_k3$variances,
                                 cluster_colors=c(1,2,3),
                                 curve_indices = c(58,67),
                                 curve_colors = 1,
                                 function_value = "unemployment rate")
### K=4 : funHDDC
funhddc_res_avg_k4 <- funHDDC(all_avg_smooth,K=4)
funhddc_res_psk_k4 <- funHDDC(psk_all_smooth,K=4)
funhddc_res_une_k4 <- funHDDC(une_all_smooth,K=4)














