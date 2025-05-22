# Multivariate Functional Data Clustering

This repository contains R code for model‐based clustering of univariate functional data. It implements smoothing, descriptive analysis, outlier detection, and two clustering approaches (Funclust and funHDDC) on three socio–economic indicators (average income, proportion of skilled population, and unemployment rate) from four different scenarios.

---
---

##  Dependencies

- **R** (≥ 4.0)  
- **R packages**  
  - Rcpp  
  - fda  
  - funHDDC  
  - robustbase (for fbplot)  

Install required packages:
install.packages(c("Rcpp", "fda", "funHDDC", "robustbase"))
Installation & Compilation
Clone the repository
git clone https://github.com/Michaelmzx/Clustering-analysis-for-cities.git


Compile C++ backends
In R or RStudio:
library(Rcpp)
sourceCpp("src/RfunclustMain.cpp")
sourceCpp("src/mfpca.cpp")

Load R code

source("R/funclust.r")


source("R/cpp_data.R")


source("R/input.R")


source("R/output.R")

 Usage
 
Load and merge data

base_avg <- read.csv("data/baseline_average_income.csv", header = TRUE)

ct_avg   <- read.csv("data/ct_average_income.csv",        header = TRUE)

tc_avg   <- read.csv("data/tc_average_income.csv",        header = TRUE)

ts_avg   <- read.csv("data/ts_average_income.csv",        header = TRUE)

all_avg  <- rbind(base_avg, ct_avg, tc_avg, ts_avg)
.....

Smoothing functional data

t          <- seq_len(ncol(all_avg))

nbasis     <- 35

basisobj   <- create.bspline.basis(rangeval = range(t), nbasis = nbasis)

all_avg_fd <- Data2fd(t, t(all_avg), basisobj = basisobj)

Descriptive plots

plot.fd(all_avg_fd, col = rep(1:4, each = 19),
        xlab = "Time", ylab = "Average Income")

Outlier detection (functional boxplot)

library(robustbase)

fbplot(eval.fd(t, all_avg_fd), x = t, method = "MBD",
       outliercol = "green", barcol = "blue", fullout = TRUE)
       
###
Clustering with Funclust

res_avg_k2 <- funclust(all_avg_fd, K = 2)

print(res_avg_k2$cls)              # cluster assignments

plot(res_avg_k2$loglikTotal, type = "l")

plot.fd(all_avg_fd, col = res_avg_k2$cls)

###
Clustering with funHDDC

library(funHDDC)

hddc_res <- funHDDC(all_avg_fd, K = 2)

table(hddc_res$class)

plot.fd(all_avg_fd, col = hddc_res$class)

For the full workflow and additional plots (P.O.S & unemployment rate; K = 2, 3, 4), see R script.



 License
 
Released under the MIT License. See LICENSE for details.

 Contact
 
For questions or bug reports, open an Issue or contact Minzhen at mmmxi@leeds.ac.uk.
