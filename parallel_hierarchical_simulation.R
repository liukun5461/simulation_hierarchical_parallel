# -----------------------------------------------------------------------------------------
# Simulation Study: Estimate mean and variance of treatment effect using Monte Carlo Method
# -----------------------------------------------------------------------------------------

# Summary:
# This project simulate students which are nested within clusters, their potential
# outcomes Y_ij(0) and Y_ij(1), and estimate the mean and variance of the average 
# treatment effect (ATE) using parallel computing.

# Steps:

# (1) Simulate students for different situations of sample size and ICC:

# mean cluster size m∈{10, 20, 30}, actual cluster size were simulated from a zero-truncated Poisson distribution
# number of clusters k∈{30, 50, 100}
# ICC∈{0.05, 0.25, 0.5}

# (2) The outcome model
# This script performs a simulation to study the variance of the treatment effect 
# (ATE). The outcome is simulated using the Multilevel Model (MLM) below:


# Multilevel Model (MLM) with Random Intercepts and Heterogeneous Treatment
# Effects: 

# The treatment outcome Y_ij for individual i in cluster j is modeled as:

# Level 1 (Individual-Level Model):
#   Y_ij = lambda_0j + zeta_01 * X_1ij + zeta_02 * X_2ij + epsilon_ij
#   where:
#   - Y_ij: Outcome for individual i in cluster j
#   - X_1ij, X_2ij: Individual-level covariates
#   - epsilon_ij ~ N(0, sigma^2): Individual-level error term

# Level 2 (Cluster-Level Model):
#   lambda_0j = gamma_00 + gamma_01 * V_1j + gamma_02 * V_2j + phi_00 * Z_j
#          + omega_01 * V_1j * Z_j + omega_02 * V_2j * Z_j + mu_0j
#   zeta_01 = eta_01 + theta_01 * Z_j
#   zeta_02 = eta_02 + theta_02 * Z_j
#   where:
#   - V_1j, V_2j: Cluster-level covariates
#   - Z_j: Treatment assignment at cluster level (Z_j = 1 for treatment, 0 for control)
#   - mu_0j ~ N(0, tau^2): Cluster-level random effect

# Combined Model:
#   Y_ij = gamma_00 + gamma_01 * V_1j + gamma_02 * V_2j +
#          Z_j * (phi_00 + omega_01 * V_1j + omega_02 * V_2j + theta_01 * X_1ij + theta_02 * X_2ij) +
#          eta_01 * X_1ij + eta_02 * X_2ij + epsilon_ij + mu_0j

# Treatment Effect:
#   - Heterogeneous treatment effect: The effect varies based on cluster (V_1j, V_2j) and individual (X_1ij, X_2ij) covariates.
#   - Each individual has two potential outcomes: Y_ij(0) and Y_ij(1).
#   - Observed outcome: Y_ij = Z_j * Y_ij(1) + (1 - Z_j) * Y_ij(0).

# 

# 
#------------set the coding environment----------------------------------

options(max.print=4000)
options(scipen=999)

rm(list = ls());gc()

if(!require( "pacman")){
  install.packages("pacman")
}

pacman::p_load(
  dplyr,
  ggplot2,
  wakefield,
  tibble,
  tidyr,
  doParallel,
  extraDistr,
  foreach
)


#------------set the parallel running environment------------------------

parallel::detectCores()
n.cores <- parallel::detectCores() 
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

foreach::getDoParWorkers()

#------------set the simulation parameters ------------------------------

## Parameters to vary: k, variance of xi1, ICC, coefficient of X on treatment effect
sim.iteration = 5000
x.sd = 1 # variance of individual level covariates
v.variance = 1 # variance of cluster level covariates
rho = 0.2 # correlation between X1 and X2
treat.int = 2 

tot.error<-1

## empty vectors to store values
ATE.vec <-numeric(sim.iteration)
sd.vec <-numeric(sim.iteration)

#------------run the simulation ------------------------------
sim.results<-foreach(ICC = c(0.05,0.25,0.5), .combine = rbind)%:%
  foreach(k = c(30, 50, 100),  .combine = rbind)%:%
  foreach(m = c(10, 20,30), .combine = rbind, 
          .packages = c("dplyr",
                               "ggplot2",
                               "wakefield",
                               "tibble",
                               "tidyr",
                               "extraDistr"))%dopar% {
                                 for (itr in 1:sim.iteration){
                                 ##simulate cluster size
                                 n = extraDistr::rtpois(k,m,0)
                                 ## ______________________________________________________________________Simulate teacher-level data
                                 TeaID <- wakefield::id(k)
                                 
                                 # cluster level variable
                                 V1 <- rnorm(k,0,v.variance)
                                 V2 <- rnorm(k,0,v.variance)
                                 
                                 # cluster level error term for outcome
                                 mu0 <- rnorm(k,0,sqrt(tot.error*ICC))
                                 # cluster level error term for treatment effect
                                 
                                 # combine data and calculate true PS score
                                 mydata <-data.frame(TeaID,n, V1,V2, mu0)

                                 # define parameters
                                 alpha1 = 0.5
                                 alpha2 = -0.2
                                 
                                 beta1 = 0.5
                                 beta2 = -0.2
                                 
                                 ##  _______________________________________________________________________Simulate student level data
                                 
                                 # student ID
                                 StuID <- wakefield::id(sum(n))
                                 stu.data <-data.frame(StuID)
                                 # merge with cluster level info
                                 tea <- data.frame(TeaID, n)
                                 rep <- rep(tea$TeaID,tea$n)
                                 
                                 stu.data <- data.frame(rep, StuID)
                                 
                                 mydata$TeaID <- as.character(mydata$TeaID)
                                 mydata2  <- dplyr::left_join(mydata, stu.data, by = c( "TeaID" = "rep"))
                                 
                                 # parameters for outcome model;
                                 ## intercept
                                 Gam00 = 0
                                 ## Coefficient of V1 and V2
                                 Gam01 = 0.2
                                 Gam02 = -0.1
                                 ## coefficient of intercept of treatment effect
                                 phi00 = treat.int
                                 ## coefficient of impact of V1 and V2 on treatment effect
                                 omega01 = 0.2
                                 omega02 = -0.1
                                 ## coefficient of impact of X on treatment effect
                                 theta1 = 0.2
                                 theta2 = -0.1
                                 
                                 ## Coefficient of X
                                 eta1 = 0.2
                                 eta2 = -0.1
                                 
                                 mydata2<-as.data.frame(mydata2)
                                 
                                 # student level covariates
                                 
                                 mydata2 <- mydata2 %>%
                                   group_by(StuID)%>%
                                   # generate covariate X
                                   mutate(z1 = rnorm(1,0, x.sd))%>%
                                   mutate(z2 = rnorm(1,0, x.sd))%>%
                                   mutate(X1 = z1)%>%
                                   mutate(X2 = rho*z1 + sqrt(1-rho^2)*z2)%>%
                                   
                                   # generate individual level error term for outcome
                                   mutate(eps = rnorm(1,0,sqrt(tot.error-tot.error*ICC)))%>%
                                   # calculate y0 and y1
                                   mutate(y0 = Gam00+Gam01*V1+Gam02*V2+eta1*X1+eta2*X2+mu0+eps)%>%
                                   mutate(y1 = y0+phi00+omega01*V1+omega02*V2+theta1*X1+theta2*X2)%>%
                                   ungroup()%>%
                                   dplyr::select(-z1,-z2)
                                 mydata2 <- as.data.frame(mydata2)
                                 # student level outcome
                                 sd.vec[itr]<-sd(mydata2$y1-mydata2$y0)
                                 ATE.vec[itr]<-mean(mydata2$y1-mydata2$y0)
                                 }
                                 ate = mean(ATE.vec)
                                 sd = mean(sd.vec)
                                 results = as.data.frame(cbind(ate,sd))%>%
                                   dplyr::mutate(mean_size = m,
                                                 Cluster_num = k,
                                                 ICC_=ICC)
                                 results
                               }

parallel::stopCluster(cl = my.cluster)

sim.results

