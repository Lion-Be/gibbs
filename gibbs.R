# -------------------------------
# Gibbs Sampling algorithm in R
# -------------------------------


# Constructing some data

  # Var/Covar matrix in population, first column is y-var
  sigma <- matrix(c(1,    0.25,  0.25,  0.25,   
                      0.25,   1,     0.25,  0.25,    
                      0.25,   0.25,  1,     0.25,     
                      0.25,   0.25,  0.25,   1),nrow=4,ncol=4)
  
  # Mean vector in population
  mu <- rep(0,4)
  
  # Sample from this multivariate normal distribution
  library(MASS)
  data <- as.data.frame(mvrnorm(n = 1000, mu, sigma, empirical=FALSE))


  
# ---------------------------------------------------------
# For Gibbs sampling, we need starting values. Let's take
# parameters of frequentist regression. 
# ---------------------------------------------------------

regression <- function(y, ...) { 
  
  # Constructing Matrix for Dependent Variable
  y.matrix <- as.matrix(y)
  
  # Constructing Matrix for Independent Variables
  predictors <- list(...) 
  x.matrix.initial <- do.call(cbind,predictors) 
  x.matrix <- cbind(x.matrix.initial, rep(1, nrow(x.matrix.initial)))
    
  # Yield ML-Estimates
  
    # Regression Coefficients
    tX.X <-solve(t(x.matrix)%*%x.matrix) 
    X <- tX.X%*%t(x.matrix) 
    coefs <<- X%*%y.matrix 
  
    # Residual Variance
    yhat <- x.matrix%*%coefs
    n = nrow(x.matrix)
    k = ncol(x.matrix)
    sigma2 <<- (sum((y.matrix - yhat)^2))/(n-k)
    
    # Predicted Values
    yhat <<- x.matrix%*%coefs
  
}

# Construct Maximum Likelihood Estimates
with(data, regression(V1,V2,V3))



# ------------------------------------------------
# Running two Gibbs samplers as separate chains
# ------------------------------------------------
  
# Define Variables, Sample Size, Number of Iterations
y <- data$V1
x1 <- data$V2
x2 <- data$V3
n <- length(data)
n.iter <- 10000

# Define whatever (informative) prior distributions
mu <- c(0.143527386, 0.003159547, 0.073368232)
tau <- c(0.005776833, 0.002544000, 0.014750001)
scale <- 0.001
shape <- 0.001

# Construct Matrices in which MCMC Results are stored

  # For First Chain
  mcmc <- matrix(0,n.iter,4)
  names <- c("Intercept","Beta 1", "Beta 2", "Sigma2")
  colnames(mcmc) <- names

  # For Second Chain
  mcmc2 <- matrix(0,n.iter,4)
  names <- c("Intercept","Beta 1", "Beta 2", "Sigma2")
  colnames(mcmc) <- names
  
  # For Posterior Predictive Check
  ppp <- rep(NA, n.iter)


# First Chain: Gibbs Sampling with ML-Estimates as Initial Values
for(i in 1:n.iter){
  
  # Sample from Conditional Posterior of Sigma2
  SSE <- sum( (y-coefs[1] - x1*coefs[2] - x2*coefs[3]) ^2)
  sigma2 <- 1/rgamma(1,n/2+scale,SSE/2+shape)
  
  # Sample from Conditional Posterior of Intercept
  numerator1 <- sum(y-x1*coefs[2]-x2*coefs[3])/sigma2 + mu[1]/tau[1]^2
  denominator1 <- n/sigma2 + 1/tau[1]^2
  coefs[1] <- rnorm(1,numerator1/denominator1,sqrt(1/denominator1))
 
  # Sample from Conditional Posterior of V1
  numerator2 <- sum(x1*(y-coefs[1] - x2*coefs[3]))/sigma2 + mu[2]/tau[2]^2
  denominator2 <- sum(x1^2)/sigma2 + 1/tau[2]^2
  coefs[2] <- rnorm(1,numerator2/denominator2,sqrt(1/denominator2))
  
  # Sample from Conditional Posterior of V2
  numerator3 <- sum(x2*(y-coefs[1] - x1*coefs[2]))/sigma2 + mu[3]/tau[3]^2
  denominator3 <- sum(x2^2)/sigma2 + 1/tau[3]^2
  coefs[3] <- rnorm(1,numerator3/denominator3,sqrt(1/denominator3))
  
  # Store Results in prepared Matrix
  mcmc[i,] <- c(coefs,sigma2)
  
}



# Second Chain: Gibbs Sampling with arbitrary Initial Values

initials <- c(0.01, 0.01, 0.01)
sigma2.2 <- 0.01

# GIBBS Sampling and PPC
for(i in 1:n.iter){
  
  # Sample from Conditional Posterior of Sigma2
  SSE <- sum( (y-initials[1] - x1*initials[2] - x2*initials[3]) ^2)
  sigma2.2 <- 1/rgamma(1,n/2+scale,SSE/2+shape)
  
  # Sample from Conditional Posterior of Intercept
  numerator1 <- sum(y-x1*initials[2]-x2*initials[3])/sigma2.2 + mu[1]/tau[1]^2
  denominator1 <- n/sigma2.2 + 1/tau[1]^2
  initials[1] <- rnorm(1,numerator1/denominator1,sqrt(1/denominator1))
  
  # Sample from Conditional Posterior of V1
  numerator2 <- sum(x1*(y-initials[1] - x2*initials[3]))/sigma2.2 + mu[2]/tau[2]^2
  denominator2 <- sum(x1^2)/sigma2.2 + 1/tau[2]^2
  initials[2] <- rnorm(1,numerator2/denominator2,sqrt(1/denominator2))
  
  # Sample from Conditional Posterior of V2
  numerator3 <- sum(x2*(y-initials[1] - x1*initials[2]))/sigma2.2 + mu[3]/tau[3]^2
  denominator3 <- sum(x2^2)/sigma2.2 + 1/tau[3]^2
  initials[3] <- rnorm(1,numerator3/denominator3,sqrt(1/denominator3))
  
  # Store Results in prepared Matrix
  mcmc2[i,] <- c(initials,sigma2.2)
  
}


# Cut off Burn-in period for all following tasks
burn <- 1000
mcmc <- mcmc[-c(1:burn),]
mcmc2 <- mcmc2[-c(1:burn),]



  

# -------------------------------------------------
# 4. Assessing Autocorrelation and History Plots
# -------------------------------------------------

# History Plots
  
  par(mfrow=c(2,2))
  
  # Intercept
  plot(mcmc[,1], type="l", col="green", main="Intercept")
  lines(mcmc2[,1], type="l", col="red")
  
  # V1
  plot(mcmc[,2], type="l", col="green", main="V1")
  lines(mcmc2[,2], type="l", col="red")
  
  # V2
  plot(mcmc[,3], type="l", col="green", main="V2")
  lines(mcmc2[,3], type="l", col="red")
  
  # Sigma2
  plot(mcmc[,4], type="l", col="green", main="Sigma2")
  lines(mcmc2[,4], type="l", col="red")
  

# Autocorrelation Plots
  
  #### Computation ####
  
  # Create empty vectors to store autocorrelation estimates
  autocorr.int1 <- 0
  autocorr.int2 <- 0
  autocorr.v1.1 <- 0
  autocorr.v1.2 <- 0
  autocorr.v2.1 <- 0
  autocorr.v2.2 <- 0
  autocorr.sigma21 <- 0
  autocorr.sigma22 <- 0
  
  # Intercept - Chain 1
  
    # Centering Values
    centered <- mcmc[,1] - mean(mcmc[,1])             
    # Calculate Variance                                    
    var <- sum(centered * centered) / (n.iter-burn)                   
    # Set first Element to 1
    autocorr.int1[1] <- 1                                 
    # Calculate Autocorrelation and Store Results
    for(i in 1:100){                              
      lag <- centered[-c(1:i)]                        
      lag.alt <- centered[1:(length(lag))]            
      autocorr.int1[i+1] <- (sum(lag.alt * lag) / (n.iter-burn)) / var  
    }
    
  
  # Intercept - Chain 2
    
    # Centering Values
    centered <- mcmc2[,1] - mean(mcmc2[,1])             
    # Calculate Variance                                    
    var <- sum(centered * centered) / n                  
    # Set first Element to 1
    autocorr.int2[1] <- 1                                 
    # Calculate Autocorrelation and Store Results
    for(i in 1:100){                              
      lag <- centered[-c(1:i)]                        
      lag.alt <- centered[1:(length(lag))]            
      autocorr.int2[i+1] <- (sum(lag.alt * lag) / n) / var  
    }
  
   
  # V1 - Chain 1
  
    # Centering Values
    centered <- mcmc[,2] - mean(mcmc[,2])             
    # Calculate Variance                                    
    var <- sum(centered * centered) / (n.iter-burn)                   
    # Set first Element to 1
    autocorr.v1.1[1] <- 1                                 
    # Calculate Autocorrelation and Store Results
    for(i in 1:100){                              
      lag <- centered[-c(1:i)]                        
      lag.alt <- centered[1:(length(lag))]            
      autocorr.v1.1[i+1] <- (sum(lag.alt * lag) / (n.iter-burn)) / var  
    }
  
  
  # V1 - Chain 2
  
    # Centering Values
    centered <- mcmc2[,2] - mean(mcmc2[,2])             
    # Calculate Variance                                    
    var <- sum(centered * centered) / (n.iter-burn)                   
    # Set first Element to 1
    autocorr.v1.2[1] <- 1                                 
    # Calculate Autocorrelation and Store Results
    for(i in 1:100){                              
      lag <- centered[-c(1:i)]                        
      lag.alt <- centered[1:(length(lag))]            
      autocorr.v1.2[i+1] <- (sum(lag.alt * lag) / (n.iter-burn)) / var  
    }
  
    
  # V2 - Chain 1
    
    # Centering Values
    centered <- mcmc[,3] - mean(mcmc[,3])             
    # Calculate Variance                                    
    var <- sum(centered * centered) / (n.iter-burn)                   
    # Set first Element to 1
    autocorr.v2.1[1] <- 1                                 
    # Calculate Autocorrelation and Store Results
    for(i in 1:100){                              
      lag <- centered[-c(1:i)]                        
      lag.alt <- centered[1:(length(lag))]            
      autocorr.v2.1[i+1] <- (sum(lag.alt * lag) / (n.iter-burn)) / var  
    }
    
    
  # V2 - Chain 2
    
    # Centering Values
    centered <- mcmc2[,3] - mean(mcmc2[,3])             
    # Calculate Variance                                    
    var <- sum(centered * centered) / (n.iter-burn)                   
    # Set first Element to 1
    autocorr.v2.2[1] <- 1                                 
    # Calculate Autocorrelation and Store Results
    for(i in 1:100){                              
      lag <- centered[-c(1:i)]                        
      lag.alt <- centered[1:(length(lag))]            
      autocorr.v2.2[i+1] <- (sum(lag.alt * lag) / (n.iter-burn)) / var  
    }
    
   
  # Sigma2 - Chain 1
    
    # Centering Values
    centered <- mcmc[,4] - mean(mcmc[,4])             
    # Calculate Variance                                    
    var <- sum(centered * centered) / (n.iter-burn)                   
    # Set first Element to 1
    autocorr.sigma21[1] <- 1                                 
    # Calculate Autocorrelation and Store Results
    for(i in 1:100){                              
      lag <- centered[-c(1:i)]                        
      lag.alt <- centered[1:(length(lag))]            
      autocorr.sigma21[i+1] <- (sum(lag.alt * lag) / (n.iter-burn)) / var  
    }
    
    
  # Sigma2 - Chain 2
    
    # Centering Values
    centered <- mcmc2[,4] - mean(mcmc2[,4])             
    # Calculate Variance                                    
    var <- sum(centered * centered) / (n.iter-burn)                   
    # Set first Element to 1
    autocorr.sigma22[1] <- 1                                 
    # Calculate Autocorrelation and Store Results
    for(i in 1:100){                              
      lag <- centered[-c(1:i)]                        
      lag.alt <- centered[1:(length(lag))]            
      autocorr.sigma22[i+1] <- (sum(lag.alt * lag) / (n.iter-burn)) / var  
    }
    
   

    #### Plotting ####
    # Plot Autocorrelation for Intercept
    par(mfrow=c(2,2))
    plot(autocorr.int1, type="h", col="red", main="Autocorrelation Intercept") 
    lines(autocorr.int2, type="h", col="green")
    
    # Plot Autocorrelation for V1
    plot(autocorr.v1.1, type="h", col="red", main="Autocorrelation V1") 
    lines(autocorr.v1.2, type="h", col="green")
    
    # Plot Autocorrelation for V2
    plot(autocorr.v2.1, type="h", col="red", main="Autocorrelation V2") 
    lines(autocorr.v2.2, type="h", col="green")
    
    # Plot Autocorrelation for Sigma2
    plot(autocorr.sigma21, type="h", col="red", main="Autocorrelation Sigma2") 
    lines(autocorr.sigma22, type="h", col="green")
  
    
    
# ------------------------
# 5. Evaluating Results
# ------------------------

# Merging both chains together
joint.mcmc <- rbind(mcmc, mcmc2)

# Construct Bayesian Regression Table
results <- matrix(nrow=4, ncol=4)

names2 <- c("Intercept", "V1", "V2", "Sigma2")
rownames(results) <- names2
names3 <- c("Mean", "SD", "2.5 CI", "97.5 CI")
colnames(results) <- names3

results[1,1] <- mean(joint.mcmc[,1])
results[1,2] <- sd(joint.mcmc[,1])
results[1,3] <- mean(joint.mcmc[,1])-1.96*sd(joint.mcmc[,1])
results[1,4] <- mean(joint.mcmc[,1])+1.96*sd(joint.mcmc[,1])

results[2,1] <- mean(joint.mcmc[,2])
results[2,2] <- sd(joint.mcmc[,2])
results[2,3] <- mean(joint.mcmc[,2])-1.96*sd(joint.mcmc[,2])
results[2,4] <- mean(joint.mcmc[,2])+1.96*sd(joint.mcmc[,2])

results[3,1] <- mean(joint.mcmc[,3])
results[3,2] <- sd(joint.mcmc[,3])
results[3,3] <- mean(joint.mcmc[,3])-1.96*sd(joint.mcmc[,3])
results[3,4] <- mean(joint.mcmc[,3])+1.96*sd(joint.mcmc[,3])

results[4,1] <- mean(joint.mcmc[,4])
results[4,2] <- sd(joint.mcmc[,4])
results[4,3] <- mean(joint.mcmc[,4])-1.96*sd(joint.mcmc[,4])
results[4,4] <- mean(joint.mcmc[,4])+1.96*sd(joint.mcmc[,4])

# Display Bayesian Regression Table
results

# Display Posterior Distributions (add titles, axes)
par(mfrow=c(2,2))
plot(density(joint.mcmc[,1]), col="blue", lwd=3, main="Posterior Distribution Intercept")
plot(density(joint.mcmc[,2]), col="blue", lwd=3, main="Posterior Distribution V1")
plot(density(joint.mcmc[,3]), col="blue", lwd=3, main="Posterior Distribution V2")
plot(density(joint.mcmc[,4]), col="blue", lwd=3, main="Posterior Distribution Sigma2")





