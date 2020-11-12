###-----------------------------------------------------###
###           Code for data generation                  ###
###-----------------------------------------------------###
set.seed(1)


## settings
m <- 50    # number of areas
G <- 5    # number of groups
DD <- c(0.2, 0.4, 0.6, 0.8, 1)    # group-wise sampling variance 
Di <- sort(rep(DD, m/G))    # area-wise sampling variance 
X <- cbind(rep(1,m), runif(m,0,1))    # covariate matrix
A <- 0.5      # random effect variance 
Beta <- c(0, 2)   # regression coefficients


## data generation
ch <- rbinom(m, 1, 0.15)
V <- (1-ch)*rnorm(m, 0, sqrt(A)) + ch*rnorm(m, 0, 10*sqrt(A))   # random effect
Th <- as.vector(X%*%Beta) + V     # true areal mean
Y <- Th + rnorm(m, 0, sqrt(Di))   # data

plot(X[,2], Y)
abline(Beta)


save(Y, X, Di, Th, file="simdata.RData")
