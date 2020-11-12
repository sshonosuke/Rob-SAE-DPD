# load functions and dataset 
source("RSAE-function.R")
load("simdata.RData")


# selection of alpha 
alpha1 <- select.alpha(Y, X, Di, IR=5)
alpha2 <- select.alpha(Y, X, Di, IR=10)

# robust method
fit1 <- RFH.DPD(Y, X, Di, alpha=alpha1)
fit2 <- RFH.DPD(Y, X, Di, alpha=alpha2)
fit3 <- RFH.DPD(Y, X, Di, alpha=0.2)

# standard FH
EB.FH <- as.vector( eblupFH(Y~X[,-1], vardir=Di)$eblup )

# comparison 
mean( (fit1$EB-Th)^2 )
mean( (fit2$EB-Th)^2 )
mean( (fit3$EB-Th)^2 )
mean( (EB.FH-Th)^2 )


# mse for robust method
mse1 <- mse.RFH.DPD(Y, X, Di, alpha=alpha1, B=200)
mse2 <- mse.RFH.DPD(Y, X, Di, alpha=alpha2, B=200)

mse1$MSE
mse2$MSE
