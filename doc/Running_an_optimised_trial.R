## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(out.width = "60%", fig.width=8, fig.height=6, fig.pos = "center",
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(OptimalTrial)
r.1      =   0.3
lambda   =   0.5
lambda.1 =   0.5
w.1      =   0.5
theta_1  =   0.3
theta_2  =   0
mu_1     =   .1
psi_1    =   sqrt(.1)
mu_2     =   0
psi_2    =   sqrt(.1)
rho0     =   .5
alpha    =   0.05
n        = 700
gamma1   = function(X_1, X_2, theta_1, theta_2, lambda, N) {lambda}
gamma2   = function(X_1, X_2, theta_1, theta_2, lambda, N) {1-lambda}
gamma12  = function(X_1, X_2, theta_1, theta_2, lambda, N) {0}
tol.zero =   0.00001
tol.lim  =   0.0001
verbose  = TRUE
use.prior.as.true = FALSE
optim_func = "hjkb"

## ------------------------------------------------------------------------
set.seed(551)
# Simulate 1st stage estimates
Sigma.1     = diag(c(4/lambda.1/n/r.1, 4/(1-lambda.1)/n/r.1))
true.mean = c(theta_1, theta_2)
theta.hat.1 = mvtnorm::rmvnorm(n = 1,
                               mean  = true.mean,
                               sigma = Sigma.1)
cat("First stage estimates are: ", theta.hat.1)

## ------------------------------------------------------------------------
lambda.2 = seq(0 + tol.lim, 1 - tol.lim, length.out = 20)
w.2      = seq(0, 1, length.out = 20)
grid  = expand.grid(lambda.2, w.2)
ngrid = length(w.2)

thetahat.1_1 = theta.hat.1[1]
thetahat.1_2 = theta.hat.1[2] 

results = sapply(X = 1:nrow(grid),
                 FUN = function(x){
                    interimOptimalTrial_EU(lambda.2 = grid[x,1],
                                          w.2 = grid[x,2],
                                          x = c(thetahat.1_1, thetahat.1_2),
                                          r = r.1,
                                          rho = 0,
                                          mu_1 = mu_1, psi_1 = psi_1,
                                          mu_2 = mu_2, psi_2 = psi_2,
                                          rho0 = rho0,
                                          lambda = lambda,
                                          lambda.1 = lambda.1,
                                          n = n,
                                          w.1 = w.1, alpha=alpha,
                                          tol.zero = 1e-10,
                                          gamma1  = function(X_1, X_2, theta_1, theta_2, lambda, N) {lambda}, 
                                          gamma2  = function(X_1, X_2, theta_1, theta_2, lambda, N) {1-lambda},  
                                          gamma12 = function(X_1, X_2, theta_1, theta_2, lambda, N) {0}, 
                                          verbose = FALSE)
                   })
results.df = data.frame(grid,
                        expected_utility = unlist(results[1,]),
                        case = unlist(results[2,]))
results.m = matrix(results.df$expected_utility, nrow = ngrid, ncol = ngrid)
colnames(results.m) = w.2
rownames(results.m) = lambda.2

par(mar = c(4, 5, 1, 2) + 0.1)
filled.contour(x = lambda.2, y = w.2, z = results.m,
  xlab =expression(r[1]^{(2)}), ylab = expression(omega[1]^{(2)}),
  zlim = c(0, .75),
  color.palette =  colorRampPalette(rev(colorspace::sequential_hcl(n = 4, power = 2))))

## ------------------------------------------------------------------------
interimDecision = interimOptimalTrial(x = theta.hat.1,
                                      r = r.1,
                                      sigma_1 = 1, sigma_2 = 1,
                                      rho  = 0,
                                      mu_1 = mu_1, psi_1 = psi_1,
                                      mu_2 = mu_2, psi_2 = psi_1,
                                      rho0 = rho0,
                                      lambda = lambda, lambda.1 = lambda.1,
                                      n   = n,
                                      w.1 = w.1, alpha = alpha,
                                      gamma1 = gamma1,
                                      gamma2 = gamma2,
                                      gamma12 = gamma12,
                                      verbose = FALSE, optim_func = optim_func)

cat("Interim decision: lambda.2", interimDecision$par[1], "w.2", interimDecision$par[2])

par(mar = c(4, 5, 1, 2) + 0.1)
filled.contour(x = lambda.2, y = w.2, z = results.m,
               xlab =expression(r[1]^{(2)}), ylab = expression(omega[1]^{(2)}),
               zlim = c(0, .75),
               color.palette =  colorRampPalette(rev(colorspace::sequential_hcl(n = 4, power = 2))),
               plot.axes = {axis(side = 1, at = seq(0,1,0.1))
                 axis(side = 2, at = seq(0,1,0.1))
                 points(x = interimDecision$par[1], y = interimDecision$par[2], col = "red", pch = 19)
                 abline(v = interimDecision$par[1], col = "red")
                 abline(h = interimDecision$par[2], col = "red")
                 mtext(text = round(interimDecision$par[1], 3), at = c(interimDecision$par[1]), adj = 0, col = "red")
                 mtext(side = 4,
                       text = round(interimDecision$par[2], 3),
                       at = c(interimDecision$par[2]), las = 1, adj = 0, col = "red")
               })

## ------------------------------------------------------------------------
# Store values of lambda.2 and w.2
lambda.2 = interimDecision$par[1]
w.2      = interimDecision$par[2]
# Calculates the variances of the second stage estimates
Sigma.2_1 = 4 / (lambda.2) / n / (1 - r.1)
Sigma.2_2 = 4 / (1 - lambda.2) / n / (1 - r.1)
# And simulate the estimates. A NA pre-allocated vector is used so that
# We only simulate the value if it is actually needed.
theta.hat.2 = c(NA, NA)
if (lambda.2 != 0) theta.hat.2[1] = rnorm(n = 1, mean  = true.mean[1], sd = sqrt(Sigma.2_1))
if (lambda.2 != 1) theta.hat.2[2] = rnorm(n = 1, mean  = true.mean[2], sd = sqrt(Sigma.2_2))
cat("Second stage estimates are: ", theta.hat.2)

## ------------------------------------------------------------------------
# Test the hypotheses given first and second stage data.
decision = decisionOptimalTrial(lambda.2 = lambda.2, w.2 = w.2,
                                x.1 = theta.hat.1,
                                x.2 = theta.hat.2,
                                r = r.1,
                                lambda = lambda, lambda.1 = lambda.1,
                                n = n,
                                w.1 = w.1,
                                alpha = alpha,
                                gamma1 = gamma1,
                                gamma2 = gamma2,
                                gamma12 = gamma12,
                                verbose = verbose)

## ------------------------------------------------------------------------
set.seed(551)
optimal1Trial(r.1      =   0.3,
              lambda   =   0.5,
              lambda.1 =   0.5,
              w.1      =   0.5,
              theta_1  =   0.3,
              theta_2  =   0,
              mu_1     =   0.1,
              psi_1    =   sqrt(.1),
              mu_2     =   0,
              psi_2    =   sqrt(.1),
              rho0     =   0,
              alpha    =   0.05,
              n        = 700,
              gamma1   = ..gamma1..,
              gamma2   = ..gamma2..,
              gamma12  = ..gamma12..,
              tol.zero =   0.00001,
              tol.lim  =   0.0001,
              verbose  = TRUE,
              use.prior.as.true = FALSE,
              optim_func = "hjkb")

