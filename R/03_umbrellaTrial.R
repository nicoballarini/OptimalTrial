#' Expected utility for values of lambda.2
#'
#' Simplified version where the utility function does not depend on the estimates or the true
#' values of the parameters and can be taken out of the first integration
#'
#' @param x        the values of the first stage estimates theta1 and theta2
#' @param lambda.2 proportion of patients in sub-population 1 in the trial 2nd stage
#' @param r        proportion of patients in 1st stage of the trial
#' @param sigma_1  standard deviation in subpopulation 1
#' @param sigma_2  standard deviation in subpopulation 2
#' @param rho      correlation between estimates in subpopulation 1 and subpopulation 2 (ignored)
#' @param mu_1     prior mean for the treatment effect in sub-population 1
#' @param psi_1    prior standard deviation for the treatment effect in sub-population 1
#' @param mu_2     prior mean for the treatment effect in sub-population 2
#' @param psi_2    prior standard deviation for the treatment effect in sub-population 2
#' @param rho_0    prior correlation between treatment effect in sub-population 1 and sub-population 2 (ignored)
#' @param lambda   proportion of patients in sub-population 1 in the overal patient population
#' @param lambda.1 proportion of patients in sub-population 1 in the trial 1st stage
#' @param n        trial sample size
#' @param w.1      factor for allocation of alpha to H01 in the testing correction
#' @param alpha    significance level
#' @param gamma1   gain when rejecting H01. a function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param gamma2   gain when rejecting H02. a function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param gamma12  gain when rejecting H01 and H02. a function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param verbose  Print intermediate values for calculations
#' @param tol.zero Absolute distance to 0 or 1 to be considered 0 or 1
#'
#' @export
umbrellaSingleStage_EU <- function(lambda.2,
                                   x = NULL,
                                   r = 0,
                                   sigma_1 = 1, sigma_2 = 1, sigma_3 = 1,
                                   rho = 0,
                                   mu_1 = 0, psi_1 = 1,
                                   mu_2 = 0, psi_2 = 1,
                                   mu_3 = 0, psi_3 = 1, rho0 = 0,
                                   lambda   = c(0.15, 0.05),
                                   lambda.1 = c(0.15, 0.05),
                                   n = 1484,
                                   alpha = c(0.05, 0.05, 0.025),
                                   gamma1  = ..umbrella_gamma1..,  #USES SIMPLE UTILITIES
                                   gamma2  = ..umbrella_gamma2..,  #USES SIMPLE UTILITIES
                                   gamma3  = ..umbrella_gamma3..,  #USES SIMPLE UTILITIES
                                   verbose = TRUE,
                                   tol.zero = 1e-10,
                                   forOpt = TRUE, sig.digits = 4){
  if(r!=0) stop("Single stage only")
  # x is a 3D vector with the estimates from 1st stage.
  # theta.1_1 <- NA
  # theta.1_2 <- NA
  # theta.1_3 <- NA
  mu = c(mu_1, mu_2, mu_3)

  # When lambda.2 and w.2 are very close to 0 or 1, we just set them to 0 or 1 respectively.
  # if (lambda.2 < tol.zero) lambda.2 = 0
  # if (lambda.2 > 1 - tol.zero) lambda.2 = 1
  # if (lambda.2 <= 0) return(0)
  # if (lambda.2 >= 1) return(0)

  # Calculate trial prevalences in 1st Stage
  # lambda.1_1 = lambda.1[1]
  # lambda.1_2 = lambda.1[2]
  # lambda.1_3 = 1 - lambda.1_1 - lambda.1_2
  # Calculate trial prevalences in 2nd Stage
  lambda.2_1 = round(lambda.2[1], sig.digits)
  lambda.2_2 = round(lambda.2[2], sig.digits)
  lambda.2_3 = abs(round(1 - lambda.2[1] - lambda.2[2], sig.digits))
  # Calculate trial prevalences in population
  lambda_1 = lambda[1]
  lambda_2 = lambda[2]
  lambda_3 = 1 - lambda_1 - lambda_2

  # 1st Arm: Conditional expectation given 1st stage ---------------------------
  psi2_1   = psi_1^2
  psi2.1_1 = psi2_1
  sigma2_1 = sigma_1^2
  sigma2_theta.2_1 = 4 * sigma2_1 / (lambda.2_1 * (1 - r) * n)
  # 2st Arm: Conditional expectation given 1st stage ---------------------------
  psi2_2   = psi_2^2
  psi2.1_2 = psi2_2
  sigma2_2 = sigma_2^2
  sigma2_theta.2_2 = 4 * sigma2_2 / (lambda.2_2 * (1 - r) * n)
  # 3rd Arm: Conditional expectation given 1st stage ---------------------------
  psi2_3   = psi_3^2
  psi2.1_3 = psi2_3
  sigma2_3 = sigma_3^2
  sigma2_theta.2_3 = 4 * sigma2_3 / (lambda.2_3 * (1 - r) * n)




  Sigma0  = matrix(c(psi2_1,               rho0 * psi_1 * psi_2, rho0 * psi_1 * psi_3,
                     rho0 * psi_1 * psi_2,               psi2_2, rho0 * psi_2 * psi_3,
                     rho0 * psi_1 * psi_2, rho0 * psi_2 * psi_3,               psi2_3),
                   nrow = 3)
  # Sigma.1_theta  = matrix(c(sigma2_theta.1_1, 0, 0,
  #                           0, sigma2_theta.1_2, 0,
  #                           0, 0, sigma2_theta.1_3),
  #                         nrow = 3, byrow = TRUE)
  Sigma.2_theta  = matrix(c(sigma2_theta.2_1, 0, 0,
                            0, sigma2_theta.2_2, 0,
                            0, 0, sigma2_theta.2_3),
                          nrow = 3)

  # Prior sigma for second stage is the same as prior as there is no first stage
  Sigma.1 = Sigma0

  # Prior Mu for second stage is the same as prior as there is no first stage
  mu.1    = mu
  mu.1_1  = mu.1[1]
  mu.1_2  = mu.1[2]
  mu.1_3  = mu.1[3]

  # Covariance matrix of the prior predictive distribution
  Sigma_ppd = Sigma.1 + Sigma.2_theta

  # Calculate 1st stage p.values
  # P.1_1 = pnorm(Z.1_1, lower.tail = FALSE)
  # P.1_2 = pnorm(Z.1_2, lower.tail = FALSE)
  # P.1_3 = pnorm(Z.1_3, lower.tail = FALSE)

  if (verbose) {
    cat(paste("", "lambda.2: ", paste0(c(lambda.2_1, lambda.2_2, lambda.2_3), collapse = " - "), "\n"))
  }

  # First calculate critical boundaries for first stage
  crit_1 = qnorm(1 - alpha[1])
  crit_2 = qnorm(1 - alpha[2])
  crit_3 = qnorm(1 - alpha[3])
  theta_crit_1 = crit_1 * sqrt(sigma2_theta.2_1)
  theta_crit_2 = crit_2 * sqrt(sigma2_theta.2_2)
  theta_crit_3 = crit_3 * sqrt(sigma2_theta.2_3)

  # Calculate the probabilities of rejection using the prior distribution
  # Now calculate probability of rejection
  Sigma_ppd[Sigma_ppd == Inf] = 1000 # Handling cases where subgroups are not selected
  pr1 <- mvtnorm::pmvnorm(lower = c(theta_crit_1, -Inf, -Inf),
                         upper = c(Inf,           Inf,  Inf),
                         mean  = mu.1,
                         sigma = Sigma_ppd,
                         algorithm = mvtnorm::GenzBretz())[1]
  pr2 <- mvtnorm::pmvnorm(lower = c(-Inf, theta_crit_2, -Inf),
                         upper = c( Inf,          Inf,  Inf),
                         mean  = mu.1,
                         sigma = Sigma_ppd,
                         algorithm = mvtnorm::GenzBretz())[1]
  pr3 <- mvtnorm::pmvnorm(lower = c(-Inf, -Inf, theta_crit_3),
                         upper = c( Inf,  Inf, Inf),
                         mean  = mu.1,
                         sigma = Sigma_ppd,
                         algorithm = mvtnorm::GenzBretz())[1]
  if(verbose){
      cat("Probability of rejection H01",pr1,"\n")
      cat("Probability of rejection H02",pr2,"\n")
      cat("Probability of rejection H03",pr3,"\n")
  }

  # Calculate expected utility. Note this should not depend on the estimates.
  # or true values in the function in this implementation.
  EU <- (gamma1(lambda = lambda_1) * pr1 +
         gamma2(lambda = lambda_2) * pr2 +
         gamma3(lambda = lambda_3) * pr3)

  if(forOpt){
    out = EU
    return(out)
  } else{
    out = c(EU = EU,
            lambda.2_1 = lambda.2_1, lambda.2_2 = lambda.2_2, lambda.2_3 = lambda.2_3,
            pr1 = pr1, pr2 = pr2, pr3 = pr3)
    return(out)
  }
}

#' Find the optimized lambda.2 for the trial
#'
#' Uses the umbrellaSingleStage_EU function to find the values that yield to the maximum utility.
#' Simplified version when the utility functino does not depend on the estimates or the true
#' values of the parameters and can be taken out of the first integration
#'
#' @param x        the values of the first stage estimates theta1 and theta2
#' @param r        proportion of patients in 1st stage of the trial
#' @param sigma_1  standard deviation in subpopulation 1
#' @param sigma_2  standard deviation in subpopulation 2
#' @param rho      correlation between estimates in subpopulation 1 and subpopulation 2 (ignored)
#' @param mu_1     prior mean for the treatment effect in sub-population 1
#' @param psi_1    prior standard deviation for the treatment effect in sub-population 1
#' @param mu_2     prior mean for the treatment effect in sub-population 2
#' @param psi_2    prior standard deviation for the treatment effect in sub-population 2
#' @param rho_0    prior correlation between treatment effect in sub-population 1 and sub-population 2 (ignored)
#' @param lambda   proportion of patients in sub-population 1 in the overal patient population
#' @param lambda.1 proportion of patients in sub-population 1 in the trial 1st stage
#' @param n        trial sample size
#' @param w.1      factor for allocation of alpha to H01 in the testing correction
#' @param alpha    significance level
#' @param gamma1   gain when rejecting H01. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param gamma2   gain when rejecting H02. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param gamma12  gain when rejecting H01 and H02. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param verbose  Print intermediate values for calculations
#' @param tol.zero Absolute distance to 0 or 1 to be considered 0 or 1
#'
#' @export
umbrellaSingleStage = function(x = NULL,
                               r = 0,
                               sigma_1 = 1, sigma_2 = 1, sigma_3 = 1,
                               rho = 0,
                               mu_1 = 0, psi_1 = 1,
                               mu_2 = 0, psi_2 = 1,
                               mu_3 = 0, psi_3 = 1, rho0 = 0,
                               lambda   = c(0.15, 0.05),
                               lambda.1 = c(0.15, 0.05),
                               n = 1484,
                               alpha = c(0.05, 0.05, 0.025),
                               gamma1  = ..umbrella_gamma1..,  #USES SIMPLE UTILITIES
                               gamma2  = ..umbrella_gamma2..,  #USES SIMPLE UTILITIES
                               gamma3  = ..umbrella_gamma3..,  #USES SIMPLE UTILITIES
                               verbose = FALSE,
                               tol.zero = 1e-10,
                               tol.lim = 0.0001,
                               optim_func  = "optim", sig.digits = 4){
  if(r != 0){
    message("r is set to 0 for this function as this is a single Stage Trial")
    r = 0
  }

  constraints = matrix(c(-1, -1,
                          1,  0,
                          0,  1),
                       nrow = 3, byrow = TRUE)

  constants = c(-1, 0, 0)
  initial   = c(0.3, 0.3)

  optimal_values = stats::constrOptim(theta = initial,
                                      f = function(y, ...) umbrellaSingleStage_EU(lambda.2 = y, ...),
                                      ui = constraints,
                                      ci = constants,
                                      control = list(fnscale = -1),
                                      method = "Nelder-Mead",
                                      x = x,
                                      r = r,
                                      sigma_1 = sigma_1, sigma_2 = sigma_2, sigma_3 = sigma_3,
                                      rho = rho,
                                      mu_1 = mu_1, psi_1 = psi_1,
                                      mu_2 = mu_2, psi_2 = psi_2,
                                      mu_3 = mu_3, psi_3 = psi_3, rho0 = rho0,
                                      lambda   = lambda,
                                      lambda.1 = lambda.1,
                                      n = n,
                                      alpha   = alpha,
                                      gamma1  = gamma1,  #USES SIMPLE UTILITIES
                                      gamma2  = gamma2,  #USES SIMPLE UTILITIES
                                      gamma3  = gamma3,  #USES SIMPLE UTILITIES
                                      verbose = verbose,
                                      tol.zero = tol.zero,
                                      forOpt = TRUE, sig.digits = sig.digits)
  optimal_values
}

#' Find the optimized lambda.2 for the trial
#'
#' Uses the umbrellaSingleStage_EU function to find the values that yield to the maximum utility.
#' Simplified version when the utility functino does not depend on the estimates or the true
#' values of the parameters and can be taken out of the first integration
#'
#' @param x        the values of the first stage estimates theta1 and theta2
#' @param r        proportion of patients in 1st stage of the trial
#' @param sigma_1  standard deviation in subpopulation 1
#' @param sigma_2  standard deviation in subpopulation 2
#' @param rho      correlation between estimates in subpopulation 1 and subpopulation 2 (ignored)
#' @param mu_1     prior mean for the treatment effect in sub-population 1
#' @param psi_1    prior standard deviation for the treatment effect in sub-population 1
#' @param mu_2     prior mean for the treatment effect in sub-population 2
#' @param psi_2    prior standard deviation for the treatment effect in sub-population 2
#' @param rho_0    prior correlation between treatment effect in sub-population 1 and sub-population 2 (ignored)
#' @param lambda   proportion of patients in sub-population 1 in the overal patient population
#' @param lambda.1 proportion of patients in sub-population 1 in the trial 1st stage
#' @param n        trial sample size
#' @param w.1      factor for allocation of alpha to H01 in the testing correction
#' @param alpha    significance level
#' @param gamma1   gain when rejecting H01. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param gamma2   gain when rejecting H02. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param gamma12  gain when rejecting H01 and H02. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param verbose  Print intermediate values for calculations
#' @param tol.zero Absolute distance to 0 or 1 to be considered 0 or 1
#'
#' @export
umbrellaSingleStage_grid = function(x = NULL,
                               r = 0,
                               sigma_1 = 1, sigma_2 = 1, sigma_3 = 1,
                               rho = 0,
                               mu_1 = 0, psi_1 = 1,
                               mu_2 = 0, psi_2 = 1,
                               mu_3 = 0, psi_3 = 1, rho0 = 0,
                               lambda   = c(0.15, 0.05),
                               lambda.1 = c(0.15, 0.05),
                               n = 1484,
                               alpha = c(0.05, 0.05, 0.025),
                               gamma1  = ..umbrella_gamma1..,  #USES SIMPLE UTILITIES
                               gamma2  = ..umbrella_gamma2..,  #USES SIMPLE UTILITIES
                               gamma3  = ..umbrella_gamma3..,  #USES SIMPLE UTILITIES
                               verbose = FALSE,
                               tol.zero = 1e-10,
                               tol.lim = 0.0001,
                               optim_func  = "optim", sig.digits = 4){
  if(r != 0){
    message("r is set to 0 for this function as this is a single Stage Trial")
    r = 0
  }
  l1 = round(seq(0,1,0.01), 2)
  l2 = round(seq(0,1,0.01), 2)
  grid = expand.grid(lambda_1 = l1, lambda_2 = l2)
  grid = grid[grid$lambda_1 + grid$lambda_2 <= 1, ]
  f.out = apply(grid, MARGIN = 1,
                function(y, ...) umbrellaSingleStage_EU(lambda.2 = y, ...),
                x = x,
                r = r,
                sigma_1 = sigma_1, sigma_2 = sigma_2, sigma_3 = sigma_3,
                rho = rho,
                mu_1 = mu_1, psi_1 = psi_1,
                mu_2 = mu_2, psi_2 = psi_2,
                mu_3 = mu_3, psi_3 = psi_3, rho0 = rho0,
                lambda   = lambda,
                lambda.1 = lambda.1,
                n = n,
                alpha   = alpha,
                gamma1  = gamma1,  #USES SIMPLE UTILITIES
                gamma2  = gamma2,  #USES SIMPLE UTILITIES
                gamma3  = gamma3,  #USES SIMPLE UTILITIES
                verbose = verbose,
                tol.zero = tol.zero,
                forOpt = TRUE, sig.digits = sig.digits)
  gridout = data.frame(grid, EU = f.out)
  gridout
}
