#' Expected utility for values of lambda.2 and w.2 after 1st stage
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
interimOptimalTrial_EU <- function(lambda.2, w.2,
                                   x, r,
                                   sigma_1 = 1, sigma_2 = 1, rho = 0,
                                   mu_1 = 0, psi_1 = 1,
                                   mu_2 = 0, psi_2 = 1, rho0 = 0,
                                   lambda = 0.5, lambda.1 = 0.5,
                                   n = 100,
                                   w.1 = 0.5,
                                   alpha = 0.025,
                                   gamma1  = ..gamma1..,  #USES SIMPLE UTILITIES
                                   gamma2  = ..gamma2..,  #USES SIMPLE UTILITIES
                                   gamma12 = ..gamma12.., #USES SIMPLE UTILITIES
                                   verbose = TRUE, show.case = TRUE, tol.zero = 1e-10){

  # x is a 2D vector with the estimates from 1st stage.
  theta.1_1 <- (x[1])
  theta.1_2 <- (x[2])
  mu = c(mu_1, mu_2)
  # When lambda.2 and w.2 are very close to 0 or 1, we just set them to 0 or 1 respectively.
  # if (lambda.2 < tol.zero) lambda.2 = 0
  # if (lambda.2 > 1 - tol.zero) lambda.2 = 1
  # if (w.2 < tol.zero)      w.2 = 0
  # if (w.2 > 1 - tol.zero)  w.2 = 1
  # if (w.2 < 0)  return(0)
  # if (w.2 > 1)  return(0)
  if (w.2 < 0)  {
    # cat("\nWarning! w.2<0 (", w.2, ") and will be set to 0\n\n")
    w.2=0
  }
  if (w.2 > 1)  {
    # cat("\nWarning! w.2>1 (", w.2, ") and will be set to 1\n\n")
    w.2 = 1
  }

  # if (lambda.2 <= 0) return(0)
  # if (lambda.2 >= 1) return(0)
  # Calculate trial prevalences
  lambda.1_1 = lambda.1
  lambda.1_2 = (1 - lambda.1)
  lambda.2_1 = lambda.2
  lambda.2_2 = (1 - lambda.2)

  # Calculate trial weights
  w.1_1 = w.1
  w.1_2 = (1-w.1)
  w.2_1 = w.2
  w.2_2 = (1-w.2)

  # Subgroup: Conditional expectation given 1st stage
  if (lambda.1_1 != 0){ # Case when subgroup 1 trial prevalence is > 0 in first stage
    psi2_1 = psi_1^2
    sigma2_1 = sigma_1^2
    psi2.1_1 = sigma2_1 * psi2_1 / (lambda.1_1 * r * n * psi2_1 / 4 + sigma2_1)
    sigma2_theta.1_1 = 4 * sigma2_1 / (lambda.1_1 * r * n)
    sigma2_theta.2_1_fixed = 4 * sigma2_1 / (lambda.1_1 * (1-r) * n)
    sigma2_theta.2_1       = 4 * sigma2_1 / (lambda.2_1 * (1-r) * n)
    sigma2_theta.p_1_fixed = 4 * sigma2_1 / n / (lambda.1_1 * r + lambda.1_1 * (1 - r))
    sigma2_theta.p_1       = 4 * sigma2_1 / n / (lambda.1_1 * r + lambda.2_1 * (1 - r))
    mu.1_1 = psi2.1_1 * (mu_1 / psi2_1 + theta.1_1 / sigma2_theta.1_1)
    Z.1_1 = theta.1_1 / (2 * sigma_1 / sqrt(n * lambda.1_1 * r))
    f.1_1 = lambda.1_1 * r / (lambda.1_1 * r  + lambda.2_1 * (1-r))
    f.2_1 = 1 - f.1_1
  } else { # If subgroup 1 trial prevalence is 0 in first stage, we only have 2nd stage data
    psi2_1 = psi_1^2
    sigma2_1 = sigma_1^2
    psi2.1_1 = psi2_1
    sigma2_theta.1_1 = 0
    sigma2_theta.2_1 = 4 * sigma2_1 / (lambda.2_1 * (1 - r) * n)
    sigma2_theta.p_1 = sigma2_theta.2_1
    mu.1_1 = mu_1
    Z.1_1 = NA
    f.1_1 = 0
    f.2_1 = 1
  }
  if (lambda.2_1 == 0){ # If subgroup 1 trial prevalence is 0 in second stage, we only have 1st stage data
    sigma2_theta.2_1 = 0
  }

  Etheta.p_1_cond = r * theta.1_1
  Vtheta.p_1_cond = (1 - r)^2 * (sigma2_theta.2_1_fixed)

  # Complement: Conditional expectation given 2nd stage
  if (lambda.1_2 != 0){# Case when subgroup 2 trial prevalence is >0 in first stage
    sigma2_2 = sigma_2^2
    psi2_2 = psi_2^2
    psi2.1_2 = sigma2_2 * psi2_2 / (lambda.1_2 * r * n * psi2_2 / 4 + sigma2_2)
    sigma2_theta.1_2 =       4 * sigma2_2 / (lambda.1_2 * r * n)
    sigma2_theta.2_2_fixed = 4 * sigma2_2 / (lambda.1_2 * (1 - r) * n)
    sigma2_theta.2_2 =       4 * sigma2_2 / (lambda.2_2 * (1 - r) * n)
    sigma2_theta.p_2_fixed = 4 * sigma_2^2 / n / (lambda.1_2 * r + lambda.1_2 * (1 - r))
    sigma2_theta.p_2 =       4 * sigma_2^2 / n / (lambda.1_2 * r + lambda.2_2 * (1 - r))
    mu.1_2 = psi2.1_2 * (mu_2 / psi2_2 + theta.1_2 / sigma2_theta.1_2)
    Z.1_2 = theta.1_2 / (2 * sigma_2 / sqrt(n * lambda.1_2 * r))
    f.1_2 = lambda.1_2 * r / (lambda.1_2 * r  + lambda.2_2 * (1-r))
    f.2_2 = 1 - f.1_2
  } else { # If subgroup 2 trial prevalence is 0 in first stage, we only have 2nd stage data
    psi2_2 = psi_2^2
    sigma2_2 = sigma_2^2
    psi2.1_2 = psi2_2
    sigma2_theta.1_2 = 0
    sigma2_theta.2_2 = 4 * sigma2_2 / (lambda.2_2 * (1 - r) * n)
    sigma2_theta.p_2 = sigma2_theta.2_2
    mu.1_2 = mu_2
    Z.1_2 = NA
    f.1_2 = 0
    f.2_2 = 1
  }
  if (lambda.2_2 == 0){# If subgroup 2 trial prevalence is 0 in second stage, we only have 1st stage data
    sigma2_theta.2_2 = 0
  }

  Etheta.p_2_cond = r * theta.1_2
  Vtheta.p_2_cond = (1 - r)^2 * (sigma2_theta.2_2_fixed)

  Sigma0  = matrix(c(psi2_1,
                     rho0 * psi_1 * psi_2,
                     rho0 * psi_1 * psi_2,
                     psi2_2),
                   nrow = 2)
  Sigma.1_theta  = matrix(c(sigma2_theta.1_1,
                           0 , 0,
                           sigma2_theta.1_2),
                         nrow = 2)
  Sigma.2_theta  = matrix(c(sigma2_theta.2_1,
                            0 , 0,
                            sigma2_theta.2_2),
                         nrow = 2)

  SigmaPlusInv = solve(Sigma0 + Sigma.1_theta)
  Sigma.1 = Sigma0 %*% SigmaPlusInv %*% Sigma.1_theta
  mu.1    = drop(Sigma0 %*% SigmaPlusInv %*% drop(x) + Sigma.1_theta %*% SigmaPlusInv %*% mu)
  mu.1_1  = mu.1[1]
  mu.1_2  = mu.1[2]

  # Covariance matrix of the prior predictive distribution
  Sigma_ppd = Sigma.1 + Sigma.2_theta

  # Calculate 1st stage p.values
  P.1_1 = pnorm(Z.1_1, lower.tail = FALSE)
  P.1_2 = pnorm(Z.1_2, lower.tail = FALSE)

  if (verbose) {
    cat(paste("Essential values", "lambda.2", lambda.2, "w.2", w.2, "\n"))
    cat(paste("First Stage estimates", theta.1_1, theta.1_2, "\n"))
    cat(paste("First Stage Z values", Z.1_1, Z.1_2, "\n"))
  }

  # Calculate A_1, A_2, B_1, B_2
  # First calculate critical boundaries for first stage
  crit_w.1_1 = qnorm(1 - alpha * w.1_1)
  crit_w.1_2 = qnorm(1 - alpha * (1-w.1_1))
  crit       = qnorm(1 - alpha)
  # Based on estiamtes from first stage we derive the distributions of
  # the conditional statistics
  # Compute conditional error rates
  A12_1 = pnorm(crit_w.1_1 * sqrt(sigma2_theta.p_1_fixed),
             mean = Etheta.p_1_cond, sd = sqrt(Vtheta.p_1_cond), lower.tail = FALSE)
  A12_2 = pnorm(crit_w.1_2 * sqrt(sigma2_theta.p_2_fixed),
             mean = Etheta.p_2_cond, sd = sqrt(Vtheta.p_2_cond), lower.tail = FALSE)
  A12 = A12_1 + A12_2 - A12_1*A12_2

  A1 = pnorm(crit * sqrt(sigma2_theta.p_1_fixed),
             mean = Etheta.p_1_cond, sd = sqrt(Vtheta.p_1_cond), lower.tail = FALSE)
  A2 = pnorm(crit * sqrt(sigma2_theta.p_2_fixed),
             mean = Etheta.p_2_cond, sd = sqrt(Vtheta.p_2_cond), lower.tail = FALSE)



  # Use predictive distribution of theta.2_1 given theta.1_1
  z_A1 = qnorm(A1, lower.tail = FALSE)
  z_A2 = qnorm(A2, lower.tail = FALSE)
  z_A12_1 = qnorm(w.2_1 * A12, lower.tail = FALSE)
  z_A12_2 = qnorm(w.2_2 * A12, lower.tail = FALSE)

  crit_A1 = z_A1 * sqrt(sigma2_theta.2_1)
  crit_A2 = z_A2 * sqrt(sigma2_theta.2_2)
  crit_A12_1 = z_A12_1 * sqrt(sigma2_theta.2_1)
  crit_A12_2 = z_A12_2 * sqrt(sigma2_theta.2_2)

  # Calculate the probabilities of rejection using the prior distribution
  case = 5
  # Check whether we have a two stage trial of one stage (r=1)
  isTwoStageTrial = (r != 1)
  if (isTwoStageTrial == FALSE){
    # If we have a one stage trials, simply test the hypotheses using 1st stage p-vaues
    if (lambda.1_1 == 0) P.1_1 = 1
    if (lambda.1_2 == 0) P.1_2 = 1
    pr1  <- 1*(P.1_1 < alpha & (P.1_1 < w.2_1 * alpha | P.1_1 < w.2_2 * alpha))
    pr2  <- 1*(P.1_2 < alpha & (P.1_1 < w.2_1 * alpha | P.1_2 < w.2_2 * alpha))
    pr12 <- pr1 * pr2
  }
  # Case: 0 < lambda.2_1 < 1 -----------------------------------------------
  if (isTwoStageTrial & (0 < lambda.2_1 & lambda.2_1 < 1)){
    # Now calculate probability of rejection
    q1 <- mvtnorm::pmvnorm(lower = c(-Inf,                   max(crit_A12_2, crit_A2)),
                           upper = c(min(crit_A1, crit_A12_1),  Inf),
                           mean  = mu.1,
                           sigma = Sigma_ppd,
                           algorithm = mvtnorm::GenzBretz())[1]
    q2 <- mvtnorm::pmvnorm(lower = c(min(crit_A1, crit_A12_1), max(crit_A12_2, crit_A2)),
                           upper = c(max(crit_A1, crit_A12_1), Inf),
                           mean  = mu.1,
                           sigma = Sigma_ppd,
                           algorithm = mvtnorm::GenzBretz())[1]
    q3 <- mvtnorm::pmvnorm(lower = c(max(crit_A1, crit_A12_1), max(crit_A12_2, crit_A2)),
                           upper = c(Inf,                   Inf),
                           mean  = mu.1,
                           sigma = Sigma_ppd,
                           algorithm = mvtnorm::GenzBretz())[1]
    q4 <- mvtnorm::pmvnorm(lower = c(max(crit_A1, crit_A12_1), min(crit_A12_2, crit_A2)),
                           upper = c(Inf,                   max(crit_A12_2, crit_A2)),
                           mean  = mu.1,
                           sigma = Sigma_ppd,
                           algorithm = mvtnorm::GenzBretz())[1]
    q5 <- mvtnorm::pmvnorm(lower = c(max(crit_A1, crit_A12_1),   -Inf),
                           upper = c(Inf,                      min(crit_A12_2, crit_A2)),
                           mean  = mu.1,
                           sigma = Sigma_ppd,
                           algorithm = mvtnorm::GenzBretz())[1]
    q6 <- mvtnorm::pmvnorm(lower = c(min(crit_A1, crit_A12_1), min(crit_A12_2, crit_A2)),
                           upper = c(max(crit_A1, crit_A12_1), max(crit_A12_2, crit_A2)),
                           mean  = mu.1,
                           sigma = Sigma_ppd,
                           algorithm = mvtnorm::GenzBretz())[1]

    if        (crit_A1 <= crit_A12_1 & crit_A2 <= crit_A12_2){
      pr1  <- q2 + q3 + q4 + q5
      pr2  <- q1 + q2 + q3 + q4
      pr12 <- q2 + q3 + q4
      case <- 1
    } else if (crit_A1 >  crit_A12_1 & crit_A2 <= crit_A12_2){
      pr1  <- q3 + q4 + q5
      pr2  <- q1 + q2 + q3 + q4 + q6
      pr12 <- q3 + q4
      case <- 2
    } else if (crit_A1 <= crit_A12_1 & crit_A2  > crit_A12_2){
      pr1  <- q2 + q3 + q4 + q5 + q6
      pr2  <- q1 + q2 + q3
      pr12 <- q2 + q3
      case <- 3
    } else if (crit_A1  > crit_A12_1 & crit_A2  > crit_A12_2){
      pr1  <- q3 + q4 + q5
      pr2  <- q1 + q2 + q3
      pr12 <- q3
      case <- 4
    }

    if(verbose){
      cat("q1:",q1,"q2:",q2,"q3:",q3,"q4:",q4,"q5:",q5,"q6:",q6, "\n")
      cat("Probability of rejection H01",pr1,"\n")
      cat("Probability of rejection H02",pr2,"\n")
      cat("Probability of rejection H01^H02",pr12,"\n")
      cat("Power = ",(pr1 + pr2 - pr12), "\n")
    }
  }

  # Calculate expected utility. Note this should not depend on the estimates.
  # or true values in the function in this implementation.
  EU <- (gamma1(lambda = lambda)  * pr1 +
         gamma2(lambda = lambda)  * pr2 +
         gamma12(lambda = lambda) * pr12)
  if(show.case) {
    out = list(EU = EU, case = case)
  } else {
    out = EU
  }
  out
}

#' Find the optimized w.2 and lambda.2 for the second stage after 1st stage
#'
#' Uses the twoStage_optimizedWeight function to find the values that yield to the maximum utility.
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
interimOptimalTrial = function(x, r,
                               sigma_1 = 1, sigma_2 = 1, rho = 0,
                               mu_1 = 0, psi_1 = 1,
                               mu_2 = 0, psi_2 = 1, rho0 = 0,
                               lambda = 0.5, lambda.1 = 0.5,
                               n = 100,
                               w.1 = 0.5,
                               alpha = 0.025,
                               gamma1  = ..gamma1..,  #USES SIMPLE UTILITIES
                               gamma2  = ..gamma2..,  #USES SIMPLE UTILITIES
                               gamma12 = ..gamma12.., #USES SIMPLE UTILITIES
                               verbose = TRUE,
                               tol.zero = 0.00001,
                               tol.lim = 0.0001,
                               optim_func  = "optim",
                               grid_points = 10, starting_points = c(0.5,0.5)){
  # optim_func = match.arg(optim_func)
  # Use the optim function to find the values lambda.2 w.2  with maximum utility
  if(optim_func == "optim"){
  out = stats::optim(par = starting_points,
                fn = function(y, ...) interimOptimalTrial_EU(lambda.2 = y[1], w.2 = y[2], ...),
                x = x, r = r,
                sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                mu_1 = mu_1, psi_1 = psi_1,
                mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                lambda = lambda, lambda.1 = lambda.1,
                n = n,
                w.1 = w.1,
                alpha = alpha,
                gamma1  = gamma1,
                gamma2  = gamma2,
                gamma12 = gamma12,
                verbose = verbose, show.case = FALSE,
                tol.zero = tol.zero,
                method = "L-BFGS-B",
                lower = c(tol.lim,     0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                upper = c(1 - tol.lim, 1),  # Constrain the search to the [0,1]X[0,1] space
                # lower = c(0, 0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                # upper = c(1, 1),  # Constrain the search to the [0,1]X[0,1] space
                control = list(fnscale = -1))
    # if (out$par[2] < tol.zero)      out$par[2] = 0
    # if (out$par[2] > 1 - tol.zero)  out$par[2] = 1
    if (out$par[2] < 0)  out$par[2] = 0
    if (out$par[2] > 1)  out$par[2] = 1
  } else if(optim_func == "hjkb"){
    out = dfoptim::hjkb(par = starting_points,
                        fn = function(y, ...) interimOptimalTrial_EU(lambda.2 = y[1], w.2 = y[2], ...),
                        x = x, r = r,
                        sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                        mu_1 = mu_1, psi_1 = psi_1,
                        mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                        lambda = lambda, lambda.1 = lambda.1,
                        n = n,
                        w.1 = w.1,
                        alpha = alpha,
                        gamma1  = gamma1,
                        gamma2  = gamma2,
                        gamma12 = gamma12,
                        verbose = verbose, show.case = FALSE, tol.zero = tol.zero,
                        lower = c(tol.lim,     0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                        upper = c(1 - tol.lim, 1),  # Constrain the search to the [0,1]X[0,1] space
                        control = list(maximize = TRUE))
  } else if(optim_func == "Rvmmin"){
    out = Rvmmin::Rvmmin(par = starting_points,
                         fn = function(y, ...) (-1)*interimOptimalTrial_EU(lambda.2 = y[1], w.2 = y[2], ...),
                         x = x, r = r,
                         sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                         mu_1 = mu_1, psi_1 = psi_1,
                         mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                         lambda = lambda, lambda.1 = lambda.1,
                         n = n,
                         w.1 = w.1,
                         alpha = alpha,
                         gamma1  = gamma1,
                         gamma2  = gamma2,
                         gamma12 = gamma12,
                         verbose = verbose, show.case = FALSE, tol.zero = tol.zero,
                         lower = c(tol.lim,     0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                         upper = c(1 - tol.lim, 1))
  } else if(optim_func == "grid"){
    out = grid_search(x = x, r = r,
                      sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                      mu_1 = mu_1, psi_1 = psi_1,
                      mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                      lambda = lambda, lambda.1 = lambda.1,
                      n = n,
                      w.1 = w.1,
                      alpha = alpha,
                      gamma1  = gamma1,
                      gamma2  = gamma2,
                      gamma12 = gamma12,
                      verbose = verbose,
                      tol.zero = tol.zero,
                      lower = c(tol.lim,     0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                      upper = c(1 - tol.lim, 1),
                      grid_points = grid_points)
  } else if(optim_func == "fast"){
    out = dfoptim::hjkb(par = c(0.5,0.5),
                        fn = function(y, ...) fast_interimOptimalTrial_EU(lambda.2_1 = y[1], w.2_1 = y[2], ...),
                        x = x,
                        theta.1_1 = x[1], theta.1_2 = x[2],
                        r = r,
                        sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                        mu = c(mu_1, mu_2),
                        mu_1 = mu_1, psi_1 = psi_1,
                        mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                        Sigma0  = matrix(c(psi_1^2,
                                           rho0 * psi_1 * psi_2,
                                           rho0 * psi_1 * psi_2,
                                           psi_2^2),
                                         nrow = 2),
                        lambda = lambda,
                        lambda.1_1 = lambda.1,
                        lambda.1_2 = 1 - lambda.1,
                        n = n,
                        w.1_1 = w.1,
                        w.1_2 = 1-w.1,
                        alpha = alpha,
                        gamma1  = gamma1,
                        gamma2  = gamma2,
                        gamma12 = gamma12,
                        verbose = verbose,
                        tol.zero = tol.zero,
                        lower = c(tol.lim,     0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                        upper = c(1 - tol.lim, 1),  # Constrain the search to the [0,1]X[0,1] space
                        control = list(maximize = TRUE))
  } else {
    out = dfoptim::nmkb(par = c(0.5,0.5),
                        fn = function(y, ...) interimOptimalTrial_EU(lambda.2 = y[1], w.2 = y[2], ...),
                        x = x, r = r,
                        sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                        mu_1 = mu_1, psi_1 = psi_1,
                        mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                        lambda = lambda, lambda.1 = lambda.1,
                        n = n,
                        w.1 = w.1,
                        alpha = alpha,
                        gamma1  = gamma1,
                        gamma2  = gamma2,
                        gamma12 = gamma12,
                        verbose = verbose, show.case = FALSE, tol.zero = tol.zero,
                        lower = c(tol.lim,     0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                        upper = c(1 - tol.lim, 1),  # Constrain the search to the [0,1]X[0,1] space
                        control = list(maximize = TRUE))
  }
  out
}

grid_search <- function(lower = c(0,0), upper = c(1,1),
                        grid_points = 20,
                        x, r,
                        sigma_1 = 1, sigma_2 = 1, rho = 0,
                        mu_1 = 0, psi_1 = 1,
                        mu_2 = 0, psi_2 = 1, rho0 = 0,
                        lambda = 0.5, lambda.1 = 0.5,
                        n = 100,
                        w.1 = 0.5,
                        alpha = 0.025,
                        gamma1  = ..gamma1..,  #USES SIMPLE UTILITIES
                        gamma2  = ..gamma2..,  #USES SIMPLE UTILITIES
                        gamma12 = ..gamma12.., #USES SIMPLE UTILITIES
                        verbose = TRUE,
                        tol.zero = 0.00001,
                        tol.lim = 0.0001
                        ){
  lambda.2 = seq(lower[1], upper[1], length.out = grid_points)
  w.2      = seq(lower[2], upper[2], length.out = grid_points)
  grid = expand.grid(lambda.2, w.2)
  ngrid_total = nrow(grid)
  # print(ngrid_total)
  # print(ngrid)
  results = sapply(X = 1:nrow(grid),
                   FUN = function(xx){
                     interimOptimalTrial_EU(lambda.2 = grid[xx,1],
                                            w.2      = grid[xx,2],
                                            x = x,
                                            r = r,
                                            rho = 0,
                                            mu_1 = mu_1, psi_1 = psi_1,
                                            mu_2 = mu_2, psi_2 = psi_2,
                                            rho0 = rho0,
                                            lambda = lambda,
                                            lambda.1 = lambda.1,
                                            n = n,
                                            w.1 = w.1,
                                            alpha = alpha,
                                            tol.zero = tol.zero,
                                            gamma1  = gamma1,
                                            gamma2  = gamma2,
                                            gamma12 = gamma12,
                                            verbose = FALSE)
                   })
  # results.df = data.frame(grid,
  #                         expected_utility = unlist(results[1,]),
  #                         case = unlist(results[2,]))
  # results.df
  grid_max = which.max(results[1,])
  list(par = unname(unlist(grid[grid_max, ])),
       value = results[1,][[grid_max]],
       convergence = 1 * (length(grid_max)!=1),
       feval = ngrid_total,
       niter = NA,
       message = if (length(grid_max)==1) "Successful convergence" else "Check!")
}

#' Performs the statistical test using first and second stage data
#'
#' See manuscript
#'
#' @param x.1      the values of the first stage estimates thetahat.1_1 and thetahat.1_2
#' @param x.2      the values of the second stage estimates  thetahat.2_1 and thetahat.2_2
#' @param r        proportion of patients in 1st stage of the trial
#' @param sigma_1  standard deviation in subpopulation 1
#' @param sigma_2  standard deviation in subpopulation 2
#' @param rho      correlation between estimates in subpopulation 1 and subpopulation 2 (ignored)
#' @param lambda   proportion of patients in sub-population 1 in the overal patient population
#' @param lambda.1 proportion of patients in sub-population 1 in the trial 1st stage
#' @param lambda.2 proportion of patients in sub-population 1 in the trial 2nd stage
#' @param n        trial sample size
#' @param w.1      factor for allocation of alpha to H01 in the testing correction (1st stage)
#' @param w.2      factor for allocation of alpha to H01 in the testing correction (final testing)
#' @param alpha    significance level
#' @param gamma1   gain when rejecting H01. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param gamma2   gain when rejecting H02. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param gamma12  gain when rejecting H01 and H02. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param verbose  Print intermediate values for calculations
#' @param tol.zero Absolute distance to 0 or 1 to be considered 0 or 1
#'
#' @export
decisionOptimalTrial <- function(lambda.2, w.2,
                              x.1,
                              x.2,
                              r,
                              sigma_1 = 1, sigma_2 = 1, rho = 0,
                              lambda = 0.5, lambda.1 = 0.5,
                              n = 100,
                              w.1 = 0.5,
                              alpha = 0.025,
                              gamma1  = ..gamma1..,  #USES SIMPLE UTILITIES
                              gamma2  = ..gamma2..,  #USES SIMPLE UTILITIES
                              gamma12 = ..gamma12.., #USES SIMPLE UTILITIES
                              fixed.weights = FALSE,
                              verbose = TRUE,
                              show.case = TRUE){
  # x.1 and x.2 are 2D vectors with the estimates from 1st stage and 2nd stage resp.
  theta.1_1 <- (x.1[1])
  theta.1_2 <- (x.1[2])
  theta.2_1 <- (x.2[1])
  theta.2_2 <- (x.2[2])

  lambda.1 = unname(lambda.1)
  lambda.2 = unname(lambda.2)
  w.1 = unname(w.1)
  w.2 = unname(w.2)
  r = unname(r)

  # Calculate trial prevalences
  lambda.1_1 = (lambda.1)
  lambda.1_2 = (1 - lambda.1)
  lambda.2_1 = (lambda.2)
  lambda.2_2 = (1 - lambda.2)

  # Calculate trial weights
  w.1_1 = w.1
  w.1_2 = (1 - w.1)
  w.2_1 = w.2
  w.2_2 = (1 - w.2)
  # Subgroup: Conditional expectation given 1st stage
  if (lambda.1_1 != 0){ # Case when subgroup 1 trial prevalence is >0 in first stage
    sigma2_1 = sigma_1^2
    sigma2_theta.1_1 =       4 * sigma2_1 / (lambda.1_1 * r * n)
    sigma2_theta.2_1_fixed = 4 * sigma2_1 / (lambda.1_1 * (1 - r) * n)
    sigma2_theta.2_1 =       4 * sigma2_1 / (lambda.2_1 * (1 - r) * n)
    sigma2_theta.p_1_fixed = 4 * sigma2_1 / n / (lambda.1_1 * r + lambda.1_1 * (1 - r))
    sigma2_theta.p_1 =       4 * sigma2_1 / n / (lambda.1_1 * r + lambda.2_1 * (1 - r))
    Z.1_1 = sqrt(n * lambda.1_1 * r) * theta.1_1 / (2 * sigma_1)
    f.1_1 = lambda.1_1 * r / (lambda.1_1 * r  + lambda.2_1 * (1-r))
    f.2_1 = 1 - f.1_1
  } else {   # If subgroup 1 trial prevalence is 0 in first stage, we only have 2nd stage data
    sigma2_1 = sigma_1^2
    sigma2_theta.1_1 = 0
    sigma2_theta.2_1 = 4 * sigma2_1 / (lambda.2_1 * (1-r) * n)
    sigma2_theta.p_1 = sigma2_theta.2_1
    Z.1_1 = NA
    f.1_1 = 0
    f.2_1 = 1
  }
  if (lambda.2_1 == 0){ # If subgroup 1 trial prevalence is 0 in second stage, we only have 1st stage data
    sigma2_theta.2_1 = 0
  }

  # Etheta.p_1_cond = f.1_1 * theta.1_1
  # Vtheta.p_1_cond = f.2_1^2 * (sigma2_theta.2_1)
  Etheta.p_1_cond = r * theta.1_1
  Vtheta.p_1_cond = (1 - r)^2 * (sigma2_theta.2_1_fixed)
  theta.p_1 = f.1_1 * theta.1_1 + f.2_1 * theta.2_1
  Z.p_1 = theta.p_1 / sqrt(sigma2_theta.p_1)
  P.p_1 = pnorm(Z.p_1, lower.tail = FALSE)

  P.1_1 = pnorm(Z.1_1, lower.tail = FALSE)
  if (lambda.2_1 == 0){# If subgroup 2 trial prevalence is 0 in second stage, we only have 1st stage data
    P.p_1 = P.1_1
  }
  # Complement: Conditional expectation given 2nd stage
  if (lambda.1_2 != 0){# Case when subgroup 2 trial prevalence is >0 in first stage
    sigma2_2 = sigma_2^2
    sigma2_theta.1_2 =       4 * sigma2_2 / (lambda.1_2 * r * n)
    sigma2_theta.2_2_fixed = 4 * sigma2_2 / (lambda.1_2 * (1-r) * n)
    sigma2_theta.2_2 =       4 * sigma2_2 / (lambda.2_2 * (1-r) * n)
    sigma2_theta.p_2_fixed = 4 * sigma2_2 / n / (lambda.1_2 * r + lambda.1_2 * (1 - r))
    sigma2_theta.p_2 =       4 * sigma2_2 / n / (lambda.1_2 * r + lambda.2_2 * (1 - r))
    Z.1_2 = sqrt(n * lambda.1_2 * r) * theta.1_2/ (2 * sigma_2)
    f.1_2 = lambda.1_2 * r / (lambda.1_2 * r  + lambda.2_2 * (1 - r))
    f.2_2 = 1 - f.1_2
  } else {# If subgroup 2 trial prevalence is 0 in first stage, we only have 2nd stage data
    sigma2_2 = sigma_2^2
    sigma2_theta.1_2 = Inf
    sigma2_theta.2_2 = 4 * sigma2_2 / (lambda.2_2 * (1 - r) * n)
    sigma2_theta.p_2 = sigma2_theta.2_2
    Z.1_2 = NA
    f.1_2 = 0
    f.2_2 = 1
  }
  if (lambda.2_2 == 0){# If subgroup 2 trial prevalence is 0 in second stage, we only have 1st stage data
    sigma2_theta.2_2 = 0
  }

  Etheta.p_2_cond = r * theta.1_2
  Vtheta.p_2_cond = (1 - r)^2 * (sigma2_theta.2_2_fixed)
  theta.p_2 = f.1_2 * theta.1_2 + f.2_2 * theta.2_2
  Z.p_2 = theta.p_2 / sqrt(sigma2_theta.p_2)
  P.p_2 = pnorm(Z.p_2, lower.tail = FALSE)

  P.1_2 = pnorm(Z.1_2, lower.tail = FALSE)

  if (lambda.2_2 == 0){# If subgroup 2 trial prevalence is 0 in second stage, we only have 1st stage data
    P.p_2 = P.1_2
  }
  # Calculate A_1, A_2, B_1, B_2
  # First calculate critical boundaries for first stage
  crit_w.1_1 = qnorm(1 - alpha * w.1_1)
  crit_w.1_2 = qnorm(1 - alpha * (1 - w.1_1))
  crit       = qnorm(1 - alpha)

  # Based on estimates from first stage we derive the distributions of
  # the conditional statistics
  # Compute conditional error rates
  A12_1 = pnorm(crit_w.1_1 * sqrt(sigma2_theta.p_1_fixed),
                mean = Etheta.p_1_cond, sd = sqrt(Vtheta.p_1_cond), lower.tail = FALSE)
  A12_2 = pnorm(crit_w.1_2 * sqrt(sigma2_theta.p_2_fixed),
                mean = Etheta.p_2_cond, sd = sqrt(Vtheta.p_2_cond), lower.tail = FALSE)
  A12 = A12_1 + A12_2 - A12_1*A12_2

  A1 = pnorm(crit * sqrt(sigma2_theta.p_1_fixed),
             mean = Etheta.p_1_cond, sd = sqrt(Vtheta.p_1_cond), lower.tail = FALSE)
  A2 = pnorm(crit * sqrt(sigma2_theta.p_2_fixed),
             mean = Etheta.p_2_cond, sd = sqrt(Vtheta.p_2_cond), lower.tail = FALSE)

  # Calculate P.values from first stage
  P.1_1 = pnorm(Z.1_1, lower.tail = FALSE)
  P.1_2 = pnorm(Z.1_2, lower.tail = FALSE)
  # Calculate Z values from 2nd stage and their respective p.values
  Z.2_1 = theta.2_1 / sqrt(sigma2_theta.2_1)
  Z.2_2 = theta.2_2 / sqrt(sigma2_theta.2_2)
  P.2_1 = pnorm(Z.2_1, lower.tail = FALSE)
  P.2_2 = pnorm(Z.2_2, lower.tail = FALSE)

  if (verbose) {
    cat("lambda.1", lambda.1, "w.1", w.1, "r.1", r, "\n")
    cat("lambda.2", lambda.2, "w.2", w.2, "\n")
    cat("1st stage error rates", "A1: ", A1, "- A2: ", A2, "-\n A12_1: ", A12_1, "- A12_2: ", A12_2, "- A12: ", A12, "\n")
    cat("Second Stage P values", P.2_1, P.2_2, "\n")
  }

  # Pre-allocate values for rejection
  RejectH_01 = RejectH_02 = RejectH_0102 = 0
  RejectH_01.p = RejectH_02.p = RejectH_0102.p = 0
  # Check whether we have a two stage trial of one stage (r=1)
  isTwoStageTrial = (r != 0 & r!= 1)
  if (isTwoStageTrial == FALSE){
    if(r == 1){
      # If we have a one stage trials, simply test the hypotheses using 1st stage p-vaues
      if (P.1_1 < alpha & (P.1_1 < w.1_1 * alpha | P.1_2 < w.1_2 * alpha)){
        RejectH_01 = 1
      }
      if (P.1_2 < alpha & (P.1_1 < w.1_1 * alpha | P.1_2 < w.1_2 * alpha)){
        RejectH_02 = 1
      }
      if (P.1_1 < w.1_1 * alpha | P.1_2 < w.1_2 * alpha){
        RejectH_0102 = 1
      }
    } else {
      # If we have a one stage trials, simply test the hypotheses using 2nd stage p-vaues
      if (P.2_1 < alpha & (P.2_1 < w.2_1 * alpha | P.2_2 < w.2_2 * alpha)){
        RejectH_01 = 1
      }
      if (P.2_2 < alpha & (P.2_1 < w.2_1 * alpha | P.2_2 < w.2_2 * alpha)){
        RejectH_02 = 1
      }
      if (P.2_1 < w.2_1 * alpha | P.2_2 < w.2_2 * alpha){
        RejectH_0102 = 1
      }
    }
  } else {
    # Case: 0 < lambda.2_1 < 1 -----------------------------------------------
    if (isTwoStageTrial & !fixed.weights & (0 < lambda.2_1 & lambda.2_1 < 1)){
      if (P.2_1 <= A1 & (P.2_1 <= w.2_1*A12 | P.2_2 <= w.2_2*A12)){
        RejectH_01 = 1
      }
      if (P.2_2 <= A2 & (P.2_1 <= w.2_1*A12 | P.2_2 <= w.2_2*A12)){
        RejectH_02 = 1
      }
      if (P.2_1 <= w.2_1*A12 | P.2_2 <= w.2_2*A12){
        RejectH_0102 = 1
      }
    }
    # Case: 0 < lambda.2_1 < 1 - Fixed weights ----------------------------------
    if (isTwoStageTrial & fixed.weights & (0 < lambda.2_1 & lambda.2_1 < 1)){
      if (P.2_1 <= A1 & (P.2_1 <= A12_1 | P.2_2 <= A12_2)){
        RejectH_01 = 1
      }
      if (P.2_2 <= A2 & (P.2_1 <= A12_1 | P.2_2 <= A12_2)){
        RejectH_02 = 1
      }
      if (P.2_1 <= A12_1 | P.2_2 <= A12_2){
        RejectH_0102 = 1
      }
    }

    # Use predictive distribution of theta.2_1 given theta.1_1
    z_A1 = qnorm(A1, lower.tail = FALSE)
    z_A2 = qnorm(A2, lower.tail = FALSE)
    z_A12_1 = qnorm(w.2_1 * A12, lower.tail = FALSE)
    z_A12_2 = qnorm(w.2_2 * A12, lower.tail = FALSE)

    crit_A1 = z_A1 * sqrt(sigma2_theta.2_1)
    crit_A2 = z_A2 * sqrt(sigma2_theta.2_2)
    crit_A12_1 = z_A12_1 * sqrt(sigma2_theta.2_1)
    crit_A12_2 = z_A12_2 * sqrt(sigma2_theta.2_2)

    if        (crit_A1 <= crit_A12_1 & crit_A2 <= crit_A12_2){
      case <- 1
    } else if (crit_A1  > crit_A12_1 & crit_A2 <= crit_A12_2){
      case <- 2
    } else if (crit_A1 <= crit_A12_1 & crit_A2  > crit_A12_2){
      case <- 3
    } else if (crit_A1  > crit_A12_1 & crit_A2  > crit_A12_2){
      case <- 4
    }

    if (!is.na(P.p_1) & !is.na(P.p_2)){
      if (P.p_1 < alpha & (P.p_1 < w.2_1 * alpha | P.p_2 <  w.2_2 * alpha)){
        RejectH_01.p = 1
      } else{
        RejectH_01.p = 0
      }
      if (P.p_2 < alpha & (P.p_1 <  w.2_1 * alpha | P.p_2 <  w.2_2 * alpha)){
        RejectH_02.p = 1
      } else{
        RejectH_02.p = 0
      }
    }
  }

  if (verbose) {
    first  = c(theta.1_1, Z.1_1, P.1_1,
               theta.2_1, Z.2_1, P.2_1,
               theta.p_1, Z.p_1, P.p_1, RejectH_01, RejectH_01.p)
    second = c(theta.1_2, Z.1_2, P.1_2,
               theta.2_2, Z.2_2, P.2_2,
               theta.p_2, Z.p_2, P.p_2, RejectH_02, RejectH_02.p)
    what   = c("First Stage estimates",  "First Stage Z-values", "First Stage p-values",
               "Second Stage estimates", "Second Stage Z-values", "Second Stage p-values",
               "Pooled estimates",       "Pooled Z-values", "Pooled p-values",
               "Decision", "Decision from pooled p-values")
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|", "Trial summary", "Subpop 1", "Subpop 2")))
    # cat("\n", paste0(sprintf("%-30s|%10s|%10s|", "------------------------------", "----------", "----------")))
    # cat("\n", paste0(sprintf("%-30s|%10s|%10s|", what, round(first, 3), round(second, 3)), "\n"), "\n")
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|", "------------------------------", "----------", "----------")))
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|", what[1:3], round(first[1:3], 3), round(second[1:3], 3)), "\n"))
    cat("", paste0(sprintf("%-30s|%10s|%10s|", "------------------------------", "----------", "----------")))
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|", what[4:6], round(first[4:6], 3), round(second[4:6], 3)), "\n"))
    cat("", paste0(sprintf("%-30s|%10s|%10s|", "------------------------------", "----------", "----------")))
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|", what[7:9], round(first[7:9], 3), round(second[7:9], 3)), "\n"))
    cat("", paste0(sprintf("%-30s|%10s|%10s|", "------------------------------", "----------", "----------")))
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|", what[10],  round(first[10], 3),  round(second[10], 3)), "\n"))
    cat("", paste0(sprintf("%-30s|%10s|%10s|", "------------------------------", "----------", "----------")))
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|", what[11],  round(first[11], 3),  round(second[11], 3)), "\n"), "\n")
  }

  # Calculate expected utility. Note this should not depend on the estimates.
  # or true values in the function.
  EU <- (gamma1(lambda = lambda) * RejectH_01 +
         gamma2(lambda = lambda) * RejectH_02 +
         gamma12(lambda = lambda) * RejectH_01*RejectH_02)

  out = c(RejectH_01 = RejectH_01, RejectH_02 = RejectH_02,
          RejectH_0102 = RejectH_0102, EU = EU,
          theta.p_1 = theta.p_1, theta.p_2 = theta.p_2,

          Z.1_1 = Z.1_1, Z.1_2 = Z.1_2,
          Z.2_1 = Z.2_1, Z.2_2 = Z.2_2,
          # Z.p_1 = Z.p_1, Z.p_2 = Z.p_2,

          P.1_1 = P.1_1, P.1_2 = P.1_2,
          P.2_1 = P.2_1, P.2_2 = P.2_2,
          P.p_1 = P.p_1, P.p_2 = P.p_2,
          # RejectH_01.p = RejectH_01.p, RejectH_02.p = RejectH_02.p,
          # case = case,
          # A1 = A1, A2 = A2, A12 = A12,
          # C = C, C. = C., C1, C2, B1 = B1, B2 = B2,
          r = r)
  out
}




#' Simulates a clinical trial with optimized sample allocation.
#'
#' Given true values for theta_1 and theta_2, the function generates first stage data,
#' then find optimal values lambda.2 and w.2, generates second stage data, and
#' finally test the hypotheses.
#'
#' @param r.1      proportion of patients in 1st stage of the trial
#' @param lambda   proportion of patients in sub-population 1 in the overal patient population
#' @param lambda.1 proportion of patients in sub-population 1 in the trial 1st stage
#' @param w.1      factor for allocation of alpha to H01 in the testing correction (1st stage)
#' @param theta_1  true parameter value for sub-population 1
#' @param theta_2  true parameter value for sub-population 2
#' @param mu_1     prior mean for the treatment effect in sub-population 1
#' @param psi_1    prior standard deviation for the treatment effect in sub-population 1
#' @param mu_2     prior mean for the treatment effect in sub-population 2
#' @param psi_2    prior standard deviation for the treatment effect in sub-population 2
#' @param rho_0    prior correlation between treatment effect in sub-population 1 and sub-population 2 (ignored)
#' @param alpha    significance level
#' @param n        trial sample size
#' @param gamma1   gain when rejecting H01. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param gamma2   gain when rejecting H02. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param gamma12  gain when rejecting H01 and H02. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param tol.zero Absolute distance to 0 or 1 to be considered 0 or 1
#' @param verbose  Print intermediate values for calculations
#' @param use.prior.as.true logical. If TRUE, then the true theta is generate using the prior distribution and theta_1 and theta_2 are not used.
#'
#' @export
optimal1Trial <- function(r.1      =   0.1,
                          lambda   =   0.5,
                          lambda.1 =   0.5,
                          w.1      =   0.5,
                          theta_1  =   0,
                          theta_2  =   0,
                          mu_1     =   0,
                          psi_1    =   1,
                          mu_2     =   0,
                          psi_2    =   1,
                          rho0     =   0,
                          alpha    =   0.05,
                          n        = 100,
                          gamma1   = ..gamma1..,  #USES SIMPLE UTILITIES
                          gamma2   = ..gamma2..,  #USES SIMPLE UTILITIES
                          gamma12  = ..gamma12.., #USES SIMPLE UTILITIES
                          tol.zero =   0.00001,
                          tol.lim  =   0.0001,
                          verbose  = TRUE,
                          use.prior.as.true = FALSE,
                          optim_func = "optim"){

  if(use.prior.as.true){
    Sigma.0 = matrix(c(psi_1^2,
                       rho0 * psi_1 * psi_2,
                       rho0 * psi_1 * psi_2,
                       psi_2^2),
                     nrow = 2)
    true.mean = mvtnorm::rmvnorm(n = 1,
                           mean = c(mu_1, mu_2),
                           sigma = Sigma.0)
  } else {
    true.mean = c(theta_1, theta_2)
  }

  if (r.1 != 0){
    # Simulate 1st stage estimates
    Sigma.1     = diag(c(4/lambda.1/n/r.1, 4/(1-lambda.1)/n/r.1))
    theta.hat.1 = mvtnorm::rmvnorm(n = 1,
                                   mean  = true.mean,
                                   sigma = Sigma.1)
  } else {
    theta.hat.1 = c(NA, NA)
  }
  # cat(theta.hat.1)
  if (verbose) cat("Simulate 1st stage: DONE \n")

  # Given First stage estimates, calculate the values of lambda.2 and w.2
  #  with the maximum expected utility
  if (r.1 != 1 & r.1 != 0){
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
  if (verbose) cat("Interim decision: DONE \n", interimDecision$feval, "\n")

  # Store values of lambda.2 and w.2
    lambda.2 = interimDecision$par[1]
    w.2      = interimDecision$par[2]
    interimEU = interimDecision$value
  } else {
    lambda.2  = lambda.1
    w.2       = w.1
    interimEU = NA
  }
  # When lambda.2 and w.2 are very close to 0 or 1, we just set them to 0 or 1 respectively.
  # if (lambda.2 < tol.zero)   lambda.2 = 0
  # if (lambda.2 > 1-tol.zero) lambda.2 = 1
  # if (w.2 < tol.zero)        w.2 = 0
  # if (w.2 > 1-tol.zero)      w.2 = 1

  # Simulates second stage estimates given the optimized lambda.2
  # Calculates first the variances
  Sigma.2_1 = 4 / (lambda.2) / n / (1 - r.1)
  Sigma.2_2 = 4 / (1 - lambda.2) / n / (1 - r.1)
  # And simulate the estimates. A NA pre-allocated vector is used so that
  # We only simulate the value if it is actually needed.
  theta.hat.2 = c(NA, NA)
  if (lambda.2 != 0) theta.hat.2[1] = rnorm(n = 1, mean  = true.mean[1], sd = sqrt(Sigma.2_1))
  if (lambda.2 != 1) theta.hat.2[2] = rnorm(n = 1, mean  = true.mean[2], sd = sqrt(Sigma.2_2))
  if (verbose) cat("Simulate 2nd stage: DONE \n")

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
  if (verbose) cat("Decision: DONE \n")
  out = c(decision,
          lambda = lambda,
          lambda.1 = lambda.1,
          lambda.2 = lambda.2,
          w.1 = w.1,
          w.2 = w.2,
          theta.hat.1_ = theta.hat.1,
          theta.hat.2_ = theta.hat.2,
          interimEU = interimEU)
  class(out) = "simOptimalTrial"
  out
}




#' Simulates M clinical trials with optimized sample allocation using simulate1trial
#'
#' @exportClass simOptimalTrial
#' @export
optimalMtrials = function(nsim   =  1000,
                           mc.cores =   1,
                           r.1      =   0.1,
                           lambda   =   0.5,
                           lambda.1 =   0.5,
                           w.1      =   0.5,
                           theta_1  =   0,
                           theta_2  =   0,
                           mu_1     =   0,
                           psi_1    =   1,
                           mu_2     =   0,
                           psi_2    =   1,
                           rho0     =   0,
                           alpha    =   0.05,
                           n        = 100,
                           gamma1   = ..gamma1..,  #USES SIMPLE UTILITIES
                           gamma2   = ..gamma2..,  #USES SIMPLE UTILITIES
                           gamma12  = ..gamma12.., #USES SIMPLE UTILITIES
                           tol.zero =   0.00001,
                           tol.lim  =   0.0001,
                           verbose  = TRUE,
                           summarize = FALSE,
                           use.prior.as.true = FALSE,
                           optim_func = "optim"){
  results = parallel::mclapply(1:nsim, FUN = function(x){
    optimal1Trial(r.1      =   r.1,
                   lambda   =   lambda,
                   lambda.1 =   lambda.1,
                   w.1      =   w.1,
                   theta_1   =   theta_1,
                   theta_2   =   theta_2,
                   mu_1  = mu_1,
                   psi_1 = psi_1,
                   mu_2  = mu_2,
                   psi_2 = psi_2,
                   rho0  = rho0,
                   alpha = alpha,
                   n        = n,
                   gamma1   = gamma1,
                   gamma2   = gamma2,
                   gamma12  = gamma12,
                   tol.zero = tol.zero,
                   tol.lim  = tol.lim,
                   verbose  = verbose,
                   use.prior.as.true = use.prior.as.true,
                  optim_func = optim_func)
  }, mc.cores = mc.cores)

  res = data.frame(do.call(rbind,results))
  out = res
  if(summarize){
    res2 = colMeans(res)
    out = t(res2)
  }
  results.df = list(simresults = out,
                    summarize = summarize)
  class(results.df) = "simOptimalTrial"
  results.df
}



#' Simulates a clinical trial with optimized sample allocation.
#'
#' Given true values for theta_1 and theta_2, the function generates first stage data,
#' then find optimal values lambda.2 and w.2, generates second stage data, and
#' finally test the hypotheses.
#'
#' @param r.1      proportion of patients in 1st stage of the trial
#' @param lambda   proportion of patients in sub-population 1 in the overal patient population
#' @param lambda.1 proportion of patients in sub-population 1 in the trial 1st stage
#' @param w.1      factor for allocation of alpha to H01 in the testing correction (1st stage)
#' @param theta_1  true parameter value for sub-population 1
#' @param theta_2  true parameter value for sub-population 2
#' @param mu_1     prior mean for the treatment effect in sub-population 1
#' @param psi_1    prior standard deviation for the treatment effect in sub-population 1
#' @param mu_2     prior mean for the treatment effect in sub-population 2
#' @param psi_2    prior standard deviation for the treatment effect in sub-population 2
#' @param rho_0    prior correlation between treatment effect in sub-population 1 and sub-population 2 (ignored)
#' @param alpha    significance level
#' @param n        trial sample size
#' @param gamma1   gain when rejecting H01. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param gamma2   gain when rejecting H02. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param gamma12  gain when rejecting H01 and H02. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param tol.zero Absolute distance to 0 or 1 to be considered 0 or 1
#' @param verbose  Print intermediate values for calculations
#'
#' @export
fullOptimization <- function(M = 1000, mc.cores = 40,
                          r.1      =   0.1,
                          lambda   =   0.5,
                          lambda.1 =   0.5,
                          w.1      =   0.5,
                          mu_1     =   0,
                          psi_1    =   1,
                          mu_2     =   0,
                          psi_2    =   1,
                          rho0     =   0,
                          alpha    =   0.05,
                          n        = 100,
                          gamma1   = ..gamma1..,  #USES SIMPLE UTILITIES
                          gamma2   = ..gamma2..,  #USES SIMPLE UTILITIES
                          gamma12  = ..gamma12.., #USES SIMPLE UTILITIES
                          tol.zero =   0.00001,
                          tol.lim  =   0.0001,
                          verbose  = TRUE,
                          summarize = TRUE,
                          optim_func = "hjkb",
                          grid_points = 10){
  # Simulate 1st stage estimates
  Sigma.1     = diag(c(4/lambda.1/n/r.1, 4/(1-lambda.1)/n/r.1))

  Sigma.0 = matrix(c(psi_1^2,
                     rho0 * psi_1 * psi_2,
                     rho0 * psi_1 * psi_2,
                     psi_2^2),
                   nrow = 2)

  theta.hat.1 = mvtnorm::rmvnorm(n = M,
                                 mean  = c(mu_1, mu_2),
                                 sigma = Sigma.1 + Sigma.0)
  if (verbose) cat("Simulate 1st stage: DONE \n")

  # Given First stage estimates, calculate the values of lambda.2 and w.2
  #  with the maximum expected utility
  res = parallel::mclapply(X = 1:M, FUN = function(i) {
    interimDecision = interimOptimalTrial(x = theta.hat.1[i, ],
                                          r = r.1,
                                          sigma_1 = 1, sigma_2 = 1,
                                          rho  = 0,
                                          mu_1 = mu_1, psi_1 = psi_1,
                                          mu_2 = mu_2, psi_2 = psi_1,
                                          rho0 = rho0,
                                          lambda = lambda,
                                          lambda.1 = lambda.1,
                                          n   = n,
                                          w.1 = w.1,
                                          alpha = alpha,
                                          gamma1 = gamma1,
                                          gamma2 = gamma2,
                                          gamma12 = gamma12,
                                          tol.zero = tol.zero,
                                          tol.lim = tol.lim,
                                          optim_func = optim_func,
                                          grid_points = grid_points,
                                          verbose = FALSE)
  if (verbose) cat("Interim decision: DONE \n")
  # Store values of lambda.2 and w.2
  c(lambda.1 = lambda.1,
    w.1 = w.1,
    r.1 = r.1,
    lambda.2 = interimDecision$par[1],
    w.2      = interimDecision$par[2],
    EU       = interimDecision$value)
  }, mc.cores = mc.cores)
  res = do.call(rbind, res)
  out = res
  if(summarize){
    res2 = colMeans(res)
    out = res2
  }
  out
}

#'@export
summary.simOptimalTrial <- function(object, verbose = FALSE, ...){
  results.df = object$simresults
  nsim = nrow(results.df)
  if (object$summarize) {
    h01   = results.df[1]
    h02   = results.df[2]
    h0102 = results.df[3]
    avgU  = results.df[4]
  } else{
    h01   = mean(results.df[, 1])
    h02   = mean(results.df[, 2])
    h0102 = mean(results.df[, 3])
    avgU  = mean(results.df$EU)
  }

  first = c(h01, h02, h0102, avgU)
  what = c("Power to reject H_01", "Power to reject H_02",
           "Power to reject H_01 int H_02",
           "Average utility")
  pp <- function(){
    cat("\n", paste0(sprintf("%-30s|%10s|", "Trial Results", "")))
    cat("\n", paste0(sprintf("%-30s|%10s|", "------------------------------", "----------")))
    cat("\n", paste0(sprintf("%-30s|%10s|", what, round(first, 5)), "\n"), "\n")
  }
  if (verbose) pp()
  data.frame(`Property` = what, `Trial Results` = first)
}

#' Expected utility for values of lambda.2 and w.2 after 1st stage
#'
#' Trying to speed up the calculation of interimOptimalTrial_EU
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
fast_interimOptimalTrial_EU <- function(lambda.2_1, w.2_1,
                                        x,
                                        theta.1_1,
                                        theta.1_2,
                                        r,
                                        sigma_1 = 1, sigma_2 = 1, rho = 0,
                                        Sigma0 = diag(1,nrow = 2),
                                        mu = c(0,0), mu_1 = 0, psi_1 = 1,
                                        mu_2 = 0, psi_2 = 1, rho0 = 0,
                                        lambda = 0.5, lambda.1_1 = 0.5, lambda.1_2 = 0.5,
                                        n = 100,
                                        w.1_1 = 0.5,
                                        w.1_2 = 0.5,
                                        alpha = 0.025,
                                        gamma1  = ..gamma1..,  #USES SIMPLE UTILITIES
                                        gamma2  = ..gamma2..,  #USES SIMPLE UTILITIES
                                        gamma12 = ..gamma12.., #USES SIMPLE UTILITIES
                                        verbose = TRUE, tol.zero = 1e-10){
  lambda.2_2 = (1 - lambda.2_1)
  w.2_2 = (1 - w.2_1)

  # Subgroup: Conditional expectation given 1st stage
  if (lambda.1_1 != 0){ # Case when subgroup 1 trial prevalence is > 0 in first stage
    psi2_1 = psi_1^2
    sigma2_1 = sigma_1^2
    psi2.1_1 = sigma2_1 * psi2_1 / (lambda.1_1 * r * n * psi2_1 / 4 + sigma2_1)
    sigma2_theta.1_1 = 4 * sigma2_1 / (lambda.1_1 * r * n)
    sigma2_theta.2_1_fixed = 4 * sigma2_1 / (lambda.1_1 * (1-r) * n)
    sigma2_theta.2_1       = 4 * sigma2_1 / (lambda.2_1 * (1-r) * n)
    sigma2_theta.p_1_fixed = 4 * sigma2_1 / n / (lambda.1_1 * r + lambda.1_1 * (1 - r))
    sigma2_theta.p_1       = 4 * sigma2_1 / n / (lambda.1_1 * r + lambda.2_1 * (1 - r))
    mu.1_1 = psi2.1_1 * (mu_1 / psi2_1 + theta.1_1 / sigma2_theta.1_1)
    Z.1_1 = theta.1_1 / (2 * sigma_1 / sqrt(n * lambda.1_1 * r))
    f.1_1 = lambda.1_1 * r / (lambda.1_1 * r  + lambda.2_1 * (1-r))
    f.2_1 = 1 - f.1_1
  } else { # If subgroup 1 trial prevalence is 0 in first stage, we only have 2nd stage data
    psi2_1 = psi_1^2
    sigma2_1 = sigma_1^2
    psi2.1_1 = psi2_1
    sigma2_theta.1_1 = 0
    sigma2_theta.2_1 = 4 * sigma2_1 / (lambda.2_1 * (1 - r) * n)
    sigma2_theta.p_1 = sigma2_theta.2_1
    mu.1_1 = mu_1
    Z.1_1 = NA
    f.1_1 = 0
    f.2_1 = 1
  }
  if (lambda.2_1 == 0){ # If subgroup 1 trial prevalence is 0 in second stage, we only have 1st stage data
    sigma2_theta.2_1 = 0
  }

  Etheta.p_1_cond = r * theta.1_1
  Vtheta.p_1_cond = (1 - r)^2 * (sigma2_theta.2_1_fixed)

  # Complement: Conditional expectation given 2nd stage
  if (lambda.1_2 != 0){# Case when subgroup 2 trial prevalence is >0 in first stage
    sigma2_2 = sigma_2^2
    psi2_2 = psi_2^2
    psi2.1_2 = sigma2_2 * psi2_2 / (lambda.1_2 * r * n * psi2_2 / 4 + sigma2_2)
    sigma2_theta.1_2 =       4 * sigma2_2 / (lambda.1_2 * r * n)
    sigma2_theta.2_2_fixed = 4 * sigma2_2 / (lambda.1_2 * (1 - r) * n)
    sigma2_theta.2_2 =       4 * sigma2_2 / (lambda.2_2 * (1 - r) * n)
    sigma2_theta.p_2_fixed = 4 * sigma_2^2 / n / (lambda.1_2 * r + lambda.1_2 * (1 - r))
    sigma2_theta.p_2 =       4 * sigma_2^2 / n / (lambda.1_2 * r + lambda.2_2 * (1 - r))
    mu.1_2 = psi2.1_2 * (mu_2 / psi2_2 + theta.1_2 / sigma2_theta.1_2)
    Z.1_2 = theta.1_2 / (2 * sigma_2 / sqrt(n * lambda.1_2 * r))
    f.1_2 = lambda.1_2 * r / (lambda.1_2 * r  + lambda.2_2 * (1-r))
    f.2_2 = 1 - f.1_2
  } else { # If subgroup 2 trial prevalence is 0 in first stage, we only have 2nd stage data
    psi2_2 = psi_2^2
    sigma2_2 = sigma_2^2
    psi2.1_2 = psi2_2
    sigma2_theta.1_2 = 0
    sigma2_theta.2_2 = 4 * sigma2_2 / (lambda.2_2 * (1 - r) * n)
    sigma2_theta.p_2 = sigma2_theta.2_2
    mu.1_2 = mu_2
    Z.1_2 = NA
    f.1_2 = 0
    f.2_2 = 1
  }
  if (lambda.2_2 == 0){# If subgroup 2 trial prevalence is 0 in second stage, we only have 1st stage data
    sigma2_theta.2_2 = 0
  }

  Etheta.p_2_cond = r * theta.1_2
  Vtheta.p_2_cond = (1 - r)^2 * (sigma2_theta.2_2_fixed)

  Sigma.1_theta  = matrix(c(sigma2_theta.1_1,
                            0 , 0,
                            sigma2_theta.1_2),
                          nrow = 2)
  Sigma.2_theta  = matrix(c(sigma2_theta.2_1,
                            0 , 0,
                            sigma2_theta.2_2),
                          nrow = 2)

  SigmaPlusInv = solve(Sigma0 + Sigma.1_theta)
  Sigma.1 = Sigma0 %*% SigmaPlusInv %*% Sigma.1_theta
  mu.1    = drop(Sigma0 %*% SigmaPlusInv %*% drop(x) + Sigma.1_theta %*% SigmaPlusInv %*% mu)
  mu.1_1  = mu.1[1]
  mu.1_2  = mu.1[2]

  # Covariance matrix of the prior predictive distribution
  Sigma_ppd = Sigma.1 + Sigma.2_theta

  # Calculate 1st stage p.values
  P.1_1 = pnorm(Z.1_1, lower.tail = FALSE)
  P.1_2 = pnorm(Z.1_2, lower.tail = FALSE)

  # Calculate A_1, A_2, B_1, B_2
  # First calculate critical boundaries for first stage
  crit_w.1_1 = qnorm(1 - alpha * w.1_1)
  crit_w.1_2 = qnorm(1 - alpha * (1-w.1_1))
  crit       = qnorm(1 - alpha)
  # Based on estiamtes from first stage we derive the distributions of
  # the conditional statistics
  # Compute conditional error rates
  A12_1 = pnorm(crit_w.1_1 * sqrt(sigma2_theta.p_1_fixed),
                mean = Etheta.p_1_cond, sd = sqrt(Vtheta.p_1_cond), lower.tail = FALSE)
  A12_2 = pnorm(crit_w.1_2 * sqrt(sigma2_theta.p_2_fixed),
                mean = Etheta.p_2_cond, sd = sqrt(Vtheta.p_2_cond), lower.tail = FALSE)
  A12 = A12_1 + A12_2 - A12_1*A12_2

  A1 = pnorm(crit * sqrt(sigma2_theta.p_1_fixed),
             mean = Etheta.p_1_cond, sd = sqrt(Vtheta.p_1_cond), lower.tail = FALSE)
  A2 = pnorm(crit * sqrt(sigma2_theta.p_2_fixed),
             mean = Etheta.p_2_cond, sd = sqrt(Vtheta.p_2_cond), lower.tail = FALSE)

  # Use predictive distribution of theta.2_1 given theta.1_1
  z_A1 = qnorm(A1, lower.tail = FALSE)
  z_A2 = qnorm(A2, lower.tail = FALSE)
  z_A12_1 = qnorm(w.2_1 * A12, lower.tail = FALSE)
  z_A12_2 = qnorm(w.2_2 * A12, lower.tail = FALSE)

  crit_A1 = z_A1 * sqrt(sigma2_theta.2_1)
  crit_A2 = z_A2 * sqrt(sigma2_theta.2_2)
  crit_A12_1 = z_A12_1 * sqrt(sigma2_theta.2_1)
  crit_A12_2 = z_A12_2 * sqrt(sigma2_theta.2_2)

  # Now calculate probability of rejection
  q1 <- mvtnorm::pmvnorm(lower = c(-Inf,                   max(crit_A12_2, crit_A2)),
                         upper = c(min(crit_A1, crit_A12_1),  Inf),
                         mean  = mu.1,
                         sigma = Sigma_ppd,
                         algorithm = mvtnorm::GenzBretz())[1]
  q2 <- mvtnorm::pmvnorm(lower = c(min(crit_A1, crit_A12_1), max(crit_A12_2, crit_A2)),
                         upper = c(max(crit_A1, crit_A12_1), Inf),
                         mean  = mu.1,
                         sigma = Sigma_ppd,
                         algorithm = mvtnorm::GenzBretz())[1]
  q3 <- mvtnorm::pmvnorm(lower = c(max(crit_A1, crit_A12_1), max(crit_A12_2, crit_A2)),
                         upper = c(Inf,                   Inf),
                         mean  = mu.1,
                         sigma = Sigma_ppd,
                         algorithm = mvtnorm::GenzBretz())[1]
  q4 <- mvtnorm::pmvnorm(lower = c(max(crit_A1, crit_A12_1), min(crit_A12_2, crit_A2)),
                         upper = c(Inf,                   max(crit_A12_2, crit_A2)),
                         mean  = mu.1,
                         sigma = Sigma_ppd,
                         algorithm = mvtnorm::GenzBretz())[1]
  q5 <- mvtnorm::pmvnorm(lower = c(max(crit_A1, crit_A12_1),   -Inf),
                         upper = c(Inf,                      min(crit_A12_2, crit_A2)),
                         mean  = mu.1,
                         sigma = Sigma_ppd,
                         algorithm = mvtnorm::GenzBretz())[1]
  q6 <- mvtnorm::pmvnorm(lower = c(min(crit_A1, crit_A12_1), min(crit_A12_2, crit_A2)),
                         upper = c(max(crit_A1, crit_A12_1), max(crit_A12_2, crit_A2)),
                         mean  = mu.1,
                         sigma = Sigma_ppd,
                         algorithm = mvtnorm::GenzBretz())[1]

  pr1  <- q2 + q3 + q4 + q5
  pr2  <- q1 + q2 + q3 + q4
  pr12 <- q2 + q3 + q4

  # Calculate expected utility. Note this should not depend on the estimates.
  # or true values in the function in this implementation.
  (gamma1(lambda = lambda)  * pr1 +
      gamma2(lambda = lambda)  * pr2 +
      gamma12(lambda = lambda) * pr12)
}







#' Find the optimized w.2 and lambda.2 for the second stage after 1st stage
#'
#' Uses the twoStage_optimizedWeight function to find the values that yield to the maximum utility.
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
fast_interimOptimalTrial = function(x, r,
                               sigma_1 = 1, sigma_2 = 1, rho = 0,
                               mu_1 = 0, psi_1 = 1,
                               mu_2 = 0, psi_2 = 1, rho0 = 0,
                               lambda = 0.5, lambda.1 = 0.5,
                               n = 100,
                               w.1 = 0.5,
                               alpha = 0.025,
                               gamma1  = ..gamma1..,  #USES SIMPLE UTILITIES
                               gamma2  = ..gamma2..,  #USES SIMPLE UTILITIES
                               gamma12 = ..gamma12.., #USES SIMPLE UTILITIES
                               verbose = TRUE,
                               tol.zero = 0.00001,
                               tol.lim = 0.0001, epsilon = 1e-3,
                               optim_func  = "optim",
                               grid_points = 10){
  # optim_func = match.arg(optim_func)




  # Use the optim function to find the values lambda.2 w.2  with maximum utility


  l1w1 = interimOptimalTrial_all(lambda.2 = 1 - tol.lim,
                                w.2 = 1,
                                x = x, r = r,
                                sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                                mu_1 = mu_1, psi_1 = psi_1,
                                mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                                lambda = lambda, lambda.1 = lambda.1,
                                n = n,
                                w.1 = w.1,
                                alpha = alpha,
                                gamma1  = gamma1,
                                gamma2  = gamma2,
                                gamma12 = gamma12,
                                verbose = verbose, show.case = TRUE,
                                tol.zero = tol.zero)
  l0w0 = interimOptimalTrial_all(lambda.2 = tol.lim,
                                w.2 = 0,
                                x = x, r = r,
                                sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                                mu_1 = mu_1, psi_1 = psi_1,
                                mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                                lambda = lambda, lambda.1 = lambda.1,
                                n = n,
                                w.1 = w.1,
                                alpha = alpha,
                                gamma1  = gamma1,
                                gamma2  = gamma2,
                                gamma12 = gamma12,
                                verbose = verbose, show.case = TRUE,
                                tol.zero = tol.zero)

  if (l1w1$p1elem < epsilon &  l0w0$p2elem < epsilon) {
    out = list(par = c(NA, NA),
               value = max(l1w1$EU, l0w0$EU),
               convergence = paste0("Utility < ", epsilon, " for all values"),
               feval = 2)
    return(out)
  }


  l5w5 = interimOptimalTrial_all(lambda.2 = .5,
                                w.2 = .5,
                                x = x, r = r,
                                sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                                mu_1 = mu_1, psi_1 = psi_1,
                                mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                                lambda = lambda, lambda.1 = lambda.1,
                                n = n,
                                w.1 = w.1,
                                alpha = alpha,
                                gamma1  = gamma1,
                                gamma2  = gamma2,
                                gamma12 = gamma12,
                                verbose = verbose, show.case = TRUE,
                                tol.zero = tol.zero)

  if (l5w5$p1elem > 1-epsilon &  l5w5$p2elem > 1-epsilon) {
    out = list(par = c(NA, NA),
               value = l5w5$EU,
               convergence = paste0("Utility > ", 1-epsilon, " for equally distributed weights and prevalences"),
               feval = 3)
    return(out)
  }

  if (l1w1$p1elem < epsilon) {
    # In this case, P(rejectH01) is 0, even if we put all the patients in subgr1
    out = stats::optim(par = c(0.5),
                fn = function(y, ...) interimOptimalTrial_EU(lambda.2 = tol.lim, w.2 = y, ...),
                x = x, r = r,
                sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                mu_1 = mu_1, psi_1 = psi_1,
                mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                lambda = lambda, lambda.1 = lambda.1,
                n = n,
                w.1 = w.1,
                alpha = alpha,
                gamma1  = gamma1,
                gamma2  = gamma2,
                gamma12 = gamma12,
                verbose = verbose, show.case = FALSE,
                tol.zero = tol.zero,
                method = "L-BFGS-B",
                lower = c(0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                upper = c(1),  # Constrain the search to the [0,1]X[0,1] space
                # lower = c(0, 0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                # upper = c(1, 1),  # Constrain the search to the [0,1]X[0,1] space
                control = list(fnscale = -1))
    out$par = c(tol.lim, out$par)
    out$OptimalTrial_message = "l1w1$p1elem < epsilon"
    return(out)
  }
  if (l0w0$p2elem < epsilon) {
    # In this case, P(rejectH02) is 0, even if we put all the patients in subgr2
    out = stats::optim(par = c(0.5),
                fn = function(y, ...)
                  interimOptimalTrial_EU(lambda.2 = 1 - tol.lim, w.2 = y, ...),
                x = x, r = r,
                sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                mu_1 = mu_1, psi_1 = psi_1,
                mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                lambda = lambda, lambda.1 = lambda.1,
                n = n,
                w.1 = w.1,
                alpha = alpha,
                gamma1  = gamma1,
                gamma2  = gamma2,
                gamma12 = gamma12,
                verbose = verbose, show.case = FALSE,
                tol.zero = tol.zero,
                method = "L-BFGS-B",
                lower = c(0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                upper = c(1),  # Constrain the search to the [0,1]X[0,1] space
                # lower = c(0, 0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                # upper = c(1, 1),  # Constrain the search to the [0,1]X[0,1] space
                control = list(fnscale = -1))
    out$par = c(1-tol.lim, out$par)
    out$OptimalTrial_message = "l0w0$p2elem < epsilon"
    return(out)
  }

  if(optim_func == "optim"){
    out = stats::optim(par = c(0.5,0.5),
                fn = function(y, ...) interimOptimalTrial_EU(lambda.2 = y[1], w.2 = y[2], ...),
                x = x, r = r,
                sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                mu_1 = mu_1, psi_1 = psi_1,
                mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                lambda = lambda, lambda.1 = lambda.1,
                n = n,
                w.1 = w.1,
                alpha = alpha,
                gamma1  = gamma1,
                gamma2  = gamma2,
                gamma12 = gamma12,
                verbose = verbose, show.case = FALSE,
                tol.zero = tol.zero,
                method = "L-BFGS-B",
                lower = c(tol.lim,     0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                upper = c(1 - tol.lim, 1),  # Constrain the search to the [0,1]X[0,1] space
                # lower = c(0, 0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                # upper = c(1, 1),  # Constrain the search to the [0,1]X[0,1] space
                control = list(fnscale = -1))
    # if (out$par[2] < tol.zero)      out$par[2] = 0
    # if (out$par[2] > 1 - tol.zero)  out$par[2] = 1
    if (out$par[2] < 0)  out$par[2] = 0
    if (out$par[2] > 1)  out$par[2] = 1
  } else if(optim_func == "hjkb"){
    out = dfoptim::hjkb(par = c(0.5,0.5),
                        fn = function(y, ...) interimOptimalTrial_EU(lambda.2 = y[1], w.2 = y[2], ...),
                        x = x, r = r,
                        sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                        mu_1 = mu_1, psi_1 = psi_1,
                        mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                        lambda = lambda, lambda.1 = lambda.1,
                        n = n,
                        w.1 = w.1,
                        alpha = alpha,
                        gamma1  = gamma1,
                        gamma2  = gamma2,
                        gamma12 = gamma12,
                        verbose = verbose, show.case = FALSE, tol.zero = tol.zero,
                        lower = c(tol.lim,     0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                        upper = c(1 - tol.lim, 1),  # Constrain the search to the [0,1]X[0,1] space
                        control = list(maximize = TRUE))
  } else if(optim_func == "Rvmmin"){
    out = Rvmmin::Rvmmin(par = c(0.5,0.5),
                         fn = function(y, ...) (-1)*interimOptimalTrial_EU(lambda.2 = y[1], w.2 = y[2], ...),
                         x = x, r = r,
                         sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                         mu_1 = mu_1, psi_1 = psi_1,
                         mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                         lambda = lambda, lambda.1 = lambda.1,
                         n = n,
                         w.1 = w.1,
                         alpha = alpha,
                         gamma1  = gamma1,
                         gamma2  = gamma2,
                         gamma12 = gamma12,
                         verbose = verbose, show.case = FALSE, tol.zero = tol.zero,
                         lower = c(tol.lim,     0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                         upper = c(1 - tol.lim, 1))
  } else if(optim_func == "grid"){
    out = grid_search(x = x, r = r,
                      sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                      mu_1 = mu_1, psi_1 = psi_1,
                      mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                      lambda = lambda, lambda.1 = lambda.1,
                      n = n,
                      w.1 = w.1,
                      alpha = alpha,
                      gamma1  = gamma1,
                      gamma2  = gamma2,
                      gamma12 = gamma12,
                      verbose = verbose,
                      tol.zero = tol.zero,
                      lower = c(tol.lim,     0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                      upper = c(1 - tol.lim, 1),
                      grid_points = grid_points)
  } else if(optim_func == "fast"){
    out = dfoptim::hjkb(par = c(0.5,0.5),
                        fn = function(y, ...) fast_interimOptimalTrial_EU(lambda.2_1 = y[1], w.2_1 = y[2], ...),
                        x = x,
                        theta.1_1 = x[1], theta.1_2 = x[2],
                        r = r,
                        sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                        mu = c(mu_1, mu_2),
                        mu_1 = mu_1, psi_1 = psi_1,
                        mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                        Sigma0  = matrix(c(psi_1^2,
                                           rho0 * psi_1 * psi_2,
                                           rho0 * psi_1 * psi_2,
                                           psi_2^2),
                                         nrow = 2),
                        lambda = lambda,
                        lambda.1_1 = lambda.1,
                        lambda.1_2 = 1 - lambda.1,
                        n = n,
                        w.1_1 = w.1,
                        w.1_2 = 1-w.1,
                        alpha = alpha,
                        gamma1  = gamma1,
                        gamma2  = gamma2,
                        gamma12 = gamma12,
                        verbose = verbose,
                        tol.zero = tol.zero,
                        lower = c(tol.lim,     0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                        upper = c(1 - tol.lim, 1),  # Constrain the search to the [0,1]X[0,1] space
                        control = list(maximize = TRUE))
  } else {
    out = dfoptim::nmkb(par = c(0.5,0.5),
                        fn = function(y, ...) interimOptimalTrial_EU(lambda.2 = y[1], w.2 = y[2], ...),
                        x = x, r = r,
                        sigma_1 = sigma_1, sigma_2 = sigma_2, rho = rho,
                        mu_1 = mu_1, psi_1 = psi_1,
                        mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                        lambda = lambda, lambda.1 = lambda.1,
                        n = n,
                        w.1 = w.1,
                        alpha = alpha,
                        gamma1  = gamma1,
                        gamma2  = gamma2,
                        gamma12 = gamma12,
                        verbose = verbose, show.case = FALSE, tol.zero = tol.zero,
                        lower = c(tol.lim,     0),  # Lambda cannot be 0 or 1, but w.2 can be 0 or 1
                        upper = c(1 - tol.lim, 1),  # Constrain the search to the [0,1]X[0,1] space
                        control = list(maximize = TRUE))
  }
  out
}






#' Expected utility for values of lambda.2 and w.2 after 1st stage
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
interimOptimalTrial_all <- function(lambda.2, w.2,
                                   x, r,
                                   sigma_1 = 1, sigma_2 = 1, rho = 0,
                                   mu_1 = 0, psi_1 = 1,
                                   mu_2 = 0, psi_2 = 1, rho0 = 0,
                                   lambda = 0.5, lambda.1 = 0.5,
                                   n = 100,
                                   w.1 = 0.5,
                                   alpha = 0.025,
                                   gamma1  = ..gamma1..,  #USES SIMPLE UTILITIES
                                   gamma2  = ..gamma2..,  #USES SIMPLE UTILITIES
                                   gamma12 = ..gamma12.., #USES SIMPLE UTILITIES
                                   verbose = TRUE, show.case = TRUE, tol.zero = 1e-10){

  # x is a 2D vector with the estimates from 1st stage.
  theta.1_1 <- (x[1])
  theta.1_2 <- (x[2])
  mu = c(mu_1, mu_2)
  # When lambda.2 and w.2 are very close to 0 or 1, we just set them to 0 or 1 respectively.
  # if (lambda.2 < tol.zero) lambda.2 = 0
  # if (lambda.2 > 1 - tol.zero) lambda.2 = 1
  # if (w.2 < tol.zero)      w.2 = 0
  # if (w.2 > 1 - tol.zero)  w.2 = 1
  # if (w.2 < 0)  return(0)
  # if (w.2 > 1)  return(0)
  if (w.2 < 0)  {
    # cat("\nWarning! w.2<0 (", w.2, ") and will be set to 0\n\n")
    w.2=0
  }
  if (w.2 > 1)  {
    # cat("\nWarning! w.2>1 (", w.2, ") and will be set to 1\n\n")
    w.2 = 1
  }

  # if (lambda.2 <= 0) return(0)
  # if (lambda.2 >= 1) return(0)
  # Calculate trial prevalences
  lambda.1_1 = lambda.1
  lambda.1_2 = (1 - lambda.1)
  lambda.2_1 = lambda.2
  lambda.2_2 = (1 - lambda.2)

  # Calculate trial weights
  w.1_1 = w.1
  w.1_2 = (1-w.1)
  w.2_1 = w.2
  w.2_2 = (1-w.2)

  # Subgroup: Conditional expectation given 1st stage
  if (lambda.1_1 != 0){ # Case when subgroup 1 trial prevalence is > 0 in first stage
    psi2_1 = psi_1^2
    sigma2_1 = sigma_1^2
    psi2.1_1 = sigma2_1 * psi2_1 / (lambda.1_1 * r * n * psi2_1 / 4 + sigma2_1)
    sigma2_theta.1_1 = 4 * sigma2_1 / (lambda.1_1 * r * n)
    sigma2_theta.2_1_fixed = 4 * sigma2_1 / (lambda.1_1 * (1-r) * n)
    sigma2_theta.2_1       = 4 * sigma2_1 / (lambda.2_1 * (1-r) * n)
    sigma2_theta.p_1_fixed = 4 * sigma2_1 / n / (lambda.1_1 * r + lambda.1_1 * (1 - r))
    sigma2_theta.p_1       = 4 * sigma2_1 / n / (lambda.1_1 * r + lambda.2_1 * (1 - r))
    mu.1_1 = psi2.1_1 * (mu_1 / psi2_1 + theta.1_1 / sigma2_theta.1_1)
    Z.1_1 = theta.1_1 / (2 * sigma_1 / sqrt(n * lambda.1_1 * r))
    f.1_1 = lambda.1_1 * r / (lambda.1_1 * r  + lambda.2_1 * (1-r))
    f.2_1 = 1 - f.1_1
  } else { # If subgroup 1 trial prevalence is 0 in first stage, we only have 2nd stage data
    psi2_1 = psi_1^2
    sigma2_1 = sigma_1^2
    psi2.1_1 = psi2_1
    sigma2_theta.1_1 = 0
    sigma2_theta.2_1 = 4 * sigma2_1 / (lambda.2_1 * (1 - r) * n)
    sigma2_theta.p_1 = sigma2_theta.2_1
    mu.1_1 = mu_1
    Z.1_1 = NA
    f.1_1 = 0
    f.2_1 = 1
  }
  if (lambda.2_1 == 0){ # If subgroup 1 trial prevalence is 0 in second stage, we only have 1st stage data
    sigma2_theta.2_1 = 0
  }

  Etheta.p_1_cond = r * theta.1_1
  Vtheta.p_1_cond = (1 - r)^2 * (sigma2_theta.2_1_fixed)

  # Complement: Conditional expectation given 2nd stage
  if (lambda.1_2 != 0){# Case when subgroup 2 trial prevalence is >0 in first stage
    sigma2_2 = sigma_2^2
    psi2_2 = psi_2^2
    psi2.1_2 = sigma2_2 * psi2_2 / (lambda.1_2 * r * n * psi2_2 / 4 + sigma2_2)
    sigma2_theta.1_2 =       4 * sigma2_2 / (lambda.1_2 * r * n)
    sigma2_theta.2_2_fixed = 4 * sigma2_2 / (lambda.1_2 * (1 - r) * n)
    sigma2_theta.2_2 =       4 * sigma2_2 / (lambda.2_2 * (1 - r) * n)
    sigma2_theta.p_2_fixed = 4 * sigma_2^2 / n / (lambda.1_2 * r + lambda.1_2 * (1 - r))
    sigma2_theta.p_2 =       4 * sigma_2^2 / n / (lambda.1_2 * r + lambda.2_2 * (1 - r))
    mu.1_2 = psi2.1_2 * (mu_2 / psi2_2 + theta.1_2 / sigma2_theta.1_2)
    Z.1_2 = theta.1_2 / (2 * sigma_2 / sqrt(n * lambda.1_2 * r))
    f.1_2 = lambda.1_2 * r / (lambda.1_2 * r  + lambda.2_2 * (1-r))
    f.2_2 = 1 - f.1_2
  } else { # If subgroup 2 trial prevalence is 0 in first stage, we only have 2nd stage data
    psi2_2 = psi_2^2
    sigma2_2 = sigma_2^2
    psi2.1_2 = psi2_2
    sigma2_theta.1_2 = 0
    sigma2_theta.2_2 = 4 * sigma2_2 / (lambda.2_2 * (1 - r) * n)
    sigma2_theta.p_2 = sigma2_theta.2_2
    mu.1_2 = mu_2
    Z.1_2 = NA
    f.1_2 = 0
    f.2_2 = 1
  }
  if (lambda.2_2 == 0){# If subgroup 2 trial prevalence is 0 in second stage, we only have 1st stage data
    sigma2_theta.2_2 = 0
  }

  Etheta.p_2_cond = r * theta.1_2
  Vtheta.p_2_cond = (1 - r)^2 * (sigma2_theta.2_2_fixed)

  Sigma0  = matrix(c(psi2_1,
                     rho0 * psi_1 * psi_2,
                     rho0 * psi_1 * psi_2,
                     psi2_2),
                   nrow = 2)
  Sigma.1_theta  = matrix(c(sigma2_theta.1_1,
                            0 , 0,
                            sigma2_theta.1_2),
                          nrow = 2)
  Sigma.2_theta  = matrix(c(sigma2_theta.2_1,
                            0 , 0,
                            sigma2_theta.2_2),
                          nrow = 2)

  SigmaPlusInv = solve(Sigma0 + Sigma.1_theta)
  Sigma.1 = Sigma0 %*% SigmaPlusInv %*% Sigma.1_theta
  mu.1    = drop(Sigma0 %*% SigmaPlusInv %*% drop(x) + Sigma.1_theta %*% SigmaPlusInv %*% mu)
  mu.1_1  = mu.1[1]
  mu.1_2  = mu.1[2]

  # Covariance matrix of the prior predictive distribution
  Sigma_ppd = Sigma.1 + Sigma.2_theta

  # Calculate 1st stage p.values
  P.1_1 = pnorm(Z.1_1, lower.tail = FALSE)
  P.1_2 = pnorm(Z.1_2, lower.tail = FALSE)

  if (verbose) {
    cat(paste("Essential values", "lambda.2", lambda.2, "w.2", w.2, "\n"))
    cat(paste("First Stage estimates", theta.1_1, theta.1_2, "\n"))
    cat(paste("First Stage Z values", Z.1_1, Z.1_2, "\n"))
  }

  # Calculate A_1, A_2, B_1, B_2
  # First calculate critical boundaries for first stage
  crit_w.1_1 = qnorm(1 - alpha * w.1_1)
  crit_w.1_2 = qnorm(1 - alpha * (1-w.1_1))
  crit       = qnorm(1 - alpha)
  # Based on estiamtes from first stage we derive the distributions of
  # the conditional statistics
  # Compute conditional error rates
  A12_1 = pnorm(crit_w.1_1 * sqrt(sigma2_theta.p_1_fixed),
                mean = Etheta.p_1_cond, sd = sqrt(Vtheta.p_1_cond), lower.tail = FALSE)
  A12_2 = pnorm(crit_w.1_2 * sqrt(sigma2_theta.p_2_fixed),
                mean = Etheta.p_2_cond, sd = sqrt(Vtheta.p_2_cond), lower.tail = FALSE)
  A12 = A12_1 + A12_2 - A12_1*A12_2

  A1 = pnorm(crit * sqrt(sigma2_theta.p_1_fixed),
             mean = Etheta.p_1_cond, sd = sqrt(Vtheta.p_1_cond), lower.tail = FALSE)
  A2 = pnorm(crit * sqrt(sigma2_theta.p_2_fixed),
             mean = Etheta.p_2_cond, sd = sqrt(Vtheta.p_2_cond), lower.tail = FALSE)



  # Use predictive distribution of theta.2_1 given theta.1_1
  z_A1 = qnorm(A1, lower.tail = FALSE)
  z_A2 = qnorm(A2, lower.tail = FALSE)
  z_A12_1 = qnorm(w.2_1 * A12, lower.tail = FALSE)
  z_A12_2 = qnorm(w.2_2 * A12, lower.tail = FALSE)

  crit_A1 = z_A1 * sqrt(sigma2_theta.2_1)
  crit_A2 = z_A2 * sqrt(sigma2_theta.2_2)
  crit_A12_1 = z_A12_1 * sqrt(sigma2_theta.2_1)
  crit_A12_2 = z_A12_2 * sqrt(sigma2_theta.2_2)

  # Calculate the probabilities of rejection using the prior distribution
  case = 5
  # Check whether we have a two stage trial of one stage (r=1)
  isTwoStageTrial = (r != 1)
  if (isTwoStageTrial == FALSE){
    # If we have a one stage trials, simply test the hypotheses using 1st stage p-vaues
    if (lambda.1_1 == 0) P.1_1 = 1
    if (lambda.1_2 == 0) P.1_2 = 1
    pr1  <- 1*(P.1_1 < alpha & (P.1_1 < w.2_1 * alpha | P.1_1 < w.2_2 * alpha))
    pr2  <- 1*(P.1_1 < alpha & (P.1_1 < w.2_1 * alpha | P.1_1 < w.2_2 * alpha))
    pr12 <- pr1 * pr2
  }
  # Case: 0 < lambda.2_1 < 1 -----------------------------------------------
  if (isTwoStageTrial & (0 < lambda.2_1 & lambda.2_1 < 1)){
    # Now calculate probability of rejection
    q1 <- mvtnorm::pmvnorm(lower = c(-Inf,                   max(crit_A12_2, crit_A2)),
                           upper = c(min(crit_A1, crit_A12_1),  Inf),
                           mean  = mu.1,
                           sigma = Sigma_ppd,
                           algorithm = mvtnorm::GenzBretz())[1]
    q2 <- mvtnorm::pmvnorm(lower = c(min(crit_A1, crit_A12_1), max(crit_A12_2, crit_A2)),
                           upper = c(max(crit_A1, crit_A12_1), Inf),
                           mean  = mu.1,
                           sigma = Sigma_ppd,
                           algorithm = mvtnorm::GenzBretz())[1]
    q3 <- mvtnorm::pmvnorm(lower = c(max(crit_A1, crit_A12_1), max(crit_A12_2, crit_A2)),
                           upper = c(Inf,                   Inf),
                           mean  = mu.1,
                           sigma = Sigma_ppd,
                           algorithm = mvtnorm::GenzBretz())[1]
    q4 <- mvtnorm::pmvnorm(lower = c(max(crit_A1, crit_A12_1), min(crit_A12_2, crit_A2)),
                           upper = c(Inf,                   max(crit_A12_2, crit_A2)),
                           mean  = mu.1,
                           sigma = Sigma_ppd,
                           algorithm = mvtnorm::GenzBretz())[1]
    q5 <- mvtnorm::pmvnorm(lower = c(max(crit_A1, crit_A12_1),   -Inf),
                           upper = c(Inf,                      min(crit_A12_2, crit_A2)),
                           mean  = mu.1,
                           sigma = Sigma_ppd,
                           algorithm = mvtnorm::GenzBretz())[1]
    q6 <- mvtnorm::pmvnorm(lower = c(min(crit_A1, crit_A12_1), min(crit_A12_2, crit_A2)),
                           upper = c(max(crit_A1, crit_A12_1), max(crit_A12_2, crit_A2)),
                           mean  = mu.1,
                           sigma = Sigma_ppd,
                           algorithm = mvtnorm::GenzBretz())[1]

    if        (crit_A1 <= crit_A12_1 & crit_A2 <= crit_A12_2){
      pr1  <- q2 + q3 + q4 + q5
      pr2  <- q1 + q2 + q3 + q4
      pr12 <- q2 + q3 + q4
      case <- 1
    } else if (crit_A1 >  crit_A12_1 & crit_A2 <= crit_A12_2){
      pr1  <- q3 + q4 + q5
      pr2  <- q1 + q2 + q3 + q4 + q6
      pr12 <- q3 + q4
      case <- 2
    } else if (crit_A1 <= crit_A12_1 & crit_A2  > crit_A12_2){
      pr1  <- q2 + q3 + q4 + q5 + q6
      pr2  <- q1 + q2 + q3
      pr12 <- q2 + q3
      case <- 3
    } else if (crit_A1  > crit_A12_1 & crit_A2  > crit_A12_2){
      pr1  <- q3 + q4 + q5
      pr2  <- q1 + q2 + q3
      pr12 <- q3
      case <- 4
    }

    if(verbose){
      cat("q1:",q1,"q2:",q2,"q3:",q3,"q4:",q4,"q5:",q5,"q6:",q6, "\n")
      cat("Probability of rejection H01",pr1,"\n")
      cat("Probability of rejection H02",pr2,"\n")
      cat("Probability of rejection H01^H02",pr12,"\n")
      cat("Power = ",(pr1 + pr2 - pr12), "\n")
    }
  }

  p1elem <- mvtnorm::pmvnorm(lower = c(crit_A1, -Inf),
                             upper = c(Inf,      Inf),
                             mean  = mu.1,
                             sigma = Sigma_ppd,
                             algorithm = mvtnorm::GenzBretz())[1]
  p2elem <- mvtnorm::pmvnorm(lower = c(-Inf, crit_A2),
                             upper = c( Inf,  Inf),
                             mean  = mu.1,
                             sigma = Sigma_ppd,
                             algorithm = mvtnorm::GenzBretz())[1]


  # Calculate expected utility. Note this should not depend on the estimates.
  # or true values in the function in this implementation.
  EU <- (gamma1(lambda = lambda)  * pr1 +
           gamma2(lambda = lambda)  * pr2 +
           gamma12(lambda = lambda) * pr12)
  if(show.case) {
    out = list(EU = EU, case = case,
               p1elem = p1elem, p2elem = p2elem,
               pr1 = pr1,
               pr2 = pr2,
               pr12 = pr12)
  } else {
    out = EU
  }
  out
}

#' Simulates a clinical trial with optimized sample allocation.
#'
#' Given true values for theta_1 and theta_2, the function generates first stage data,
#' then find optimal values lambda.2 and w.2, generates second stage data, and
#' finally test the hypotheses.
#'
#' @param r.1      proportion of patients in 1st stage of the trial
#' @param lambda   proportion of patients in sub-population 1 in the overal patient population
#' @param lambda.1 proportion of patients in sub-population 1 in the trial 1st stage
#' @param w.1      factor for allocation of alpha to H01 in the testing correction (1st stage)
#' @param theta_1  true parameter value for sub-population 1
#' @param theta_2  true parameter value for sub-population 2
#' @param mu_1     prior mean for the treatment effect in sub-population 1
#' @param psi_1    prior standard deviation for the treatment effect in sub-population 1
#' @param mu_2     prior mean for the treatment effect in sub-population 2
#' @param psi_2    prior standard deviation for the treatment effect in sub-population 2
#' @param rho_0    prior correlation between treatment effect in sub-population 1 and sub-population 2 (ignored)
#' @param alpha    significance level
#' @param n        trial sample size
#' @param gamma1   gain when rejecting H01. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param gamma2   gain when rejecting H02. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param gamma12  gain when rejecting H01 and H02. A function with parameters (mu1, mu2, theta1, theta2, lambda)
#' @param tol.zero Absolute distance to 0 or 1 to be considered 0 or 1
#' @param verbose  Print intermediate values for calculations
#'
#' @export
fast_fullOptimization <- function(M = 1000, mc.cores = 40,
                             r.1      =   0.1,
                             lambda   =   0.5,
                             lambda.1 =   0.5,
                             w.1      =   0.5,
                             mu_1     =   0,
                             psi_1    =   1,
                             mu_2     =   0,
                             psi_2    =   1,
                             rho0     =   0,
                             alpha    =   0.05,
                             n        = 100,
                             gamma1   = ..gamma1..,  #USES SIMPLE UTILITIES
                             gamma2   = ..gamma2..,  #USES SIMPLE UTILITIES
                             gamma12  = ..gamma12.., #USES SIMPLE UTILITIES
                             tol.zero =   0.00001,
                             tol.lim  =   0.0001,
                             verbose  = TRUE,
                             summarize = TRUE,
                             optim_func = "hjkb",
                             grid_points = 10){
  # Simulate 1st stage estimates
  Sigma.1     = diag(c(4/lambda.1/n/r.1, 4/(1-lambda.1)/n/r.1))

  Sigma.0 = matrix(c(psi_1^2,
                     rho0 * psi_1 * psi_2,
                     rho0 * psi_1 * psi_2,
                     psi_2^2),
                   nrow = 2)

  theta.hat.1 = mvtnorm::rmvnorm(n = M,
                                 mean  = c(mu_1, mu_2),
                                 sigma = Sigma.1 + Sigma.0)
  if (verbose) cat("Simulate 1st stage: DONE \n")

  # Given First stage estimates, calculate the values of lambda.2 and w.2
  #  with the maximum expected utility
  res = parallel::mclapply(X = 1:M, FUN = function(i) {
    interimDecision = fast_interimOptimalTrial(x = theta.hat.1[i, ],
                                          r = r.1,
                                          sigma_1 = 1, sigma_2 = 1,
                                          rho  = 0,
                                          mu_1 = mu_1, psi_1 = psi_1,
                                          mu_2 = mu_2, psi_2 = psi_1,
                                          rho0 = rho0,
                                          lambda = lambda,
                                          lambda.1 = lambda.1,
                                          n   = n,
                                          w.1 = w.1,
                                          alpha = alpha,
                                          gamma1 = gamma1,
                                          gamma2 = gamma2,
                                          gamma12 = gamma12,
                                          tol.zero = tol.zero,
                                          tol.lim = tol.lim,
                                          optim_func = optim_func,
                                          grid_points = grid_points,
                                          verbose = FALSE)
    if (verbose) cat("Interim decision: DONE \n")
    # Store values of lambda.2 and w.2
    c(lambda.1 = lambda.1,
      w.1 = w.1,
      r.1 = r.1,
      lambda.2 = interimDecision$par[1],
      w.2      = interimDecision$par[2],
      EU       = interimDecision$value)
  }, mc.cores = mc.cores)
  res = do.call(rbind, res)
  out = res
  if(summarize){
    res2 = colMeans(res)
    out = res2
  }
  out
}
