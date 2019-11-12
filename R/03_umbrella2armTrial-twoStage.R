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
umbrella2arm1Trial <- function(r = 0.5,
                           theta_1 = 0, theta_2 = 0,
                           sigma_1 = 1, sigma_2 = 1,
                           rho = 0,
                           mu_1 = 0, psi_1 = 1,
                           mu_2 = 0, psi_2 = 1,
                           rho0 = 0,
                           lambda   = c(0.15),
                           lambda.1 = c(0.15),
                           n = 1484,
                           alpha = c(0.05, 0.05),
                           gamma1  = ..umbrella_gamma1..,  #USES SIMPLE UTILITIES
                           gamma2  = ..umbrella_gamma2..,  #USES SIMPLE UTILITIES
                           verbose = TRUE,
                           tol.zero = 1e-10,
                           use.prior.as.true = FALSE,
                           sig.digits = 4){

  if(use.prior.as.true){
    psi2_1 = psi_1^2
    psi2_2 = psi_2^2
    Sigma0  = matrix(c(psi2_1,               rho0 * psi_1 * psi_2,
                       rho0 * psi_1 * psi_2,               psi2_2),
                     nrow = 2)
    true.mean = mvtnorm::rmvnorm(n = 1,
                                 mean = c(mu_1, mu_2),
                                 sigma = Sigma0)
  } else {
    true.mean = c(theta_1, theta_2)
  }

  lambda.1_1 = round(lambda.1[1], sig.digits)
  lambda.1_2 = 1-lambda.1_1

  if (r != 0){
    # Simulate 1st stage estimates
    Sigma.1     = diag(c(4/lambda.1_1/n/r,
                         4/lambda.1_2/n/r))
    theta.hat.1 = mvtnorm::rmvnorm(n = 1,
                                   mean  = true.mean,
                                   sigma = Sigma.1)
  } else {
    theta.hat.1 = c(NA, NA)
  }
  # Simulate 1st stage: DONE ---------------------------------------------------
  if (verbose) cat("Simulate 1st stage: DONE \n")

  # Given First stage estimates, calculate the values of lambda.2 and w.2
  #  with the maximum expected utility
  if (r != 1 & r != 0){
    interimDecision = interimUmbrella2armTrial(x = theta.hat.1,
                                           r          = r,
                                           lambda     = lambda,
                                           lambda.1   = lambda.1,
                                           mu_1       = mu_1,
                                           mu_2       = mu_2,
                                           psi_1      = psi_1,
                                           psi_2      = psi_2,
                                           alpha      = alpha,
                                           rho0       = rho0,
                                           n          = n,
                                           gamma1  = gamma1,
                                           gamma2  = gamma2,
                                           verbose    = FALSE,
                                           sig.digits = sig.digits)

    # Store values of lambda.2 and w.2
    lambda.2   = c(interimDecision$par[1], 1-interimDecision$par[1])
    lambda.2_1 = round(interimDecision$par[1], sig.digits)
    lambda.2_2 = round(1-interimDecision$par[1], sig.digits)
    interimEU  = interimDecision$value
  } else {
    lambda.2   = c(lambda.1_1, lambda.1_2)
    lambda.2_1 = lambda.1_1
    lambda.2_2 = lambda.1_2
    interimEU = NA
  }

  # When lambda.2 and w.2 are very close to 0 or 1, we just set them to 0 or 1 respectively.
  # if (lambda.2 < tol.zero)   lambda.2 = 0
  # if (lambda.2 > 1-tol.zero) lambda.2 = 1
  # if (w.2 < tol.zero)        w.2 = 0
  # if (w.2 > 1-tol.zero)      w.2 = 1

  # Simulates second stage estimates given the optimized lambda.2
  # Calculates first the variances
  Sigma.2_1 = 4 / (lambda.2_1) / n / (1 - r)
  Sigma.2_2 = 4 / (lambda.2_2) / n / (1 - r)
  # And simulate the estimates. A NA pre-allocated vector is used so that
  # We only simulate the value if it is actually needed.
  theta.hat.2 = c(NA, NA)
  if(r < 1){
    if (lambda.2_1 > 0) theta.hat.2[1] = rnorm(n = 1, mean  = true.mean[1], sd = sqrt(Sigma.2_1))
    if (lambda.2_2 > 0) theta.hat.2[2] = rnorm(n = 1, mean  = true.mean[2], sd = sqrt(Sigma.2_2))
  }
  if (verbose) cat("Simulate 2nd stage: DONE \n")
  # Simulate 2nd stage: DONE ---------------------------------------------------
  if (verbose) {
    cat(paste("lambda.1: ", paste0(round(c(lambda.1_1, lambda.1_2), 4), collapse = " - "), "\n"))
    cat(paste("theta_hat.1: ", paste0(round(c(theta.hat.1), 4), collapse = " - "), "\n"))
    cat(paste("lambda.2: ", paste0(c(lambda.2_1, lambda.2_2), collapse = " - "), "\n"))
    cat(paste("lambda.2: ", paste0(c(lambda.2), collapse = " - "), "\n"))
    cat(paste("theta_hat.2: ", paste0(round(c(theta.hat.2), 4), collapse = " - "), "\n"))
  }
  # Test the hypotheses given first and second stage data.
  decision = decisionUmbrella2armTrial(lambda.2 = lambda.2_1,
                                        x.1      = theta.hat.1,
                                        x.2      = theta.hat.2,
                                        r        = r,
                                        lambda   = lambda,
                                        lambda.1 = lambda.1,
                                        n        = n,
                                        alpha    = alpha,
                                        gamma1   = gamma1,
                                        gamma2   = gamma2,
                                        verbose  = verbose,
                                        sig.digits = sig.digits)
  if (verbose) cat("Decision: DONE \n")
  out = c(decision,
          lambda       = lambda,
          lambda.1_1    = lambda.1_1,
          lambda.2_1    = lambda.2_1,
          theta.hat.1_ = theta.hat.1,
          theta.hat.2_ = theta.hat.2,
          alpha_       = alpha,
          interimEU    = interimEU)
  class(out) = "simOptimalTrial"
  out
}
#' Simulates M clinical trials with optimized sample allocation using simulate1trial
#'
#' @exportClass simOptimalTrial
#' @export
umbrella2armMTrials = function(nsim   =  1000,
                          mc.cores =   1,
                          r = 0.5,
                          theta_1 = 0, theta_2 = 0,
                          sigma_1 = 1, sigma_2 = 1,
                          rho = 0,
                          mu_1 = 0, psi_1 = 1,
                          mu_2 = 0, psi_2 = 1,
                          rho0 = 0,
                          lambda   = c(0.15),
                          lambda.1 = c(0.15),
                          n = 1484,
                          alpha = c(0.05, 0.05),
                          gamma1  = ..umbrella_gamma1..,  #USES SIMPLE UTILITIES
                          gamma2  = ..umbrella_gamma2..,  #USES SIMPLE UTILITIES
                          verbose = FALSE,
                          tol.zero = 1e-10,
                          summarize = TRUE,
                          use.prior.as.true = FALSE,
                          sig.digits = 4){
  results = parallel::mclapply(1:nsim, FUN = function(x){
    umbrella2arm1Trial(r = r,
                   theta_1 = theta_1, theta_2 = theta_2,
                   sigma_1 = sigma_1, sigma_2 = sigma_2,
                   rho = rho,
                   mu_1 = mu_1, psi_1 = psi_1,
                   mu_2 = mu_2, psi_2 = psi_2,rho0 = rho0,
                   lambda   = lambda,
                   lambda.1 = lambda.1,
                   n = n,
                   alpha = alpha,
                   gamma1  = gamma1,
                   gamma2  = gamma2,
                   verbose = verbose,
                   tol.zero = tol.zero,
                   use.prior.as.true = use.prior.as.true,
                   sig.digits = sig.digits)
  }, mc.cores = mc.cores)

  res = data.frame(do.call(rbind,results))
  out = res
  if(summarize){
    res2 = colMeans(res, na.rm = TRUE)
    out = t(res2)
  }
  results.df = list(simresults = out,
                    summarize = summarize)
  class(results.df) = "simOptimalTrial"
  results.df
}

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
interimUmbrella2armTrial_EU <- function(lambda.2,
                                   x = NULL,
                                   r = 0.5,
                                   sigma_1 = 1, sigma_2 = 1,
                                   rho = 0,
                                   mu_1 = 0, psi_1 = 1,
                                   mu_2 = 0, psi_2 = 1,
                                   rho0 = 0,
                                   lambda   = c(0.15),
                                   lambda.1 = c(0.15),
                                   n = 1484,
                                   alpha = c(0.05, 0.05),
                                   gamma1  = ..umbrella_gamma1..,  #USES SIMPLE UTILITIES
                                   gamma2  = ..umbrella_gamma2..,  #USES SIMPLE UTILITIES
                                   verbose = TRUE,
                                   tol.zero = 1e-10,
                                   forOpt = TRUE, sig.digits = 4){
  if(r==0) stop("Multiple stage only")
  if(is.null(x)) stop("Please provide first stage data x")

  # x is a 3D vector with the estimates from 1st stage.
  theta.1_1 <- (x[1])
  theta.1_2 <- (x[2])
  mu = c(mu_1, mu_2)

  # When lambda.2 and w.2 are very close to 0 or 1, we just set them to 0 or 1 respectively.
  # if (lambda.2 < tol.zero) lambda.2 = 0
  # if (lambda.2 > 1 - tol.zero) lambda.2 = 1
  # if (lambda.2 <= 0) return(0)
  # if (lambda.2 >= 1) return(0)

  # Calculate trial prevalences in 1st Stage
  lambda.1_1 = abs(round(lambda.1[1], sig.digits))
  lambda.1_2 = 1 - lambda.1_1
  # Calculate trial prevalences in 2nd Stage
  lambda.2_1 = abs(round(lambda.2[1], sig.digits))
  lambda.2_2 = 1 - lambda.2_1

  # Calculate trial prevalences in population
  lambda_1 = lambda[1]
  lambda_2 = 1 - lambda_1

  # 1st Arm: Conditional expectation given 1st stage ---------------------------
  psi2_1   = psi_1^2
  sigma2_1 = sigma_1^2
  if (lambda.1_1 != 0){ # Case when subgroup 1 trial prevalence is > 0 in first stage
    psi2.1_1               = sigma2_1 * psi2_1 / (lambda.1_1 * r * n * psi2_1 / 4 + sigma2_1)
    sigma2_theta.1_1       = 4 * sigma2_1 / (lambda.1_1 * r * n)
    sigma2_theta.2_1_fixed = 4 * sigma2_1 / (lambda.1_1 * (1-r) * n)
    sigma2_theta.2_1       = 4 * sigma2_1 / (lambda.2_1 * (1-r) * n)
    sigma2_theta.p_1_fixed = 4 * sigma2_1 / n / (lambda.1_1 * r + lambda.1_1 * (1 - r))
    sigma2_theta.p_1       = 4 * sigma2_1 / n / (lambda.1_1 * r + lambda.2_1 * (1 - r))
    mu.1_1 = psi2.1_1 * (mu_1 / psi2_1 + theta.1_1 / sigma2_theta.1_1)
    Z.1_1  = theta.1_1 / (2 * sigma_1 / sqrt(n * lambda.1_1 * r))
  } else { # If subgroup 1 trial prevalence is 0 in first stage, we only have 2nd stage data
    psi2.1_1 = psi2_1
    sigma2_theta.1_1 = 0
    sigma2_theta.2_1 = 4 * sigma2_1 / (lambda.2_1 * (1 - r) * n)
    sigma2_theta.p_1 = sigma2_theta.2_1
    mu.1_1 = mu_1
    Z.1_1 = NA
  }
  Etheta.p_1_cond = r * theta.1_1
  Vtheta.p_1_cond = (1 - r)^2 * (sigma2_theta.2_1_fixed)

  # 2st Arm: Conditional expectation given 1st stage ---------------------------
  psi2_2   = psi_2^2
  sigma2_2 = sigma_2^2
  if (lambda.1_2 != 0){ # Case when subgroup 1 trial prevalence is > 0 in first stage
    psi2.1_2               = sigma2_2 * psi2_2 / (lambda.1_2 * r * n * psi2_2 / 4 + sigma2_2)
    sigma2_theta.1_2       = 4 * sigma2_2 / (lambda.1_2 * r * n)
    sigma2_theta.2_2_fixed = 4 * sigma2_2 / (lambda.1_2 * (1-r) * n)
    sigma2_theta.2_2       = 4 * sigma2_2 / (lambda.2_2 * (1-r) * n)
    sigma2_theta.p_2_fixed = 4 * sigma2_2 / n / (lambda.1_2 * r + lambda.1_2 * (1 - r))
    sigma2_theta.p_2       = 4 * sigma2_2 / n / (lambda.1_2 * r + lambda.2_2 * (1 - r))
    mu.1_2 = psi2.1_2 * (mu_2 / psi2_2 + theta.1_2 / sigma2_theta.1_2)
    Z.1_2  = theta.1_2 / (2 * sigma_2 / sqrt(n * lambda.1_2 * r))
  } else { # If subgroup 1 trial prevalence is 0 in first stage, we only have 2nd stage data
    psi2.1_2 = psi2_2
    sigma2_theta.1_2 = 0
    sigma2_theta.2_2 = 4 * sigma2_2 / (lambda.2_2 * (1 - r) * n)
    sigma2_theta.p_2 = sigma2_theta.2_2
    mu.1_2 = mu_2
    Z.1_2 = NA
  }
  Etheta.p_2_cond = r * theta.1_2
  Vtheta.p_2_cond = (1 - r)^2 * (sigma2_theta.2_2_fixed)



  Sigma0  = matrix(c(psi2_1,               rho0 * psi_1 * psi_2,
                     rho0 * psi_1 * psi_2,               psi2_2),
                   nrow = 2)
  Sigma.1_theta  = matrix(c(sigma2_theta.1_1, 0,
                            0, sigma2_theta.1_2),
                          nrow = 2)
  Sigma.2_theta  = matrix(c(sigma2_theta.2_1, 0,
                            0, sigma2_theta.2_2),
                          nrow = 2)

  # Prior sigma for second stage is the same as prior as there is no first stage
  SigmaPlusInv = solve(Sigma0 + Sigma.1_theta)
  Sigma.1 = Sigma0 %*% SigmaPlusInv %*% Sigma.1_theta

  # Prior Mu for second stage is the updated mu using 1st stage
  mu.1    = drop(Sigma0 %*% SigmaPlusInv %*% drop(x) + Sigma.1_theta %*% SigmaPlusInv %*% mu)
  mu.1_1 = mu.1[1]
  mu.1_2 = mu.1[2]
  # Covariance matrix of the prior predictive distribution
  Sigma_ppd = Sigma.1 + Sigma.2_theta

  # Calculate 1st stage p.values
  P.1_1 = pnorm(Z.1_1, lower.tail = FALSE)
  P.1_2 = pnorm(Z.1_2, lower.tail = FALSE)

  if (verbose) {
    cat(paste("Essential values", "lambda.2: ", paste0(c(lambda.2_1, lambda.2_2), collapse = " - "), "\n"))
    cat(paste("Essential values", "theta_hat.1: ", paste0(c(theta.1_1, theta.1_2), collapse = " - "), "\n"))
    cat(paste("Essential values", "theta_hat.1: ", paste0(c(theta.1_1, theta.1_2), collapse = " - "), "\n"))
  }

  # Calculate A_1, A_2, B_1, B_2
  # First calculate critical boundaries for first stage
  crit_1 = qnorm(1 - alpha[1])
  crit_2 = qnorm(1 - alpha[2])
  # theta_crit_1 = crit_1 * sqrt(sigma2_theta.2_1)
  # theta_crit_2 = crit_2 * sqrt(sigma2_theta.2_2)
  # theta_crit_3 = crit_3 * sqrt(sigma2_theta.2_3)

  # First calculate critical boundaries for first stage
  A1 = pnorm(crit_1 * sqrt(sigma2_theta.p_1_fixed),
             mean = Etheta.p_1_cond, sd = sqrt(Vtheta.p_1_cond), lower.tail = FALSE)
  A2 = pnorm(crit_2 * sqrt(sigma2_theta.p_2_fixed),
             mean = Etheta.p_2_cond, sd = sqrt(Vtheta.p_2_cond), lower.tail = FALSE)

  # Use predictive distribution of theta.2_1 given theta.1_1
  z_A1 = qnorm(A1, lower.tail = FALSE)
  z_A2 = qnorm(A2, lower.tail = FALSE)
  crit_A1 = z_A1 * sqrt(sigma2_theta.2_1)
  crit_A2 = z_A2 * sqrt(sigma2_theta.2_2)

  # Check whether we have a two stage trial of one stage (r=1)
  isTwoStageTrial = (r != 1)
  if (isTwoStageTrial == FALSE){
    # If we have a one stage trials, simply test the hypotheses using 1st stage p-vaues
    if (lambda.1_1 == 0) P.1_1 = 1
    if (lambda.1_2 == 0) P.1_2 = 1
    pr1  <- 1*(P.1_1 < alpha[1])
    pr2  <- 1*(P.1_2 < alpha[2])
  } else{
    # Calculate the probabilities of rejection using the prior distribution
    # Now calculate probability of rejection
    Sigma_ppd[Sigma_ppd == Inf] = 1000 # Handling cases where subgroups are not selected
    # We just use a large variance to avoid error of pmvnorm.
    # The corresponding pr will be =0 if crit_A# in cases that the prevalence is set to 0
    pr1 <- mvtnorm::pmvnorm(lower = c(crit_A1, -Inf),
                            upper = c(Inf,      Inf),
                            mean  = mu.1,
                            sigma = Sigma_ppd,
                            algorithm = mvtnorm::GenzBretz())[1]
    pr2 <- mvtnorm::pmvnorm(lower = c(-Inf, crit_A2),
                            upper = c( Inf,     Inf),
                            mean  = mu.1,
                            sigma = Sigma_ppd,
                            algorithm = mvtnorm::GenzBretz())[1]
  }
  if(verbose){
    cat("Probability of rejection H01",pr1,"\n")
    cat("Probability of rejection H02",pr2,"\n")
  }

  # Calculate expected utility. Note this should not depend on the estimates.
  # or true values in the function in this implementation.
  EU <- (gamma1(lambda = lambda_1) * pr1 +
         gamma2(lambda = lambda_2) * pr2)

  if(forOpt){
    out = EU
    return(out)
  } else{
    out = c(EU = EU,
            lambda.2_1 = lambda.2_1, lambda.2_2 = lambda.2_2,
            pr1 = pr1, pr2 = pr2)
    return(out)
  }
}

#' Find the optimized lambda.2 for the trial
#'
#' Uses the interimUmbrella2armTrial_EU function to find the values that yield to the maximum utility.
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
interimUmbrella2armTrial = function(x = NULL,
                               r = 0.5,
                               sigma_1 = 1, sigma_2 = 1,
                               rho = 0,
                               mu_1 = 0, psi_1 = 1,
                               mu_2 = 0, psi_2 = 1,
                               rho0 = 0,
                               lambda   = c(0.15),
                               lambda.1 = c(0.15),
                               n = 1484,
                               alpha = c(0.05, 0.05),
                               gamma1  = ..umbrella_gamma1..,  #USES SIMPLE UTILITIES
                               gamma2  = ..umbrella_gamma2..,  #USES SIMPLE UTILITIES
                               verbose = FALSE,
                               tol.zero = 1e-10,
                               tol.lim = 0.0001,
                               optim_func  = "optimise", sig.digits = 4){
  if(r == 0){
    message("If r is 0 then use the single Stage Trial function")
  }

  if(optim_func=="optimise"){
    opt = stats::optimize(f = function(y, ...) interimUmbrella2armTrial_EU(lambda.2 = y, ...),
                          interval = c(0, 1),
                          maximum = TRUE,
                          x = x,
                          r = r,
                          sigma_1 = sigma_1, sigma_2 = sigma_2,
                          rho = rho,
                          mu_1 = mu_1, psi_1 = psi_1,
                          mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                          lambda   = lambda,
                          lambda.1 = lambda.1,
                          n = n,
                          alpha   = alpha,
                          gamma1  = gamma1,  #USES SIMPLE UTILITIES
                          gamma2  = gamma2,  #USES SIMPLE UTILITIES
                          verbose = verbose,
                          tol.zero = tol.zero,
                          forOpt = TRUE, sig.digits = sig.digits)
    optimal_values = list(par = opt$maximum, value = opt$objective)
  } else if(optim_func=="optim"){
    optimal_values = stats::optim(par = lambda.1,
                                  f = function(y, ...) interimUmbrella2armTrial_EU(lambda.2 = y, ...),
                                  lower = 0,
                                  upper = 1,
                                  method = "L-BFGS-B",
                                  control = list(fnscale = -1),
                                  x = x,
                                  r = r,
                                  sigma_1 = sigma_1, sigma_2 = sigma_2,
                                  rho = rho,
                                  mu_1 = mu_1, psi_1 = psi_1,
                                  mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                                  lambda   = lambda,
                                  lambda.1 = lambda.1,
                                  n = n,
                                  alpha   = alpha,
                                  gamma1  = gamma1,  #USES SIMPLE UTILITIES
                                  gamma2  = gamma2,  #USES SIMPLE UTILITIES
                                  verbose = verbose,
                                  tol.zero = tol.zero,
                                  forOpt = TRUE, sig.digits = sig.digits)
  }
  optimal_values
}


#' Find the optimized lambda.2 for the trial
#'
#' Uses the interimUmbrella2armTrial_EU function to find the values that yield to the maximum utility.
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
interimUmbrella2armTrial_grid = function(x = NULL,
                                    r = 0.5,
                                    sigma_1 = 1, sigma_2 = 1,
                                    rho = 0,
                                    mu_1 = 0, psi_1 = 1,
                                    mu_2 = 0, psi_2 = 1,
                                    rho0 = 0,
                                    lambda   = c(0.15),
                                    lambda.1 = c(0.15),
                                    n = 1484,
                                    alpha = c(0.05, 0.05),
                                    gamma1  = ..umbrella_gamma1..,  #USES SIMPLE UTILITIES
                                    gamma2  = ..umbrella_gamma2..,  #USES SIMPLE UTILITIES
                                    verbose = FALSE,
                                    tol.zero = 1e-10,
                                    tol.lim = 0.0001,
                                    optim_func  = "optim", sig.digits = 4){
  if(r == 0){
    message("If r is 0 then use the single Stage Trial function")
  }
  l1 = round(seq(0,1,0.01), 2)
  grid = data.frame(lambda_1 = l1)
  f.out = apply(grid, MARGIN = 1,
                function(y, ...) interimUmbrella2armTrial_EU(lambda.2 = y, ...),
                x = x,
                r = r,
                sigma_1 = sigma_1, sigma_2 = sigma_2,
                rho = rho,
                mu_1 = mu_1, psi_1 = psi_1,
                mu_2 = mu_2, psi_2 = psi_2, rho0 = rho0,
                lambda   = lambda,
                lambda.1 = lambda.1,
                n = n,
                alpha   = alpha,
                gamma1  = gamma1,  #USES SIMPLE UTILITIES
                gamma2  = gamma2,  #USES SIMPLE UTILITIES
                verbose = verbose,
                tol.zero = tol.zero,
                forOpt = TRUE, sig.digits = sig.digits)
  gridout = data.frame(grid, EU = f.out)
  gridout
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
decisionUmbrella2armTrial <- function(lambda.2,
                                 x.1,
                                 x.2,
                                 r,
                                 sigma_1 = 1, sigma_2 = 1,
                                 rho = 0,
                                 lambda   = c(0.15),
                                 lambda.1 = c(0.15),
                                 n = 1484,
                                 alpha = c(0.05, 0.05),
                                 gamma1  = ..umbrella_gamma1..,  #USES SIMPLE UTILITIES
                                 gamma2  = ..umbrella_gamma2..,  #USES SIMPLE UTILITIES
                                 fixed.weights = FALSE,
                                 verbose = TRUE,
                                 show.case = TRUE,
                                 sig.digits = 4){
  # x.1 and x.2 are 2D vectors with the estimates from 1st stage and 2nd stage resp.
  theta.1_1 <- (x.1[1])
  theta.1_2 <- (x.1[2])
  theta.2_1 <- (x.2[1])
  theta.2_2 <- (x.2[2])

  r = unname(r)
  lambda.1 = unname(lambda.1)
  lambda.2 = unname(lambda.2)
  # Calculate trial prevalences in 1st Stage
  lambda.1_1 = round(lambda.1[1], sig.digits)
  lambda.1_2 = 1-lambda.1_1
  # Calculate trial prevalences in 2nd Stage
  lambda.2_1 = round(lambda.2[1], sig.digits)
  lambda.2_2 = 1-lambda.2_1
  # Calculate trial prevalences in population
  lambda_1 = lambda[1]
  lambda_2 = 1 - lambda_1


  # 1st Arm: Conditional expectation given 1st stage ---------------------------
  sigma2_1 = sigma_1^2
  if (lambda.1_1 != 0){ # Case when subgroup 1 trial prevalence is > 0 in first stage
    sigma2_theta.1_1       = 4 * sigma2_1 / (lambda.1_1 * r * n)
    sigma2_theta.2_1_fixed = 4 * sigma2_1 / (lambda.1_1 * (1-r) * n)
    sigma2_theta.2_1       = 4 * sigma2_1 / (lambda.2_1 * (1-r) * n)
    sigma2_theta.p_1_fixed = 4 * sigma2_1 / n / (lambda.1_1 * r + lambda.1_1 * (1 - r))
    sigma2_theta.p_1       = 4 * sigma2_1 / n / (lambda.1_1 * r + lambda.2_1 * (1 - r))
    Z.1_1  = theta.1_1 / (2 * sigma_1 / sqrt(n * lambda.1_1 * r))
  } else { # If subgroup 1 trial prevalence is 0 in first stage, we only have 2nd stage data
    sigma2_theta.1_1 = 0
    sigma2_theta.2_1 = 4 * sigma2_1 / (lambda.2_1 * (1 - r) * n)
    sigma2_theta.p_1 = sigma2_theta.2_1
    Z.1_1 = NA
  }

  Etheta.p_1_cond = r * theta.1_1
  Vtheta.p_1_cond = (1 - r)^2 * (sigma2_theta.2_1_fixed)
  f.1_1 = lambda.1_1 * r / (lambda.1_1 * r  + lambda.2_1 * (1 - r))
  f.2_1 = 1 - f.1_1
  theta.p_1 = f.1_1 * theta.1_1 + f.2_1 * theta.2_1
  Z.p_1 = theta.p_1 / sqrt(sigma2_theta.p_1)
  P.p_1 = pnorm(Z.p_1, lower.tail = FALSE)

  P.1_1 = pnorm(Z.1_1, lower.tail = FALSE)
  if (lambda.2_1 == 0){# If Arm prevalence is 0 in second stage, we only have 1st stage data
    P.p_1 = P.1_1
  }

  # 2st Arm: Conditional expectation given 1st stage ---------------------------
  sigma2_2 = sigma_2^2
  if (lambda.1_2 != 0){ # Case when subgroup 1 trial prevalence is > 0 in first stage
    sigma2_theta.1_2       = 4 * sigma2_2 / (lambda.1_2 * r * n)
    sigma2_theta.2_2_fixed = 4 * sigma2_2 / (lambda.1_2 * (1-r) * n)
    sigma2_theta.2_2       = 4 * sigma2_2 / (lambda.2_2 * (1-r) * n)
    sigma2_theta.p_2_fixed = 4 * sigma2_2 / n / (lambda.1_2 * r + lambda.1_2 * (1 - r))
    sigma2_theta.p_2       = 4 * sigma2_2 / n / (lambda.1_2 * r + lambda.2_2 * (1 - r))
    Z.1_2  = theta.1_2 / (2 * sigma_2 / sqrt(n * lambda.1_2 * r))
  } else { # If subgroup 1 trial prevalence is 0 in first stage, we only have 2nd stage data
    sigma2_theta.1_2 = 0
    sigma2_theta.2_2 = 4 * sigma2_2 / (lambda.2_2 * (1 - r) * n)
    sigma2_theta.p_2 = sigma2_theta.2_2
    Z.1_2 = NA
  }
  Etheta.p_2_cond = r * theta.1_2
  Vtheta.p_2_cond = (1 - r)^2 * (sigma2_theta.2_2_fixed)
  f.1_2 = lambda.1_2 * r / (lambda.1_2 * r  + lambda.2_2 * (1 - r))
  f.2_2 = 1 - f.1_2
  theta.p_2 = f.1_2 * theta.1_2 + f.2_2 * theta.2_2
  Z.p_2 = theta.p_2 / sqrt(sigma2_theta.p_2)
  P.p_2 = pnorm(Z.p_2, lower.tail = FALSE)

  P.1_2 = pnorm(Z.1_2, lower.tail = FALSE)
  if (lambda.2_2 == 0){# If Arm prevalence is 0 in second stage, we only have 1st stage data
    P.p_2 = P.1_2
  }



  # Calculate A_1, A_2, B_1, B_2
  # First calculate critical boundaries for first stage
  crit_1 = qnorm(1 - alpha[1])
  crit_2 = qnorm(1 - alpha[2])
  # theta_crit_1 = crit_1 * sqrt(sigma2_theta.2_1)
  # theta_crit_2 = crit_2 * sqrt(sigma2_theta.2_2)
  # theta_crit_3 = crit_3 * sqrt(sigma2_theta.2_3)

  # First calculate critical boundaries for first stage
  A1 = pnorm(crit_1 * sqrt(sigma2_theta.p_1_fixed),
             mean = Etheta.p_1_cond, sd = sqrt(Vtheta.p_1_cond), lower.tail = FALSE)
  A2 = pnorm(crit_2 * sqrt(sigma2_theta.p_2_fixed),
             mean = Etheta.p_2_cond, sd = sqrt(Vtheta.p_2_cond), lower.tail = FALSE)

  # Use predictive distribution of theta.2_1 given theta.1_1
  z_A1 = qnorm(A1, lower.tail = FALSE)
  z_A2 = qnorm(A2, lower.tail = FALSE)
  crit_A1 = z_A1 * sqrt(sigma2_theta.2_1)
  crit_A2 = z_A2 * sqrt(sigma2_theta.2_2)

  # Calculate P.values from first stage
  P.1_1 = pnorm(Z.1_1, lower.tail = FALSE)
  P.1_2 = pnorm(Z.1_2, lower.tail = FALSE)
  # Calculate Z values from 2nd stage and their respective p.values
  Z.2_1 = theta.2_1 / sqrt(sigma2_theta.2_1)
  Z.2_2 = theta.2_2 / sqrt(sigma2_theta.2_2)
  P.2_1 = P.2_2 = 1 # Initialize p-values in 1.
  # If we recruit patients from the arm, then update them
  if(lambda.2_1 > 0) {
    P.2_1 = pnorm(Z.2_1, lower.tail = FALSE)
  }
  if(lambda.2_2 > 0) {
    P.2_2 = pnorm(Z.2_2, lower.tail = FALSE)
  }

  if (verbose) {
    cat(paste("lambda.1: ", paste0(round(c(lambda.1_1, lambda.1_2), 4), collapse = " - "), "\n"))
    cat(paste("theta_hat.1: ", paste0(round(c(x.1), 4), collapse = " - "), "\n"))
    cat(paste("lambda.2: ", paste0(round(c(lambda.2_1, lambda.2_2), 4), collapse = " - "), "\n"))
    cat(paste("theta_hat.2: ", paste0(round(c(x.2), 4), collapse = " - "), "\n"))
    cat("1st stage error rates", "A1: ", A1, "- A2: ", A2, "\n")
  }

  # Pre-allocate values for rejection
  RejectH_01   = RejectH_02  = 0
  RejectH_01.p = RejectH_02.p  = 0
  # Check whether we have a two stage trial of one stage (r=1)
  isTwoStageTrial = (r != 0 & r!= 1)
  if (isTwoStageTrial == FALSE){
    if(r == 1){
      P.p_1 = P.1_1
      P.p_2 = P.1_2
      # If we have a one stage trials, simply test the hypotheses using 1st stage p-vaues
      if (P.1_1 <= alpha[1]){
        RejectH_01 = 1
      }
      if (P.1_2 <= alpha[2]){
        RejectH_02 = 1
      }
    } else {
      P.p_1 = P.2_1
      P.p_2 = P.2_2
      # If we have a one stage trials, simply test the hypotheses using 2nd stage p-vaues
      if (P.2_1 <= alpha[1]){
        RejectH_01 = 1
      }
      if (P.2_2 <= alpha[2]){
        RejectH_02 = 1
      }
    }
  } else {
    if (P.2_1 <= A1){
      RejectH_01 = 1
    }
    if (P.2_2 <= A2){
      RejectH_02 = 1
    }
  }
  if (P.p_1 <= alpha[1]){
    RejectH_01.p = 1
  }
  if (P.p_2 <= alpha[2]){
    RejectH_02.p = 1
  }

  if (verbose) {
    first  = c(theta.1_1, Z.1_1, P.1_1,
               theta.2_1, Z.2_1, P.2_1,
               theta.p_1, Z.p_1, P.p_1, RejectH_01, RejectH_01.p, A1)
    second = c(theta.1_2, Z.1_2, P.1_2,
               theta.2_2, Z.2_2, P.2_2,
               theta.p_2, Z.p_2, P.p_2, RejectH_02, RejectH_02.p, A2)
    what   = c("First Stage estimates",  "First Stage Z-values", "First Stage p-values",
               "Second Stage estimates", "Second Stage Z-values", "Second Stage p-values",
               "Pooled estimates",       "Pooled Z-values", "Pooled p-values",
               "Decision", "Decision from pooled p-values", "1st Stage Error Rates")
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|%10s|", "Trial summary", "Arm 1", "Arm 2")))
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|%10s|", "------------------------------", "----------", "----------", "----------")))
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|%10s|", what[1:3], round(first[1:3], 3), round(second[1:3], 3)), "\n"))
    cat("",   paste0(sprintf("%-30s|%10s|%10s|%10s|", "------------------------------", "----------", "----------", "----------")))
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|%10s|", what[12], round(first[12], 3), round(second[12], 3)), "\n"))
    cat("",   paste0(sprintf("%-30s|%10s|%10s|%10s|", "------------------------------", "----------", "----------", "----------")))
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|%10s|", what[4:6], round(first[4:6], 3), round(second[4:6], 3)), "\n"))
    cat("",   paste0(sprintf("%-30s|%10s|%10s|%10s|", "------------------------------", "----------", "----------", "----------")))
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|%10s|", what[7:9], round(first[7:9], 3), round(second[7:9], 3)), "\n"))
    cat("",   paste0(sprintf("%-30s|%10s|%10s|%10s|", "------------------------------", "----------", "----------", "----------")))
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|%10s|", what[10],  round(first[10], 3),  round(second[10], 3)), "\n"))
    cat("",   paste0(sprintf("%-30s|%10s|%10s|%10s|", "------------------------------", "----------", "----------", "----------")))
    cat("\n", paste0(sprintf("%-30s|%10s|%10s|%10s|", what[11],  round(first[11], 3),  round(second[11], 3)), "\n"), "\n")
  }

  # Calculate expected utility. Note this should not depend on the estimates.
  # or true values in the function.
  EU <- (gamma1(lambda = lambda_1) * RejectH_01 +
          gamma2(lambda = lambda_2) * RejectH_02)

  out = c(RejectH_01 = RejectH_01,
          RejectH_02 = RejectH_02,EU = EU,
          RejectH_01.p = RejectH_01.p,
          RejectH_02.p = RejectH_02.p,
          EU = EU,
          theta.p_1 = theta.p_1,
          theta.p_2 = theta.p_2,
          Z.1_1 = Z.1_1, Z.1_2 = Z.1_2,
          Z.2_1 = Z.2_1, Z.2_2 = Z.2_2,
          P.1_1 = P.1_1, P.1_2 = P.1_2,
          P.2_1 = P.2_1, P.2_2 = P.2_2,
          P.p_1 = P.p_1, P.p_2 = P.p_2,
          r = r)
  out
}



#' Helps find optimal values of first stage prevalences
#'
#' We calculate expected utility of the trial given that we will optimize at interim.
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
umbrella2armFullOptimization <- function(M = 1000, mc.cores = 40,
                                     r = 0.5,
                                     sigma_1 = 1, sigma_2 = 1,
                                     rho = 0,
                                     mu_1 = 0, psi_1 = 1,
                                     mu_2 = 0, psi_2 = 1,
                                     rho0 = 0,
                                     lambda   = c(0.15),
                                     lambda.1 = c(0.15),
                                     n = 1484,
                                     alpha = c(0.05, 0.05),
                                     gamma1  = ..umbrella_gamma1..,  #USES SIMPLE UTILITIES
                                     gamma2  = ..umbrella_gamma2..,  #USES SIMPLE UTILITIES
                                     verbose = TRUE,
                                     tol.zero = 1e-10,
                                     sig.digits = 4,
                                     summarize = TRUE){
  # Simulate 1st stage estimates
  lambda.1_1 = round(lambda.1[1], sig.digits)
  lambda.1_2 = 1 - lambda.1_1

  # Simulate 1st stage estimates
  Sigma.1     = diag(c(4/lambda.1_1/n/r,
                       4/lambda.1_2/n/r))

  psi2_1 = psi_1^2
  psi2_2 = psi_2^2
  Sigma0  = matrix(c(psi2_1,               rho0 * psi_1 * psi_2,
                     rho0 * psi_1 * psi_2,               psi2_2),
                   nrow = 2)

  theta.hat.1 = mvtnorm::rmvnorm(n = M,
                                 mean  = c(mu_1, mu_2),
                                 sigma = Sigma.1 + Sigma0)
  if (verbose) cat("Simulate 1st stage: DONE \n")

  # Given First stage estimates, calculate the values of lambda.2 and w.2
  #  with the maximum expected utility
  res = parallel::mclapply(X = 1:M, function(i) {
    interimDecision = interimUmbrella2armTrial(x = theta.hat.1[i, ],
                                           r          = r,
                                           lambda     = lambda,
                                           lambda.1   = lambda.1,
                                           mu_1       = mu_1,
                                           mu_2       = mu_2,
                                           psi_1      = psi_1,
                                           psi_2      = psi_2,
                                           alpha      = alpha,
                                           rho0       = rho0,
                                           n          = n,
                                           gamma1  = gamma1,
                                           gamma2  = gamma2,
                                           verbose    = verbose,
                                           sig.digits = sig.digits)
    if (verbose) cat("Interim decision: DONE \n")
    # Store values of lambda.2 and w.2
    c(lambda.1_1 = lambda.1,
      lambda.1_2 = 1-lambda.1,
      theta_hat.1_ = theta.hat.1[i, ],
      r = r,
      lambda.2_1 = round(interimDecision$par[1], 4),
      lambda.2_2 = 1-round(interimDecision$par[1], 4),
      EU         = interimDecision$value)
  }, mc.cores = mc.cores)
  res = do.call(rbind, res)
  out = res
  if(summarize){
    res2 = colMeans(res)
    out = res2
  }
  out
}
