#' Simulates a clinical trial with adaptive enrichment
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
#' @details The function performs an adaptive enrichment design.
#' At interim analysis, there are three posibilities:
#'
#' 1. We continue with the same trial prevalence and testing weight.
#' 2. We enrich for subgroup and therefore lambda.1 = 1 and w.1 = 1
#' 3. We enrich for complement and therefore lambda.1 = 0 and w.1 = 0
#'
#' @importFrom stats rnorm pnorm qnorm
#'
#' @export
AdaptiveEnrichmentTrial <- function(r.1      =   0.1,
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
                                     use.prior.as.true = FALSE){

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

  # Simulate 1st stage estimates ----------------------------------------------
  Sigma.1     = diag(c(1/lambda.1/n/r.1, 1/(1-lambda.1)/n/r.1))
  theta.hat.1 = mvtnorm::rmvnorm(n = 1,
                                 mean = true.mean,
                                 sigma = Sigma.1)
  # cat(theta.hat.1)
  if (verbose) cat("Simulate 1st stage: DONE \n")

  # Given First stage estimates, calculate the values of lambda.2 and w.2
  #  with the maximum expected utility ----------------------------------------------
  interim3 = interimAdaptiveEnrichment(lambda.2 = lambda.1, w.2 = w.1,
                                       x = theta.hat.1,
                                       r = r.1,
                                       sigma_1 = 1, sigma_2 = 1,
                                       rho  = 0,
                                       mu_1 = mu_1, psi_1 = psi_1,
                                       mu_2 = mu_2, psi_2 = psi_1,
                                       rho0 = rho0,
                                       lambda = lambda, lambda.1 = lambda.1,
                                       n   = n,
                                       w.1 = w.1, alpha = alpha,
                                       adj = "bonferroni-holm",
                                       gamma1 = gamma1,
                                       gamma2 = gamma2,
                                       gamma12 = gamma12,
                                       verbose = FALSE)
  interim2 = interimAdaptiveEnrichment(lambda.2 = 1, w.2 = 1,
                                       x = theta.hat.1,
                                       r = r.1,
                                       sigma_1 = 1, sigma_2 = 1,
                                       rho  = 0,
                                       mu_1 = mu_1, psi_1 = psi_1,
                                       mu_2 = mu_2, psi_2 = psi_1,
                                       rho0 = rho0,
                                       lambda = lambda, lambda.1 = lambda.1,
                                       n   = n,
                                       w.1 = w.1, alpha = alpha,
                                       adj = "bonferroni-holm",
                                       gamma1 = gamma1,
                                       gamma2 = gamma2,
                                       gamma12 = gamma12,
                                       verbose = FALSE)
  interim1 = interimAdaptiveEnrichment(lambda.2 = 0, w.2 = 0,
                                       x = theta.hat.1,
                                       r = r.1,
                                       sigma_1 = 1, sigma_2 = 1,
                                       rho  = 0,
                                       mu_1 = mu_1, psi_1 = psi_1,
                                       mu_2 = mu_2, psi_2 = psi_1,
                                       rho0 = rho0,
                                       lambda = lambda, lambda.1 = lambda.1,
                                       n   = n,
                                       w.1 = w.1, alpha = alpha,
                                       adj = "bonferroni-holm",
                                       gamma1 = gamma1,
                                       gamma2 = gamma2,
                                       gamma12 = gamma12,
                                       verbose = FALSE)

  EUs = c(interim1$EU, interim2$EU, interim3$EU)
  interimDecision = which(EUs == max(EUs))
  if (verbose) cat("Interim decision: DONE \n")

  if (interimDecision == 1){
    lambda.2 = 0
    w.2      = 0
  } else if (interimDecision == 2){
    lambda.2 = 1
    w.2      = 1
  } else {
    lambda.2 = lambda.1
    w.2      = w.1
  }
  # Store values of lambda.2 and w.2

  # When lambda.2 and w.2 are very close to 0 or 1, we just set them to 0 or 1 respectively.
  if (lambda.2 < tol.zero)   lambda.2 = 0
  if (lambda.2 > 1-tol.zero) lambda.2 = 1
  if (w.2 < tol.zero)        w.2 = 0
  if (w.2 > 1-tol.zero)      w.2 = 1

  # Simulates second stage estimates given the optimized lambda.2 --------------
  # Calculates first the variances
  Sigma.2_1 = 1/(lambda.2)/n/(1-r.1)
  Sigma.2_2 = 1/(1-lambda.2)/n/(1-r.1)
  # And simulate the estimates. A NA pre-allocated vector is used so that
  # We only simulate the value if it actually is needed.
  theta.hat.2 = c(NA, NA)
  if (lambda.2 != 0) theta.hat.2[1] = rnorm(n = 1, mean  = true.mean[1], sd = sqrt(Sigma.2_1))
  if (lambda.2 != 1) theta.hat.2[2] = rnorm(n = 1, mean  = true.mean[2], sd = sqrt(Sigma.2_2))
  if (verbose) cat("Simulate 2nd stage: DONE \n")

  # Test the hypotheses given first and second stage data. ---------------------
  decision = decisionAdaptiveEnrichment(lambda.2 = lambda.2, w.2 = w.2,
                                        x.1 = theta.hat.1,
                                        x.2 = theta.hat.2,
                                        r = r.1,
                                        lambda = lambda, lambda.1 = lambda.1,
                                        n = n,
                                        w.1 = w.1,
                                        alpha = alpha,
                                        adj = "bonferroni-holm",
                                        gamma1 = gamma1,
                                        gamma2 = gamma2,
                                        gamma12 = gamma12,
                                        verbose = verbose)
  if (verbose) cat("Decision: DONE \n")
  out = c(decision,
          lambda = lambda,
          lambda.1 = lambda.1, lambda.2 = lambda.2,
          w.1 = w.1, w.2 = w.2,
          theta.hat.1_ = theta.hat.1, theta.hat.2_ = theta.hat.2)
  class(out) = "simOptimalTrial"
  out
}



#' Expected utility for values of lambda.2 and w.2 after 1st stage
#'
#' Simplified version where the utility function does not depend on the
#' estimates or the true values of the parameters and can be taken out of
#' the first integration
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
#' @param adj      adjustment for multiple testing. Either "bonferroni", "bonferroni-holm", or "simes"
#' @param verbose  Print intermediate values for calculations
#' @param tol.zero Absolute distance to 0 or 1 to be considered 0 or 1
#'
#' @export
interimAdaptiveEnrichment <- function(lambda.2, w.2,
                                       x, r,
                                       sigma_1 = 1, sigma_2 = 1, rho = 0,
                                       mu_1 = 0, psi_1 = 1,
                                       mu_2 = 0, psi_2 = 1, rho0 = 0,
                                       lambda = 0.5, lambda.1 = 0.5,
                                       n = 100,
                                       w.1 = 0.5,
                                       alpha = 0.025,
                                       adj = c("bonferroni", "bonf", "bonferroni-holm", "simes"),
                                       gamma1  = ..gamma1..,  #USES SIMPLE UTILITIES
                                       gamma2  = ..gamma2..,  #USES SIMPLE UTILITIES
                                       gamma12 = ..gamma12.., #USES SIMPLE UTILITIES
                                       verbose = TRUE, show.case = TRUE, tol.zero = 0.001){


  adj = match.arg(adj)
  # x is a 2D vector with the estimates from 1st stage.
  theta.1_1 <- (x[1])
  theta.1_2 <- (x[2])
  mu = c(mu_1, mu_2)

  # When lambda.2 and w.2 are very close to 0 or 1, we just set them to 0 or 1 respectively.
  if (lambda.2 < tol.zero) lambda.2 = 0
  if (w.2 < tol.zero)      w.2 = 0
  if (lambda.2 > 1-tol.zero) lambda.2 = 1
  if (w.2 > 1-tol.zero)      w.2 = 1

  # Calculate trial prevalences
  lambda.1_1 = lambda.1
  lambda.1_2 = (1-lambda.1)
  lambda.2_1 = lambda.2
  lambda.2_2 = (1-lambda.2)

  # Calculate trial weights
  w.1_1 = w.1
  w.1_2 = (1-w.1)
  w.2_1 = w.2
  w.2_2 = (1-w.2)


  # Subgroup: Conditional expectation given 1st stage
  if (lambda.1_1 != 0){ # Case when subgroup 1 trial prevalence is >0 in first stage
    psi2_1 = psi_1^2
    sigma2_1 = sigma_1^2
    psi2.1_1 = sigma2_1 * psi2_1 / (lambda.1_1 * r * n *psi2_1 + sigma2_1)
    sigma2_theta.1_1 = sigma2_1 / (lambda.1_1 * r * n)
    sigma2_theta.2_1_fixed = sigma2_1 / (lambda.1_1 * (1-r) * n)
    sigma2_theta.2_1 = sigma2_1 / (lambda.2_1 * (1-r) * n)
    sigma2_theta.p_1_fixed = sigma2_1 / n / (lambda.1_1 * r + lambda.1_1 * (1 - r))
    sigma2_theta.p_1 = sigma2_1 / n / (lambda.1_1 * r + lambda.2_1 * (1 - r))
    mu.1_1 = psi2.1_1 * (mu_1 / psi2_1 + theta.1_1 / sigma2_theta.1_1)
    Z.1_1 = theta.1_1/(sigma_1/sqrt(n * lambda.1_1 * r))
    f.1_1 = lambda.1_1 * r / (lambda.1_1 * r  + lambda.2_1 * (1-r))
    f.2_1 = 1 - f.1_1
  } else { # If subgroup 1 trial prevalence is 0 in first stage, we only have 2nd stage data
    psi2_1 = psi_1^2
    sigma2_1 = sigma_1^2
    psi2.1_1 = psi2_1
    sigma2_theta.1_1 = 0
    sigma2_theta.2_1 = sigma2_1 / (lambda.2_1 * (1-r) * n)
    sigma2_theta.p_1 = sigma2_1 / (lambda.2_1 * (1-r) * n)
    mu.1_1 = mu_1
    Z.1_1 = NA
    f.1_1 = 0
    f.2_1 = 1
  }
  if (lambda.2_1 == 0){ # If subgroup 1 trial prevalence is 0 in second stage, we only have 1st stage data
    sigma2_theta.2_1 = 0
  }

  Etheta.p_1_cond = r * theta.1_1
  Vtheta.p_1_cond = (1-r)^2 * (sigma2_theta.2_1_fixed)

  # Complement: Conditional expectation given 2nd stage
  if (lambda.1_2 != 0){# Case when subgroup 2 trial prevalence is >0 in first stage
    sigma2_2 = sigma_2^2
    psi2_2 = psi_2^2
    psi2.1_2 = sigma2_2 * psi2_2 / (lambda.1_2 * r * n *psi2_2 + sigma2_2)
    sigma2_theta.1_2 = sigma2_2 / (lambda.1_2 * r * n)
    sigma2_theta.2_2_fixed = sigma2_2 / (lambda.1_2 * (1-r) * n)
    sigma2_theta.2_2 = sigma2_2 / (lambda.2_2 * (1-r) * n)
    sigma2_theta.p_2_fixed = sigma_2^2 / n / (lambda.1_2 * r + lambda.1_2 * (1 - r))
    sigma2_theta.p_2 = sigma_2^2 / n / (lambda.1_2 * r + lambda.2_2 * (1 - r))
    mu.1_2 = psi2.1_2 * (mu_2 / psi2_2 + theta.1_2 / sigma2_theta.1_2)
    Z.1_2 = theta.1_2/(sigma_2/sqrt(n * lambda.1_2 * r))
    f.1_2 = lambda.1_2 * r / (lambda.1_2 * r  + lambda.2_2 * (1-r))
    f.2_2 = 1 - f.1_2
  } else { # If subgroup 2 trial prevalence is 0 in first stage, we only have 2nd stage data
    psi2_2 = psi_2^2
    sigma2_2 = sigma_2^2
    psi2.1_2 = psi2_2
    sigma2_theta.1_2 = 0
    sigma2_theta.2_2 = sigma2_2 / (lambda.2_2 * (1-r) * n)
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
  Vtheta.p_2_cond = (1-r)^2 * (sigma2_theta.2_2_fixed)

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
  Sigma = Sigma.1 + Sigma.2_theta

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
  A1 = pnorm(crit_w.1_1 * sqrt(sigma2_theta.p_1_fixed),
             mean = Etheta.p_1_cond, sd = sqrt(Vtheta.p_1_cond), lower.tail = FALSE)
  A2 = pnorm(crit_w.1_2 * sqrt(sigma2_theta.p_2_fixed),
             mean = Etheta.p_2_cond, sd = sqrt(Vtheta.p_2_cond), lower.tail = FALSE)
  A12 = A1 + A2 - A1*A2

  B1 = pnorm(crit * sqrt(sigma2_theta.p_1_fixed),
             mean = Etheta.p_1_cond, sd = sqrt(Vtheta.p_1_cond), lower.tail = FALSE)
  B2 = pnorm(crit * sqrt(sigma2_theta.p_2_fixed),
             mean = Etheta.p_2_cond, sd = sqrt(Vtheta.p_2_cond), lower.tail = FALSE)


  # Use predictive distribution of theta.2_1 given theta.1_1
  z_A1 = qnorm(A1, lower.tail = FALSE)
  z_A2 = qnorm(A2, lower.tail = FALSE)
  z_A12 = qnorm(A12, lower.tail = FALSE)
  z_B1 = qnorm(B1, lower.tail = FALSE)
  z_B2 = qnorm(B2, lower.tail = FALSE)

  crit_A1 = z_A1 * sqrt(sigma2_theta.2_1)
  crit_A2 = z_A2 * sqrt(sigma2_theta.2_2)
  crit_B1 = z_B1 * sqrt(sigma2_theta.2_1)
  crit_B2 = z_B2 * sqrt(sigma2_theta.2_2)

  # Calculate the probabilities of rejection using the prior distribution
  case = 5
  # Check whether we have a two stage trial or one stage (r=1)
  isTwoStageTrial = (r != 1)
  # Case 1: Enrich with complement Lambda.2_1 = 0 ------------------------------
  if (isTwoStageTrial & lambda.2_1 == 0 & w.2_1 == 0){# Case w.2 = 0
    crit_A12 = z_A12 * sqrt(sigma2_theta.2_2)
    pr2 <-  pnorm(q = max(crit_B2, crit_A12),
                  mean = mu.1_2, sd = sqrt(Sigma[2,2]), lower.tail = FALSE)
    pr1 <-  0
    pr12 <- pr1 * pr2
    # cat("Probability of rejection H01",pr1,"\n")
    # cat("Probability of rejection H02",pr2,"\n")
    # cat("Probability of rejection H01^H02",pr12,"\n")
    # cat("Power = ",(pr1 + pr2 - pr12), "\n")
  }
  # Case 2: Enrich the subgroup Lambda.2_1 = 1 ---------------------------------
  if (isTwoStageTrial & lambda.2_1 == 1 & w.2_1 == 1){# Case w.2 = 1
    crit_A12 = z_A12 * sqrt(sigma2_theta.2_1)
    pr1 <-  pnorm(q = max(crit_B1, crit_A12),
                  mean = mu.1_1, sd = sqrt(Sigma[1,1]), lower.tail = FALSE)
    pr2 <-  0
    pr12 <- pr1 * pr2
    # cat("Probability of rejection H01",pr1,"\n")
    # cat("Probability of rejection H02",pr2,"\n")
    # cat("Probability of rejection H01^H02",pr12,"\n")
    # cat("Power = ",(pr1 + pr2 - pr12), "\n")
  }
  # Case 3: Continue with the same as in first stage 0 < lambda.2_1 < 1 --------
  if (isTwoStageTrial & (0 < lambda.2_1 & lambda.2_1 < 1)){
    # Now calculate probability of rejection
    q1 <- mvtnorm::pmvnorm(lower = c(-Inf, max(crit_A2, crit_B2)),
                           upper = c(min(crit_B1, crit_A1),  Inf),
                           mean  = c(mu.1_1, mu.1_2),
                           sigma = Sigma,
                           algorithm = mvtnorm::GenzBretz())[1]
    q2 <- mvtnorm::pmvnorm(lower = c(min(crit_B1, crit_A1), max(crit_A2, crit_B2)),
                           upper = c(max(crit_B1, crit_A1), Inf),
                           mean  = c(mu.1_1, mu.1_2),
                           sigma = Sigma,
                           algorithm = mvtnorm::GenzBretz())[1]
    q3 <- mvtnorm::pmvnorm(lower = c(max(crit_B1, crit_A1), max(crit_A2, crit_B2)),
                           upper = c(Inf,   Inf),
                           mean  = c(mu.1_1, mu.1_2),
                           sigma = Sigma,
                           algorithm = mvtnorm::GenzBretz())[1]
    q4 <- mvtnorm::pmvnorm(lower = c(max(crit_B1, crit_A1), min(crit_A2, crit_B2)),
                           upper = c(Inf,     max(crit_A2, crit_B2)),
                           mean  = c(mu.1_1, mu.1_2),
                           sigma = Sigma,
                           algorithm = mvtnorm::GenzBretz())[1]
    q5 <- mvtnorm::pmvnorm(lower = c(max(crit_B1, crit_A1),   -Inf),
                           upper = c(Inf,    min(crit_A2, crit_B2)),
                           mean  = c(mu.1_1, mu.1_2),
                           sigma = Sigma,
                           algorithm = mvtnorm::GenzBretz())[1]
    q6 <- mvtnorm::pmvnorm(lower = c(min(crit_B1, crit_A1), min(crit_A2, crit_B2)),
                           upper = c(max(crit_B1, crit_A1), max(crit_A2, crit_B2)),
                           mean  = c(mu.1_1, mu.1_2),
                           sigma = Sigma,
                           algorithm = mvtnorm::GenzBretz())[1]

    if (crit_B1 <= crit_A1 & crit_B2 <= crit_A2){
      pr1  <- q2 + q3 + q4 + q5
      pr2  <- q1 + q2 + q3 + q4
      pr12 <- q2 + q3 + q4
      case <- 1
    } else {
      stop("Check code")
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
  EU <- (gamma1(lambda = lambda) * pr1 +
           gamma2(lambda = lambda) * pr2 +
           gamma12(lambda = lambda) * pr12)
  if(show.case) {
    out = list(EU = EU, case = case)
  } else {
    out = EU
  }
  out
}

#' Expected utility for values of lambda.2 and w.2 after 1st stage
#'
#' Simplified version where the utility function does not depend on the
#' estimates or the true values of the parameters and can be taken out of
#' the first integration
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
#' @param adj      adjustment for multiple testing. Either "bonferroni", "bonferroni-holm", or "simes"
#' @param verbose  Print intermediate values for calculations
#' @param tol.zero Absolute distance to 0 or 1 to be considered 0 or 1
#'
#' @export
interimAdaptiveEnrichment_val <- function(x, r,
                                      sigma_1 = 1, sigma_2 = 1, rho = 0,
                                      mu_1 = 0, psi_1 = 1,
                                      mu_2 = 0, psi_2 = 1, rho0 = 0,
                                      lambda = 0.5, lambda.1 = 0.5,
                                      n = 100,
                                      w.1 = 0.5,
                                      alpha = 0.025,
                                      adj = c("bonferroni", "bonf", "bonferroni-holm", "simes"),
                                      gamma1  = ..gamma1..,  #USES SIMPLE UTILITIES
                                      gamma2  = ..gamma2..,  #USES SIMPLE UTILITIES
                                      gamma12 = ..gamma12.., #USES SIMPLE UTILITIES
                                      verbose = TRUE, show.case = TRUE, tol.zero = 0.001){


  # Given First stage estimates, calculate the values of lambda.2 and w.2
  #  with the maximum expected utility ----------------------------------------------
  interim3 = interimAdaptiveEnrichment(lambda.2 = lambda.1, w.2 = w.1,
                                       x = x,
                                       r = r,
                                       sigma_1 = 1, sigma_2 = 1,
                                       rho  = 0,
                                       mu_1 = mu_1, psi_1 = psi_1,
                                       mu_2 = mu_2, psi_2 = psi_1,
                                       rho0 = rho0,
                                       lambda = lambda, lambda.1 = lambda.1,
                                       n   = n,
                                       w.1 = w.1, alpha = alpha,
                                       adj = "bonferroni-holm",
                                       gamma1 = gamma1,
                                       gamma2 = gamma2,
                                       gamma12 = gamma12,
                                       verbose = TRUE)
  interim2 = interimAdaptiveEnrichment(lambda.2 = 1, w.2 = 1,
                                       x = x,
                                       r = r,
                                       sigma_1 = 1, sigma_2 = 1,
                                       rho  = 0,
                                       mu_1 = mu_1, psi_1 = psi_1,
                                       mu_2 = mu_2, psi_2 = psi_1,
                                       rho0 = rho0,
                                       lambda = lambda, lambda.1 = lambda.1,
                                       n   = n,
                                       w.1 = w.1, alpha = alpha,
                                       adj = "bonferroni-holm",
                                       gamma1 = gamma1,
                                       gamma2 = gamma2,
                                       gamma12 = gamma12,
                                       verbose = TRUE)
  interim1 = interimAdaptiveEnrichment(lambda.2 = 0, w.2 = 0,
                                       x = x,
                                       r = r,
                                       sigma_1 = 1, sigma_2 = 1,
                                       rho  = 0,
                                       mu_1 = mu_1, psi_1 = psi_1,
                                       mu_2 = mu_2, psi_2 = psi_1,
                                       rho0 = rho0,
                                       lambda = lambda, lambda.1 = lambda.1,
                                       n   = n,
                                       w.1 = w.1, alpha = alpha,
                                       adj = "bonferroni-holm",
                                       gamma1 = gamma1,
                                       gamma2 = gamma2,
                                       gamma12 = gamma12,
                                       verbose = FALSE)

  EUs = c(Compl = interim1$EU, Subgr = interim2$EU, Both = interim3$EU)
  interimDecision = which(EUs == max(EUs))
  if (verbose) cat("Interim decision: DONE \n")

  if (interimDecision == 1){
    lambda.2 = 0
    w.2      = 0
  } else if (interimDecision == 2){
    lambda.2 = 1
    w.2      = 1
  } else {
    lambda.2 = lambda.1
    w.2      = w.1
  }
  # Store values of lambda.2 and w.2
  list(par = c(lambda.2, w.2),
       value =  max(EUs),
       utilities = EUs)
}

#' Performs hypothesis test for a clinical trial with adaptive enrichment
#'
#' Given values for first and second stage estimates, the function generates
#' test the hypotheses of no effect in the subgroup
#'
#' @param r.1      proportion of patients in 1st stage of the trial
#' @param x.1      1st stage data
#' @param x.2      2nd stage data
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
#'
#' @importFrom stats rnorm pnorm qnorm
#'
#' @export
decisionAdaptiveEnrichment <- function(lambda.2, w.2,
                              x.1,
                              x.2,
                              r,
                              sigma_1 = 1, sigma_2 = 1, rho = 0,
                              lambda = 0.5, lambda.1 = 0.5,
                              n = 100,
                              w.1 = 0.5,
                              alpha = 0.025,
                              adj = c("bonferroni", "bonf", "bonferroni-holm", "simes"),
                              gamma1  = ..gamma1..,  #USES SIMPLE UTILITIES
                              gamma2  = ..gamma2..,  #USES SIMPLE UTILITIES
                              gamma12 = ..gamma12.., #USES SIMPLE UTILITIES
                              fixed.weights = FALSE,
                              verbose = TRUE,
                              show.case = TRUE){
  adj = match.arg(adj)
  # x.1 and x.2 are 2D vectors with the estimates from 1st stage and 2nd stage resp.
  theta.1_1 <- (x.1[1])
  theta.1_2 <- (x.1[2])
  theta.2_1 <- (x.2[1])
  theta.2_2 <- (x.2[2])

  # Calculate trial prevalences
  lambda.1_1 = lambda.1
  lambda.1_2 = (1 - lambda.1)
  lambda.2_1 = lambda.2
  lambda.2_2 = (1 - lambda.2)

  # Calculate trial weights
  w.1_1 = w.1
  w.1_2 = (1 - w.1)
  w.2_1 = w.2
  w.2_2 = (1 - w.2)

  # Subgroup: Conditional expectation given 1st stage
  if (lambda.1_1 != 0){ # Case when subgroup 1 trial prevalence is >0 in first stage
    sigma2_1 = sigma_1^2
    sigma2_theta.1_1 = sigma2_1 / (lambda.1_1 * r * n)
    sigma2_theta.2_1_fixed = sigma2_1 / (lambda.1_1 * (1-r) * n)
    sigma2_theta.2_1 = sigma2_1 / (lambda.2_1 * (1-r) * n)
    sigma2_theta.p_1_fixed = sigma2_1 / n / (lambda.1_1 * r + lambda.1_1 * (1 - r))
    sigma2_theta.p_1 = sigma2_1 / n / (lambda.1_1 * r + lambda.2_1 * (1 - r))
    Z.1_1 = sqrt(n * lambda.1_1 * r) * theta.1_1 / (sigma_1)
    f.1_1 = lambda.1_1 * r / (lambda.1_1 * r  + lambda.2_1 * (1-r))
    f.2_1 = 1 - f.1_1
  } else {   # If subgroup 1 trial prevalence is 0 in first stage, we only have 2nd stage data
    sigma2_1 = sigma_1^2
    sigma2_theta.1_1 = 0
    sigma2_theta.2_1 = sigma2_1 / (lambda.2_1 * (1-r) * n)
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
  Vtheta.p_1_cond = (1-r)^2 * (sigma2_theta.2_1_fixed)
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
    sigma2_theta.1_2 = sigma2_2 / (lambda.1_2 * r * n)
    sigma2_theta.2_2_fixed = sigma2_2 / (lambda.1_2 * (1-r) * n)
    sigma2_theta.2_2 = sigma2_2 / (lambda.2_2 * (1-r) * n)
    sigma2_theta.p_2_fixed = sigma2_2 / n / (lambda.1_2 * r + lambda.1_2 * (1 - r))
    sigma2_theta.p_2 = sigma2_2 / n / (lambda.1_2 * r + lambda.2_2 * (1 - r))
    Z.1_2 = sqrt(n * lambda.1_2 * r) * theta.1_2/ (sigma_2)
    f.1_2 = lambda.1_2 * r / (lambda.1_2 * r  + lambda.2_2 * (1-r))
    f.2_2 = 1 - f.1_2
  } else {# If subgroup 2 trial prevalence is 0 in first stage, we only have 2nd stage data
    sigma2_2 = sigma_2^2
    sigma2_theta.1_2 = Inf
    sigma2_theta.2_2 = sigma2_2 / (lambda.2_2 * (1-r) * n)
    sigma2_theta.p_2 = sigma2_theta.2_2
    Z.1_2 = NA
    f.1_2 = 0
    f.2_2 = 1
  }
  if (lambda.2_2 == 0){# If subgroup 2 trial prevalence is 0 in second stage, we only have 1st stage data
    sigma2_theta.2_2 = 0
  }

  Etheta.p_2_cond = r * theta.1_2
  Vtheta.p_2_cond = (1-r)^2 * (sigma2_theta.2_2_fixed)
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
  crit_w.1_2 = qnorm(1 - alpha * (1-w.1_1))
  crit       = qnorm(1 - alpha)

  # Based on estimates from first stage we derive the distributions of
  # the conditional statistics
  # Compute conditional error rates
  A1 = pnorm(crit_w.1_1 * sqrt(sigma2_theta.p_1_fixed),
             mean = Etheta.p_1_cond, sd = sqrt(Vtheta.p_1_cond), lower.tail = FALSE)
  A2 = pnorm(crit_w.1_2 * sqrt(sigma2_theta.p_2_fixed),
             mean = Etheta.p_2_cond, sd = sqrt(Vtheta.p_2_cond), lower.tail = FALSE)
  A12 = A1 + A2 - A1*A2
  B1 = pnorm(crit * sqrt(sigma2_theta.p_1_fixed),
             mean = Etheta.p_1_cond, sd = sqrt(Vtheta.p_1_cond), lower.tail = FALSE)
  B2 = pnorm(crit * sqrt(sigma2_theta.p_2_fixed),
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
    cat("1st stage error rates", A1, A2, B1, B2, "\n")
    cat("Second Stage P values", P.2_1, P.2_2, "\n")
  }

  # Pre-allocate values for rejection
  RejectH_01 = RejectH_02 = RejectH_0102 = 0
  # Check whether we have a two stage trial of one stage (r=1)
  isTwoStageTrial = (r != 1)
  if (isTwoStageTrial == FALSE){
    stop("Not implemented")
  }


  # Case 1: Enrich with complement Lambda.2_1 = 0 ------------------------------
  if (isTwoStageTrial & lambda.2_1 == 0 & w.2_1 == 0){# Case w.2 = 0
    if (P.2_2 <= B2 & P.2_2 <= A12){
      RejectH_02 = 1
      RejectH_0102 = 1
    }
  }
  # Case 2: Enrich the subgroup Lambda.2_1 = 1 ---------------------------------
  if (isTwoStageTrial & lambda.2_1 == 1 & w.2_1 == 1){# Case w.2 = 1
    if (P.2_1 < B1 & P.2_1 < A12){
      RejectH_01 = 1
      RejectH_0102 = 1
    }
  }
  # Case 3: Continue with the same as in first stage 0 < lambda.2_1 < 1 --------
  if (isTwoStageTrial & (0 < lambda.2_1 & lambda.2_1 < 1)){
    if (P.2_1 <= B1 & (P.2_1 <= A1 | P.2_2 <= A2)){
      RejectH_01 = 1
    }
    if (P.2_2 <= B2 & (P.2_1 <= A1 | P.2_2 <= A2)){
      RejectH_02 = 1
    }
    if (P.2_1 <= A1 | P.2_2 <= A2){
      RejectH_0102 = 1
    }
  }



  # Use predictive distribution of theta.2_1 given theta.1_1
  z_A1 = qnorm(A1, lower.tail = FALSE)
  z_A2 = qnorm(A2, lower.tail = FALSE)
  z_B1 = qnorm(B1, lower.tail = FALSE)
  z_B2 = qnorm(B2, lower.tail = FALSE)

  crit_A1 = z_A1 * sqrt(sigma2_theta.2_1)
  crit_A2 = z_A2 * sqrt(sigma2_theta.2_2)
  crit_B1 = z_B1 * sqrt(sigma2_theta.2_1)
  crit_B2 = z_B2 * sqrt(sigma2_theta.2_2)


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

#' Simulates M clinical trial with adaptive enrichment
#'
#' Uses the AdaptiveEnrichmentTrial to simulate nsim trials.
#' Just a wrapper for using mclapply.
#'
#'
#' @export
AdaptiveEnrichmentMtrials = function(nsim   =  1000,
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
         use.prior.as.true = FALSE){
  results = parallel::mclapply(1:nsim, FUN = function(x){
    AdaptiveEnrichmentTrial(r.1      =   r.1,
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
                            use.prior.as.true = use.prior.as.true)
  }, mc.cores = mc.cores)
  res = data.frame(do.call(rbind,results))
  out = res
  if(summarize){
    res2 = colMeans(res)
    out = res2
  }
  results.df = list(simresults = out,
                    summarize = summarize)
  class(results.df) = "simOptimalTrial"
  results.df
}
