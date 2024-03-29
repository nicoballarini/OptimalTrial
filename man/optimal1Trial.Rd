% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_twostage_withrho.R
\name{optimal1Trial}
\alias{optimal1Trial}
\title{Simulates a clinical trial with optimized sample allocation.}
\usage{
optimal1Trial(r.1 = 0.1, lambda = 0.5, lambda.1 = 0.5, w.1 = 0.5,
  theta_1 = 0, theta_2 = 0, mu_1 = 0, psi_1 = 1, mu_2 = 0,
  psi_2 = 1, rho0 = 0, alpha = 0.05, n = 100,
  gamma1 = ..gamma1.., gamma2 = ..gamma2.., gamma12 = ..gamma12..,
  tol.zero = 1e-05, tol.lim = 1e-04, verbose = TRUE,
  use.prior.as.true = FALSE, optim_func = "optim")
}
\arguments{
\item{r.1}{proportion of patients in 1st stage of the trial}

\item{lambda}{proportion of patients in sub-population 1 in the overal patient population}

\item{lambda.1}{proportion of patients in sub-population 1 in the trial 1st stage}

\item{w.1}{factor for allocation of alpha to H01 in the testing correction (1st stage)}

\item{theta_1}{true parameter value for sub-population 1}

\item{theta_2}{true parameter value for sub-population 2}

\item{mu_1}{prior mean for the treatment effect in sub-population 1}

\item{psi_1}{prior standard deviation for the treatment effect in sub-population 1}

\item{mu_2}{prior mean for the treatment effect in sub-population 2}

\item{psi_2}{prior standard deviation for the treatment effect in sub-population 2}

\item{alpha}{significance level}

\item{n}{trial sample size}

\item{gamma1}{gain when rejecting H01. A function with parameters (mu1, mu2, theta1, theta2, lambda)}

\item{gamma2}{gain when rejecting H02. A function with parameters (mu1, mu2, theta1, theta2, lambda)}

\item{gamma12}{gain when rejecting H01 and H02. A function with parameters (mu1, mu2, theta1, theta2, lambda)}

\item{tol.zero}{Absolute distance to 0 or 1 to be considered 0 or 1}

\item{verbose}{Print intermediate values for calculations}

\item{use.prior.as.true}{logical. If TRUE, then the true theta is generate using the prior distribution and theta_1 and theta_2 are not used.}

\item{rho_0}{prior correlation between treatment effect in sub-population 1 and sub-population 2 (ignored)}
}
\description{
Given true values for theta_1 and theta_2, the function generates first stage data,
then find optimal values lambda.2 and w.2, generates second stage data, and
finally test the hypotheses.
}
