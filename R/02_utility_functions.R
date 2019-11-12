# Default utility functions to use in other functions
# The utility may depend on the estimates from the trial X1 and X2
# The true parameters theta1, theta2,
# The pervalence of the subgroup lambda
# And the Market size N
#' @export
..gamma1..  <- function(X_1, X_2, theta_1, theta_2, lambda, N) {lambda}
#' @export
..gamma2..  <- function(X_1, X_2, theta_1, theta_2, lambda, N) {1-lambda}
#' @export
..gamma12.. <- function(X_1, X_2, theta_1, theta_2, lambda, N) {0}


#' @export
..umbrella_gamma1..  <- function(X_1, theta_1, lambda, N) {lambda}
#' @export
..umbrella_gamma2..  <- function(X_2, theta_2, lambda, N) {lambda}
#' @export
..umbrella_gamma3..  <- function(X_3, theta_3, lambda, N) {lambda}
