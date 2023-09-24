#' Calculate the Gaussian kernel.
#'
#' Calculate the Gaussian kernel.
#' @param a A vector containing the train/test points.
#' @param b A vector containing the train/test points.
#' @param sigma_sq Signal variance parameter.
#' @param lambda Length scale parameter.
#' @return A matrix containing the kernel between each element of \code{a} and \code{b}.
#' @export
kern_gaussian <- function(a, b, sigma_sq, lambda) {

  k_ab <- matrix(rep(0, length(a) * length(b)), nrow = length(a))
  for (i in 1:length(a)) {
    for (j in 1:length(b)) {
      k_ab[i, j] <- sigma_sq * exp(-0.5 * t((a[i] - b[j])) * (a[i] - b[j]) / (lambda^2)) # guassian kernel
    }
  }
  return(k_ab)
}

#' Constant prior mean structure of the Gaussian process
#'
#' @param input_points A vector containing the input points.
#' @return A matrix containing the constant prior mean structure used for the Gaussian process regression.
#' @export
const_prior <- function(input_points){
  h <- rep(1, length(input_points))
  h_t <- t(t(h))
  return(h_t)
}

#' Linear prior mean structure of the Gaussian process
#'
#' @param input_points A vector containing the input points.
#' @return A matrix containing the linear prior mean structure used for the Gaussian process regression.
#' @export
linear_prior <- function(input_points){
  h <- rbind(1, input_points)
  rownames(h) <- NULL
  h_t <- t(h)
  return(h_t)
}

#' Quadratic prior mean structure of the Gaussian process
#'
#' @param input_points A vector containing the input points.
#' @return A matrix containing the quadratic prior mean structure used for the Gaussian process regression.
#' @export
quad_prior <- function(input_points){
  h <- rbind(1, input_points, input_points^2)
  rownames(h) <- NULL
  h_t <- t(h)
  return(h_t)
}

#' Compute the single level Gaussian process posterior.
#'
#' Compute the mean vector and the covariance matrix of the two level Gaussian process posterior.
#' @param train_points A vector containing the training points.
#' @param test_points A vector containing the test points.
#' @param sigma_sq Signal variance parameter of the Gaussian process.
#' @param lambda Length scale parameter of the Gaussian process.
#' @param l_train A vector containing the function approximation at the training points.
#' @param p The nugget.
#' @param prior_mean_str A function giving the structure of the prior mean.
#' @return A list \code{posterior} containing the posterior mean vector \code{mean} and
#' the posterior covariance matrix \code{covariance}.
#' @export
gp_posterior <- function(train_points, test_points, sigma_sq, lambda, l_train, p, prior_mean_str = NULL) {

  ker_trtr <- kern_gaussian(train_points, train_points, sigma_sq, lambda)
  ker_trte <- kern_gaussian(train_points, test_points, sigma_sq, lambda)
  ker_tetr <- kern_gaussian(test_points, train_points, sigma_sq, lambda)
  ker_tete <- kern_gaussian(test_points, test_points, sigma_sq, lambda)

  posterior <- list()

  k_trtr_inv <- solve(ker_trtr + p * diag(1, nrow = nrow(ker_trtr)))

  H_train <- prior_mean_str(train_points)
  H_test <- prior_mean_str(test_points)

  R_T <- H_test - t(ker_trte) %*% k_trtr_inv %*% H_train

  common_inv <- solve(t(H_train) %*% k_trtr_inv %*% H_train)

  beta_bar <- common_inv %*% t(H_train) %*% k_trtr_inv %*% (l_train)

  posterior[["mean"]] <- ker_tetr %*% k_trtr_inv %*% l_train + R_T %*% beta_bar
  posterior[["covariance"]] <- ker_tete - ker_tetr %*% k_trtr_inv %*% ker_trte + R_T %*% common_inv %*% t(R_T)
  return(posterior)
}

# single level hyper estimation
cov_V_one <- function(train_points, hypers) {

  sigma_sq_1 <- hypers[1]; lambda_1 <- hypers[2]
  kernel_A1_d1 <- kern_gaussian(train_points, train_points, 1, lambda_1)
  V <- sigma_sq_1 * kernel_A1_d1
  return(V)
}

#' Compute the estimated intercept parameters
#'
#' Compute the estimated intercept parameters
#' of the low-level and high-level.
#' @param train_points A vector containing the training points.
#' @param l_train A vector containing the function approximation at the training points.
#' @param hypers A vector containing the hyperparameters of the Gaussian porcess: signal variance and length scale parameters.
#' @param p The nugget.
#' @param prior_mean_str A function giving the structure of the prior mean.
#' @return A list \code{beta_estimates} containing the estimated intercept parameters.
beta_estim_one <- function(train_points, l_train, hypers, p, prior_mean_str){

  z <- matrix(c(l_train), nrow = (length(l_train)), ncol = 1)

  H <- cbind(prior_mean_str(train_points))

  V_z <- cov_V_one(train_points, hypers)
  V_inv <- solve(V_z + p * diag(1, nrow = nrow(V_z)))

  beta_estimates <- solve(t(H) %*% V_inv %*% H) %*% t(H) %*% V_inv %*% (z)
  return(beta_estimates)
}

# function to compute the mean vector - prior
mean_fun_one <- function(hypers, train_points, l_train, prior_mean_str){
  p <- 1e-4
  H <- cbind(prior_mean_str(train_points))

  beta_values <- beta_estim_one(train_points,
                                l_train, hypers, p, prior_mean_str)
  mean_vec <- H %*% beta_values

  return(mean_vec)
}
# function to compute the covarinace - prior
cov_fun_one <- function(hypers, train_points){

  sigma_sq_1 <- hypers[1]; lambda_1 <- hypers[2]
  K11 <- sigma_sq_1*kern_gaussian(train_points, train_points, 1, lambda_1)

  return(K11)
}

#' Compute the profile log-likelihood of the intercept coefficients.
#'
#' Compute the profile log-likelihood of the intercept coefficients.
#' @param hypers A vector containing the hyperparameters of the Gaussian porcess: signal variance and length scale parameters.
#' @param train_points A vector containing the training points.
#' @param l_train A vector containing the function approximation at the training points.
#' @param prior_mean_str A function giving the structure of the prior mean.
#' @return A vector \code{b_prof_log_lik_one} containing the profile log-likelihood
#' at the given values of hyperparameters and training points.
b_profile_one <- function(hypers, train_points,
                          l_train, prior_mean_str){

  mean_vec <- mean_fun_one(hypers, train_points,
                           l_train,
                           prior_mean_str)
  cov_matrix <- cov_fun_one(hypers, train_points)
  y_train_all <- c(l_train)
  b_prof_log_lik_one <- mvtnorm::dmvnorm(y_train_all, mean = mean_vec, sigma = cov_matrix, log = TRUE)

  return(b_prof_log_lik_one)
}

#' Compute the block covariance matrix V.
#'
#' Compute the block covariance matrix V of the two level Gaussian process.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @param hypers A vector containing the hyperparameters of the Gaussian porcess of each level:
#' signal variance, length scale and autoregressive parameters.
#' @return A block matrix \code{V} containing the kernel matrices for both levels of
#' approximation.
cov_V <- function(train_points_low, train_points_high, hypers) {

  sigma_sq_1 <- hypers[1]; sigma_sq_2 <- hypers[3];
  lambda_1 <- hypers[2]; lambda_2 <- hypers[4]; r_1 <- hypers[5]

  kernel_A1_d1 <- kern_gaussian(train_points_low, train_points_low, 1, lambda_1)
  kernel_A2_d2 <- kern_gaussian(train_points_high, train_points_high, 1, lambda_2)
  kernel_A1_d1d2 <- kern_gaussian(train_points_low, train_points_high, 1, lambda_1)
  kernel_A1_d2d1 <- kern_gaussian(train_points_high, train_points_low, 1, lambda_1)
  kernel_A1_d2 <- kern_gaussian(train_points_high, train_points_high, 1, lambda_1)

  V_11 <- sigma_sq_1 * kernel_A1_d1
  V_12 <- r_1 * sigma_sq_1 * kernel_A1_d1d2
  V_21 <- r_1 * sigma_sq_1 * kernel_A1_d2d1
  V_22 <- r_1^2 * sigma_sq_1 * kernel_A1_d2 + sigma_sq_2 * kernel_A2_d2

  V <- rbind(cbind(V_11, V_12), cbind(V_21, V_22))

  return(V)
}

#' Compute the estimated intercept parameters
#'
#' Compute the estimated intercept parameters
#' of the low-level and high-level.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @param l_train_low A vector containing the low-level function approximation at the low-level training points.
#' @param l_train_high A vector containing the high-level function approximation at the high-level training points.
#' @param hypers A vector containing the hyperparameters of the Gaussian process of each level:
#' signal variance, length scale and autoregressive parameters.
#' @param prior_mean_str A function giving the structure of the prior mean.
#' @return A list \code{beta_estimates} containing the estimated intercept parameters
#' for the low and high-levels.
beta_estim <- function(train_points_low, train_points_high,
                       l_train_low, l_train_high,
                       hypers, prior_mean_str){

  p <- 1e-4
  r_1 <- hypers[5]
  z <- matrix(c(l_train_low, l_train_high), nrow = (length(l_train_low) + length(l_train_high)), ncol = 1)

  H_bottom <- cbind(r_1 * prior_mean_str(train_points_high), prior_mean_str(train_points_high))
  H_top <- cbind(prior_mean_str(train_points_low))
  for (k in (dim(H_top)[2]+1):dim(H_bottom)[2]) {
    H_top <- cbind(H_top, 0)
  }
  H <- rbind(H_top, H_bottom)

  V_z <- cov_V(train_points_low, train_points_high, hypers)
  V_inv <- solve(V_z + p * diag(1, nrow = nrow(V_z)))

  beta_estimates <- solve(t(H) %*% V_inv %*% H) %*% t(H) %*% V_inv %*% (z)

  return(beta_estimates)
}

#' Compute the mean vector of the Gaussian process posterior.
#'
#' Compute the mean vector of the two level Gaussian process posterior.
#' @param hypers A vector containing the hyperparameters of the Gaussian process of each level:
#' signal variance, length scale and autoregressive parameters.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @param l_train_low A vector containing the low-level function approximation at the low-level training points.
#' @param l_train_high A vector containing the high-level function approximation at the high-level training points.
#' @param prior_mean_str A function giving the structure of the prior mean.
#' @return A vector \code{mean_vec} containing the posterior mean vector of the
#' GP posterior for the two-level approximation.
mean_fun_two <- function(hypers, train_points_low, train_points_high,
                         l_train_low, l_train_high,
                         prior_mean_str){

  r_1 <- hypers[5]

  H_bottom <- cbind(r_1 * prior_mean_str(train_points_high), prior_mean_str(train_points_high))
  H_top <- cbind(prior_mean_str(train_points_low))
  for (k in (dim(H_top)[2]+1):dim(H_bottom)[2]) {
    H_top <- cbind(H_top, 0)
  }
  H <- rbind(H_top, H_bottom)

  beta_values <- beta_estim(train_points_low, train_points_high,
                            l_train_low, l_train_high,
                            hypers, prior_mean_str)
  mean_vec <- H %*% beta_values

  return(mean_vec)
}

#' Compute the covariance matrix of the Gaussian process posterior.
#'
#' Compute the ovariance matrix of the two-level Gaussian process posterior.
#' @param hypers A vector containing the hyperparameters of the Gaussian process of each level:
#' signal variance, length scale and autoregressive parameters.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @return A block matrix \code{cov_matr_block} containing the posterior covariance of the
#' GP posterior for the two-level approximation.
cov_fun_two <- function(hypers, train_points_low, train_points_high){

  sigma_sq_1 <- hypers[1]; sigma_sq_2 <- hypers[3]
  lambda_1 <- hypers[2]; lambda_2 <- hypers[4]
  r_1 <- hypers[5]

  # covariance
  K11 <- sigma_sq_1*kern_gaussian(train_points_low, train_points_low, 1, lambda_1)
  K12 <- r_1*sigma_sq_1*kern_gaussian(train_points_low, train_points_high, 1, lambda_1)
  K21 <- r_1*sigma_sq_1*kern_gaussian(train_points_high, train_points_low, 1, lambda_1)
  K22 <- r_1^2*sigma_sq_1*kern_gaussian(train_points_high, train_points_high, 1, lambda_1) +
    sigma_sq_2*kern_gaussian(train_points_high, train_points_high, 1, lambda_2)

  cov_matr_block <- rbind(cbind(K11, K12), cbind(K21, K22))

  return(cov_matr_block)
}


#' Compute the profile log-likelihood of the intercept coefficients.
#'
#' Compute the profile log-likelihood of the intercept coefficients.
#' @param hypers A vector containing the hyperparameters of the Gaussian process of each level:
#' signal variance, length scale and autoregressive parameters.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @param l_train_low A vector containing the low-level function approximation at the low-level training points.
#' @param l_train_high A vector containing the high-level function approximation at the high-level training points.
#' @param prior_mean_str A function giving the structure of the prior mean.
#' @return A vector \code{b_prof_log_lik_two} containing the profile log-likelihood
#' at the given values of hyperparameters and training points for the two-level approximation.
b_profile <- function(hypers, train_points_low, train_points_high,
                      l_train_low, l_train_high, prior_mean_str){

  mean_vec <- mean_fun_two(hypers, train_points_low, train_points_high,
                           l_train_low, l_train_high,
                           prior_mean_str)
  cov_matrix <- cov_fun_two(hypers, train_points_low, train_points_high)

  # given beta we have the beta profile log marginal likelihood
  # y values for both
  y_train_all <- c(l_train_low, l_train_high)
  b_prof_log_lik_two <- mvtnorm::dmvnorm(y_train_all, mean = mean_vec, sigma = cov_matrix, log = TRUE)

  return(b_prof_log_lik_two)
}


#' Hyperparameter estimation of the two-level Gaussian process
#'
#' Compute the estimated hyperparameters for each level: signal variance, length scale parameters
#' and the autoregressive parameter of the two-level Gaussian process using the Nelder-Mead method of optimisation.
#' @param hypers A vector containing the initial values of the hyperparameters,
#' used for optimisation, of the Gaussian process of each level:
#' signal variance, length scale and autoregressive parameters.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @param l_train_low A vector containing the low-level function approximation at the low-level training points.
#' @param l_train_high A vector containing the high-level function approximation at the high-level training points.
#' @param prior_mean_str A function giving the structure of the prior mean.
#' @return A list containing the estimated hyperparameters, \code{optim_hypers},
#' and the optimisation result \code{optim_result}.
#' @export
hyper_optim <- function(hypers, train_points_low, train_points_high,
                        l_train_low, l_train_high, prior_mean_str){

  sigma_sq_1 <- hypers[1]; sigma_sq_2 <- hypers[3]
  lambda_1 <- hypers[2]; lambda_2 <- hypers[4]; r_1 <- hypers[5]

  optim_NM <- stats::optim(c(sigma_sq_1, lambda_1, sigma_sq_2, lambda_2, r_1), fn = b_profile, method = "Nelder-Mead",
                    train_points_low = train_points_low, train_points_high = train_points_high,
                    l_train_low = l_train_low,
                    l_train_high = l_train_high, prior_mean_str = prior_mean_str,
                    control = list(fnscale = -1, maxit = 1000))
  optim_hypers <- optim_NM$par
  names(optim_hypers) <- c("ssquare1", "l1","ssquare2", "l2", "r1")
  print(paste("The signal variance parameter of the low-level is", round(optim_hypers["ssquare1"], 5)))
  print(paste("The length-scale parameter of the low-level is", round(optim_hypers["l1"], 5)))
  print(paste("The signal variance parameter of the high-level is", round(optim_hypers["ssquare2"], 5)))
  print(paste("The length-scale parameter of the high-level is", round(optim_hypers["l2"], 5)))
  print(paste("The autoregressive parameter is", round(optim_hypers["r1"], 5)))
  return(list(optim_hypers = optim_hypers, optim_result = optim_NM))
}

#' Gaussian process posterior plot for the two-level approximation
#'
#' Plot the Gaussian process posterior mean and credible intervals for the two-level approximation.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @param test_points A vector containing the test points used for prediction.
#' @param l_train_low A vector containing the low-level function approximation at the low-level training points.
#' @param l_train_high A vector containing the high-level function approximation at the high-level training points.
#' @param gp_mean A vector containing Gaussian process posterior mean - result of the \code{gp_post()} function.
#' @param gp_cov A matrix containing Gaussian process posterior covariance - result of the \code{gp_post()} function.
#' @param accurate Exact - Accurate approximation of the high-level if available. Default value is \code{NULL}
#' @param train_points_mid A vector containing the training points of the middle-level.
#' Default value is \code{NULL}. Used when using three levels of approximation.
#' @param l_train_mid A vector containing the middle-level function approximation at the middle-level training points.
#' Default value is \code{NULL} Used when using three levels of approximation.
#' @return Plot of the GP posterior distribution.
#' @export
gp_posterior_plot <- function(train_points_low, train_points_high, test_points,
                              l_train_low, l_train_high, gp_mean, gp_cov, accurate = NULL,
                              train_points_mid  = NULL, l_train_mid = NULL){

  gp_plot <- ggplot2::ggplot()+
    ggplot2::geom_line(data = data.frame(test_points, gp_mean),
                       ggplot2::aes(test_points, gp_mean), color = "red", lwd = 1)+
    ggplot2::geom_ribbon(data = data.frame(test_points, y_up = gp_mean + 1.96*sqrt(diag(gp_cov)),
                                  y_low = gp_mean - 1.96*sqrt(diag(gp_cov))),
                         ggplot2::aes(test_points, ymax = y_up, ymin = y_low), fill = "skyblue", alpha = 0.5)+
    ggplot2::geom_point(data = data.frame(train_points_low, l_train_low),
                        ggplot2::aes(train_points_low, l_train_low), size = 2.5, col = "green")+
    ggplot2::geom_point(data = data.frame(train_points_high, l_train_high),
                        ggplot2::aes(train_points_high, l_train_high), size = 2.5, col = "blue")+
    ggplot2::theme_minimal()+
    ggplot2::labs(title = "",
         subtitle = "",
         y = "", x = "x") +
    ggplot2::theme(axis.text=ggplot2::element_text(size=14), axis.title=ggplot2::element_text(size=16))

  if (is.null(accurate) == F) {
    gp_plot <- gp_plot + ggplot2::geom_line(data = data.frame(test_points, accurate),
                                            ggplot2::aes(test_points, accurate), color = "orange", linetype = "dashed", lwd = 1)
  }
  if (is.null(train_points_mid) == F) {
    gp_plot <- gp_plot + ggplot2::geom_point(data = data.frame(train_points_mid, l_train_mid),
                                             ggplot2::aes(train_points_mid, l_train_mid), size = 2.5, col = "pink")
  }

  return(gp_plot)
}

#' Expected gain in utility for two-level Gaussian processes
#'
#' Compute the experimental design of each level of approximation using expected gain in utility.
#' @param result_gp A list containing the posterior mean vector and the posterior covariance matrix of the
#' multi-level GP posterior distribution.
#' @param hypers A vector containing the initial values of the hyperparameters,
#' used for optimisation, of the Gaussian process of each level:
#' signal variance, length scale and autoregressive parameters.
#' @param hyper_initial A vector containing the initial values of the hyperparameters for estimation.
#' @param egu_min Expected gain in utility minimum tolerance as a stopping rule.
#' @param max_points Maximum points to be added to the experimental design as a stopping rule.
#' @param cand_points A vector containing the candidate points.
#' @param log_lik_high Function of the high-level points giving the high-level approximation.
#' @param log_lik_low Function of the low-level points giving the low-level approximation.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @param l_train_low A vector containing the low-level function approximation at the low-level training points.
#' @param l_train_high A vector containing the high-level function approximation at the high-level training points.
#' @param test_points A vector containing the test points.
#' @param cost_h Constant indicating the cost of the high-level.
#' @param cost_l Constant indicating the cost of the low-level.
#' @param prior_mean_str A function giving the structure of the prior mean.
#' @return A list containing the estimated hyperparameters of each level, \code{ssquare1}, \code{ssquare2}, \code{l1} and \code{l2},
#' the autoregressive parameter \code{r1},
#' the training points of each level, \code{train_points_low} and \code{train_points_high},
#' and the function evaluations \code{l_train_low} and \code{l_train_high},
#' the EGU per cost for each level \code{egu_low} and \code{egu_high}
#' after the addition of new points and the new points added to the design, \code{new_point}.
#' @export
egu <- function(result_gp, hypers, hyper_initial, egu_min, max_points,
                cand_points, log_lik_high, log_lik_low,
                train_points_low, train_points_high,
                l_train_low, l_train_high, test_points, cost_h, cost_l, prior_mean_str){

  p <- 1e-4
  # current utility - we consider the maximum of the posterior mean - this is our target
  exp_util_current <- max(result_gp$gp_mean)
  # hyperparameters initial values
  ssquare1 <- hyper_initial[1]; l1 <- hyper_initial[2];
  ssquare2 <- hyper_initial[3]; l2 <- hyper_initial[4]; r1 <- hyper_initial[5]

  # initial values
  y_quad <- NULL; max_mean_post_high <- NULL; max_mean_post_low <- NULL
  sum_vmu_high <- NULL; approx_eu_high <- NULL; egu_high <- NULL; sum_vmu_low <- NULL; approx_eu_low <- NULL
  egu_low <- NULL; high_p <- NULL; low_p <- NULL; new_p <- NULL; new_log <- NULL
  egu_cost_h <- NULL; egu_cost_l <- NULL
  # store to plot later
  egu_list <- list()
  egu_high_list <- list(); egu_low_list <- list(); post_mean_list <- list(); post_covar_list <- list()
  post_exp_util_list <- list(); train_points_low_list <- list()
  train_points_high_list <- list(); l_train_low_list <- list(); l_train_high_list <- list()
  # to store them
  ssquare1_egu <- c(); l1_egu <- c(); ssquare2_egu <- c(); l2_egu <- c(); r_1_new <- c()
  ssquare1_egu[1] <- hypers[["ssquare1"]]; l1_egu[1] <- hypers[["l1"]];
  ssquare2_egu[1] <- hypers[["ssquare2"]]; l2_egu[1] <- hypers[["l2"]];
  r_1_new[1] <- hypers[["r1"]];
  re_hyper_estim <- hypers

  new_points <- 1; egu_stop <- 10
  n_q <- 50 # how many quadrature points to use for the approximation - default
  cost <- 0 # total cost
  cand_p <- length(cand_points)
  while(egu_stop >= egu_min & new_points < max_points){
    print(paste("new point", new_points))
    for (cand_i in 1:cand_p) {
      print(paste("candidate point", cand_i))
      # Step 1. Find the mean and sd of posterior predictive distribution at the candidate point of the GP approximation
      post_results_cand <- gp_post(train_points_low, train_points_high, l_train_low, l_train_high,
                                   cand_points[cand_i], re_hyper_estim,
                                   prior_mean_str, 2)

      mean_cand <- post_results_cand$gp_mean
      sd_cand <- sqrt(diag(post_results_cand$gp_cov))

      # Step 2. Find quadrature points
      quadr_two <- fastGHQuad::gaussHermiteData(n_q)
      quad_points_two <- quadr_two$x; quad_weights_two <- quadr_two$w
      for (i in 1:n_q) {
        y_quad[i] <- mean_cand + sqrt(2) * sd_cand * quad_points_two[i]
      }

      # Step 3. Find the maximum of the new posterior mean given the additional candidate point in each level of
      # approximation separately and response from quadrature at the candidate point for each quadrature point from Step 2.
      train_p_high_new <- c(train_points_high, cand_points[cand_i]) # high
      train_p_low_new <- c(train_points_low, cand_points[cand_i]) # low

      for (i in 1:n_q) {
        log_train_high_new <- c(l_train_high, y_quad[i]) # high
        log_train_low_new <- c(l_train_low, y_quad[i])  # low

        # Fit a new GP with the new point (cand, y_quad_two) added to the training set of each level separately
        # keeping the hyperparameters fixed
        hyper_est_high <- re_hyper_estim; hyper_est_low <- re_hyper_estim

        # added to the high-level
        post_results_high <-  gp_post(train_points_low, train_p_high_new, l_train_low, log_train_high_new,
                                      test_points, hyper_est_high,
                                      prior_mean_str, 2)
        # added to the low-level
        post_results_low <- gp_post(train_p_low_new, train_points_high, log_train_low_new, l_train_high,
                                    test_points, hyper_est_low, prior_mean_str, 2)

        # Find the maximum of the new posterior mean for each new design (added quad_points)
        max_mean_post_high[i] <- max(post_results_high$gp_mean)
        max_mean_post_low[i] <- max(post_results_low$gp_mean)
      }

      # Step 4. calculate the sum of ... as an approximation to the expected utility
      v_quad_two <- quad_weights_two/sqrt(pi)
      sum_vmu_high[1] <- v_quad_two[1]*max_mean_post_high[1] # high
      sum_vmu_low[1] <- v_quad_two[1]*max_mean_post_low[1] # low
      for (j in 2:n_q) {
        sum_vmu_high[j] <- sum_vmu_high[j-1] + v_quad_two[j]*max_mean_post_high[j] # high
        sum_vmu_low[j] <- sum_vmu_low[j-1] + v_quad_two[j]*max_mean_post_low[j]# low
      }
      # approximate EU for each candidate point
      approx_eu_high[cand_i] <- sum_vmu_high[n_q]
      approx_eu_low[cand_i] <- sum_vmu_low[n_q]
      # compute expected gain in utility for each candidate point
      egu_high[cand_i] <- approx_eu_high[cand_i] - exp_util_current # high
      egu_low[cand_i] <- approx_eu_low[cand_i] - exp_util_current # low
    }
    # save EGU
    # EGU per cost for all candidate points
    egu_high_list[[new_points]] <- egu_high/cost_h
    egu_low_list[[new_points]] <- egu_low/cost_l

    # computing the EGU per COST for each level
    egu_cost_h[new_points] <- max(egu_high, na.rm = T)/cost_h
    egu_cost_l[new_points] <- max(egu_low, na.rm = T)/cost_l

    # compare EGU per cost
    if (egu_cost_h[new_points] > egu_cost_l[new_points]) {
      # new point at the high-level
      high_p[new_points] <- T
      low_p[new_points] <- F
      new_p[new_points] <- cand_points[which.max(egu_high)]
      train_points_high <- c(train_points_high, new_p[new_points])
      # what is the approximation method used - evaluate at the new points
      l_train_high <- c(l_train_high, log_lik_high(new_p[new_points]))
      new_log[new_points] <- l_train_high[length(l_train_high)]
      cost <- cost + cost_h # total cost
    }else{
      # new point at the low-level
      high_p[new_points] <- F
      low_p[new_points] <- T
      new_p[new_points] <- cand_points[which.max(egu_low)]
      train_points_low <- c(train_points_low, new_p[new_points])
      # computing the log-likelihood at the new point
      l_train_low <- c(l_train_low, log_lik_low(new_p[new_points]))
      new_log[new_points] <- l_train_low[length(l_train_low)]
      cost <- cost + cost_l # total cost
    }
    # re-estimate hyperparameters
    optim_hyper_EGU <- invisible(hyper_optim(c(ssquare1, l1, ssquare2, l2, r1), train_points_low, train_points_high,
                                             l_train_low, l_train_high, prior_mean_str))
    re_hyper_estim <- optim_hyper_EGU$optim_hypers

    ssquare1_egu[new_points+1] <- re_hyper_estim[["ssquare1"]]; l1_egu[new_points+1] <- re_hyper_estim[["l1"]]
    ssquare2_egu[new_points+1] <- re_hyper_estim[["ssquare2"]]; l2_egu[new_points+1] <- re_hyper_estim[["l2"]]
    r_1_new[new_points+1] <- re_hyper_estim[["r1"]]

    # new posterior - new point added
    post_n_point <- gp_post(train_points_low, train_points_high, l_train_low, l_train_high,
                            test_points, re_hyper_estim,
                            prior_mean_str, 2)

    # store new training points - for plotting purposes
    train_points_low_list[[new_points]] <- train_points_low; train_points_high_list[[new_points]] <- train_points_high
    # store new log-likelihoods - for plotting purposes
    l_train_low_list[[new_points]] <- l_train_low; l_train_high_list[[new_points]] <- l_train_high

    # new expected untility
    exp_util_current <- max(post_n_point$gp_mean)

    # save posterior
    post_mean_list[[new_points]] <- post_n_point$gp_mean; post_covar_list[[new_points]] <- post_n_point$gp_cov
    post_exp_util_list[[new_points]] <- exp_util_current

    new_points <- new_points + 1 # moving on to the next point
    egu_stop <- max(egu_low, egu_high, na.rm = T)
  }
  print("The training points of the low-level are")
  print(round(train_points_low_list[[new_points-1]], 3))
  print("The training points of the high-level are")
  print(round(train_points_high_list[[new_points-1]], 3))
  return(list(ssquare1 = ssquare1_egu, ssquare2 = ssquare2_egu, l1 = l1_egu, l2 = l2_egu, r1 = r_1_new,
              train_points_low = train_points_low_list, train_points_high = train_points_high_list,
              l_train_low = l_train_low_list, l_train_high = l_train_high_list,
              high_p = high_p, low_p = low_p,
              egu_high = egu_high_list, egu_low = egu_low_list, post_mean = post_mean_list,
              post_covar = post_covar_list, total_cost = cost, points_n = new_points-1,
              new_point = new_p, new_lik = new_log))
}


cov_V_three <- function(train_points_low, train_points_mid, train_points_high, hypers) {

  r_1 <- hypers[7]; r_2 <- hypers[8]
  lambda_1 <- hypers[2]; lambda_2 <- hypers[4]; lambda_3 <- hypers[6]
  sigma_sq_1 <- hypers[1]; sigma_sq_2 <- hypers[3]; sigma_sq_3 <- hypers[5]

  # calculate the matrices of correlation
  # Note : we compute correlation so the s_2 in the kern_gaussian function is 1
  kernel_A1_d1 <- kern_gaussian(train_points_low, train_points_low, 1, lambda_1)
  kernel_A2_d2 <- kern_gaussian(train_points_mid, train_points_mid, 1, lambda_2)
  kernel_A3_d3 <- kern_gaussian(train_points_high, train_points_high, 1, lambda_3)

  kernel_A1_d1d2 <- kern_gaussian(train_points_low, train_points_mid, 1, lambda_1)
  kernel_A1_d1d3 <- kern_gaussian(train_points_low, train_points_high, 1, lambda_1)
  kernel_A1_d2d3 <- kern_gaussian(train_points_mid, train_points_high, 1, lambda_1)

  kernel_A1_d2d1 <- kern_gaussian(train_points_mid, train_points_low, 1, lambda_1)
  kernel_A1_d3d1 <- kern_gaussian(train_points_high, train_points_low, 1, lambda_1)
  kernel_A1_d3d2 <- kern_gaussian(train_points_high, train_points_mid, 1, lambda_1)

  kernel_A1_d2 <- kern_gaussian(train_points_mid, train_points_mid, 1, lambda_1)
  kernel_A1_d3 <- kern_gaussian(train_points_high, train_points_high, 1, lambda_1)

  kernel_A2_d2d3 <- kern_gaussian(train_points_mid, train_points_high, 1, lambda_2)
  kernel_A2_d3d2 <- kern_gaussian(train_points_high, train_points_mid, 1, lambda_2)
  kernel_A2_d3 <- kern_gaussian(train_points_high, train_points_high, 1, lambda_2)

  # compute the entries of the block matrix
  V_11 <- sigma_sq_1 * kernel_A1_d1
  V_22 <- sigma_sq_2 * kernel_A2_d2 + r_1^2 * sigma_sq_1 * kernel_A1_d2
  V_33 <- sigma_sq_3 * kernel_A3_d3 + r_2^2 * sigma_sq_2 * kernel_A2_d3 + r_1^2 * r_2^2 * sigma_sq_1 * kernel_A1_d3

  V_12 <- r_1 * sigma_sq_1 * kernel_A1_d1d2
  V_21 <- r_1 * sigma_sq_1 * kernel_A1_d2d1
  V_13 <- r_1 * r_2 * sigma_sq_1 * kernel_A1_d1d3
  V_31 <- r_1 * r_2 * sigma_sq_1 * kernel_A1_d3d1
  V_23 <- r_2 * sigma_sq_2 * kernel_A2_d2d3 + r_1^2 * r_2 * sigma_sq_1 * kernel_A1_d2d3
  V_32 <- r_2 * sigma_sq_2 * kernel_A2_d3d2 + r_1^2 * r_2 * sigma_sq_1 * kernel_A1_d3d2

  V <- rbind(cbind(V_11, V_12, V_13), cbind(V_21, V_22, V_23), cbind(V_31, V_32, V_33))
  return(V)
}

#' Compute the estimated intercept parameters
#'
#' Compute the estimated intercept parameters
#' of the low-level, middle-level and high-level.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_mid A vector containing the training points of the middle-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @param l_train_low A vector containing the low-level function approximation at the low-level training points.
#' @param l_train_mid A vector containing the middle-level function approximation at the middle-level training points.
#' @param l_train_high A vector containing the high-level function approximation at the high-level training points.
#' @param hypers A vector containing the hyperparameters of the Gaussian process of each level:
#' signal variance, length scale and autoregressive parameters.
#' @param prior_mean_str A function giving the structure of the prior mean.
#' @return A list \code{beta_estimates} containing the estimated intercept parameters
#' for the low, middle and high levels.
beta_estim_three <- function(train_points_low, train_points_mid, train_points_high,
                             l_train_low, l_train_mid, l_train_high,
                             hypers, prior_mean_str){

  r_1 <- hypers[7]; r_2 <- hypers[8]; p <- 1e-4

  z <- matrix(c(l_train_low, l_train_mid, l_train_high), nrow = (length(l_train_low) + length(l_train_mid) +
                                                                   + length(l_train_high)), ncol = 1)

  H_bottom <- cbind(r_1 * r_2 * prior_mean_str(train_points_high), r_2 * prior_mean_str(train_points_high), prior_mean_str(train_points_high))
  H_mid <- cbind(r_1 * prior_mean_str(train_points_mid), prior_mean_str(train_points_mid))
  H_top <- cbind(prior_mean_str(train_points_low))
  for (k in (dim(H_top)[2]+1):dim(H_bottom)[2]) {
    H_top <- cbind(H_top, 0)
  }
  for (k in (dim(H_mid)[2]+1):dim(H_bottom)[2]) {
    H_mid <- cbind(H_mid, 0)
  }
  H <- rbind(H_top, H_mid, H_bottom)

  V_z <- cov_V_three(train_points_low, train_points_mid, train_points_high, hypers)
  V_inv <- solve(V_z + p * diag(1, nrow = nrow(V_z)))

  beta_estimates <- solve(t(H) %*% V_inv %*% H) %*% t(H) %*% V_inv %*% (z)
  return(beta_estimates)
}

#' Compute the mean vector of the Gaussian process posterior.
#'
#' Compute the mean vector of the three level Gaussian process posterior.
#' @param hypers A vector containing the hyperparameters of the Gaussian process of each level:
#' signal variance, length scale and autoregressive parameters.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_mid A vector containing the training points of the middle-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @param l_train_low A vector containing the low-level function approximation at the low-level training points.
#' @param l_train_mid A vector containing the middle-level function approximation at the middle-level training points.
#' @param l_train_high A vector containing the high-level function approximation at the high-level training points.
#' @param prior_mean_str A function giving the structure of the prior mean.
#' @return A vector \code{mean_vec} containing the posterior mean vector of the
#' GP posterior for the three-level approximation.
mean_fun_three <- function(hypers, train_points_low, train_points_mid, train_points_high,
                           l_train_low, l_train_mid, l_train_high,
                           prior_mean_str){

  r_1 <- hypers[7]; r_2 <- hypers[8]

  H_bottom <- cbind(r_1 * r_2 * prior_mean_str(train_points_high), r_2 * prior_mean_str(train_points_high), prior_mean_str(train_points_high))
  H_mid <- cbind(r_1 * prior_mean_str(train_points_mid), prior_mean_str(train_points_mid))
  H_top <- cbind(prior_mean_str(train_points_low))
  for (k in (dim(H_top)[2]+1):dim(H_bottom)[2]) {
    H_top <- cbind(H_top, 0)
  }
  for (k in (dim(H_mid)[2]+1):dim(H_bottom)[2]) {
    H_mid <- cbind(H_mid, 0)
  }
  H <- rbind(H_top, H_mid, H_bottom)

  beta_values <- beta_estim_three(train_points_low, train_points_mid, train_points_high,
                                  l_train_low, l_train_mid, l_train_high,
                                  hypers, prior_mean_str)
  mean_vec <- H %*% beta_values

  return(mean_vec)
}

#' Compute the covariance matrix of the Gaussian process posterior.
#'
#' Compute the ovariance matrix of the three-level Gaussian process posterior.
#' @param hypers A vector containing the hyperparameters of the Gaussian process of each level:
#' signal variance, length scale and autoregressive parameters.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_mid A vector containing the training points of the middle-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @return A block matrix \code{cov_matr_block} containing the posterior covariance of the
#' GP posterior for the three-level approximation.
cov_fun_three <- function(hypers, train_points_low, train_points_mid, train_points_high){

  sigma_sq_1 <- hypers[1]; sigma_sq_2 <- hypers[3]; sigma_sq_3 <- hypers[5]
  lambda_1 <- hypers[2]; lambda_2 <- hypers[4]; lambda_3 <- hypers[6]
  r_1 <- hypers[7]; r_2 <- hypers[8]

  # covariance
  K11 <- sigma_sq_1*kern_gaussian(train_points_low, train_points_low, 1, lambda_1)
  K12 <- r_1*sigma_sq_1*kern_gaussian(train_points_low, train_points_mid, 1, lambda_1)
  K13 <- r_1*r_2*sigma_sq_1*kern_gaussian(train_points_low, train_points_high, 1, lambda_1)

  K21 <- r_1*sigma_sq_1*kern_gaussian(train_points_mid, train_points_low, 1, lambda_1)
  K22 <- r_1^2*sigma_sq_1*kern_gaussian(train_points_mid, train_points_mid, 1, lambda_1) +
    sigma_sq_2*kern_gaussian(train_points_mid, train_points_mid, 1, lambda_2)
  K23 <- r_1^2*r_2*sigma_sq_1*kern_gaussian(train_points_mid, train_points_high, 1, lambda_1) +
    r_2*sigma_sq_2*kern_gaussian(train_points_mid, train_points_high, 1, lambda_1)

  K31 <- r_1*r_2*sigma_sq_1*kern_gaussian(train_points_high, train_points_low, 1, lambda_1)
  K32 <- r_1^2*r_2*sigma_sq_1*kern_gaussian(train_points_high, train_points_mid, 1, lambda_1) +
    r_2*sigma_sq_2*kern_gaussian(train_points_high, train_points_mid, 1, lambda_1)
  K33 <- r_1^2*r_2^2*sigma_sq_1*kern_gaussian(train_points_high, train_points_high, 1, lambda_1) +
    r_2^2*sigma_sq_2*kern_gaussian(train_points_high, train_points_high, 1, lambda_2) +
    sigma_sq_3*kern_gaussian(train_points_high, train_points_high, 1, lambda_3)

  cov_matr_block <- rbind(cbind(K11, K12, K13), cbind(K21, K22, K23), cbind(K31, K32, K33))

  return(cov_matr_block)
}


#' Compute the profile log-likelihood of the intercept coefficients.
#'
#' Compute the profile log-likelihood of the intercept coefficients.
#' @param hypers A vector containing the hyperparameters of the Gaussian process of each level:
#' signal variance, length scale and autoregressive parameters.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_mid A vector containing the training points of the middle-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @param l_train_low A vector containing the low-level function approximation at the low-level training points.
#' @param l_train_mid A vector containing the middle-level function approximation at the middle-level training points.
#' @param l_train_high A vector containing the high-level function approximation at the high-level training points.
#' @param prior_mean_str A function giving the structure of the prior mean.
#' @return A vector \code{b_prof_log_lik_three} containing the profile log-likelihood
#' at the given values of hyperparameters and training points for the three-level approximation.
b_profile_three <- function(hypers, train_points_low, train_points_mid, train_points_high,
                            l_train_low, l_train_mid, l_train_high,
                            prior_mean_str){

  sigma_sq_1 <- hypers[1]; sigma_sq_2 <- hypers[3]; sigma_sq_3 <- hypers[5]
  lambda_1 <- hypers[2]; lambda_2 <- hypers[4]; lambda_3 <- hypers[6]
  r_1 <- hypers[7]; r_2 <- hypers[8]

  mean_vec <- mean_fun_three(hypers, train_points_low, train_points_mid, train_points_high,
                             l_train_low, l_train_mid, l_train_high,
                             prior_mean_str)
  cov_matrix <- cov_fun_three(hypers, train_points_low, train_points_mid, train_points_high)

  # given beta we have the beta profile log marginal likelihood
  # y values for all the levels
  y_train_all <- c(l_train_low, l_train_mid, l_train_high)
  b_prof_log_lik_three <- mvtnorm::dmvnorm(y_train_all, mean = mean_vec, sigma = cov_matrix, log = TRUE)

  return(b_prof_log_lik_three)
}

#' Hyperparameter estimation of the three-level Gaussian process
#'
#' Compute the estimated hyperparameters for each level: signal variance, length scale parameters
#' and the autoregressive parameter of the three-level Gaussian process using the Nelder-Mead method of optimisation.
#' @param hypers A vector containing the initial values of the hyperparameters,
#' used for optimisation, of the Gaussian process of each level:
#' signal variance, length scale and autoregressive parameters.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_mid A vector containing the training points of the middle-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @param l_train_low A vector containing the low-level function approximation at the low-level training points.
#' @param l_train_mid A vector containing the middle-level function approximation at the middle-level training points.
#' @param l_train_high A vector containing the high-level function approximation at the high-level training points.
#' @param prior_mean_str A function giving the structure of the prior mean.
#' @return A list containing the posterior estimated hyperparameters, \code{optim_hypers},
#' and the optimisation result \code{optim_result}.
hyper_optim_three <- function(hypers, train_points_low, train_points_mid, train_points_high,
                              l_train_low, l_train_mid, l_train_high, prior_mean_str){

  sigma_sq_1 <- hypers[1]; sigma_sq_2 <- hypers[3]; sigma_sq_3 <- hypers[5]
  lambda_1 <- hypers[2]; lambda_2 <- hypers[4]; lambda_3 <- hypers[6];
  r_1 <- hypers[7]; r_2 <- hypers[8]

  optim_NM <- stats::optim(c(sigma_sq_1, lambda_1, sigma_sq_2, lambda_2, sigma_sq_3, lambda_3, r_1, r_2),
                    fn = b_profile_three, method = "Nelder-Mead",
                    train_points_low = train_points_low, train_points_mid = train_points_mid,
                    train_points_high = train_points_high, l_train_low = l_train_low, l_train_mid = l_train_mid,
                    l_train_high = l_train_high,
                    prior_mean_str = prior_mean_str, control = list(fnscale = -1, maxit = 500))
  optim_hypers <- optim_NM$par
  names(optim_hypers) <- c("ssquare1", "l1","ssquare2", "l2", "ssquare3", "l3", "r1", "r2")
  print(paste("The signal variance parameter of the low-level is", round(optim_hypers["ssquare1"], 5)))
  print(paste("The length-scale parameter of the low-level is", round(optim_hypers["l1"], 5)))
  print(paste("The signal variance parameter of the middle-level is", round(optim_hypers["ssquare2"], 5)))
  print(paste("The length-scale parameter of the middle-level is", round(optim_hypers["l2"], 5)))
  print(paste("The signal variance parameter of the high-level is", round(optim_hypers["ssquare3"], 5)))
  print(paste("The length-scale parameter of the high-level is", round(optim_hypers["l3"], 5)))
  print(paste("The autoregressive parameter between the low and middle level is", round(optim_hypers["r1"], 5)))
  print(paste("The autoregressive parameter between the middle and high level is", round(optim_hypers["r2"], 5)))
  return(list(optim_hypers = optim_hypers, optim_result = optim_NM))
}

#' Expected gain in utility for thee-level Gaussian processes
#'
#' Compute the experimental design of each level of approximation using expected gain in utility.
#' @param result_gp A list containing the posterior mean vector and the posterior covariance matrix of the
#' multi-level GP posterior distribution.
#' @param hypers A vector containing the initial values of the hyperparameters,
#' used for optimisation, of the Gaussian process of each level:
#' signal variance, length scale and autoregressive parameters.
#' @param hyper_initial A vector containing the initial values of the hyperparameters for estimation.
#' @param egu_min Expected gain in utility minimum tolerance as a stopping rule.
#' @param max_points Maximum points to be added to the experimental design as a stopping rule.
#' @param cand_points A vector containing the candidate points.
#' @param log_lik_high Function of the high-level points giving the high-level approximation.
#' @param log_lik_mid Function of the middle-level points giving the middle-level approximation.
#' @param log_lik_low Function of the low-level points giving the low-level approximation.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_mid A vector containing the training points of the middle-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @param l_train_low A vector containing the low-level function approximation at the low-level training points.
#' @param l_train_mid A vector containing the middle-level function approximation at the middle-level training points.
#' @param l_train_high A vector containing the high-level function approximation at the high-level training points.
#' @param test_points A vector containing the test points.
#' @param cost_h Constant indicating the cost of the high-level.
#' @param cost_m Constant indicating the cost of the middle-level.
#' @param cost_l Constant indicating the cost of the low-level.
#' @param prior_mean_str A function giving the structure of the prior mean.
#' @return A list containing the estimated hyperparameters of each level, \code{ssquare1}, \code{ssquare2}, \code{ssquare3},
#' \code{l1}, \code{l2}, \code{l3}, the autoregressive parameters \code{r1} and \code{r2},
#' the training points of each level, \code{train_points_low}, \code{train_points_mid} and \code{train_points_high},
#' and the function evaluations \code{l_train_low}, \code{l_train_mid} and \code{l_train_high},
#' the EGU per cost for each level \code{egu_low}, \code{egu_mid} and \code{egu_high}
#' after the addition of new points and the new points added to the design, \code{new_point}.
#' @export
egu_three <- function(result_gp, hypers, hyper_initial, egu_min, max_points,
                      cand_points, log_lik_high, log_lik_mid, log_lik_low,
                      train_points_low, train_points_mid, train_points_high,
                      l_train_low, l_train_mid, l_train_high, test_points, cost_h, cost_m, cost_l, prior_mean_str){

  p <- 1e-4
  # current utility - we consider the maximum of the posterior mean - this is our target
  exp_util_current <- max(result_gp$gp_mean)

  # hyperparameters initial values
  ssquare1 <- hyper_initial[1]; l1 <- hyper_initial[2];
  ssquare2 <- hyper_initial[3]; l2 <- hyper_initial[4]; r1 <- hyper_initial[7]
  ssquare3 <- hyper_initial[5]; l3 <- hyper_initial[6]; r2 <- hyper_initial[8]

  # initial values
  y_quad <- NULL; max_mean_post_high <- NULL; max_mean_post_mid <- NULL; max_mean_post_low <- NULL
  sum_vmu_high <- NULL; approx_eu_high <- NULL; egu_high <- NULL;
  sum_vmu_mid <- NULL; approx_eu_mid <- NULL; egu_mid <- NULL;
  sum_vmu_low <- NULL; approx_eu_low <- NULL; egu_low <- NULL;
  high_p <- NULL; mid_p <- NULL; low_p <- NULL; new_p <- NULL; new_log <- NULL
  egu_cost_h <- NULL;  egu_cost_m <- NULL; egu_cost_l <- NULL
  # store to plot later
  egu_list <- list()
  egu_high_list <- list(); egu_mid_list <- list(); egu_low_list <- list();
  post_mean_list <- list(); post_covar_list <- list()
  post_exp_util_list <- list();
  train_points_low_list <- list(); train_points_mid_list <- list(); train_points_high_list <- list();
  l_train_low_list <- list(); l_train_mid_list <- list(); l_train_high_list <- list()
  # to store them
  ssquare1_egu <- c(); l1_egu <- c(); ssquare2_egu <- c(); l2_egu <- c();
  ssquare3_egu <- c(); l3_egu <- c(); r_1_new <- c();  r_2_new <- c()
  ssquare1_egu[1] <- hypers[["ssquare1"]]; l1_egu[1] <- hypers[["l1"]];
  ssquare2_egu[1] <- hypers[["ssquare2"]]; l2_egu[1] <- hypers[["l2"]];
  ssquare3_egu[1] <- hypers[["ssquare3"]]; l3_egu[1] <- hypers[["l3"]];
  r_1_new[1] <- hypers[["r1"]]; r_2_new[1] <- hypers[["r2"]];
  re_hyper_estim <- hypers

  new_points <- 1; egu_stop <- 10
  n_q <- 50 # how many quadrature points to use for the approximation - default
  cost <- 0 # total cost
  cand_p <- length(cand_points)
  while(egu_stop >= egu_min & new_points < max_points){
    print(paste("new point", new_points))
    for (cand_i in 1:cand_p) {
      print(paste("candidate point", cand_i))
      # Step 1. Find the mean and sd of posterior predictive distribution at the candidate point of the GP approximation
      post_results_cand <- gp_post(train_points_low, train_points_high, l_train_low, l_train_high,
                                   cand_points[cand_i], re_hyper_estim,
                                   prior_mean_str, 3, train_points_mid, l_train_mid)

      mean_cand <- post_results_cand$gp_mean
      sd_cand <- sqrt(abs(diag(post_results_cand$gp_cov)))

      # Step 2. Find quadrature points
      quadr_two <- fastGHQuad::gaussHermiteData(n_q)
      quad_points_two <- quadr_two$x; quad_weights_two <- quadr_two$w
      for (i in 1:n_q) {
        y_quad[i] <- mean_cand + sqrt(2) * sd_cand * quad_points_two[i]
      }

      # Step 3. Find the maximum of the new posterior mean given the additional candidate point in each level of
      # approximation separately and response from quadrature at the candidate point for each quadrature point from Step 2.
      train_p_high_new <- c(train_points_high, cand_points[cand_i]) # high
      train_p_mid_new <- c(train_points_mid, cand_points[cand_i]) # middle
      train_p_low_new <- c(train_points_low, cand_points[cand_i]) # low

      for (i in 1:n_q) {
        log_train_high_new <- c(l_train_high, y_quad[i]) # high
        log_train_mid_new <- c(l_train_mid, y_quad[i])  # middle
        log_train_low_new <- c(l_train_low, y_quad[i])  # low

        # Fit a new GP with the new point (cand, y_quad_two) added to the training set of each level separately
        # keeping the hyperparameters fixed
        hyper_est_high <- re_hyper_estim; hyper_est_low <- re_hyper_estim

        # added to the high-level
        post_results_high <-  gp_post(train_points_low, train_p_high_new, l_train_low, log_train_high_new,
                                      test_points, hyper_est_high,
                                      prior_mean_str, 3, l_train_mid, l_train_mid)
        # added to the middle-level
        post_results_mid <-  gp_post(train_points_low, train_points_high, l_train_low, l_train_high,
                                     test_points, hyper_est_high,
                                     prior_mean_str, 3, train_p_mid_new, log_train_mid_new)
        # added to the low-level
        post_results_low <- gp_post(train_p_low_new, train_points_high, log_train_low_new, l_train_high,
                                    test_points, hyper_est_low, prior_mean_str, 3, l_train_mid, l_train_mid)

        # Find the maximum of the new posterior mean for each new design (added quad_points)
        max_mean_post_high[i] <- max(post_results_high$gp_mean)
        max_mean_post_mid[i] <- max(post_results_mid$gp_mean)
        max_mean_post_low[i] <- max(post_results_low$gp_mean)
      }

      # Step 4. calculate the sum of ... as an approximation to the expected utility
      v_quad_two <- quad_weights_two/sqrt(pi)
      sum_vmu_high[1] <- v_quad_two[1]*max_mean_post_high[1] # high
      sum_vmu_mid[1] <- v_quad_two[1]*max_mean_post_mid[1] # middle
      sum_vmu_low[1] <- v_quad_two[1]*max_mean_post_low[1] # low
      for (j in 2:n_q) {
        sum_vmu_high[j] <- sum_vmu_high[j-1] + v_quad_two[j]*max_mean_post_high[j] # high
        sum_vmu_mid[j] <- sum_vmu_mid[j-1] + v_quad_two[j]*max_mean_post_mid[j] # middle
        sum_vmu_low[j] <- sum_vmu_low[j-1] + v_quad_two[j]*max_mean_post_low[j]# low
      }
      # approximate EU for each candidate point
      approx_eu_high[cand_i] <- sum_vmu_high[n_q]
      approx_eu_mid[cand_i] <- sum_vmu_mid[n_q]
      approx_eu_low[cand_i] <- sum_vmu_low[n_q]
      # compute expected gain in utility for each candidate point
      egu_high[cand_i] <- approx_eu_high[cand_i] - exp_util_current # high
      egu_mid[cand_i] <- approx_eu_mid[cand_i] - exp_util_current # middle
      egu_low[cand_i] <- approx_eu_low[cand_i] - exp_util_current # low
    }
    # save EGU
    # EGU per cost for all candidate points
    egu_high_list[[new_points]] <- egu_high/cost_h
    egu_mid_list[[new_points]] <- egu_mid/cost_m
    egu_low_list[[new_points]] <- egu_low/cost_l

    # computing the EGU per COST for each level
    egu_cost_h[new_points] <- max(egu_high, na.rm = T)/cost_h
    egu_cost_m[new_points] <- max(egu_mid, na.rm = T)/cost_m
    egu_cost_l[new_points] <- max(egu_low, na.rm = T)/cost_l

    # compare EGU per cost
    egu_max <- max(egu_cost_l[new_points], egu_cost_m[new_points], egu_cost_h[new_points])
    if (egu_max == egu_cost_l[new_points]) {
      # new point at the low-level
      high_p[new_points] <- F; mid_p[new_points] <- F; low_p[new_points] <- T
      new_p[new_points] <- cand_points[which.max(egu_low)]
      train_points_low <- c(train_points_low, new_p[new_points])
      # what is the approximation method used - evaluate at the new point
      l_train_low <- c(l_train_low, log_lik_low(new_p[new_points]))
      new_log[new_points] <- l_train_low[length(l_train_low)]
      cost <- cost + cost_l # total cost
    }
    if (egu_max == egu_cost_m[new_points]) {
      # new point at the middle-level
      high_p[new_points] <- F; mid_p[new_points] <- T; low_p[new_points] <- F
      new_p[new_points] <- cand_points[which.max(egu_mid)]
      train_points_mid <- c(train_points_mid, new_p[new_points])
      # what is the approximation method used - evaluate at the new point
      l_train_mid <- c(l_train_mid, log_lik_mid(new_p[new_points]))
      new_log[new_points] <- l_train_mid[length(l_train_mid)]
      cost <- cost + cost_m # total cost
    }
    if (egu_max == egu_cost_h[new_points]) {
      # new point at the high-level
      high_p[new_points] <- T; mid_p[new_points] <- F; low_p[new_points] <- F
      new_p[new_points] <- cand_points[which.max(egu_high)]
      train_points_high <- c(train_points_high, new_p[new_points])
      # what is the approximation method used - evaluate at the new point
      l_train_high <- c(l_train_high, log_lik_high(new_p[new_points]))
      new_log[new_points] <- l_train_high[length(l_train_high)]
      cost <- cost + cost_h # total cost
    }

    # re-estimate hyperparameters
    optim_hyper_EGU <- invisible(hyper_optim_three(c(ssquare1, l1, ssquare2, l2, ssquare3, l3, r1, r2),
                                                   train_points_low, train_points_mid, train_points_high,
                                                   l_train_low, l_train_mid, l_train_high, prior_mean_str))
    re_hyper_estim <- optim_hyper_EGU$optim_hypers

    ssquare1_egu[new_points+1] <- re_hyper_estim[["ssquare1"]]; l1_egu[new_points+1] <- re_hyper_estim[["l1"]]
    ssquare2_egu[new_points+1] <- re_hyper_estim[["ssquare2"]]; l2_egu[new_points+1] <- re_hyper_estim[["l2"]]
    ssquare3_egu[new_points+1] <- re_hyper_estim[["ssquare3"]]; l3_egu[new_points+1] <- re_hyper_estim[["l3"]]
    r_1_new[new_points+1] <- re_hyper_estim[["r1"]]; r_2_new[new_points+1] <- re_hyper_estim[["r2"]]

    # new posterior - new point added
    post_n_point <- gp_post(train_points_low, train_points_high, l_train_low, l_train_high,
                            test_points, re_hyper_estim,
                            prior_mean_str, 3, train_points_mid, l_train_mid)

    # store new training points - for plotting purposes
    train_points_low_list[[new_points]] <- train_points_low; train_points_mid_list[[new_points]] <- train_points_mid;
    train_points_high_list[[new_points]] <- train_points_high
    # store new log-likelihoods - for plotting purposes
    l_train_low_list[[new_points]] <- l_train_low; l_train_mid_list[[new_points]] <- l_train_mid; l_train_high_list[[new_points]] <- l_train_high

    # new expected untility
    exp_util_current <- max(post_n_point$gp_mean)

    # save posterior
    post_mean_list[[new_points]] <- post_n_point$gp_mean; post_covar_list[[new_points]] <- post_n_point$gp_cov
    post_exp_util_list[[new_points]] <- exp_util_current

    new_points <- new_points + 1 # moving on to the next point
    egu_stop <- max(egu_low, egu_mid, egu_high, na.rm = T)
  }
  print("The training points of the low-level are")
  print(round(train_points_low_list[[new_points-1]], 3))
  print("The training points of the middle-level are")
  print(round(train_points_mid_list[[new_points-1]], 3))
  print("The training points of the high-level are")
  print(round(train_points_high_list[[new_points-1]], 3))
  return(list(ssquare1 = ssquare1_egu, ssquare2 = ssquare2_egu, ssquare3 = ssquare3_egu,
              l1 = l1_egu, l2 = l2_egu, l3 = l3_egu, r1 = r_1_new, r2 = r_2_new,
              train_points_low = train_points_low_list, train_points_mid = train_points_mid_list, train_points_high = train_points_high_list,
              l_train_low = l_train_low_list, l_train_mid = l_train_mid_list, l_train_high = l_train_high_list,
              high_p = high_p, mid_p = mid_p, low_p = low_p,
              egu_high = egu_high_list, egu_mid = egu_mid_list, egu_low = egu_low_list, post_mean = post_mean_list,
              post_covar = post_covar_list, total_cost = cost, points_n = new_points-1,
              new_point = new_p, new_lik = new_log))
}



#' Compute the multi-level Gaussian process posterior.
#'
#' Compute the mean vector and the covariance matrix of the two level Gaussian process posterior.
#' @param train_points_low A vector containing the training points of the low-level.
#' @param train_points_high A vector containing the training points of the high-level.
#' @param l_train_low A vector containing the low-level function approximation at the low-level training points.
#' @param l_train_high A vector containing the high-level function approximation at the high-level training points.
#' @param test_points A vector containing the test points used for prediction.
#' @param hypers A vector containing the hyperparameters of the Gaussian process of each level:
#' signal variance, length scale and autoregressive parameters.
#' @param prior_mean_str A function giving the structure of the prior mean.
#' @param approx_levels Constant indicating how many levels of approximation to use. Choose between 2 or 3.
#' @param train_points_mid A vector containing the training points of the middle-level for a three-level approximation.
#' Default value is \code{NULL}.
#' @param l_train_mid A vector containing the middle-level function approximation at the middle-level training points.
#' for a three-level approximation. Default value is \code{NULL}.
#' @return A list containing the posterior mean vector \code{gp_mean} and
#' the posterior covariance matrix \code{gp_cov}.
#' @export
gp_post <- function(train_points_low, train_points_high, l_train_low, l_train_high, test_points,
                    hypers, prior_mean_str, approx_levels,
                    train_points_mid = NULL, l_train_mid = NULL){

  p <- 1e-4

  if (approx_levels == 2) {
    train_n <- c(length(train_points_low), length(train_points_high))
    train <- 1:sum(train_n)
    train_p <- list(train_points_low, train_points_high)
    l_train <- list(l_train_low, l_train_high)
    V_train_test <- cov_V(train_p[[1]], c(train_p[[2]], test_points), hypers)
    y_train <- matrix(c(l_train[[1]], l_train[[2]]), nrow = sum(train_n), ncol = 1)
    h_prime_T <- cbind(hypers["r1"] * prior_mean_str(test_points), prior_mean_str(test_points))
    H_bottom <- cbind(hypers["r1"] * prior_mean_str(train_p[[2]]), prior_mean_str(train_p[[2]]))
    H_top <- cbind(prior_mean_str(train_p[[1]]))
    for (k in (dim(H_top)[2]+1):dim(H_bottom)[2]) {
      H_top <- cbind(H_top, 0)
    }
    H <- rbind(H_top, H_bottom)
  }
  if (approx_levels == 3) {
    train_n <- c(length(train_points_low), length(train_points_mid), length(train_points_high))
    train <- 1:sum(train_n)
    train_p <- list(train_points_low, train_points_mid, train_points_high)
    l_train <- list(l_train_low, l_train_mid, l_train_high)
    V_train_test <- cov_V_three(train_p[[1]], train_p[[2]], c(train_p[[3]], test_points), hypers)
    y_train <- matrix(c(l_train[[1]], l_train[[2]], l_train[[3]]), nrow = sum(train_n), ncol = 1)
    h_prime_T <- cbind(hypers["r1"] * hypers["r2"] * prior_mean_str(test_points), hypers["r2"] * prior_mean_str(test_points), prior_mean_str(test_points))
    H_bottom <- cbind(hypers["r1"] * hypers["r2"] * prior_mean_str(train_p[[3]]), hypers["r2"] * prior_mean_str(train_p[[3]]), prior_mean_str(train_p[[3]]))
    H_mid <- cbind(hypers["r1"] * prior_mean_str(train_p[[2]]), prior_mean_str(train_p[[2]]))
    H_top <- cbind(prior_mean_str(train_p[[1]]))
    for (k in (dim(H_top)[2]+1):dim(H_bottom)[2]) {
      H_top <- cbind(H_top, 0)}
    for (k in (dim(H_mid)[2]+1):dim(H_bottom)[2]) {
      H_mid <- cbind(H_mid, 0)}
    H <- rbind(H_top, H_mid, H_bottom)
  }

  V_22 <- matrix(V_train_test[-train, -train], nrow = length(test_points))
  V_12 <- V_train_test[train, -train]
  V_21 <- V_train_test[-train, train]
  V_11 <- V_train_test[train, train]

  V_11_inv <- solve(V_11 + p * diag(1, nrow = nrow(V_11)))
  # posterior mean
  beta_est <- solve(t(H) %*% V_11_inv %*% H) %*% t(H) %*% V_11_inv %*% (y_train)
  mu_1 <- H %*% beta_est
  mu_2 <- h_prime_T %*% beta_est
  gp_mean <- mu_2 + V_21 %*% V_11_inv %*% (y_train - mu_1)
  # posterior covariance
  gp_cov <- V_22 - V_21 %*% V_11_inv %*% V_12

  return(list(gp_mean = gp_mean, gp_cov = gp_cov))
}
