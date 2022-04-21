### A simulation to compare the accuracy of happi under different options for f
### Amy Willis & Pauline Trinh, April 20 2022
library(happi)
library(parallel)
library(tidyverse)

logit_data <- readRDS("data/logit_DRR102664.RDS")
gam_model <- mgcv::gam(presence ~ s(coverage), family = binomial, data = logit_data, method = "REML")

epsilon <- 0

univar_likelihood <- function(par, my_data){
  e = 0 # e is parameter argument 1
  vBeta = c(par[1],par[2]) # my beta vector is parameter arguments 2 and 3
  vGamma = par[3] # my gamma vector is parameter arguments 4
  Y0 <- my_data %>% dplyr::filter(Y == 0)
  Y1 <- my_data %>% dplyr::filter(Y == 1)

  Mi0 = as.matrix(Y0['coverage'])
  mX0 = as.matrix(Y0['X1'])
  mX0 = cbind(1, mX0)

  Mi1 = as.matrix(Y1['coverage'])
  mX1 = as.matrix(Y1['X1'])
  mX1 = cbind(1, mX1)

  LL <- sum(log((1-e)*exp(-(mX0 %*% vBeta))/(1+exp(-(mX0 %*% vBeta))) + # when Y = 0
                  (1-(1-exp(-(Mi0*vGamma))))*(1/(1+exp(-(mX0 %*% vBeta)))))) +
    sum(log(e*exp(-(mX1 %*% vBeta))/(1+exp(-(mX1 %*% vBeta))) +  # when Y = 1
              ((1-exp(-(Mi1*vGamma)))/(1+exp(-(mX1 %*% vBeta))))))
  return(LL)
}
quality_function_logit <- function(quality_vector, gamma_vector){
  quality_function_value <- exp(quality_vector %*% gamma_vector)/(1+exp(quality_vector %*% gamma_vector))
  return(quality_function_value)
}
expit_likelihood <- function(par, my_data){
  e = 0 # e is parameter argument 1
  vBeta = c(par[1],par[2]) # my beta vector is parameter arguments 2 and 3
  vGamma = c(par[3],par[4]) # my gamma vector is parameter arguments 4
  Y0 <- my_data %>% dplyr::filter(Y == 0)
  Y1 <- my_data %>% dplyr::filter(Y == 1)

  Mi0 = as.matrix(Y0['coverage']) # I'm going to model the regular coverage
  Mi0 = cbind(1, Mi0)
  mX0 = as.matrix(Y0['X1'])
  mX0 = cbind(1, mX0)

  Mi1 = as.matrix(Y1['coverage'])
  Mi1 = cbind(1, Mi1)
  mX1 = as.matrix(Y1['X1'])
  mX1 = cbind(1, mX1)

  LL <- sum(log((1-e)*exp(-(mX0 %*% vBeta))/(1+exp(-(mX0 %*% vBeta))) + # when Y = 0
                  (1-(quality_function_logit(quality_vector = Mi0, gamma_vector = vGamma)))*(1/(1+exp(-(mX0 %*% vBeta)))))) +
    sum(log(e*exp(-(mX1 %*% vBeta))/(1+exp(-(mX1 %*% vBeta))) +  # when Y = 1
              (1/(1+exp(-(mX1 %*% vBeta))))*(quality_function_logit(quality_vector = Mi1, gamma_vector = vGamma))))
  return(LL)
}

compare_f <- function(sim) {
  results_out <- matrix(NA, nrow = niter, ncol = 8)
  colnames(results_out) <- c("sim", "nn", "xx_sd", "beta0", "beta1", "monotone", "expit", "univar")
  results_out[, "sim"] <- sim
  results_out[, "nn"] <- combos$ns[sim]
  results_out[, "xx_sd"] <- combos$xx_sds[sim]
  results_out[, "beta0"] <- combos$beta0[sim]
  results_out[, "beta1"] <- combos$beta1[sim]

  for (bb in 1:niter) {


    nn <- combos$ns[sim]
    xx_sd <- combos$xx_sds[sim]
    beta <- c(combos$beta0[sim], combos$beta1[sim])

    mm <- seq(10, 40, length.out = nn)
    true_f <- predict(gam_model, type = "response", newdata = data.frame("coverage" = mm))
    x1 <- rnorm(nn, seq(0, 1, length.out = nn), sd=xx_sd)
    xx <- cbind(1, x1)
    true_lambda_probs <- expit(c(xx %*% beta))
    lambdas <- rbinom(n = nn, size = 1, prob = true_lambda_probs)
    true_prob_y_equal_1 <- true_f # for lambdas == 1
    true_prob_y_equal_1[lambdas == 0] <- epsilon # for lambdas == 0
    ys <- rbinom(n = nn, size = 1, prob = true_prob_y_equal_1)

    try({
      happi_monotone <- happi(outcome = ys,
                              covariate = xx,
                              quality_var = mm,
                              max_iterations = 50,
                              nstarts = 1,
                              change_threshold = 0.1,
                              epsilon = 0)

      draw <- tibble(Y = ys, X1 = x1, coverage = mm)
      happi_univar <- optimx::optimx(c(0, 0, 20),
                                     univar_likelihood,
                                     my_data = draw,
                                     method = 'L-BFGS-B',
                                     lower= c(-Inf,-Inf, 0),
                                     upper=c(Inf, Inf, Inf),
                                     control = list(trace = 0, fnscale= -1))

      univar_likelihood(c(-1, 2, 0.001), draw)
      happi_expit <- optimx::optimx(c(0, 0, 0, 0),
                                    expit_likelihood,
                                    my_data = draw,
                                    method = 'BFGS',
                                    control = list(trace = 0, fnscale= -1))

      results_out[bb, "monotone"] <- happi_monotone$beta[!is.na(happi_monotone$beta[,2]),2] %>% tail(1)
      results_out[bb, "univar"] <- happi_univar[1,2]
      results_out[bb, "expit"] <- happi_expit[1, 2]

    })
  }
  results_out
}


####################
#### 1st run ######
####################
ns <- c(30, 50)
xx_sds <- c(0.25, 0.5)
beta0 <- c(0)
beta1 <- c(0.5, 1, 2)
combos <- crossing(xx_sds, ns, beta0, beta1)
combos <- bind_cols("sim" = 1:nrow(combos), combos)
combos

####################
#### 2nd run #######
####################
niter <- 50
set.seed(104)
compare_f_run1 <- mclapply(1:12, compare_f, mc.cores = 6)
niter <- 150
set.seed(105)
compare_f_run2 <- mclapply(1:12, compare_f, mc.cores = 6)

errors <- c(compare_f_run1, compare_f_run2) %>%
  do.call(rbind, .) %>%
  as_tibble %>%
  pivot_longer(6:8) %>%
  mutate(error = (value - beta1)^2)

errors %>%
  group_by(nn, xx_sd, name, beta1) %>%
  summarise(mean(error), median(error))
## monotone is best

# relative improvements
errors %>%
  group_by(nn, xx_sd, beta1, name) %>%
  summarise(med = median(error)) %>%
  pivot_wider(names_from="name", values_from="med") %>%
  mutate(100*(expit - monotone)/expit,
         100*(univar - monotone)/univar)


errors %>% ggplot(aes(x = as.character(nn), col = name, y = error)) +
  facet_grid(xx_sd ~ beta1) +
  geom_boxplot() +
  scale_y_log10()
