### A simulation to compare the accuracy of happi under different options for f
### Amy Willis & Pauline Trinh, April 25 2022
library(happi)
library(parallel)
library(tidyverse)

logit_data <- readRDS("data/logit_DRR102664.RDS")
gam_model <- mgcv::gam(presence ~ s(coverage), family = binomial, data = logit_data, method = "REML")

epsilon <- 0

compare_f <- function(sim) {
  results_out <- matrix(NA, nrow = niter, ncol = 7)
  colnames(results_out) <- c("sim", "nn", "xx_sd", "beta0", "beta1", "monotone", "spline")
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
                              max_iterations = 100,
                              firth=T, method="isotone",
                              nstarts = 1,
                              change_threshold = 0.1,
                              epsilon = 0)

      happi_spline <- happi(outcome = ys,
                            covariate = xx,
                            quality_var = mm,
                            max_iterations = 100,
                            firth=T, method="spline",
                            nstarts = 1,
                            change_threshold = 0.1,
                            epsilon = 0)

      results_out[bb, "monotone"] <- happi_monotone$beta[!is.na(happi_monotone$beta[,2]),2] %>% tail(1)
      results_out[bb, "spline"] <- happi_spline$beta[!is.na(happi_spline$beta[,2]),2] %>% tail(1)

    })
  }
  results_out
}


####################
#### run sims ######
####################
ns <- c(30, 50, 100)
xx_sds <- c(0.25, 0.5)
beta0 <- 0
beta1 <- c(0.5, 1, 2)
combos <- crossing(xx_sds, ns, beta0, beta1)
combos <- bind_cols("sim" = 1:nrow(combos), combos)
combos

niter <- 250
set.seed(104)
compare_f_run1 <- mclapply(18:1, compare_f, mc.cores = 6)

#saveRDS(errors, "sims_results/accuracy_f_nonparametric.RDS")

# To reproduce this analysis from the happi_supplementary repo
# remove the hash to read in dataset for object errors

# errors <- readRDS("data/accuracy_f_nonparametric.RDS")

####################
#### analysis ######
####################

errors <- compare_f_run1 %>%
  do.call(rbind, .) %>%
  as_tibble %>%
  pivot_longer(6:7) %>%
  mutate(error = (value - beta1)^2)

errors %>%
  group_by(name) %>%
  summarise(mean(error), median_error = median(error)) %>%
  arrange(median_error)
## spline is best

errors %>%
  group_by(nn, xx_sd, name, beta1) %>%
  summarise(mean(error), median_error = median(error)) %>%
  group_by(nn, xx_sd, beta1) %>%
  arrange(median_error, .by_group=T) %>% print(n=Inf)

errors %>%
  filter(name %in% c("monotone", "spline")) %>%
  group_by(nn, xx_sd, name, beta1) %>%
  summarise(median_error = median(error)) %>%
  group_by(nn, xx_sd, beta1) %>%
  arrange(median_error, .by_group=T) %>%
  pivot_wider(names_from=name, values_from=median_error) %>%
  mutate(spline_winner = spline < monotone) %>%
  ungroup %>%
  summarise(sum(spline_winner), n(), 100*mean((monotone-spline)/spline))
# `sum(spline_winner)` `n()` `100 * mean((monotone - spline)/spline)`
# <int> <int>                                    <dbl>
#   1                   17    18                                     53.9

errors %>%
  filter(name %in% c("monotone", "spline")) %>%
  group_by(nn, xx_sd, name, beta1) %>%
  summarise(mean_error = mean(error)) %>%
  group_by(nn, xx_sd, beta1) %>%
  arrange(mean_error, .by_group=T) %>%
  pivot_wider(names_from=name, values_from=mean_error) %>%
  mutate(spline_winner = spline < monotone) %>%
  ungroup %>%
  summarise(sum(spline_winner), n(), 100*mean((monotone-spline)/spline))
# `sum(spline_winner)` `n()` `100 * mean((monotone - spline)/spline)`
# <int> <int>                                    <dbl>
#   1                   18    18                                     50.8

