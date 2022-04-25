### A simulation to compare the Type 2 error rates of happi
### Amy Willis, April 20 2022
library(happi)
library(parallel)
library(tidyverse)

logit_data <- readRDS("data/logit_DRR102664.RDS")
gam_model <- mgcv::gam(presence ~ s(coverage), family = binomial, data = logit_data, method = "REML")

epsilon <- 0

run_power_splines <- function(sim) {
  results_out <- matrix(NA, nrow = niter, ncol = 6)
  colnames(results_out) <- c("sim", "nn", "xx_sd", "beta", "logreg", "happi")
  results_out[, "sim"] <- sim
  results_out[, "nn"] <- combos$ns[sim]
  results_out[, "xx_sd"] <- combos$xx_sds[sim]
  results_out[, "beta"] <- combos$beta1[sim]

  for (bb in 1:niter) {

    nn <- combos$ns[sim]
    xx_sd <- combos$xx_sds[sim]
    # beta <- c(-1, combos$beta1[sim])
    beta <- c(0, combos$beta1[sim])

    mm <- seq(10, 40, length.out = nn)
    true_f <- predict(gam_model, type = "response", newdata = data.frame("coverage" = mm))

    x1 <- rnorm(nn, seq(0, 1, length.out = nn), sd=xx_sd)
    xx <- cbind(1, x1)

    true_lambda_probs <- expit(c(xx %*% beta))

    lambdas <- rbinom(n = nn, size = 1, prob = true_lambda_probs)
    true_prob_y_equal_1 <- true_f # for lambdas == 1
    true_prob_y_equal_1[lambdas == 0] <- epsilon # for lambdas == 0
    ys <- rbinom(n = nn, size = 1, prob = true_prob_y_equal_1)

    ## logistic reg
    results_out[bb, "logreg"] <- anova(glm(ys ~ xx - 1, family="binomial"),
                                       glm(ys ~ 1, family="binomial"),
                                       test = "Rao")[2, "Pr(>Chi)"]

    happi_out <- happi(outcome = ys,
          covariate = xx,
          quality_var = mm,
          max_iterations = 50,
          method="splines", firth=T, spline_df=4,
          nstarts = 1,
          change_threshold = 0.1,
          epsilon = 0)

    results_out[bb, "happi"] <- happi_out$loglik$pvalue[!is.na(happi_out$loglik$pvalue)] %>% tail(1)
  }
  results_out
}


####################
#### 1st half ######
####################
ns <- c(30, 50, 100)
xx_sds <- c(0.25, 0.5)
beta1 <- seq(from = 0.5, to = 3, by = 0.5)
combos <- crossing(xx_sds, ns, beta1)
combos <- bind_cols("sim" = 1:nrow(combos), combos)
combos
niter <- 2
run_power_splines(1)
run_power_splines(19)
run_power_splines(36)
niter <- 250
set.seed(200)
type2_splines_1 <- mclapply(36:1, run_power_splines, mc.cores = 6)

####################
#### 2nd half ######
####################
set.seed(201)
niter <- 250
type2_splines_2 <- mclapply(36:1, run_power_splines, mc.cores = 6)

type2_splines <- c(type2_splines_1, type2_splines_2) %>%
  do.call(rbind, .) %>%
  as_tibble



####################
#### more ######
####################
set.seed(202)
niter <- 250
type2_splines_3 <- mclapply(36:1, run_power_splines, mc.cores = 6)

set.seed(203)
niter <- 250
type2_splines_4 <- mclapply(36:1, run_power_splines, mc.cores = 6)

type2_splines <- c(type2_splines_1, type2_splines_2, type2_splines_3, type2_splines_4) %>%
  do.call(rbind, .) %>%
  as_tibble

# saveRDS(type2_splines, "type2_results.RDS")
