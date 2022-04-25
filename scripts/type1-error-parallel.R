### A simulation to compare the Type 1 error rates of GLM-Rao and happi
### Amy Willis, April 20 2022
library(happi)
library(parallel)
library(tidyverse)

logit_data <- readRDS("data/logit_DRR102664.RDS")
gam_model <- mgcv::gam(presence ~ s(coverage), family = binomial, data = logit_data, method = "REML")
epsilon <- 0

run_type1_splines <- function(sim) {
  results_out <- matrix(NA, nrow = niter, ncol = 7)
  colnames(results_out) <- c("sim", "nn", "xx_sd", "beta", "logreg_rao", "logreg_lrt", "happi")
  results_out[, "sim"] <- sim
  results_out[, "nn"] <- combos$ns[sim]
  results_out[, "xx_sd"] <- combos$xx_sds[sim]
  results_out[, "beta"] <- combos$beta1[sim]

  for (bb in 1:niter) {

    nn <- combos$ns[sim]
    xx_sd <- combos$xx_sds[sim]
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
    results_out[bb, "logreg_rao"] <- anova(glm(ys ~ xx - 1, family="binomial"),
                                           glm(ys ~ 1, family="binomial"),
                                           test = "Rao")[2, "Pr(>Chi)"]


    results_out[bb, "logreg_lrt"] <- anova(glm(ys ~ xx - 1, family="binomial"),
                                           glm(ys ~ 1, family="binomial"),
                                           test = "LRT")[2, "Pr(>Chi)"]

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
beta1 <- 0
combos <- crossing(xx_sds, ns, beta1)
combos <- bind_cols("sim" = 1:nrow(combos), combos)
combos
niter <- 200
set.seed(101)
type1_splines <- mclapply(1:6, run_type1_splines, mc.cores = 6)

####################
#### 2nd half ######
####################
niter <- 300
set.seed(103)
type1_splines_part_b <- mclapply(1:6, run_type1_splines, mc.cores = 6)



type1_splines_df <- c(type1_splines, type1_splines_part_b) %>%
  do.call(rbind, .) %>%
  as_tibble

# saveRDS(type1_splines_df, "type1_results.RDS")

####################
#### Statistics ####
####################
any(is.na(type1_all$happi))
mean((type1_all$logreg_rao - type1_all$logreg_lrt)^2)
cor(type1_all$logreg_rao, type1_all$logreg_lrt)

type1_splines_df %>%
  mutate("iter" = 1:nrow(.)) %>%
  pivot_longer(5:7, values_to="pvalue", names_to="Method") %>%
  select(-sim) %>%
  filter(Method != "logreg_lrt") %>%
  group_by(xx_sd, nn, Method) %>%
  summarise("T1" = mean(pvalue<0.05),
            "n" = n()) %>%
  ungroup %>%
  arrange(T1) %>%
  mutate("lower" = 100*(T1 - 1.96*sqrt(T1 * (1-T1)/n)),
         "upper" = 100*(T1 + 1.96*sqrt(T1 * (1-T1)/n)))

type1_splines_df %>%
  mutate("iter" = 1:nrow(.)) %>%
  pivot_longer(1:6, values_to="pvalue") %>%
  mutate(sim = as.numeric(str_remove(name, "sim"))) %>%
  full_join(combos) %>%
  select(-seed, -name, -sim) %>%
  group_by(xx_sds, ns) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  mutate(theoretical = seq(1/niter, 1, length.out = niter),
         Method = "GLM-Rao")
