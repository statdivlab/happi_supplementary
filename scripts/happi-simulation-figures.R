library(ggpubr)

#### f plot

logit_data <- readRDS("data/logit_DRR102664.RDS")
gam_model <- mgcv::gam(presence ~ s(coverage), family = binomial, data = logit_data, method = "REML")

f_used <- logit_data %>%
  filter(ID == "DRR102664") %>%
  group_by(coverage) %>%
  summarise(prop = mean(presence)) %>%
  ggplot(aes(y = prop, x = coverage)) +
  geom_point(cex=0.7) +
  geom_line(aes(y = gam), data = tibble("coverage" = seq(0, 30, length.out = 1000)) %>%
              mutate("gam" = predict(gam_model, type = "response", newdata = .))) +
  theme_bw() +
  theme(legend.text=element_text(size=5),
        legend.position = "none",
        legend.background = element_blank()) +
  ylab("Proportion of subsamples containing gene") +
  xlab("Subsampled Sequencing Coverage")  +
  NULL
# ggsave("../../research/pauline/pangenomics-hypothesis-testing/manuscript/figures/f_DRR102664.pdf",
#        height = 3.5, width = 5)


### type 1
type1_all <- readRDS("sims_results/type1_all.RDS")
type1_fig <- type1_all %>%
  mutate("iter" = 1:nrow(.)) %>%
  pivot_longer(5:7, values_to="pvalue", names_to="Method") %>%
  select(-sim) %>%
  group_by(xx_sd, nn, Method) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  mutate(theoretical = seq(1/500, 1, length.out = 500)) %>%
  mutate(n = as.factor(nn)) %>%
  mutate(xx_sd = recode_factor(xx_sd, `0.25` = "sigma[x] == 0.25", `0.5` = "sigma[x] == 0.5")) %>%
  filter(Method != "logreg_lrt") %>%
  ggplot(aes(x=theoretical, y = pvalue, color = Method, linetype = as.factor(n))) +
  facet_grid(~xx_sd, label = "label_parsed") +
  geom_line(size = 0.65) +
  geom_abline(intercept = 0, slope = 1) +
  ylab("Sample Quantile") + xlab("Theoretical Quantile") +
  theme_bw() +
  theme(legend.position="bottom",
        legend.text=element_text(size =10),
        plot.title = element_text(hjust = 0.5)) +
  labs(color = "", linetype = "Sample Size") +
  scale_color_manual(values = c("logreg_rao" = "orange", "happi" = "darkblue"), labels = c("GLM-Rao","happi")) +
  # scale_color_manual(values = c("orange", "darkblue", "seagreen"), labels = c("GLM-Rao","GLM-LRT","happi")) +
  NULL
type1_fig

### type 2

t2 <- readRDS("sims_results/type2_all.RDS")
type2_fig <- t2 %>%
  group_by(nn, xx_sd, beta) %>%
  summarise("GLM-Rao" = mean(logreg < 0.05), "happi" = mean(happi < 0.05)) %>%
  ungroup %>%
  pivot_longer(4:5, names_to="Method", values_to="power") %>%
  filter(Method == "happi") %>% # GLM doesn't control Type 1 error
  rename(n = nn) %>%
  mutate(xx_sd = recode_factor(xx_sd, `0.25` = "sigma[x] == 0.25", `0.5` = "sigma[x] == 0.5")) %>%
  ggplot(aes(x=beta, y = power, color = Method, linetype = as.factor(n))) +
  geom_point(aes(color = Method)) +
  geom_line(aes(linetype = as.factor(n), color = Method), size = 0.8) +
  geom_hline(yintercept=0.05,linetype=3, color = "black") +
  xlab(~paste("True value of ", beta[1])) +
  ylab("Power") +
  labs(color = "", linetype = "Sample Size") +
  theme_bw() +
  xlim(0.75, 3.25)+
  theme(legend.position="bottom",
        legend.text=element_text(size =10),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(~xx_sd, label = "label_parsed") +
  scale_color_manual(values = c("darkblue"), labels = c("happi")) +
  NULL
type2_fig

ggarrange(type1_fig, type2_fig, nrow = 1, common.legend=T, legend="right")
# ggsave("../../research/pauline/pangenomics-hypothesis-testing/manuscript/figures/simulation_04192022.pdf",
#        height = 3.5, width = 12)
