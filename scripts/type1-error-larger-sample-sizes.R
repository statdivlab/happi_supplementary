# Pauline Trinh, April 15, 2023
# Type 1 error rate control of happi-a at larger sample sizes 

df_list <- list.files(path="/Users/paulinetrinh/Documents/GitHub/pangenomics-hypothesis-testing/R/cluster/BRAIN_happi_sims/largesamplesize_sims", full.names = TRUE) %>% 
  lapply(readRDS)
my_list <- unlist(df_list, recursive = FALSE)
my_dataframe <- my_list %>% 
  do.call("rbind", .) %>%
  as_tibble

#saveRDS(my_data,"type1_pblrt_full_results_Mar7.RDS")
#my_data <- readRDS("type1_pblrt_full_results.RDS")
T1_stats  <- my_dataframe %>%
  as_tibble() %>% 
  mutate("iter" = 1:nrow(.)) %>%
  pivot_longer(5:7, values_to="pvalue", names_to="Method") %>%
  select(-sim) %>%
  group_by(xx_sd, nn, Method) %>%
  summarise("T1" = mean(pvalue<0.05),
            "n" = n()) %>%
  ungroup %>%
  arrange(Method) %>%
  mutate("lower" = 100*(T1 - 1.96*sqrt(T1 * (1-T1)/n)),
         "upper" = 100*(T1 + 1.96*sqrt(T1 * (1-T1)/n))) 


my_type1_plot_larger <- my_dataframe %>%
  mutate("iter" = 1:nrow(.)) %>%
  pivot_longer(5:7, values_to="pvalue", names_to="Method") %>%
  select(-sim) %>%
  filter(Method != "logreg_lrt") %>% 
  filter(Method != "happi_pbLRT") %>% 
  mutate(pvalue = as.numeric(pvalue)) %>% 
  group_by(xx_sd, nn, Method) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  mutate(theoretical = seq(1/530, 1, length.out = 530)) %>%
  mutate(n = as.factor(nn)) %>%
  mutate(xx_sd = recode_factor(xx_sd, `0.25` = "sigma[x] == 0.25", `0.5` = "sigma[x] == 0.5")) %>%
  ggplot(aes(x=theoretical, y = pvalue, color = Method, linetype = as.factor(n))) +
  facet_grid(~xx_sd, label = "label_parsed") +
  geom_line(size = 0.65) +
  geom_abline(intercept = 0, slope = 1) +
  ylab("Sample Quantile") + xlab("Theoretical Quantile") +
  theme_bw() +
  theme(legend.position="right",
        legend.text=element_text(size =10),
        plot.title = element_text(hjust = 0.5)) +
  labs(color = "", linetype = "Sample Size") +
  scale_color_manual(values = c("logreg_rao" = "orange",  "happi" = "darkblue"), labels = c("happi-a","GLM-Rao")) +
  # scale_color_manual(values = c("logreg_rao" = "orange", "happi_pbLRT" = "red", "happi_npLRT" = "darkblue"), labels = c("happi nonparametric LRT","happi parametric bootstrap LRT", "GLM-Rao")) +
  # scale_color_manual(values = c("orange", "darkblue", "seagreen"), labels = c("GLM-Rao","GLM-LRT","happi")) +
  NULL
dev.off()
my_type1_plot_larger
#saveRDS(my_type1_plot_larger,"my_type1_plot_larger.RDS")
ggsave("my_type1_plot_larger.png",
       height = 4, width = 8)
