# Pauline Trinh, April 15, 2023 
# Main Manuscript Type 1 Error Results
library(tidyverse)
df_list <- list.files(path="/Users/paulinetrinh/Documents/GitHub/pangenomics-hypothesis-testing/R/cluster/BRAIN_happi_sims/type1_results_mar7/n30n50", full.names = TRUE) %>% 
  lapply(readRDS)
my_list <- unlist(df_list, recursive = FALSE)
my_dataframe <- my_list %>% 
  do.call("rbind", .) %>%
  as_tibble %>% 
  select(-happi_pbLRT)

df_list2 <- list.files(path="/Users/paulinetrinh/Documents/GitHub/pangenomics-hypothesis-testing/R/cluster/BRAIN_happi_sims/type1_results_mar7/n100", full.names = TRUE) %>% 
  lapply(readRDS)
my_list2 <- unlist(df_list2, recursive = FALSE)
my_dataframe2 <- my_list2 %>% 
  do.call("rbind", .) %>%
  as_tibble

df_list3 <- list.files(path="/Users/paulinetrinh/Documents/GitHub/pangenomics-hypothesis-testing/R/cluster/BRAIN_happi_sims/type1_results_mar7/n100_part2", full.names = TRUE) %>% 
  lapply(readRDS)
my_list3 <- unlist(df_list3, recursive = FALSE)
my_dataframe3 <- my_list3 %>% 
  do.call("rbind", .) %>%
  as_tibble %>% 
  select(-happi_pbLRT)

#my_final_dataframe <- rbind(my_dataframe, my_dataframe2, my_dataframe3)

my_final_dataframe <- readRDS("data/type1_results_03302023.RDS")
T1_stats <- my_final_dataframe %>%
  as_tibble() %>% 
  mutate("iter" = 1:nrow(.)) %>%
  pivot_longer(5:8, values_to="pvalue", names_to="Method") %>%
  select(-sim) %>%
  group_by(xx_sd, nn, Method) %>%
  summarise("T1" = mean(pvalue<0.05),
            "n" = n()) %>%
  ungroup %>%
  arrange(Method) %>%
  mutate("lower" = 100*(T1 - 1.96*sqrt(T1 * (1-T1)/n)),
         "upper" = 100*(T1 + 1.96*sqrt(T1 * (1-T1)/n))) %>%
  mutate(Method_Names = ifelse(Method == "happi","happi-a",
                               ifelse(Method == "happi_npLRT","happi-np",Method))) %>%
  select(-Method) %>% 
  rename(Method = Method_Names)
#saveRDS(my_final_dataframe,"type1_results_03302023.RDS")

my_type1_plot <- my_final_dataframe %>%
  mutate("iter" = 1:nrow(.)) %>%
  pivot_longer(5:8, values_to="pvalue", names_to="Method") %>%
  select(-sim) %>%
  filter(Method != "logreg_lrt") %>% 
  filter(Method != "happi_pbLRT") %>% 
  mutate(pvalue = as.numeric(pvalue)) %>% 
  group_by(xx_sd, nn, Method) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  mutate(theoretical = seq(1/500, 1, length.out = 500)) %>%
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
  scale_color_manual(values = c("logreg_rao" = "orange",  "happi" = "red", "happi_npLRT" = "darkblue"), labels = c("happi-a","happi-np","GLM-Rao")) +
  NULL
dev.off()
saveRDS(my_type1_plot,"my_type1_plot.RDS")
my_type1_plot <- readRDS("my_type1_plot.RDS")
ggsave("proposed_limited_sizes_simulation_03132023.pdf",
       height = 3.5, width = 8)