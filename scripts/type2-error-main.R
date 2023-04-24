# Pauline Trinh, April 15, 2023 
# Power Simulations Results for Main Manuscript
library(tidyverse)

df_list_power <- list.files(path = "/Users/paulinetrinh/Documents/GitHub/pangenomics-hypothesis-testing/R/cluster/BRAIN_happi_sims/power_sims/n30_n50", full.names = TRUE) %>% 
  lapply(readRDS)
my_list_power <- unlist(df_list_power, recursive = FALSE)
my_dataframe_power <- my_list_power %>% 
  do.call("rbind", .) %>%
  as_tibble %>% 
  arrange(sim) %>% 
  mutate(nn = as.numeric(nn), 
         xx_sd = as.numeric(xx_sd), 
         beta = as.numeric(beta), 
         logreg_rao = as.numeric(logreg_rao), 
         happi = as.numeric(happi), 
         happi_npLRT = as.numeric(happi_npLRT)) 

df_list_power2 <- list.files(path = "/Users/paulinetrinh/Documents/GitHub/pangenomics-hypothesis-testing/R/cluster/BRAIN_happi_sims/power_sims/n30_n50_part2", full.names = TRUE) %>% 
  lapply(readRDS)
my_list_power2 <- unlist(df_list_power2, recursive = FALSE)
my_dataframe_power2 <- my_list_power2 %>% 
  do.call("rbind", .) %>%
  as_tibble %>% 
  arrange(sim) %>% 
  mutate(nn = as.numeric(nn), 
         xx_sd = as.numeric(xx_sd), 
         beta = as.numeric(beta), 
         logreg_rao = as.numeric(logreg_rao), 
         happi = as.numeric(happi), 
         happi_npLRT = as.numeric(happi_npLRT)) 

df_list_power3 <- list.files(path = "/Users/paulinetrinh/Documents/GitHub/pangenomics-hypothesis-testing/R/cluster/BRAIN_happi_sims/power_main_n100", full.names = TRUE) %>% 
  lapply(readRDS)
my_list_power3 <- unlist(df_list_power3, recursive = FALSE)
my_dataframe_power3 <- my_list_power3 %>% 
  do.call("rbind", .) %>%
  as_tibble  %>% 
  select(-logreg_lrt) %>% 
  arrange(sim) %>% 
  select(-sim) %>% 
  mutate(nn = as.numeric(nn), 
         xx_sd = as.numeric(xx_sd), 
         beta = as.numeric(beta), 
         logreg_rao = as.numeric(logreg_rao), 
         happi = as.numeric(happi), 
         happi_npLRT = as.numeric(happi_npLRT))

df_list_power4 <- list.files(path = "/Users/paulinetrinh/Documents/GitHub/pangenomics-hypothesis-testing/R/cluster/BRAIN_happi_sims/power_main_n100_part2", full.names = TRUE) %>% 
  lapply(readRDS)
my_list_power4 <- unlist(df_list_power4, recursive = FALSE)
my_dataframe_power4 <- my_list_power4 %>% 
  do.call("rbind", .) %>%
  as_tibble %>% arrange(sim) %>% 
  select(-sim) %>% 
  mutate(nn = as.numeric(nn), 
         xx_sd = as.numeric(xx_sd), 
         beta = as.numeric(beta), 
         logreg_rao = as.numeric(logreg_rao), 
         happi = as.numeric(happi), 
         happi_npLRT = as.numeric(happi_npLRT))

my_power_results <- my_dataframe_power %>% full_join(my_dataframe_power2) %>% full_join(my_dataframe_power3) %>% full_join(my_dataframe_power4)

power_summary <- my_power_results %>% 
  as_tibble() %>% 
  mutate("iter" = 1:nrow(.)) %>%
  pivot_longer(5:7, values_to="pvalue", names_to="Method") %>%
  select(-sim) %>%
  group_by(xx_sd, nn, Method,beta) %>%
  summarise("Power" = mean(pvalue<0.05),
            "n" = n()) %>%
  ungroup 

#saveRDS(my_power_results,"power_results_03302023.RDS")

type2_fig <- my_power_results %>%
  group_by(nn, xx_sd, beta) %>%
  summarise("happi-np" = mean(happi_npLRT < 0.05), "happi-a" = mean(happi < 0.05)) %>%
  ungroup %>%
  pivot_longer(4:5, names_to="Method", values_to="power") %>%
  filter(Method %in% c("happi-np","happi-a")) %>%
  filter(!(Method == "happi-a" & nn %in% c(30,50))) %>% 
  rename(n = nn) %>%
  mutate(xx_sd = recode_factor(xx_sd, `0.25` = "sigma[x] == 0.25", `0.5` = "sigma[x] == 0.5")) %>%
  ggplot(aes(x=as.numeric(beta), y = power, color = Method, linetype = as.factor(n))) +
  geom_point(aes(color = Method)) +
  geom_line(aes(linetype = as.factor(n), color = Method), size = 0.8) +
  geom_hline(yintercept=0.05,linetype=3, color = "black") +
  xlab(~paste("True value of ", beta[1])) +
  ylab("Power") +
  labs(color = "", linetype = "Sample Size") +
  theme_bw() +
  xlim(0.75, 3.25)  +
  theme(legend.position="bottom",
        legend.text=element_text(size =10),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(~xx_sd, label = "label_parsed") +
  scale_color_manual(values = c("happi-np"="darkblue","happi-a"="red"), labels = c("happi-a","happi-np")) +
  NULL
type2_fig

#saveRDS(type2_fig, "type2_fig_mar292023.RDS")

ggpubr::ggarrange(my_type1_plot, type2_fig, nrow = 1, common.legend=T, legend="right")
ggsave("simulation_03302023.pdf",
       height = 3.5, width = 12)


power <- my_dataframe_power %>% 
  as_tibble() %>% 
  mutate("iter" = 1:nrow(.)) %>%
  pivot_longer(5:7, values_to="pvalue", names_to="Method") %>%
  select(-sim) %>%
  group_by(xx_sd, nn, Method,beta) %>%
  summarise("Power" = mean(pvalue<0.05),
            "n" = n()) %>%
  ungroup 
