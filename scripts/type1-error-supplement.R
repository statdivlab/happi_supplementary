# Pauline Trinh, April 15, 2023 
# CAMISIM Type 1 error rate results and Figure 
library(tidyverse)
library(ggplot2)

df_list_mm100 <- list.files(path="/Users/paulinetrinh/Documents/GitHub/CAMISIM-DATA/type1_apr8/mm100", full.names = TRUE) %>% 
  lapply(readRDS)
df_list_mm200 <- list.files(path="/Users/paulinetrinh/Documents/GitHub/CAMISIM-DATA/type1_apr8/mm200", full.names = TRUE) %>% 
  lapply(readRDS)
df_list_mm300 <- list.files(path="/Users/paulinetrinh/Documents/GitHub/CAMISIM-DATA/type1_apr8/mm300", full.names = TRUE) %>% 
  lapply(readRDS)

my_list_mm100 <- unlist(df_list_mm100, recursive = FALSE)
my_dataframe_mm100 <- my_list_mm100 %>% 
  do.call("rbind", .) %>%
  as_tibble %>% 
  mutate(max_mm = 100)

my_list_mm200 <- unlist(df_list_mm200, recursive = FALSE)
my_dataframe_mm200 <- my_list_mm200 %>% 
  do.call("rbind", .) %>%
  as_tibble %>% 
  mutate(max_mm = 200)

my_list_mm300 <- unlist(df_list_mm300, recursive = FALSE)
my_dataframe_mm300 <- my_list_mm300 %>% 
  do.call("rbind", .) %>%
  as_tibble %>% 
  mutate(max_mm = 300)

my_dataframe <- my_dataframe_mm100 %>% full_join(my_dataframe_mm200) %>% 
  full_join(my_dataframe_mm300)
saveRDS(my_dataframe, "type1_results_camisim.RDS")

camisim_type1 <- my_dataframe %>% 
  as_tibble() %>% 
  mutate("iter" = 1:nrow(.)) %>%
  pivot_longer(5:8, values_to="pvalue", names_to="Method") %>%
  select(-c(sim,iter)) %>%
  group_by(xx_sd, nn, Method,max_mm) %>%
  summarise("T1" = mean(pvalue<0.05),
            "n" = n()) %>%
  ungroup %>%
  arrange(max_mm,Method) %>%
  mutate("lower" = 100*(T1 - 1.96*sqrt(T1 * (1-T1)/n)),
         "upper" = 100*(T1 + 1.96*sqrt(T1 * (1-T1)/n))) %>%
  mutate(Method_Names = ifelse(Method == "happi","happi-a",
                               ifelse(Method == "happi_npLRT","happi-np",Method))) %>%
  select(-Method) %>% 
  rename(Method = Method_Names) 
write.csv(camisim_type1,"camisim_type1.csv")
saveRDS(my_dataframe, "my_dataframe_camisim.RDS")

camisim_type1 %>% group_by(Method) %>% summarise(mean = mean(T1)) 
mean((my_dataframe$logreg_rao - my_dataframe$logreg_lrt)^2)
cor(my_dataframe$logreg_rao, my_dataframe$logreg_lrt)

my_type1_plot <- my_dataframe %>%
  mutate("iter" = 1:nrow(.)) %>%
  pivot_longer(5:8, values_to="pvalue", names_to="Method") %>%
  select(-c(sim,iter)) %>%
  filter(Method != "logreg_lrt") %>%
  mutate(pvalue = as.numeric(pvalue)) %>% 
  group_by(xx_sd, nn, Method, max_mm) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  mutate(theoretical = seq(1/500, 1, length.out = 500)) %>%
  mutate(n = as.factor(nn)) %>%
  mutate(xx_sd = recode_factor(xx_sd, `0.25` = "sigma[x] == 0.25", `0.5` = "sigma[x] == 0.5")) %>%
  mutate(max_mm = recode_factor(max_mm, `100` = "coverage[max] == 100", `200` = "coverage[max] == 200",`300` = "coverage[max] == 300")) %>%
  ggplot(aes(x=theoretical, y = pvalue, color = Method, linetype = as.factor(n))) +
  facet_grid(cols = vars(xx_sd), label = "label_parsed", 
             rows = vars(max_mm)) +
  geom_line(linewidth = 0.65) +
  geom_abline(intercept = 0, slope = 1) +
  ylab("Sample Quantile") + xlab("Theoretical Quantile") +
  theme_bw() +
  theme(legend.position="right",
        legend.text=element_text(size =10),
        plot.title = element_text(hjust = 0.5)) +
  labs(color = "", linetype = "Sample Size") +
  scale_color_manual(values = c("logreg_rao" = "orange",  "happi" = "red", "happi_npLRT" = "darkblue"), labels = c("happi-a","happi-np","GLM-Rao")) +
  # scale_color_manual(values = c("logreg_rao" = "orange", "happi_pbLRT" = "red", "happi_npLRT" = "darkblue"), labels = c("happi nonparametric LRT","happi parametric bootstrap LRT", "GLM-Rao")) +
  # scale_color_manual(values = c("orange", "darkblue", "seagreen"), labels = c("GLM-Rao","GLM-LRT","happi")) +
  NULL
#dev.off()

saveRDS(my_type1_plot,"my_type1_plot_camisim_apr12.RDS")
ggsave("type1_camisim_results_apr12.pdf",
       height = 8, width = 10)
