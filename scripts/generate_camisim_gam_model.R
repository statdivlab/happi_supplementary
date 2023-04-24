# Pauline Trinh, Apr 15, 2023 
# Code for taking CAMISIM data and creating GAM model for supplementary simulation study 
library(tidyverse)
coverage_sweep <- read_tsv("/Users/paulinetrinh/Documents/GitHub/CAMISIM-DATA/coverage_sweep/Coverage_Sweep_results.tsv")

genome_fract <- coverage_sweep %>% 
  filter(metric == "Genome_fraction_(%)") %>%
  filter(Assembler == "MEGAHIT_1.0.3") %>%
  filter(`Error rate` %in% c("2%","ART CAMI")) %>% 
  arrange(Assembler,`Error rate`,Coverage) 

genome_fract_error2 <- coverage_sweep %>% 
  filter(metric == "Genome_fraction_(%)") %>%
  filter(Assembler == "MEGAHIT_1.0.3") %>%
  filter(`Error rate` == "2%") %>% 
  arrange(Assembler,`Error rate`,Coverage) 


set.seed(8)
N = 100000
cov2 <- cbind(2,rbinom(N,1,0.20793))
cov4 <- cbind(4,rbinom(N,1,0.68923))
cov8 <- cbind(8,rbinom(N,1,0.96961))
cov16 <- cbind(16,rbinom(N,1,0.98234))
cov32 <- cbind(32,rbinom(N,1,0.98239))
cov64 <- cbind(64,rbinom(N,1,0.98014))
cov128 <- cbind(128,rbinom(N,1,0.96845))
cov256 <- cbind(256,rbinom(N,1,0.53845))
cov512 <- cbind(512,rbinom(N,1,0.06056))

my_camisim_data <- rbind(cov2,cov4,cov8,cov16,cov32,cov64,cov128,cov256,cov512) %>% 
  as_tibble() %>%
  rename(coverage = V1, 
         presence = V2)

my_camisim_data %>% 
  group_by(coverage) %>%
  summarise(sum(presence)/100000)

saveRDS(my_camisim_data,"my_camisim_data.RDS")
library(mgcv)

gam_model <- mgcv::gam(presence ~ s(coverage, k = 9), 
                       family = binomial, 
                       data = my_camisim_data, method = "REML")
plot_smooths(model = gam_model, 
             series = coverage)
#saveRDS(gam_model,"camisim_gam_model.RDS")

## Generate Supplementary Figure 4 
expit <- function(input) {
  exp(input)/ (1+exp(input))
} 


plot_smooths(model = gam_model, 
             series = coverage, 
             transform = expit) + 
  theme_bw() + 
  ylab("Probability of gene detection") + 
  xlab("Coverage")
ggsave("camisim_model_curve.png",
       height = 3.5, width = 8)