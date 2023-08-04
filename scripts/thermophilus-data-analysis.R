# Pauline Trinh, April 15, 2023 
# Thermophilus Data Analysis 

# Data Analysis of thermophillus MAGs from Sweden and Spain 
library(tidyverse)
library(happi)
library(parallel)
library(dplyr)
thermophillus <- read_tsv("/Users/paulinetrinh/Downloads/gene_presence_absence-2.Rtab")
mgnify_metadata <- read_tsv("/Users/paulinetrinh/Downloads/genomes-all_metadata.tsv")

thermophillus_ids <- colnames(thermophillus) %>% 
  as_tibble() %>% 
  filter(value != "Gene") %>% 
  rename(Genome = value)

dim(thermophillus_ids)

my_thermophillus_metadata <- thermophillus_ids %>% left_join(mgnify_metadata, by = "Genome")
dim(my_thermophillus_metadata)
thermophilus_core <- read.table("/Users/paulinetrinh/Downloads/core_genes-thermophilus.txt", header = FALSE) %>%
  mutate(Core = 1) %>% 
  rename(Gene = V1)
# 813   2

Core_genes_list <- thermophilus_core %>% full_join(thermophillus, by = "Gene")

my_thermophillus_metadata %>%
  group_by(Country) %>% 
  summarise(n = n())
# Figuring out which ones are the single-copy core genes 
# n = 157 MAGs of Streptococcus thermophilus used in fermentation 
# where 82 are from Spain and 75 are from Sweden 
my_thermophillus_data <- my_thermophillus_metadata %>% 
  filter(Country %in% c("Spain","Sweden"))
transpose_thermophillus <- thermophillus %>% column_to_rownames("Gene") %>% t() %>%
  as.data.frame() %>% 
  rownames_to_column("Genome")
min(my_thermophillus_metadata$Completeness)
max(my_thermophillus_metadata$Completeness)

min(my_thermophillus_metadata$Contamination)
max(my_thermophillus_metadata$Contamination)

sums_of_things <- colSums(transpose_thermophillus[2:3417]) %>% as_tibble() %>% 
  arrange(value)


thermophillus_data <- my_thermophillus_data %>% 
  left_join(transpose_thermophillus, by = "Genome") 
sums_of_things <- colSums(thermophillus_data[21:3437]) %>% as_tibble() %>% 
  arrange(desc(value)) #618 columns that sum to 0 


colMax <- sapply(thermophillus_data, max, na.rm=TRUE)
colMin <- sapply(thermophillus_data, min, na.rm=TRUE)                                                             
nr_thermophillus_data <- thermophillus_data[, colMin != 0 | colMax != 0]
dim(nr_thermophillus_data)
sums_of_things <- colSums(nr_thermophillus_data[21:2819]) %>% as_tibble() %>% 
  arrange(desc(value)) # got rid of all columns that are all 1 or 0 

#157 2819

saveRDS(nr_thermophillus_data,"thermophilus_data.RDS")
nr_thermophillus_data <- readRDS("data/thermophilus_data.RDS")
thermophilus_df <- nr_thermophillus_data %>% 
  select(Completeness,Country,`pbpX~~~pbpX_1~~~pbpX_2`:ncol(.) ) %>% 
  mutate(Spain = ifelse(Country == "Spain", 1,0)) %>% 
  select(Completeness, Spain, `pbpX~~~pbpX_1~~~pbpX_2`:`smc_2`)

#saveRDS(thermophilus_df,"thermophilus_df.RDS")
thermophilus_df <- readRDS("data/thermophilus_df.RDS")
thermophilus_df <- thermophilus_df %>% 
  mutate(log_completeness = log(Completeness))
x_matrix <- model.matrix(~Spain, data = thermophilus_df)

run_happi_thermophilus <- function(colnum) {
  happi_results <- happi::happi(outcome=unlist(thermophilus_df[,colnum]), 
                                     covariate=x_matrix, 
                                     quality_var=thermophilus_df$log_completeness,
                                     method="splines", 
                                     firth=T, 
                                     spline_df=4,
                                     max_iterations=1000, 
                                     change_threshold=0.01, 
                                     epsilon=0)
  
  return(happi_results)
}

run_happi_thermophilus_e01 <- function(colnum) {
  happi_results <- happi::happi(outcome=unlist(thermophilus_df[,colnum]), 
                                     covariate=x_matrix, 
                                     quality_var=thermophilus_df$log_completeness,
                                     method="splines", 
                                     firth=T, 
                                     spline_df=4,
                                     max_iterations=1000, 
                                     change_threshold=0.01, 
                                     epsilon=0.01)
  return(happi_results)
}
run_happi_thermophilus_e05 <- function(colnum) {
  happi_results <- happi::happi(outcome=unlist(thermophilus_df[,colnum]), 
                                     covariate=x_matrix, 
                                     quality_var=thermophilus_df$log_completeness,
                                     method="splines", 
                                     firth=T, 
                                     spline_df=4,
                                     max_iterations=1000, 
                                     change_threshold=0.01, 
                                     epsilon=0.05)
  
  return(happi_results)
}

run_happi_thermophilus_e1 <- function(colnum) {
  happi_results <- happi::happi(outcome=unlist(thermophilus_df[,colnum]), 
                                     covariate=x_matrix, 
                                     quality_var=thermophilus_df$log_completeness,
                                     method="splines", 
                                     firth=T, 
                                     spline_df=4,
                                     max_iterations=1000, 
                                     change_threshold=0.01, 
                                     epsilon=0.1)
  
  return(happi_results)
}

set.seed(18)
thermophilus_results_part1 <- mclapply(3:1708, run_happi_thermophilus, mc.cores=6)
#saveRDS(thermophilus_results_part1,"thermophilus_results_part1_mar23.RDS")
thermophilus_results_part2 <- mclapply(1709:2801, run_happi_thermophilus, mc.cores=6)
#saveRDS(thermophilus_results_part2,"thermophilus_results_part2_mar23.RDS")

set.seed(10)
### e 0.01 
thermophilus_results_sensitivity_e01_part1 <- mclapply(3:1708, run_happi_thermophilus_e01, mc.cores=6)
saveRDS(thermophilus_results_sensitivity_e01_part1,"thermophilus_results_sensitivity_e01_part1_mar23.RDS")
thermophilus_results_pensitivity_e01_part2 <- mclapply(1709:2801, run_happi_thermophilus_e01, mc.cores=6)
saveRDS(thermophilus_results_pensitivity_e01_part2,"thermophilus_results_sensitivity_e01_part2_mar23.RDS")
set.seed(11)
### e 0.05 
thermophilus_results_sensitivity_e05_part1 <- mclapply(3:1708, run_happi_thermophilus_e05, mc.cores=6)
saveRDS(thermophilus_results_sensitivity_e05_part1,"thermophilus_results_sensitivity_e05_part1.RDS")
thermophilus_results_sensitivity_e05_part2 <- mclapply(1709:2801, run_happi_thermophilus_e05, mc.cores=6)
saveRDS(thermophilus_results_sensitivity_e05_part2,"thermophilus_results_sensitivity_e05_part2.RDS")
set.seed(12)
### e 0.1
thermophilus_results_sensitivity_e1_part1 <- mclapply(3:1708, run_happi_thermophilus_e1, mc.cores=6)
saveRDS(thermophilus_results_sensitivity_e1_part1,"thermophilus_results_sensitivity_e1_part1.RDS")
thermophilus_results_sensitivity_e1_part2 <- mclapply(1709:2801, run_happi_thermophilus_e1, mc.cores=6)
saveRDS(thermophilus_results_sensitivity_e1_part2,"thermophilus_results_sensitivity_e1_part2.RDS")


# need to split the results into two parts 
#saveRDS(thermophilus_results_part1,"thermophilus_results_part1.RDS")
saveRDS(thermophilus_results_part2,"thermophilus_results_part2.RDS")
#saveRDS(thermophilus_results_part3,"thermophilus_results_part3.RDS")
thermophilus_results_part1 <- readRDS("thermophilus_results_part1.RDS")
thermophilus_results_part2 <- readRDS("thermophilus_results_part2.RDS")
#thermophilus_results_part3 <- readRDS("thermophilus_results_part3.RDS")
thermophilus_results_sensitivity_e01_part1 <- readRDS("thermophilus_results_sensitivity_e01_part1_mar23.RDS")
thermophilus_results_pensitivity_e01_part2 <- readRDS("thermophilus_results_sensitivity_e01_part2_mar23.RDS")

thermophilus_results_sensitivity_e05_part1 <- readRDS("thermophilus_results_sensitivity_e05_part1.RDS")
thermophilus_results_sensitivity_e05_part2 <- readRDS("thermophilus_results_sensitivity_e05_part2.RDS")

thermophilus_results_sensitivity_e1_part1 <- readRDS("thermophilus_results_sensitivity_e1_part1.RDS")
thermophilus_results_sensitivity_e1_part2 <- readRDS("thermophilus_results_sensitivity_e1_part2.RDS")


pvalue_happi_part1 <- lapply(thermophilus_results_part1, function(x) tail(x$loglik$pvalue_nopenalty[!is.na(x$loglik$pvalue_nopenalty)], 1)) %>% unlist
pvalue_happi_part2 <- lapply(thermophilus_results_part2, function(x) tail(x$loglik$pvalue_nopenalty[!is.na(x$loglik$pvalue_nopenalty)], 1)) %>% unlist
#saveRDS(pvalue_happi_part1,"pvalue_happi_part1.RDS")
#saveRDS(pvalue_happi_part2,"pvalue_happi_part2.RDS")

pvalue_happi_part3 <- lapply(thermophilus_results_sensitivity_e01_part1, function(x) tail(x$loglik$pvalue_nopenalty[!is.na(x$loglik$pvalue_nopenalty)], 1)) %>% unlist
pvalue_happi_part4 <- lapply(thermophilus_results_pensitivity_e01_part2, function(x) tail(x$loglik$pvalue_nopenalty[!is.na(x$loglik$pvalue_nopenalty)], 1)) %>% unlist
saveRDS(pvalue_happi_part3,"pvalue_happi_part3.RDS")
saveRDS(pvalue_happi_part4,"pvalue_happi_part4.RDS")

pvalue_happi_part5 <- lapply(thermophilus_results_sensitivity_e05_part1, function(x) tail(x$loglik$pvalue_nopenalty[!is.na(x$loglik$pvalue_nopenalty)], 1)) %>% unlist
pvalue_happi_part6 <- lapply(thermophilus_results_sensitivity_e05_part2, function(x) tail(x$loglik$pvalue_nopenalty[!is.na(x$loglik$pvalue_nopenalty)], 1)) %>% unlist
saveRDS(pvalue_happi_part5,"pvalue_happi_part5.RDS")
saveRDS(pvalue_happi_part6,"pvalue_happi_part6.RDS")

pvalue_happi_part7 <- lapply(thermophilus_results_sensitivity_e1_part1, function(x) tail(x$loglik$pvalue_nopenalty[!is.na(x$loglik$pvalue_nopenalty)], 1)) %>% unlist
pvalue_happi_part8 <- lapply(thermophilus_results_sensitivity_e1_part2, function(x) tail(x$loglik$pvalue_nopenalty[!is.na(x$loglik$pvalue_nopenalty)], 1)) %>% unlist
saveRDS(pvalue_happi_part7,"pvalue_happi_part7.RDS")
saveRDS(pvalue_happi_part8,"pvalue_happi_part8.RDS")


#saveRDS(pvalue_happi_part3,"pvalue_happi_part3.RDS")
pvalue_happi_part1 <- readRDS("pvalue_happi_part1.RDS")
pvalue_happi_part2 <- readRDS("pvalue_happi_part2.RDS")
pvalue_happi_part3 <- readRDS("pvalue_happi_part3.RDS")
pvalue_happi_part4 <- readRDS("pvalue_happi_part4.RDS")
pvalue_happi_part5 <- readRDS("pvalue_happi_part5.RDS")
pvalue_happi_part6 <- readRDS("pvalue_happi_part6.RDS")
pvalue_happi_part7 <- readRDS("pvalue_happi_part7.RDS")
pvalue_happi_part8 <- readRDS("pvalue_happi_part8.RDS")

pvalue_happi_e0 <- c(pvalue_happi_part1,
                     pvalue_happi_part2) %>% unlist

dim(thermophilus_results)
pvalue_happi_e01 <- c(pvalue_happi_part3,
                      pvalue_happi_part4) %>% unlist
pvalue_happi_e05 <- c(pvalue_happi_part5,
                      pvalue_happi_part6) %>% unlist
pvalue_happi_e1 <- c(pvalue_happi_part7,
                     pvalue_happi_part8) %>% unlist

sensitivity_results_BH <- tibble("Gene" = colnames(thermophilus_df)[3:2801],  
                                 pvalue_happi_e01,
                                 pvalue_happi_e05, 
                                 pvalue_happi_e1) %>%
  pivot_longer(cols = 2:5, names_to="test", values_to="pvalue") %>%
  group_by(test) %>%
  mutate("qvalue" = p.adjust(pvalue, "BH")) %>%
  mutate(test = str_remove(test, "pvalue_")) %>% 
  pivot_wider(names_from="test", values_from="pvalue":"qvalue") %>%
  left_join(thermophilus_core, by = "Gene")

run_glm_thermophilus <- function(colnum) {
  ha <- glm(unlist(thermophilus_df[,colnum]) ~ thermophilus_df$Spain, family="binomial")
  h0 <- glm(unlist(thermophilus_df[,colnum]) ~ 1, family="binomial")
  c(anova(ha, h0, test = "Rao")[2, "Pr(>Chi)"],
    anova(ha, h0, test = "LRT")[2, "Pr(>Chi)"])
}

thermophilus_df <- thermophilus_df %>% arrange(Spain)

thermophilus_results_glm <- lapply(3:2801, run_glm_thermophilus)
thermophilus_results_glm %<>% do.call(rbind, .)
#saveRDS(thermophilus_results_glm,"thermophilus_results_glm.RDS")
thermophilus_results_glm <- readRDS("thermophilus_results_glm.RDS")
thermophilus_results_glm <- readRDS("thermophilus_results_glm.RDS")
thermophilus_results_fisher <- lapply(3:2801, run_fisher_thermophilus)
thermophilus_results_fisher %<>% do.call(rbind, .)
#saveRDS(thermophilus_results_fisher,"thermophilus_results_fisher.RDS")
thermophilus_results_fisher <- readRDS("thermophilus_results_fisher.RDS")
dim(thermophilus_results_fisher)
thermophilus_results_fisher %<>% as_tibble() 
pvalue_happi <- pvalue_happi_e05 %>% unlist 

hyp_results_BH <- tibble("Gene" = colnames(thermophilus_df)[3:2801], 
                         pvalue_happi, 
                         pvalue_rao = thermophilus_results_glm[,1], 
                         pvalue_lrt = thermophilus_results_glm[,2]) %>%
  pivot_longer(cols = 2:5, names_to="test", values_to="pvalue") %>%
  group_by(test) %>%
  mutate("qvalue" = p.adjust(pvalue, "BH")) %>%
  mutate(test = str_remove(test, "pvalue_")) %>% 
  pivot_wider(names_from="test", values_from="pvalue":"qvalue")

hyp_results <- tibble("gene" = colnames(thermophilus_df)[3:2801], 
                      pvalue_happi, 
                      # "pvalue_fisher" = thermophilus_results_fisher,
                      pvalue_rao = thermophilus_results_glm[, 1], 
                      pvalue_lrt = thermophilus_results_glm[, 2]) %>%
  pivot_longer(cols = 2:4, names_to="test", values_to="pvalue") %>%
  group_by(test) %>%
  mutate("qvalue" = qvalue::qvalue(pvalue)$qvalue) %>%
  mutate(test = str_remove(test, "pvalue_")) %>% 
  pivot_wider(names_from="test", values_from="pvalue":"qvalue")

#thermophilus_core

hyp_results_BH <- hyp_results_BH %>% left_join(thermophilus_core, by = "Gene")
#saveRDS(hyp_results_BH,"hyp_results_BH_04022023.RDS")
hyp_results_BH<- readRDS("data/hyp_results_BH_04022023.RDS")


hyp_results_BH %>%
  select("Gene", starts_with("qvalue")) %>%
  pivot_longer(2:5) %>% 
  filter(value < 0.05) %>%
  group_by(name) %>%
  summarise(n = n()) 
# qvalue_happi    219
# qvalue_lrt      311
# qvalue_rao      254

hyp_results_BH %>% 
  filter(qvalue_happi < 0.05 & qvalue_lrt < 0.05) %>%
  dim()

hyp_results_BH %>% 
  filter(qvalue_happi < 0.05 & qvalue_rao < 0.05) %>%
  dim()


scgs <- hyp_results_BH %>% 
  filter(Core==1)

sig_core_list_thermophilus <- scgs %>% 
  select("Gene", starts_with("qvalue")) %>%
  pivot_longer(2:5) %>% 
  filter(value < 0.05)

scgs %>% 
  select("Gene", starts_with("qvalue")) %>%
  pivot_longer(2:5) %>% 
  filter(value < 0.05) %>%
  group_by(name) %>%
  summarise(n = n()) 
# qvalue_happi      3
# qvalue_lrt       27
# qvalue_rao        6
scgs %>% 
  filter(qvalue_happi < 0.05 & qvalue_lrt < 0.05) %>%
  dim()

scgs %>% 
  filter(qvalue_happi < 0.05 & qvalue_rao < 0.05) %>%
  dim()

scgs %>% 
  filter(qvalue_happi < 0.05 & qvalue_fisher < 0.05) %>%
  dim()


my_therm_plot <- hyp_results_BH %>% 
  mutate(log10_pvalue_happi = -log10(pvalue_happi), 
         log10_pvalue_rao = -log10(pvalue_rao), 
         log10_pvalue_lrt = -log10(pvalue_lrt), 
         log10_pvalue_fisher = -log10(pvalue_fisher)) %>%
  mutate(Genome = ifelse(is.na(Core), "Accessory","Core")) %>% 
  arrange(Genome)

library(forcats)
happi_rao <-  my_therm_plot %>% 
  ggplot(aes(x = log10_pvalue_happi, y = log10_pvalue_rao,  group = Genome, color = Genome)) +
  geom_point(size = 2.5, alpha = 0.8) + 
  geom_abline(intercept = 0, slope = 1) + 
  xlim(0,16) + ylim(0,16) + 
  ylab("-log10(pvalue Rao)") + 
  xlab("-log10(pvalue happi-a, epsilon 0.05)") + 
  scale_color_discrete(name="Genome") + 
  theme_bw() +
  scale_color_manual(values = c("Accessory" = "darkgrey", "Core" = "blue"))

happi_lrt <-  my_therm_plot  %>% 
  ggplot(aes(x = log10_pvalue_happi, y = log10_pvalue_lrt,  group = Genome, color = Genome)) +
  geom_point(size = 2.5, alpha = 0.8) + 
  geom_abline(intercept = 0, slope = 1) + 
  xlim(0,16) + ylim(0,16) + 
  ylab("-log10(pvalue LRT)") + 
  xlab("-log10(pvalue happi-a, epislon 0.05)") + 
  scale_color_discrete(name="Genome") + 
  theme_bw() +
  scale_color_manual(values = c("Accessory" = "darkgrey", "Core" = "blue"))

sensitivity_results_BH <- readRDS("data/sensitivity_results_BH.RDS")
sensitivity_results_BH %>%
  select("Gene", starts_with("qvalue")) %>%
  pivot_longer(2:5) %>% 
  filter(value < 0.05) %>%
  group_by(name) %>%
  summarise(n = n()) 

### sensitivity comparisons
sensitivity_plot <- sensitivity_results_BH %>% 
  mutate(log10_pvalue_happi_e0 = -log10(pvalue_happi_e0), 
         log10_pvalue_happi_e01 = -log10(pvalue_happi_e01), 
         log10_pvalue_happi_e05 = -log10(pvalue_happi_e05), 
         log10_pvalue_happi_e1 = -log10(pvalue_happi_e1)) %>%
  mutate(Genome = ifelse(is.na(Core), "Accessory","Core")) %>% 
  arrange(Genome)
sensitivity_plot <- readRDS("data/sensitivity_plot.RDS")
#saveRDS(sensitivity_results_BH,"sensitivity_results_BH.RDS")
#saveRDS(sensitivity_plot,"sensitivity_plot.RDS")
describe_sensitivity_results <- sensitivity_results_BH %>% 
  mutate(e1_change = 100*(pvalue_happi_e1-pvalue_happi_e05)/pvalue_happi_e05, 
         e01_change = 100*(pvalue_happi_e01 - pvalue_happi_e05)/pvalue_happi_e05, 
         e0_change = 100*(pvalue_happi_e0 - pvalue_happi_e05)/pvalue_happi_e05, 
         intersection = ifelse(qvalue_happi_e01 < 0.05 & qvalue_happi_e05 < 0.05 & qvalue_happi_e1 <0.05, 1,0)) %>%
  mutate(Genome = ifelse(is.na(Core), "Accessory","Core")) %>% 
  arrange(Genome)

describe_sensitivity_results %>%
  summarise(sum(intersection))

describe_sensitivity_results %>% 
  mutate(increase5_e1 = ifelse(e1_change > 5, 1,0),
         decrease5_e1 = ifelse(e1_change < -5, 1,0), 
         increase5_e01 = ifelse(e01_change > 5, 1, 0),
         decrease5_e01 = ifelse(e01_change < -5, 1, 0),
         increase5_e0 = ifelse(e0_change > 5,1,0), 
         decrease5_e0 = ifelse(e0_change < -5, 1,0)) %>%
  pivot_longer(15:20, values_to="percentchange", names_to="sensitivity") %>%
  group_by(sensitivity) %>%
  summarise("prop_change" = mean(percentchange),
            "n" = n()) %>%
  arrange(prop_change)

happi_e01 <-  sensitivity_plot %>% 
  ggplot(aes(x = log10_pvalue_happi_e05, y = log10_pvalue_happi_e01)) +
  geom_point(size = 2.5, alpha = 0.9, colour = "darkgrey") + 
  geom_abline(intercept = 0, slope = 1) + 
  xlim(0,16) + ylim(0,16) + 
  ylab("-log10(pvalue happi-a, epsilon 0.01)") + 
  xlab("-log10(pvalue happi-a, epsilon 0.05)") + 
  #scale_color_discrete(name="Genome") + 
  theme_bw() 
# scale_color_manual(values = c("Accessory" = "darkgrey", "Core" = "blue"))

happi_e0 <-  sensitivity_plot  %>% 
  ggplot(aes(x = log10_pvalue_happi_e05, y = log10_pvalue_happi_e0,  group = Genome, color = Genome)) +
  geom_point(size = 2.5, alpha = 0.8) + 
  geom_abline(intercept = 0, slope = 1) + 
  xlim(0,16) + ylim(0,16) + 
  ylab("-log10(pvalue happi-a, epsilon 0)") + 
  xlab("-log10(pvalue happi-a, epsilon 0.05)") + 
  scale_color_discrete(name="Genome") + 
  theme_bw() +
  scale_color_manual(values = c("Accessory" = "darkgrey", "Core" = "blue"))

happi_e1<-  sensitivity_plot  %>% 
  ggplot(aes(x = log10_pvalue_happi_e05, y = log10_pvalue_happi_e1)) +
  geom_point(size = 2.5, alpha = 0.9, colour = "darkgrey") + 
  geom_abline(intercept = 0, slope = 1) + 
  xlim(0,16) + ylim(0,16) + 
  ylab("-log10(pvalue happi-a, epsilon 0.1)") + 
  xlab("-log10(pvalue happi-a, epsilon 0.05)") + 
  #scale_color_discrete(name="Genome") + 
  theme_bw()
#scale_color_manual(values = c("Accessory" = "darkgrey", "Core" = "blue"))
ggpubr::ggarrange(happi_e01,  happi_e1, nrow = 1, common.legend=T, legend="right")
ggsave("sensitivity_comparisons_thermophilus_04062023.pdf", height = 3.5, width = 12)


