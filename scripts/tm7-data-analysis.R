# happi data analysis of TM7 update 
# Pauline Trinh 
# February 27, 2023 

library(tidyverse)
library(magrittr)
library(happi)
library(parallel)

COG_data <- read_tsv('/Users/paulinetrinh/Documents/pangenomics-hypothesis-testing/ANVIO_TM7_METAPANGENOME/TM7-Summary/TM7_gene_clusters_summary.txt')
COG_manual <- read_csv("/Users/paulinetrinh/Documents/GitHub/pangenomics-hypothesis-testing/R/cluster/SCG_COG.csv")

COG_annotated <- COG_data %>% filter(!is.na(COG_FUNCTION)) %>% 
  filter(bin_name == "Core") %>%
  select(bin_name,COG_FUNCTION) %>% 
  distinct() %>% 
  mutate(gene = COG_FUNCTION) %>% 
  select(-COG_FUNCTION)  %>% full_join(COG_manual) 

TM7_data2 <- readRDS("/Users/paulinetrinh/Documents/GitHub/happi_supplementary/data/TM7_presence_absence_04032023.RDS") %>% as_tibble 
TM7_data <- readRDS("/Users/paulinetrinh/Documents/GitHub/happi_supplementary/data/TM7_presence_absence.RDS") %>% as_tibble 

tm7_df <- TM7_data %>% 
  select(site, mean_coverage, `Cellulase/cellobiase CelA1`:ncol(.)) %>%
  mutate(tongue = ifelse(site == "tongue", 1, 0)) %>% 
  select(2, ncol(.), 3:(ncol(.) - 1))

tm7_df2 <- TM7_data2 %>% 
  select(site, mean_coverage,log_completion,log_completion_prop,percent_completion, `Cellulase/cellobiase CelA1`:ncol(.)) %>%
  mutate(tongue = ifelse(site == "tongue", 1, 0)) %>% 
  select(2, ncol(.), 3:(ncol(.) - 1))

x_matrix <- model.matrix(~tongue, data = tm7_df)

run_happi_tm7 <- function(colnum) {
  happi_results <- happi(outcome=unlist(tm7_df[,colnum]), 
                              covariate=x_matrix, 
                              quality_var=tm7_df$mean_coverage,
                              method="splines", 
                              firth=T, 
                              spline_df=4,
                              max_iterations=1000, 
                              change_threshold=0.1, 
                              epsilon=0)
  
  npLRT_results <- happi::npLRT(happi_results, 
                                method="splines", 
                                firth=T, 
                                spline_df=4,
                                max_iterations=1000, 
                                change_threshold=0.1, 
                                epsilon=0, 
                                P = 1000)
  
  return(npLRT_results)
}

run_happi_tm7_meancoverage <- function(colnum) {
  happi(outcome=unlist(tm7_df[,colnum]), 
             covariate=x_matrix, 
             quality_var=tm7_df$mean_coverage,
             method="splines", 
             firth=T, 
             spline_df=4,
             max_iterations=1000, 
             change_threshold=0.1, 
             epsilon=0)
}

x_matrix2 <- model.matrix(~tongue, data = tm7_df2)
run_happi_tm7_completion_beta <- function(colnum) {
  happi(outcome=unlist(tm7_df2[,colnum]), 
             covariate=x_matrix2, 
             quality_var=tm7_df2$log_completion,
             method="splines", 
             firth=T, 
             spline_df=4,
             max_iterations=1000, 
             change_threshold=0.1, 
             epsilon=0)
}

run_glm_tm7 <- function(colnum) {
  ha <- glm(unlist(tm7_df[,colnum]) ~ tm7_df$tongue, family="binomial")
  h0 <- glm(unlist(tm7_df[,colnum]) ~ 1, family="binomial")
  c(anova(ha, h0, test = "Rao")[2, "Pr(>Chi)"],
    anova(ha, h0, test = "LRT")[2, "Pr(>Chi)"])
}


x_matrix <- model.matrix(~tongue, data = tm7_df)


## takes ages! can speed up by lowering max_iterations (from 1000) or increasing change_threshold (from 0.01)
set.seed(1)
start <- Sys.time()
#tm7_results <- mclapply(3:715, run_happi_tm7, mc.cores=6)
tm7_results_seed1 <- mclapply(3:715, run_happi_tm7, mc.cores=6)
stop <- Sys.time()
stop-start
#Time difference of 46.84393 mins (3:13) 11
start <- Sys.time()
#tm7_results <- mclapply(3:715, run_happi_tm7, mc.cores=6)
run_happi_tm7(3)
stop <- Sys.time()
stop-start
start2 <- Sys.time()
run_happi_tm7_meancoverage(3)
stop2 <- Sys.time()
stop2-start2
#saveRDS(tm7_results, "tm7_happi_splines_Feb27.RDS")
#saveRDS(tm7_results, "tm7_happi_splines_Mar8.RDS")
#saveRDS(tm7_results_seed1,"tm7_results_seed1_Mar10.RDS")
set.seed(2)
tm7_results_completion_beta <- mclapply(6:718, run_happi_tm7_completion_beta, mc.cores=6)
saveRDS(tm7_results_completion_beta,"tm7_results_completion_beta.RDS")
tm7_results_meancoverage_beta <- mclapply(3:715, run_happi_tm7_meancoverage, mc.cores=6)
saveRDS(tm7_results_meancoverage_beta,"tm7_results_meancoverage_beta.RDS")
beta1_hat_completion <-lapply(tm7_results_completion_beta, function(x) tail(x$beta[!is.na(x$beta)], 1)) %>% unlist
beta1_hat_coverage <-lapply(tm7_results_meancoverage_beta, function(x) tail(x$beta[!is.na(x$beta)], 1)) %>% unlist

beta1_hat_comparison <- tibble("gene" = colnames(tm7_df)[3:715], 
                               beta1_hat_completion, 
                               beta1_hat_coverage) 
library(latex2exp)
beta1_hat_comparison %>% 
  ggplot(aes(x = beta1_hat_coverage, y = beta1_hat_completion)) + 
  geom_point(cex = 0.5) + 
  theme_bw() + 
  xlab(TeX("$\\hat{\\beta}_{1}$ using mean coverage")) + 
  #xlab(~paste(beta[1], " using mean coverage")) +   
  ylab(TeX("$\\hat{\\beta}_{1}$ using completion %")) + 
  theme(legend.position="none") + 
  geom_abline(intercept = 0, slope = 1)
ggsave("beta1_comparison_quality_var.png", width = 8, height = 4)


tm7_results_seed1 <- readRDS("data/tm7_results_seed1_Mar10.RDS")

tm7_results_glm <- lapply(3:715, run_glm_tm7)
tm7_results_glm %<>% do.call(rbind, .)

#tm7_results <- readRDS("tm7_happi_splines_Feb27.RDS")
#pvalue_happi <- lapply(tm7_results, function(x) tail(x$loglik$pvalue_nopenalty[!is.na(x$loglik$pvalue)], 1)) %>% unlist
pvalue_happi <- tm7_results_seed1 %>% unlist
#hyp_results <- tibble("gene" = colnames(tm7_df)[3:715], 
#                      pvalue_happi, 
#                      "pvalue_fisher" = tm7_results_fisher,
#                      pvalue_rao = tm7_results_glm[, 1], 
#                      pvalue_lrt = tm7_results_glm[, 2]) %>%
#  pivot_longer(cols = 2:5, names_to="test", values_to="pvalue") %>%
#  group_by(test) %>%
#  mutate("qvalue" = qvalue::qvalue(pvalue)$qvalue) %>%
#  mutate(test = str_remove(test, "pvalue_")) %>% 
#  pivot_wider(names_from="test", values_from="pvalue":"qvalue")

hyp_results_BH <- tibble("gene" = colnames(tm7_df)[3:715], 
                         pvalue_happi, 
                         pvalue_rao = tm7_results_glm[, 1], 
                         pvalue_lrt = tm7_results_glm[, 2]) %>%
  pivot_longer(cols = 2:4, names_to="test", values_to="pvalue") %>%
  group_by(test) %>%
  mutate("qvalue" = p.adjust(pvalue, "BH"))%>%
  mutate(test = str_remove(test, "pvalue_")) %>% 
  pivot_wider(names_from="test", values_from="pvalue":"qvalue")
saveRDS(hyp_results_BH,"hyp_results_BH.RDS")

##################
### make figure 
##################

TM7_data %>%
  select(MAGs_Id, site, mean_coverage, 
         `Ribosomal protein L27`, 
         `Ribosomal protein S18`,
         `Ribosomal protein L33`, 
         `3-phosphoglycerate kinase`,
         `Preprotein translocase subunit SecE`, 
         `Predicted amidohydrolase`) %>%
  pivot_longer(4:ncol(.), names_to="gene") %>%
  inner_join(hyp_results_BH)  %>% 
  mutate(gene = factor(gene, ordered=T, 
                       levels=c("Ribosomal protein L27", 
                                "Ribosomal protein S18",
                                "Ribosomal protein L33",
                                "3-phosphoglycerate kinase",
                                "Preprotein translocase subunit SecE", 
                                "Predicted amidohydrolase"))) %>%
  mutate(happi = paste("happi p =", ifelse(pvalue_happi > 0.001, round(pvalue_happi, 3), round(pvalue_happi, 6)))) %>%
  mutate(rao = paste("GLM-Rao p =", ifelse(pvalue_rao > 0.001, round(pvalue_rao, 3), round(pvalue_rao, 6)))) %>%
  mutate(happi = ifelse(MAGs_Id == "ORAL_P_A_F_Bin_00012", happi, NA)) %>%
  mutate(rao = ifelse(MAGs_Id == "ORAL_P_A_F_Bin_00012", rao, NA)) %>%
  ggplot() +
  facet_wrap(~gene) +
  geom_jitter(aes(x = mean_coverage, y = value, col = site, pch = site), height=0.08, width=0.00) +
  xlab("Mean coverage") + ylab("") +
  theme_bw() + 
  geom_text(x = 10, y = -.14, aes(label = happi), hjust = "left") + 
  geom_text(x = 10, y = -0.3, aes(label = rao), hjust = "left") + 
  scale_colour_manual(values= c("mediumseagreen", "dodgerblue")) + 
  theme(legend.position="right") +
  scale_y_continuous(breaks = c(0,1),
                     label = c("Not detected", "Detected"), limits=c(-0.32, 1.1)) +
  NULL

ggsave("TM7_specific_comparisons_03272023_coregenesoption.pdf", width = 10.5, height = 4.5)

####################################
### get results for text
####################################

hyp_results_BH %>%
  select("gene", starts_with("qvalue")) %>%
  pivot_longer(2:4) %>% 
  filter(value < 0.05) %>%
  group_by(name) %>%
  summarise(n = n()) 
### Results for intersection of differentially present genes 
hyp_results_BH %>%
  select("gene", starts_with("qvalue")) %>%
  pivot_longer(2:4) %>% 
  filter(value < 0.05) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  filter(!is.na(qvalue_happi) & !is.na(qvalue_rao)) %>% dim()
# 166 intersecting with Rao 
hyp_results_BH %>%
  select("gene", starts_with("qvalue")) %>%
  pivot_longer(2:4) %>% 
  filter(value < 0.05) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  filter(!is.na(qvalue_happi) & !is.na(qvalue_lrt)) %>% dim()
# 176 intersection with GLM-LRT 
#Glm-Rao identified 9 different genes 

# percentage results for main text 
hyp_results_BH %>%
  mutate(pvalue_happi < pvalue_lrt, 
         pvalue_happi < pvalue_rao) %>%
  summarise(across(starts_with("pvalue_happi <"), mean))

hyp_results_BH_scgs <- hyp_results_BH %>% 
  left_join(COG_annotated)

scgs_BH <- hyp_results_BH_scgs %>% 
  filter(bin_name == "Core")

# how many core genes identified
scgs_BH %>%
  select("gene", starts_with("qvalue")) %>%
  pivot_longer(2:4) %>% 
  filter(value < 0.05) %>%
  group_by(name) %>%
  summarise(n = n()) 
# numbers for text 

scgs_BH %>% 
  filter(qvalue_happi < 0.05 & qvalue_lrt < 0.05) %>%
  dim()
scgs_BH %>% 
  filter(qvalue_happi < 0.05 & qvalue_rao < 0.05) %>%
  dim()

