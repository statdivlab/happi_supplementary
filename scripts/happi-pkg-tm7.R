### A simulation to compare the results of happi to competitor methods on the Saccharibacteria dataset
### Amy Willis, April 20 2022

library(tidyverse)
library(magrittr)
library(happi)
library(parallel)
TM7_data <- readRDS("data/TM7_presence_absence.RDS") %>% as_tibble 

tm7_df <- TM7_data %>% 
  select(site, mean_coverage, `Cellulase/cellobiase CelA1`:ncol(.)) %>%
  mutate(tongue = ifelse(site == "tongue", 1, 0)) %>% 
  select(2, ncol(.), 3:(ncol(.) - 1))

x_matrix <- model.matrix(~tongue, data = tm7_df)
run_happi_tm7 <- function(colnum) {
  happi(outcome=unlist(tm7_df[,colnum]), 
        covariate=x_matrix, 
        quality_var=tm7_df$mean_coverage,
        method="splines", 
        firth=T, 
        spline_df=4,
        max_iterations=1000, 
        change_threshold=0.01, 
        epsilon=0)
}

run_glm_tm7 <- function(colnum) {
  ha <- glm(unlist(tm7_df[,colnum]) ~ tm7_df$tongue, family="binomial")
  h0 <- glm(unlist(tm7_df[,colnum]) ~ 1, family="binomial")
  c(anova(ha, h0, test = "Rao")[2, "Pr(>Chi)"],
    anova(ha, h0, test = "LRT")[2, "Pr(>Chi)"])
}
run_fisher_tm7 <- function(colnum) {
  fisher.test(table(unlist(tm7_df[,colnum]), tm7_df$tongue))$p.value
}

## takes ages! can speed up by lowering max_iterations (from 1000) or increasing change_threshold (from 0.01)
tm7_results <- mclapply(3:715, run_happi_tm7, mc.cores=6)
saveRDS(tm7_results, "tm7_happi_splines.RDS")

tm7_results_glm <- lapply(3:715, run_glm_tm7)
tm7_results_glm %<>% do.call(rbind, .)
tm7_results_fisher <- lapply(3:715, run_fisher_tm7)
tm7_results_fisher %<>% unlist


pvalue_happi <- lapply(tm7_results, function(x) tail(x$loglik$pvalue[!is.na(x$loglik$pvalue)], 1)) %>% unlist
hyp_results <- tibble("gene" = colnames(tm7_df)[3:715], 
                      pvalue_happi, 
                      "pvalue_fisher" = tm7_results_fisher,
                      pvalue_rao = tm7_results_glm[, 1], 
                      pvalue_lrt = tm7_results_glm[, 2]) %>%
    pivot_longer(cols = 2:5, names_to="test", values_to="pvalue") %>%
    group_by(test) %>%
    mutate("qvalue" = p.adjust(pvalue, "BH")) %>%
    mutate(test = str_remove(test, "pvalue_")) %>% 
    pivot_wider(names_from="test", values_from="pvalue":"qvalue")
#saveRDS(hyp_results, "tm7_hyp_results_summary.RDS")
hyp_results <- readRDS("data/tm7_hyp_results_summary.RDS")

hyp_results %>%
  filter(gene %in% c("Ribosomal protein L27", 
                     "Membrane protein insertase Oxa1/YidC/SpoIIIJ, required for the localization of integral membrane proteins",
                     "DNA topoisomerase IA",
                     "3-phosphoglycerate kinase",
                     "Preprotein translocase subunit SecE", 
                     "Surface antigen"))

##################
### make figure 
##################
TM7_data %>%
  select(MAGs_Id, site, mean_coverage, 
         `Ribosomal protein L27`, 
         `Membrane protein insertase Oxa1/YidC/SpoIIIJ, required for the localization of integral membrane proteins`,
         `DNA topoisomerase IA`, 
         `3-phosphoglycerate kinase`,
         `Preprotein translocase subunit SecE`, 
         `Surface antigen`) %>%
  pivot_longer(4:ncol(.), names_to="gene") %>%
  inner_join(hyp_results)  %>% 
  mutate(gene = ifelse(gene == "Membrane protein insertase Oxa1/YidC/SpoIIIJ, required for the localization of integral membrane proteins", 
                       "Membrane protein insertase Oxa1/YidC/SpoIIIJ", gene)) %>%
  mutate(gene = factor(gene, ordered=T, 
                       levels=c("Ribosomal protein L27", 
                                "Membrane protein insertase Oxa1/YidC/SpoIIIJ",
                                "DNA topoisomerase IA",
                                "3-phosphoglycerate kinase",
                                "Preprotein translocase subunit SecE", 
                                "Surface antigen"))) %>%
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
# ggsave("../manuscript/figures/TM7_specific_comparisons_04222022.pdf", width = 10.5, height = 4.5)

####################################
### get results for text
####################################

hyp_results %>% 
  ggplot(aes(pvalue_rao, pvalue_lrt)) +
  geom_point()

hyp_results %>%
  select("gene", starts_with("qvalue")) %>%
  pivot_longer(2:5) %>% 
  filter(value < 0.05) %>%
  group_by(name) %>%
  summarise(n = n())

hyp_results %>%
  mutate(pvalue_happi < pvalue_lrt, 
         pvalue_happi < pvalue_rao, 
         pvalue_happi < pvalue_fisher) %>%
  summarise(across(starts_with("pvalue_happi <"), mean))
