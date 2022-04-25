## R script to create TM7_presence_absence.RDS
## All necessary data files can be found in the folder:
# happi_supplementary/data/TM7_data_cleaning_files/

library(tidyverse)
#######################################################################################################################
#### TM7-functions-occurrence-frequency.txt was created using anvi'o version 6 with the following python/conda script:
# anvi-compute-functional-enrichment -p TM7-PAN.db \                                                                    
# -g TM7-GENOMES.db \
# --annotation-source COG_FUNCTION \
# --category-variable Clade \
# -o TM7-enriched-functions-clade.txt \
# --functional-occurrence-table-output TM7-functions-occurrence-frequency.txt
#######################################################################################################################
functions <- read_tsv("data/TM7_data_cleaning_files/TM7-functions-occurrence-frequency.txt")
functions_subset <- subset(functions, select = c(`...1`,`ORAL_P_A_F_Bin_00012`:`ORAL_T_D_F_MAG_00010`))
colnames(functions_subset)[1] <- "X1"
coverages <- readxl::read_excel("data/TM7_data_cleaning_files/Table_S02_Genomic_DB.xlsx",sheet="(c) NR_genomes_mean_coverage")
# Need to replace - with _ in column 1 of coverages file MAGs_Id
coverages$MAGs_Id <- gsub('-', '_', coverages$MAGs_Id)
MAG_tibble <- read_tsv("data/TM7_data_cleaning_files/43_MAGS_list.txt") %>% as_tibble()
# Want to subset the columns where we have the combined mean coverages 
coverages_sub <- coverages[,grep("_com|_Id|_length",colnames(coverages))]
MAG_coverages <- MAG_tibble %>% inner_join(coverages_sub)

## The mean coverage of the MAG relates to the value of _com that corresponds to the same host
#######################################################################################################################
# NOTE: I manually assigned the appropriate mean coverages for each MAG based on their host                           #
# For ex: MAG ORAL_P_A_F_Bin_00012 the mean coverage for this MAG within subject P_A_F was 7.99 found under P_A_F_com #
# For ORAL_P_A_M_Bin_00015 the mean coverage was 20.984678752 located under P_A_M_com                                 #
#######################################################################################################################
# write.csv(MAG_coverages, file = "MAG_coverages.csv")
read.csv("data/TM7_data_cleaning_files/MAG_coverages_MOD.csv") -> MAG_coverages_MOD
summary(MAG_coverages_MOD$mean_coverage)

functions_transpose <- as.data.frame(t(as.matrix(functions_subset)))
colnames(functions_transpose) <- as.character(unlist(functions_transpose[1,]))
functions_transpose = functions_transpose[-1, ]
functions_transpose[] <- lapply(functions_transpose, function(x) as.numeric(as.character(x)))
functions_transpose <- tibble::rownames_to_column(functions_transpose, "MAGs_Id")
functions_bin <- functions_transpose
functions_bin[,2:756] <- ifelse(functions_bin[,2:756] > 0,1,0)

MAG_data <- MAG_coverages_MOD %>% full_join(functions_bin, by = "MAGs_Id")
MAG_data %>% 
  mutate(tongue = ifelse(site == "tongue",1,0)) -> TM7_data
saveRDS(TM7_data, file = "TM7_presence_absence.RDS")
