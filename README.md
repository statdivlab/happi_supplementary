
# `happi` supplementary materials 
This repository contains the necessary code and data to reproduce all results in the `happi` manuscript. [preprint available here](https://www.biorxiv.org/content/10.1101/2022.04.26.489591v1). The code and data contained in this repository can also be found at: [![DOI](https://zenodo.org/badge/483091425.svg)](https://zenodo.org/badge/latestdoi/483091425). 

The most recent `happi` `R` package is available [here](https://github.com/statdivlab/happi). 

`happi` version 0.8.7 can be used to produce the manuscript results and a static version has been made available: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8216120.svg)](https://doi.org/10.5281/zenodo.8216120). 

##  How to  use this repository 
Below you'll find a description of the `R` scripts and what analyses can be reproduced with each script. 
Please refer to the comments within the `R` scripts for the associated data files needed to run or generate the figures/results. 


|              **Script**             | **Description**                                                                                                                |                     **Data Files Needed**                     |
|:-----------------------------------:|--------------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------:|
|        `TM7_data_cleaning.R`        | Data cleaning script that uses Saccharibacteria data from public repositories and creates data file `TM7_presence_absence.RDS` | `TM7_data_cleaning_files/`                                    |
| `compare-nonparametric-f-options.R` | Script for evaluating estimators for _f_                                                                                       | `logit_DRR102664.RDS` `accuracy_f_nonparametric.RDS`          |
|       `type1-error-parallel.R`      | Type 1 error simulations in parallel                                                                                           | `logit_DRR102664.RDS` `type1_results.RDS`                     |
|       `type2-error-parallel.R`      | Type 2 error simulations in parallel                                                                                           | `logit_DRR102664.RDS` `type2_results.RDS`                     |
|          `happi-pkg-tm7.R`          | Data analysis of Saccharibacteria MAGs and creation of Figure 1                                                                | `TM7_presence_absence.RDS` `tm7_hyp_results_summary.RDS`      |
|     `happi-simulation-figures.R`    | Creation of Figures 2 & 3                                                                                                      | `logit_DRR102664.RDS` `type1_results.RDS` `type2_results.RDS` |

