rm(list=ls())
# MR analysis of measures of adiposity and metabolites 

# environment ====
## library ====
#remotes::install_github("MRCIEU/TwoSampleMR")
#remotes::install_github("mrcieu/ieugwasr")

library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(RadialMR)
library(dplyr)

### methods
methods <- mr_method_list()
methods_heterogeneity <- subset(methods, heterogeneity_test == TRUE)$obj
methods_heterogeneity <- methods_heterogeneity[c(1,2,3,5)]
methods <- methods[c(1,2,3,6,10,13),1]

### colours
#install.packages("wesanderson")
library(wesanderson)
d1 <- wes_palette("Royal1", type = "discrete")
d2 <- wes_palette("GrandBudapest2", type = "discrete")
d3 <- wes_palette("Cavalcanti1", type = "discrete")
d4 <- wes_palette("Rushmore1", type = "discrete")
discrete_wes_pal <- c(d1, d2, d3, d4)
rm(d1,d2,d3,d4)

exposure_data <- extract_instruments(outcomes=c('met-d-LDL_TG', 'met-d-HDL_TG', 'met-d-LDL_PL', 'met-d-LDL_CE', 'met-d-non_HDL_C', 'met-d-Total_FC', 'met-d-LDL_FC',
'met-d-LDL_L', 'met-d-LDL_P', 'met-d-VLDL_size', 'met-d-Sphingomyelins', 'met-d-ApoB', 'met-d-Unsaturation', 'met-d-Omega_3', 'met-d-Clinical_LDL_C', 
'met-d-DHA', 'met-d-Omega_3_pct', 'met-d-LDL_C', 'met-d-Omega_6_by_Omega_3', 'met-d-Gly', 'met-d-Phe', 'met-d-Citrate', 'met-d-Total_TG', 'met-d-Albumin', 
'met-d-XXL_VLDL_P', 'met-d-XXL_VLDL_L', 'met-d-XXL_VLDL_PL', 'met-d-XXL_VLDL_CE', 'met-d-XXL_VLDL_FC', 'met-d-XXL_VLDL_TG', 'met-d-XXL_VLDL_TG', 'met-d-XL_VLDL_TG',
'met-d-L_VLDL_L', 'met-d-XS_VLDL_P', 'met-d-XS_VLDL_L', 'met-d-L_VLDL_FC', 'met-d-XS_VLDL_PL', 'met-d-XS_VLDL_C', 'met-d-XS_VLDL_CE', 'met-d-XS_VLDL_TG', 
'met-d-IDL_P', 'met-d-IDL_L', 'met-d-IDL_PL', 'met-d-IDL_C', 'met-d-IDL_CE', 'met-d-IDL_FC', 'met-d-IDL_TG', 'met-d-L_LDL_P', 'met-d-L_LDL_L', 'met-d-L_LDL_C', 'met-d-L_LDL_FC',
'met-d-L_LDL_TG', 'met-d-M_LDL_L', 'met-d-M_LDL_PL', 'met-d-M_LDL_FC', 'met-d-S_LDL_L', 'met-d-S_LDL_C', 'met-d-S_LDL_CE', 'met-d-S_LDL_FC', 'met-d-S_LDL_TG',
'met-d-XL_HDL_P', 'met-d-XL_HDL_PL', 'met-d-XL_HDL_TG', 'met-d-L_HDL_TG', 'met-d-M_VLDL_C', 'met-d-M_HDL_TG', 'met-d-M_VLDL_CE', 'met-d-XXL_VLDL_PL_pct', 
'met-d-XXL_VLDL_C_pct', 'met-d-XXL_VLDL_TG_pct', 'met-d-XL_VLDL_PL_pct', 'met-d-XL_VLDL_CE_pct', 'met-d-XL_VLDL_TG_pct', 'met-d-L_VLDL_PL_pct', 'met-d-L_VLDL_C_pct',
'met-d-L_VLDL_CE_pct', 'met-d-L_VLDL_TG_pct', 'met-d-S_VLDL_PL_pct', 'met-d-XS_VLDL_TG_pct', 'met-d-IDL_PL_pct', 'met-d-IDL_C_pct', 'met-d-IDL_CE_pct', 'met-d-IDL_TG_pct',
'met-d-L_LDL_PL_pct', 'met-d-L_LDL_C_pct', 'met-d-S_VLDL_C_pct', 'met-d-L_LDL_TG_pct', 'met-d-M_LDL_FC_pct', 'met-d-M_LDL_TG_pct', 'met-d-S_LDL_C_pct', 
'met-d-S_LDL_FC_pct', 'met-d-S_LDL_TG_pct', 'met-d-XL_HDL_C_pct', 'met-d-XL_HDL_CE_pct', 'met-d-XL_HDL_FC_pct', 'met-d-S_VLDL_FC_pct', 'met-d-M_HDL_CE_pct',
'met-d-M_HDL_TG_pct', 'met-d-S_HDL_PL_pct', 'met-d-S_HDL_C_pct', 'met-d-S_HDL_CE_pct', 'met-d-S_HDL_FC_pct', 'met-d-S_VLDL_TG_pct', 'met-d-XS_VLDL_C_pct','met-d-VLDL_TG'), p1 = 5e-8, clump = T)



## extract outcome data ====
outcome_data <- extract_outcome_data(exposure_data$SNP, c('ieu-a-835', 'ieu-a-79', 'ieu-a-27', 'ieu-a-89'), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

## harmonize data ====
harmonise_data <- harmonise_data(exposure_data, outcome_data, action = 2)

## MR ====
mr_results <- mr(harmonise_data, method_list = methods)

## Sensitivity analysis ====
mr_singlesnp <- mr_singlesnp(harmonise_data)
mr_hetrogeneity <- mr_heterogeneity(harmonise_data,
                                    method_list = methods_heterogeneity)
mr_pleiotropy <- mr_pleiotropy_test(harmonise_data)
mr_leaveoneout <- mr_leaveoneout(harmonise_data)

## Plots ====
plot_mr_scatter <- mr_scatter_plot(mr_results, harmonise_data)
plot_singlesnp_forest <- mr_forest_plot(mr_singlesnp)
plot_leaveoneout_forest <- mr_leaveoneout_plot(mr_leaveoneout)
plot_mr_funnel <- mr_funnel_plot(mr_singlesnp)

### save plots ====
pdf("~/001_projects/006_MR/plot_mr_scatter_metab_ex.pdf")
for (i in 1:length(plot_mr_scatter)) {
  print(plot_mr_scatter[[i]])
}
dev.off()

pdf("~/001_projects/006_MR/plot_singlesnp_forest_metab_ex.pdf")
for (i in 1:length(plot_singlesnp_forest)) {
  print(plot_singlesnp_forest[[i]])
}
dev.off()

pdf("~/001_projects/006_MR/plot_leaveoneout_forest_metab_ex.pdf")
for (i in 1:length(plot_leaveoneout_forest)) {
  print(plot_leaveoneout_forest[[i]])
}
dev.off()

pdf("/~/001_projects/006_MR/plot_mr_funnel_metab_ex.pdf")
for (i in 1:length(plot_mr_funnel)) {
  print(plot_mr_funnel[[i]])
}
dev.off()

## Save output ====
write.table(exposure_data, "~/001_projects/006_MR/exposure_data_metab_ex.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(harmonise_data, "~/001_projects/006_MR/harmonise_data_metab_ex.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_results, "~/001_projects/006_MR/mr_results_metab_ex.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_singlesnp, "~/001_projects/006_MR/mr_singlesnp_metab_ex.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_hetrogeneity, "~/001_projects/006_MR/mr_hetrogeneity_metab_ex.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_pleiotropy, "~/001_projects/006_MR/mr_pleiotropy_metab_ex.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_leaveoneout, "~/001_projects/006_MR/mr_leaveoneout_metab_ex.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")



