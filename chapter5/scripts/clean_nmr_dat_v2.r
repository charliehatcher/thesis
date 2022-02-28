############################################################################
#
#                           Prepare NMR data                               #
#
############################################################################

# Changes from v1:
# Added in 11/03/2021 save NMR/markers data for all females combined
# Fixed error in INT transforming markers within stratum


############################################################################
#                                Set - UP                                  #
############################################################################

#module add languages/R-3.6.2-gcc9.1.0

##
## Clear the work environment
rm(list = ls())


##
## R version 
print(sessionInfo()) # R version 3.6.2 (2019-12-12)


##
## Set working directory for UKBB NMR data (app 30418, release 2020-07-01)
setwd(paste0(Sys.getenv('UKB_NMR'), "/data/NGH_UKB_phase1_metabolomics/"))


##
## Directory for reproductive health project
repro_dir <- paste0(Sys.getenv('UKB_all'), "/phenotypic/applications/30418/dev/release_candidate/data/repro_health/")


##
## Setting repository (UoB)
options(repos = c(CRAN ="http://www.stats.bris.ac.uk/R/"))


##
## Library directory
.libPaths()


##
## Updating packages 
#update.packages(ask = "FALSE")


##
## Install packages
#install.packages(c("data.table", "purrr", "dplyr", "skimr", "xlsx", "RNOmni"))


##
## Load libraries
library(data.table)
library(dplyr)
library(purrr)


############################################################################
#                   Import and describe 'raw' data                         #
############################################################################


##
## Import data

nmr_dat <- fread("./NGH_phase1_final_delivery/UKB-phase1__results.tsv")

##
## Describe 'raw data'

# Number of observations:
nrow(nmr_dat) # 126846

# Variables distributions
sumdat <- skimr::skim(nmr_dat) 
sumdat

# Missingness distribution:
summary(sumdat$n_missing) # Range: 0 - 6012

# Save summary NMR data
write.csv(sumdat, file = paste0(repro_dir, "data/sum_nmrdat.csv"))


############################################################################
#                           Sample-level QC                                #
############################################################################


##
## Remove sample failures

# List of failed samples IDs to be removed
fail <- fread("./sample_failures/ukb_phase1_samples_to_be_removed.tsv")

# Number of failed observations:
nrow(fail) # 40

# Number of failed IDs:
length(unique(fail$sample_id)) # 39

# Duplicated ID
fail[duplicated(fail$sample_id)] 

# Remove failed samples
nmr_dat <- filter(nmr_dat, ! sample_id %in% fail$sample_id)

# N after removing failed observations:
nrow(nmr_dat) # 126806


##
## Check for poor quality samples

# Nightingale: 
# "I recommend performing the sensitivity analysis by excluding the samples tagged with "Low Protein" and biomarker values tagged with "Below_limit_of_quantification". 

# Sample level tags (indicate samples with poor quality)
stag <- fread("./NGH_phase1_final_delivery/UKB-phase1__sample_tags_and_notes.tsv") 

# Frequency of "Low protein" sample tag:
table(stag$Low_protein) # 909

# Distribution of sample notes:
table(stag$Notes)
# P.S. where samples have any of the notes defined above, the biomarker values have been set to missing (except for one sample with a "Technical error" note)


##
## Remove blind duplicates

# Blind duplicates file
bdup <- fread("./Blind_duplicates/ukb_phase1_blind_duplicates.tsv") 

# Keep one sample at random within blind duplicates 
set.seed(3723)
bdup_keep <- bdup %>%
			group_by(bsd_pair, delivery) %>%
			sample_n(1) %>%
			ungroup(bsd_pair, delivery) %>%
			mutate(sel = "yes") %>%
			merge(., bdup, by = c("sample_id", "plate_id", "plate_position", "visit", "bsd_pair", "delivery"), all.y = T) %>%
			mutate(keep = (ifelse(is.na(sel), F, T)))

# Merge blind duplicates file with NMR data ad remove blind duplicates
nmr_dat <- nmr_dat %>%
			merge(., bdup_keep, by = c("sample_id", "plate_id", "plate_position", "visit"), all.x = T) %>%
			filter(keep == T | is.na(keep)) %>%
			select(-bsd_pair, -delivery, -sel, -keep)
			
# N after removing blind duplicates:
nrow(nmr_dat) # 123670


##
## Remove sample duplicates 

# Nigthingale: "please select plate 0490000006108 for samples 1219449 and 4478891 as this plate has been measured in a single spectrometer"
nmr_dat <- filter(nmr_dat, ! (sample_id %in% c(1219449, 4478891) & plate_id == 490000006107)) 

# N after removing sample duplicates:
nrow(nmr_dat) # 123668

# Nithingale: "Phase 1 Nightingale data includes ~5,000 repeat visit samples. Among those, 1,437 samples have also corresponding baseline visit samples in Phase 1"
# Number of remaining duplicated IDs:  
table(duplicated(nmr_dat$sample_id)) # 1402 duplicated sample IDs

# Keep a single record from remaining duplicated samples (pick at random)
set.seed(5683)
nmr_dat <- nmr_dat %>%
            group_by(sample_id) %>%
            sample_n(1) %>%
			ungroup(sample_id)

# N after removing sample duplicates at random:
nrow(nmr_dat) # 122266


##
## Exclude individuals who withdrew consent 

# List of individuals who withdrew consent (up to July 2020)
woc <- fread("./withdrawals/ukb_phase1_samples_to_be_removed.tsv") 
		
# Merge list with NMR data and exclude WoC individuals
nmr_dat <- nmr_dat %>%
			filter(! sample_id %in% woc$sample_id)

# N after removing WoC participants:
nrow(nmr_dat) # 122266


############################################################################
#                           Biomarker-level QC                             #
############################################################################


# Nightingale: 
# "In general, if a biomarker has a tag, but biomarker value is still provided, it means that during the manual QC step, 
# we have evaluated that the presence of the interfering substance is low, and it is not interfering with the quantification 
# of the biomarker (the value can be trusted)"
# "I recommend performing the sensitivity analysis by excluding the samples tagged with "Low Protein" and biomarker values tagged 
# with "Below_limit_of_quantification". 
# Unfortunately, we don't yet have a clear explanation why some biomarkers tagged with "technical error" have not been set to missing 
# but the number of such cases should be very low so you could set the values to missing. 

##
## Number of observations with bad tags per metabolite

# Metabolite level tags (indicate metabolites with poor quality)
mtag <- fread("./NGH_phase1_final_delivery/UKB-phase1__biomarker_tags_and_notes.tsv")

# Count number of "Technical_error" per metabolite			  
error_count <- mtag %>%
				mutate_at(., vars(Total_C:S_HDL_TG_pct), list(~ ifelse(. == "Technical_error", 1, 0)))  %>%
				summarise_at(.,  vars(Total_C:S_HDL_TG_pct), list(~ sum(.))) %>%
				summarise(rowmax = max(., na.rm = T))				

error_count # Maximum of 135 of "Technical_error" per metabolite

# Count maximum number of "Below_limit_of_quantification" per metabolite
low_count <- mtag %>%
             mutate_at(., vars(Total_C:S_HDL_TG_pct), list(~ ifelse(. == "Below_limit_of_quantification", 1, 0))) %>%
             summarise_at(.,  vars(Total_C:S_HDL_TG_pct), list(~ sum(.))) %>%
             summarise(rowmax = max(., na.rm = T))

low_count # Maximum of 110736 of "Below_limit_of_quantification" per metabolite


##
## Set to missing observations with "Technical error" tag

# Merge NMR data to biomarker tag data
nmr_dat <- mtag %>%
		rename_at(vars(Total_C:S_HDL_TG_pct), ~ paste0(., "_TAG")) %>%
		merge(nmr_dat, ., by = c("sample_id", "plate_id", "plate_position", "visit", "spectrometer")) 

# Set to missing if tag == "Technical_error"	
for(x in names(nmr_dat[6:254])) {

		print(paste0("Var: ", x))
		
		# Position of metabolite variable
		i <- grep(paste0("^", x, "$"), names(nmr_dat))
		
		# Position of metabolite tag variable 
		itag <- grep(paste0("^", x, "_TAG$"), names(nmr_dat))
		
		# "Technical_error" logical vector 
		is_error <- grepl("Technical_error", nmr_dat[[itag]])
		print(table(is_error))
		
		# Set to missing if "Technical_error" == TRUE
		print("Number of observations set to missing:")
		print(select(nmr_dat, i) %>% filter(is_error == T & !is.na(.)) %>% count %>% pull(n))
		nmr_dat <- mutate_at(nmr_dat, vars(i), list(~ ifelse(is_error == TRUE, NA, .)))		
}  
						
# Remove tag variables
nmr_dat <- select(nmr_dat, -ends_with("_TAG"))

nrow(nmr_dat) # 122266


############################################################################
#   			 Merge to linker and covariable files                      #
############################################################################

# Bridging file between app30418 and NMR IDs
pheno_id_map <- fread("./Linking_code/Nightingale_Phase1_bridge.tsv")

# Bridging file between app30418 and IEU genetic IDs
geno_id_map <- fread(paste0(Sys.getenv('UKB_NMR'), "/data/", "linker.csv"))

# Covariables file 1: fasting time, sex and ever pregnant (generated by Gemma Clayton)
covar <- fread(paste0(repro_dir, "mrbase_variables.csv"))

# Covariables file 2: genotyping array file
array <- fread(paste0(Sys.getenv('UKB_all'), "/software/gwas_pipeline/dev/release_candidate/data/covariates/data.covariates.bolt.txt")) %>%
		select(-sex)

# Merge linker/covariables files to NMR data
nmr_dat_linked <- nmr_dat %>%
					merge(pheno_id_map, ., by = c("sample_id", "plate_id", "plate_position")) %>%
					merge(., covar, by.x = "eid_30418", by.y = "n_eid") %>%
					merge(geno_id_map, ., by.x = "app", by.y = "eid_30418") %>%
					merge(array, ., by.x = "IID", by.y = "ieu") 					

nrow(nmr_dat_linked) # 121667

# Any duplicated IDs left?
table(duplicated(nmr_dat_linked$IID)) # 90

# Remove any remaining duplicated IDs
set.seed(8420)
nmr_dat_linked <- nmr_dat_linked %>%
					group_by(IID) %>%
					sample_n(1) %>%
					ungroup(IID)
table(duplicated(nmr_dat_linked$IID)) # 0
							
# Create variable to stratify according to sex/ever pregnant
nmr_dat_linked <- nmr_dat_linked %>%
					mutate(sex_preg = case_when(
										  sex == "Male" ~ "male",
										  sex == "Female" & ever_pregnant == 0 ~ "female_nopreg",
										  sex == "Female" & ever_pregnant == 1 ~ "female_preg"
										  )) 
				
					
############################################################################
#                            Transform NMR data                            #
############################################################################


# Function to apply rank-based inverse normal transformation (INT) to metabolites
int <- function(x) { 
  
  y <- qnorm((rank(x, na.last="keep") - 0.375) / (sum(!is.na(x))+0.25))

}


##
## Whole dataset

# INT transform all metabolites for whole sample
nmr_dat_all <- nmr_dat_linked %>%
					mutate_at(vars(Total_C:S_HDL_TG_pct), int) %>%
					rename_at(vars(Total_C:S_HDL_TG_pct), ~ paste0(., "_int")) %>%
					select(FID, IID, ends_with("_int"), chip, fasting_t, sex)

# N after merging files:
nrow(nmr_dat_all) # 121577


##
## Stratified dataset (male, never preg female, ever preg female)

# INT transform all metabolites within each stratum (defined by sex and ever pregnant status)
nmr_dat_strat <- nmr_dat_linked %>%
					group_by(sex_preg) %>%
					mutate_at(vars(Total_C:S_HDL_TG_pct), int) %>%
					rename_at(vars(Total_C:S_HDL_TG_pct), ~ paste0(., "_int")) %>%
					ungroup(sex_preg) %>%
					select(FID, IID, ends_with("_int"), chip, fasting_t, sex_preg)
			
# N after merging files:
nrow(nmr_dat_strat) # 121577


##
## Stratified dataset (male vs female) - Added in 11/03/2021 to derive all females data

# INT transform all metabolites within each stratum (defined by sex)
nmr_dat_sexstrat <- nmr_dat_linked %>%
					group_by(sex) %>%
					mutate_at(vars(Total_C:S_HDL_TG_pct), int) %>%
					rename_at(vars(Total_C:S_HDL_TG_pct), ~ paste0(., "_int")) %>%
					ungroup(sex) %>%
					select(FID, IID, ends_with("_int"), chip, fasting_t, sex)
			
# N after merging files:
nrow(nmr_dat_sexstrat) # 121577


############################################################################
#                         Save clean NMR data                              #
############################################################################

##
## Whole dataset
					
write.table(nmr_dat_all, file = paste0(repro_dir, "nmr_dat_all.txt"), quote = F, sep = " ", row.names = F, col.names = T)

