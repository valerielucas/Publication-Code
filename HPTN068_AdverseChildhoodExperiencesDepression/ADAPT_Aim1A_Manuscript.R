
#### Setup ####

# load in packages

## data reading/writing
library(haven)
library(readr)

## data cleaning
library(tidyverse)
library(psych)
library(data.table)
library(stats)

## data visualization
library(table1)
library(ggplot2)
library(grid)
library(gridExtra)

# file path
file_path <- "~/Documents/R/ADAPT Aim 1A/"

# custom analysis functions
source(paste0(file_path, "Manuscript/analysis_funcs.R"))

# suppress scientific notation
options(scipen=999)

# read in data
hptn068_complete <- read_dta(paste0(file_path, "Input Data/Full Dataset/hptn068.complete_201801.dta"))

#### Methods: Exposure Cleaning ####

## new var for double orphan
hptn068_complete$double_orphan <- ifelse(hptn068_complete$g1b8 == 2 & hptn068_complete$g1b5 == 2, 1, 0)

## new var for parent/guardian not caring
hptn068_complete$not_care <- ifelse(hptn068_complete$g1j13 == 1, 1, 0)

## new var for food insecurity 
hptn068_complete$food_insecure <- hptn068_complete$g1g11

## new var for any physical ipv
hptn068_complete$ipv_physical <- ifelse(hptn068_complete$g1d1a == 1 | hptn068_complete$g1d2a == 1 | hptn068_complete$g1d3a == 1 | hptn068_complete$g1d4a == 1 | hptn068_complete$g1d5a == 1 | hptn068_complete$g1d6a == 1, 1, 0)

## new var for any sexual violence
hptn068_complete$sex_violence <- ifelse(hptn068_complete$g1d7a == 1 | hptn068_complete$g1d8a == 1, 1, 0)

## new var for any school violence
hptn068_complete$school_violence <- ifelse(hptn068_complete$g1c17h == 1 | hptn068_complete$g1c17i == 1 | hptn068_complete$g1c17j == 1 | hptn068_complete$g1c17k == 1 | hptn068_complete$g1c17m  == 1, 1, 0)

## new var for sum of ACEs
hptn068_complete$ace_sum <- hptn068_complete$school_violence + hptn068_complete$sex_violence + hptn068_complete$ipv_physical + hptn068_complete$food_insecure + hptn068_complete$not_care + hptn068_complete$double_orphan

## max category is 4+ ACES
hptn068_complete$ace_sum <- ifelse(hptn068_complete$ace_sum > 4, 4, hptn068_complete$ace_sum)



#### Methods: Outcome Cleaning ####

# new var for CDI score

## mark CDI questions in dataframe
cdi_original <- data.frame(select(hptn068_complete, g1o15:g1o24))

## indicate same-coded items with 1 and reverse-coded items with -1
keys <- c(1,-1, 1, -1, -1, -1, 1, 1, 1, -1)

## use `psych` package to reverse-code items
hptn068_complete <- data.frame(hptn068_complete, reverse.code(keys, cdi_original, mini=1, maxi=3))

## sum all CDI questions
hptn068_complete$cdi_sum <- hptn068_complete$g1o15.1 + hptn068_complete$g1o16. + hptn068_complete$g1o17.1 + hptn068_complete$g1o18. + hptn068_complete$g1o19. + hptn068_complete$g1o20. + hptn068_complete$g1o21.1 + hptn068_complete$g1o22.1 + hptn068_complete$g1o23.1 + hptn068_complete$g1o24. - 10

## use typical cutoff for CDI short version
hptn068_complete$cdi_depressed7 <- ifelse(hptn068_complete$cdi_sum < 7, 0, 1)

# new var for CESD score

## mark CESD questions in dataframe
cesd_original <- data.frame(select(hptn068_complete, g1o25:g1o44))

## indicate same-coded items with 1 and reverse-coded items with -1
cesd_keys <- c(1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,1)

## use `psych` package to reverse-code items
hptn068_complete <- data.frame(hptn068_complete, reverse.code(cesd_keys, cesd_original, mini=1, maxi=4))

## sum all CESD questions
hptn068_complete$cesd_sum <- hptn068_complete$g1o25.1 + hptn068_complete$g1o26.1 + hptn068_complete$g1o27.1 + hptn068_complete$g1o28. + hptn068_complete$g1o29.1 + hptn068_complete$g1o30.1 + hptn068_complete$g1o31.1 + hptn068_complete$g1o32. + hptn068_complete$g1o33.1 + hptn068_complete$g1o34.1 + hptn068_complete$g1o35.1 + hptn068_complete$g1o36. + hptn068_complete$g1o37.1 + hptn068_complete$g1o38.1 + hptn068_complete$g1o39.1 + hptn068_complete$g1o40. + hptn068_complete$g1o41.1 + hptn068_complete$g1o42.1 + hptn068_complete$g1o43.1 + hptn068_complete$g1o44.1 - 20

# depression determination CESD ≥ 16
hptn068_complete$depression <- ifelse(hptn068_complete$cesd_sum < 16, 0, 1)

# depression determination CESD ≥ 20
hptn068_complete$depression20 <- ifelse(hptn068_complete$cesd_sum < 20, 0, 1)

#### Methods: Covariate Cleaning ####

## create arm a name variable
hptn068_complete$arm_name <- ifelse(hptn068_complete$arm == 1234, "control", "intervention")

## create age categories 
hptn068_complete$age_cat <- case_when(
  hptn068_complete$yw_age <= 15 ~ "13-15",
  hptn068_complete$yw_age >= 16 & hptn068_complete$yw_age <= 17 ~ "16-17",
  hptn068_complete$yw_age >= 18 ~ "18+",
  .default = NA
)

## rename consumption variables
hptn068_complete$consumption <- hptn068_complete$logpcexptotal

## create consumption categories
hptn068_complete$consumption_cat <- ifelse(hptn068_complete$consumption >= median(hptn068_complete$consumption, na.rm = TRUE), "Above median", "Below median")


## create knots with for spline w/ consumption
knots_consumption <- quantile(hptn068_complete$consumption, c(0.10, 0.33, 0.67, 0.90), na.rm = TRUE)

# making knots in trial population
hptn068_complete <- hptn068_complete |>
  mutate(
    # age spline
    consumption_sq1 = ifelse(consumption > knots_consumption[1], (consumption - knots_consumption[1])^2, 0),
    consumption_sq2 = ifelse(consumption > knots_consumption[2], (consumption - knots_consumption[2])^2, 0),
    consumption_sq3 = ifelse(consumption > knots_consumption[3], (consumption - knots_consumption[3])^2, 0),
    consumption_sq4 = ifelse(consumption > knots_consumption[4], (consumption - knots_consumption[4])^2, 0),
    consumption_rs1 = ifelse(consumption > knots_consumption[1], consumption_sq4 - consumption_sq1, 0),
    consumption_rs2 = ifelse(consumption > knots_consumption[2], consumption_sq4 - consumption_sq2, 0),
    consumption_rs3 = ifelse(consumption > knots_consumption[3], consumption_sq4 - consumption_sq3, 0)
  )



#### Methods: Dataset Cleaning ####

# filter out extraneous biological data (Visit 901) and remove those 18+ at enrollment
hptn068_filtered <- hptn068_complete |>
  filter(visit != 901) |>
  group_by(uid) |>
  mutate(age_enrollment = min(yw_age, na.rm = TRUE)) |>
  filter(age_enrollment <= 17)

# make new dataframe with select variables
adapt_visit201 <- hptn068_filtered  |>
  filter(visit == 201) |>
  rename(grade_level = g1c5) |>
  select(uid, visit, arm_name, double_orphan, not_care, food_insecure, ipv_physical, sex_violence, school_violence, someprimaryorless_mother, consumption, consumption_rs1, consumption_rs2, consumption_rs3, yw_age, grade_level, cdi_depressed7, ace_sum, age_cat, consumption_cat)

# complete case analysis at baseline
adapt_visit201_nona <- adapt_visit201 |>
  drop_na(-cdi_depressed7)

# save percentage missing
percent_missing_exposure <- 1 - (nrow(adapt_visit201_nona) / nrow(adapt_visit201))

# select only depression information
adapt_depression <- hptn068_complete |>
  filter(visit != 901) |>
  select(uid, visit, depression, depression20) 


#### Methods: Figure 1, HPTN 068 follow-up visits ####

# print number of follow-up visits
visit_followup <- data.frame(table(hptn068_filtered[hptn068_filtered$uid %in% adapt_visit201_nona$uid,]$visit))

# output the code
write.csv(visit_followup , paste0(file_path, "Output Data/fig01_visit_followup.csv"))

#### Methods: G-computation model definitons ####

## Define g-computation models
gcomp_model_double_orphan <- depression ~ arm_name + double_orphan + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + yw_age*double_orphan + as.factor(visit)*yw_age + as.factor(visit)*double_orphan
gcomp_model_not_care <- depression ~ arm_name + not_care + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + yw_age*not_care + as.factor(visit)*yw_age + as.factor(visit)*not_care 
gcomp_model_food_insecure <- depression ~ arm_name + food_insecure + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + yw_age*food_insecure + as.factor(visit)*yw_age + as.factor(visit)*food_insecure
gcomp_model_ipv_physical <- depression ~ arm_name + ipv_physical + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + yw_age*ipv_physical + as.factor(visit)*yw_age + as.factor(visit)*ipv_physical
gcomp_model_sex_violence <- depression ~ arm_name + sex_violence + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + yw_age*sex_violence + as.factor(visit)*yw_age + as.factor(visit)*sex_violence 
gcomp_model_school_violence <- depression ~ arm_name + school_violence + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + yw_age*school_violence + as.factor(visit)*yw_age + as.factor(visit)*school_violence
gcomp_model_ipv_sexviolence <- depression ~ arm_name + sex_violence + ipv_physical + ipv_physical*sex_violence + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + yw_age*sex_violence + yw_age*ipv_physical + as.factor(visit)*yw_age + as.factor(visit)*sex_violence + as.factor(visit)*ipv_physical

## Create list of g-computation models
gcomp_model_list <- list(gcomp_model_double_orphan, gcomp_model_not_care, gcomp_model_food_insecure, gcomp_model_ipv_physical, gcomp_model_sex_violence, gcomp_model_school_violence, gcomp_model_ipv_sexviolence)


#### Results: Table 1, Sample description ####

# display table 1 for people who are ≤17 at enrollment
table1(~ factor(double_orphan) + factor(not_care) + factor(food_insecure) + factor(ipv_physical) + factor(sex_violence) + factor(school_violence) + factor(ace_sum) + factor(someprimaryorless_mother) + consumption + yw_age + factor(yw_age) + factor(arm_name) + factor(grade_level), data = adapt_visit201)





#### Results: Figure 2, Depression prevalence curves ####

# name ACEs
aces <- c("double_orphan", "not_care", "food_insecure", "ipv_physical", "sex_violence", "school_violence")
ace_names <- c("double orphan", "low parental care", "food insecurity", "physical IPV", "sexual violence", "school violence")

adapt_results <- data.frame(ace = character(0), visit_number = numeric(0), prev_diff = numeric(0))

# Construct dataset
adapt_analytic <- dataset_construction(uid_list = unique(adapt_visit201_nona$uid), adapt_visit201 = adapt_visit201_nona, adapt_depression = adapt_depression)

# remove visit 201
adapt_analytic <- adapt_analytic |>
  filter(visit != 201)

# imputation model
impute_missing <- depression ~ arm_name + school_violence + sex_violence + ipv_physical + food_insecure + not_care + double_orphan + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + as.factor(visit)*yw_age + yw_age*(school_violence + sex_violence + ipv_physical + food_insecure + not_care + double_orphan) + as.factor(visit)*(school_violence + sex_violence + ipv_physical + food_insecure + not_care + double_orphan)


# fit imputation model
impute_model <- glm(impute_missing, family = binomial(link = "logit"), data = adapt_analytic)

# make predictions for depression probability for each individual
adapt_analytic$probability <- predict(impute_model, type = "response", newdata = adapt_analytic)

# reframe to calculate depression prevalence, using imputed values for participants who are missing depression information
impute_results <- adapt_analytic |>
  mutate(depression_impute = ifelse(is.na(depression), probability, depression)) |>
  group_by(visit_number) |>
  reframe(prev_impute = mean(depression_impute))

# set number of bootstrap reps
reps <- 500

# create shell for bootstrap results
impute_results_boot <- data.frame(visit_number = numeric(0), prev_impute = numeric(0))

# set seed for bootstrap
set.seed(1939)

# bootstrap
for(i in 1:reps){
  
  # sample UIDs with replacement
  uid_list <- sample(x = adapt_visit201_nona$uid, size = length(adapt_visit201_nona$uid), replace = TRUE)
  
  # create dataset with sampled UIDs
  adapt_analytic_boot <- dataset_construction(uid_list = uid_list, adapt_visit201 = adapt_visit201_nona, adapt_depression = adapt_depression)
  
  # remove visit 201
  adapt_analytic_boot <- adapt_analytic_boot |>
    filter(visit != 201)
  
  # fit imputation model
  impute_model_boot <- glm(impute_missing, family = binomial(link = "logit"), data = adapt_analytic_boot)
  
  # make predictions for depression probability for each individual
  adapt_analytic_boot$probability <- predict(impute_model, type = "response", newdata = adapt_analytic_boot)
  
  # reframe to calculate depression prevalence, using imputed values for participants who are missing depression information
  impute_results_boot_loading <- adapt_analytic_boot |>
    mutate(depression_impute = ifelse(is.na(depression), probability, depression)) |>
    group_by(visit_number) |>
    reframe(prev_impute = mean(depression_impute))
  
  # bind the results for this ACE to the overall results
  impute_results_boot <- rbind(impute_results_boot, impute_results_boot_loading)
}


# take the standard deviation of the bootstrap estimates
impute_results_boot_se <- impute_results_boot |>
  group_by(visit_number) |>
  reframe(se = sd(prev_impute))

# calculate lower and upper confidence limits based on bootstraped
prevalence_impute <- merge(impute_results, impute_results_boot_se, by = "visit_number") |>
  mutate(lower_prevalence = prev_impute - 1.96*se,
         upper_prevalence = prev_impute + 1.96*se)

# graph crude prevalence
prevalence_crude <- adapt_analytic |>
  group_by(visit_number) |>
  reframe(n_nonmissing = sum(!is.na(depression)),
          n_missing = sum(is.na(depression)),
          n_depressed = sum(depression, na.rm = TRUE),
          n_notdepressed = sum((1 - depression), na.rm = TRUE),
          prevalence = mean(depression, na.rm = TRUE),
          se = sqrt(prevalence*(1-prevalence)/n_nonmissing),
          lower_prevalence = prevalence - 1.96*se,
          upper_prevalence = prevalence + 1.96*se)

# plot crude and imputed depression prevalence with 95% confidence intervals
adapt_depression_plot <- ggplot() +
  geom_line(aes(x = visit_number, y = prevalence, color = "Crude"), data = prevalence_crude) +
  geom_point(aes(x = visit_number, y = prevalence, color = "Crude"), data = prevalence_crude) +
  geom_ribbon(aes(x = visit_number, ymin = lower_prevalence, ymax = upper_prevalence), fill = "#DDCC77", alpha = 0.25, data = prevalence_crude) +
  geom_line(aes(x = visit_number, y = prev_impute, color = "Adjusted"), data = prevalence_impute) +
  geom_point(aes(x = visit_number, y = prev_impute, color = "Adjusted"), data = prevalence_impute) +
  geom_ribbon(aes(x = visit_number, ymin = lower_prevalence, ymax = upper_prevalence), fill = "#88CCEE", alpha = 0.25, data = prevalence_impute) +
  xlab(NULL) +
  ylab("Depression prevalence") +
  scale_x_continuous(breaks = c(1:3,5,7), labels = c("Year 1", "Year 2", "Year 3", "Year 5", "Year 7")) +
  scale_color_manual("Key", values = c("Crude"="#DDCC77", "Adjusted"="#88CCEE")) +
  ylim(0, 0.5) +
  theme_classic()

# save image to PDF
ggsave(file = paste0(file_path, "Charts/fig02_adapt_depression_prev.pdf"), adapt_depression_plot, width = 7, height = 3.5, units = "in")


#### Results: Table 3 and Figure 2 + Appendix 1: Figure 12, Primary prevalence, prevalence difference, and prevalence ratio results ####

# reset analytic dataset
adapt_analytic_main <- dataset_construction(uid_list = adapt_visit201_nona$uid, adapt_visit201_nona, adapt_depression)

# remove baseline visit (201), no outcomes reported
adapt_analytic_main <- adapt_analytic_main |>
  filter(visit != 201)

# calculate point estimates
adapt_results_main <- point_estimates_gcomp(adapt_no201 = adapt_analytic_main, aces = aces, gcomp_model_list = gcomp_model_list)

# set seed for bootstrap
set.seed(1939)

# bootstrap analysis
adapt_results_ci_main <- bootstrap_analysis(reps = reps, adapt_visit201 = adapt_visit201_nona, adapt_depression = adapt_depression, aces, analysis = point_estimates_gcomp, point_estimates = adapt_results_main, gcomp_model_list = gcomp_model_list)

# save chart and CSV file
results_report(adapt_results_ci = adapt_results_ci_main, file_path = file_path, label = "main")

#### Appendix 1: Table 1, Composite ACEs ####

# extract IDs final
final_uids <- unique(adapt_visit201_nona$uid)

# make table of breakdown of composite ACEs
hptn068_complete_composite <- hptn068_complete |>
  filter(visit == 201) |> 
  filter(uid %in% final_uids) |>
  reframe(ipv_physical = sum(ipv_physical),
          slapped_thrown = sum(g1d1a, na.rm = TRUE), 
          pushed_shoved = sum(g1d2a, na.rm = TRUE),
          hit = sum(g1d3a, na.rm = TRUE),
          kicked_dragged = sum(g1d4a, na.rm = TRUE),
          choked_burnt = sum(g1d5a, na.rm = TRUE),
          threatened_weapon = sum(g1d6a, na.rm = TRUE),
          sex_violence = sum(sex_violence),
          not_wanted = sum(g1d7a, na.rm = TRUE),
          fear = sum(g1d8a, na.rm = TRUE),
          school_violence = sum(school_violence),
          drug_dealing = sum(g1c17h, na.rm = TRUE),
          unsafe = sum(g1c17i, na.rm = TRUE),
          sexual_harrassment_student = sum(g1c17j, na.rm = TRUE),
          sexual_harrassment_staff = sum(g1c17k, na.rm = TRUE),
          violent_teachers = sum(g1c17m, na.rm = TRUE)) |>
  data.table::transpose(keep.names = "ACEs")

# save to CSV file
write_csv(hptn068_complete_composite, paste0(file_path, "Output Data/hptn068_complete_composite.csv"))

#### Appendix 1: Table 2, Correlation matrix ####

## filter to include only participants age ≤17 at enrollment then select only ACE variables
adapt_cor <- adapt_visit201_nona |>
  filter(yw_age <= 17) |>
  select(double_orphan, not_care, food_insecure, ipv_physical, sex_violence, school_violence)

## generate correlation matrix
cor_matrix <- cor(adapt_cor, method = "pearson")

## round correlation coeffs, then convert matrix to data frame
correlation_ace <- data.frame(round(cor_matrix, 3))

## write correlation matrix to .csv file
write_csv(correlation_ace, paste0(file_path, "Output Data/ace_covariance.csv"))

#### Appendix 1: Table 3, Crude and imputed depression prevalence ####

# merge together crude data and imputed data
adapt_depression_table_draft <- merge(x = prevalence_crude, y = prevalence_impute, by = "visit_number", suffixes = c("_crude","_adjusted"))

# create output table 
adapt_depression_table <- adapt_depression_table_draft |>
  mutate(prev_crude = paste0(round(prevalence*100, digits = 1), "% (", round(lower_prevalence_crude*100, digits = 1), "%, ", round(upper_prevalence_crude*100, digits = 1), "%)"),
         prev_adjusted = paste0(round(prev_impute*100, digits = 1), "% (", round(lower_prevalence_adjusted*100, digits = 1), "%, ", round(upper_prevalence_adjusted*100, digits = 1), "%)")) |>
  select(visit_number, n_notdepressed,  n_depressed, n_missing, prev_crude, prev_adjusted)

# write to csv
write_csv(adapt_depression_table, paste0(file_path, "Output Data/adapt_depression_prev.csv"))

#### Appendix 1: Table 4, Attrition by ACE ####

# construct analytic dataset
adapt_analytic <- dataset_construction(uid_list = unique(adapt_visit201_nona$uid), adapt_visit201 = adapt_visit201_nona, adapt_depression = adapt_depression)

# make attrition table for each visit 
adapt_attrition <- adapt_analytic |>
  group_by(visit_number) |>
  reframe(double_orphan_na = "",
          double_orphan_yes = sum(double_orphan*(!is.na(depression) | visit_number == 0), na.rm = TRUE),
          double_orphan_no = sum((1-double_orphan)*!is.na(depression | visit_number == 0), na.rm = TRUE),
          not_care_na = "",
          not_care_yes = sum(not_care*(!is.na(depression) | visit_number == 0), na.rm = TRUE),
          not_care_no = sum((1-not_care)*!is.na(depression | visit_number == 0), na.rm = TRUE),
          food_insecure_na = "",
          food_insecure_yes = sum(food_insecure*(!is.na(depression) | visit_number == 0), na.rm = TRUE),
          food_insecure_no = sum((1-food_insecure)*!is.na(depression | visit_number == 0), na.rm = TRUE),
          ipv_physical_na = "",
          ipv_physical_yes = sum(ipv_physical*(!is.na(depression) | visit_number == 0), na.rm = TRUE),
          ipv_physical_no = sum((1-ipv_physical)*!is.na(depression | visit_number == 0), na.rm = TRUE),
          sex_violence_na = "",
          sex_violence_yes = sum(sex_violence*(!is.na(depression) | visit_number == 0), na.rm = TRUE),
          sex_violence_no = sum((1-sex_violence)*!is.na(depression | visit_number == 0), na.rm = TRUE),
          school_violence_na = "",
          school_violence_yes = sum(school_violence*(!is.na(depression) | visit_number == 0), na.rm = TRUE),
          school_violence_no = sum((1-school_violence)*!is.na(depression | visit_number == 0), na.rm = TRUE)
  ) |>
  mutate(double_orphan_yes = paste0(double_orphan_yes, " (", round(double_orphan_yes / max(double_orphan_yes)*100), "%)"),
         double_orphan_no = paste0(double_orphan_no, " (", round(double_orphan_no / max(double_orphan_no)*100), "%)"),
         not_care_yes = paste0(not_care_yes, " (", round(not_care_yes / max(not_care_yes)*100), "%)"),
         not_care_no = paste0(not_care_no, " (", round(not_care_no / max(not_care_no)*100), "%)"),
         food_insecure_yes = paste0(food_insecure_yes, " (", round(food_insecure_yes / max(food_insecure_yes)*100), "%)"),
         food_insecure_no = paste0(food_insecure_no, " (", round(food_insecure_no / max(food_insecure_no)*100), "%)"),
         ipv_physical_yes = paste0(ipv_physical_yes, " (", round(ipv_physical_yes / max(ipv_physical_yes)*100), "%)"),
         ipv_physical_no = paste0(ipv_physical_no, " (", round(ipv_physical_no / max(ipv_physical_no)*100), "%)"),
         sex_violence_yes = paste0(sex_violence_yes, " (", round(sex_violence_yes / max(sex_violence_yes)*100), "%)"),
         sex_violence_no = paste0(sex_violence_no, " (", round(sex_violence_no / max(sex_violence_no)*100), "%)"),
         school_violence_yes = paste0(school_violence_yes, " (", round(school_violence_yes / max(school_violence_yes)*100), "%)"),
         school_violence_no = paste0(school_violence_no, " (", round(school_violence_no / max(school_violence_no)*100), "%)")) |>
  data.table::transpose(keep.names = "ACE", make.names = "visit_number")

# save to .csv file
write.csv(adapt_attrition, paste0(file_path, "Output Data/adapt_attrition.csv"))

#### Appendix 1: Table 5 and Figure 1: Results limited to those aged 13-15 at baseline ####

# subset data to include only those 13-15 at baseline
adapt_visit201_age1315 <- adapt_visit201_nona[adapt_visit201_nona$yw_age <= 15,]

# construct new dataset 
adapt_analytic_age1315 <- dataset_construction(uid_list = adapt_visit201_age1315$uid, adapt_visit201_age1315, adapt_depression)

# remove baseline visit (201), no outcomes reported
adapt_analytic_age1315 <- adapt_analytic_age1315 |>
  filter(visit != 201)

# calculate point estimates
adapt_results_age1315 <- point_estimates_gcomp(adapt_no201 = adapt_analytic_age1315, aces = aces, gcomp_model_list = gcomp_model_list)

# set seed for bootstrap
set.seed(1939)

# bootstrap analysis
adapt_results_ci_age1315 <- bootstrap_analysis(reps = reps, adapt_visit201 = adapt_visit201_age1315, adapt_depression = adapt_depression, aces, analysis = point_estimates_gcomp, point_estimates = adapt_results_age1315, gcomp_model_list = gcomp_model_list)

# save chart and CSV file
results_report(adapt_results_ci = adapt_results_ci_age1315, file_path = file_path, label = "age1315")



#### Appendix 1: Table 6 and Figure 2: Results limited to those aged 16-17 at baseline ####

# subset data to include only those 16-17 at baseline
adapt_visit201_age1617 <- adapt_visit201_nona[adapt_visit201_nona$yw_age >= 16,]

# construct new dataset 
adapt_analytic_age1617 <- dataset_construction(uid_list = adapt_visit201_age1617$uid, adapt_visit201_age1617, adapt_depression)

# remove baseline visit (201), no outcomes reported
adapt_analytic_age1617 <- adapt_analytic_age1617 |>
  filter(visit != 201)

# calculate point estimates
adapt_results_age1617 <- point_estimates_gcomp(adapt_no201 = adapt_analytic_age1617, aces = aces, gcomp_model_list = gcomp_model_list)

# set seed for bootstrap
set.seed(1939)

# bootstrap analysis
adapt_results_ci_age1617 <- bootstrap_analysis(reps = reps, adapt_visit201 = adapt_visit201_age1617, adapt_depression = adapt_depression, aces, analysis = point_estimates_gcomp, point_estimates = adapt_results_age1617, gcomp_model_list = gcomp_model_list)

# save chart and CSV file
results_report(adapt_results_ci = adapt_results_ci_age1617, file_path = file_path, label = "age1617")



#### Appendix 1: Table 7 and Figure 3: Results with those 18+ at baseline included ####

# filter out extraneous biological data (Visit 901) and remove those 18+ at enrollment
adapt_18included <- hptn068_complete |>
  filter(visit != 901) |>
  group_by(uid) |>
  mutate(age_enrollment = min(yw_age, na.rm = TRUE))

# make new dataframe with select variables
adapt_visit201_18included <- adapt_18included |>
  filter(visit == 201) |>
  select(uid, visit, arm_name, double_orphan, not_care, food_insecure, ipv_physical, sex_violence, school_violence, someprimaryorless_mother, consumption, consumption_rs1, consumption_rs2, consumption_rs3, yw_age,  cdi_depressed7, ace_sum)

# complete case analysis at baseline
adapt_visit201_18included_nona <- na.omit(adapt_visit201_18included)

# construct new dataset 
adapt_analytic_18included <- dataset_construction(uid_list = adapt_visit201_18included_nona$uid, adapt_visit201_18included_nona, adapt_depression)

# remove baseline visit (201), no outcomes reported
adapt_analytic_18included <- adapt_analytic_18included |>
  filter(visit != 201)

# calculate point estimates
adapt_results_18included <- point_estimates_gcomp(adapt_no201 = adapt_analytic_18included, aces = aces, gcomp_model_list = gcomp_model_list)

# set seed for bootstrap
set.seed(1939)

# bootstrap analysis
adapt_results_ci_18included <- bootstrap_analysis(reps = reps, adapt_visit201 = adapt_visit201_18included_nona, adapt_depression = adapt_depression, aces, analysis = point_estimates_gcomp, point_estimates = adapt_results_18included, gcomp_model_list = gcomp_model_list)

# save chart and CSV file
results_report(adapt_results_ci = adapt_results_ci_18included, file_path = file_path, label = "included18")



#### Appendix 1: Table 8 and Figure 4: Results limited to those not depressed at baseline ####

# subset data to include only those not depressed at baseline
adapt_visit201_notdepressed <- adapt_visit201_nona |>
  filter(cdi_depressed7 == 0)

# construct new dataset
adapt_analytic_notdepressed <- dataset_construction(uid_list = adapt_visit201_notdepressed$uid, adapt_visit201_notdepressed, adapt_depression)

# remove baseline visit (201), no outcomes reported
adapt_analytic_notdepressed <- adapt_analytic_notdepressed |>
  filter(visit != 201)

# calculate point estimates
adapt_results_notdepressed <- point_estimates_gcomp(adapt_no201 = adapt_analytic_notdepressed, aces = aces, gcomp_model_list = gcomp_model_list)

# set seed for bootstrap
set.seed(1939)

# bootstrap analysis
adapt_results_ci_notdepressed <- bootstrap_analysis(reps = reps, adapt_visit201 = adapt_visit201_notdepressed, adapt_depression = adapt_depression, aces, analysis = point_estimates_gcomp, point_estimates = adapt_results_notdepressed, gcomp_model_list = gcomp_model_list)

# save chart and CSV file
results_report(adapt_results_ci = adapt_results_ci_notdepressed, file_path = file_path, label = "notdepressed")



#### Appendix 1: Table 9 and Figure 5: Results for CESD ≥ 20 outcome definition ####

# reset depression variable to equal the CESD ≥ 20 cutoff depression variable
adapt_depression_cesd20 <- adapt_depression |>
  mutate(depression = depression20)

# reset analytic dataset
adapt_analytic_cesd20 <- dataset_construction(uid_list = adapt_visit201_nona$uid, adapt_visit201_nona, adapt_depression_cesd20)

# remove baseline visit (201), no outcomes reported
adapt_analytic_cesd20 <- adapt_analytic_cesd20 |>
  filter(visit != 201)

# calculate point estimates
adapt_results_cesd20 <- point_estimates_gcomp(adapt_no201 = adapt_analytic_cesd20, aces = aces, gcomp_model_list = gcomp_model_list)

# set seed for bootstrap
set.seed(1939)

# bootstrap analysis
adapt_results_ci_cesd20 <- bootstrap_analysis(reps = reps, adapt_visit201 = adapt_visit201_nona, adapt_depression = adapt_depression_cesd20, aces, analysis = point_estimates_gcomp, point_estimates = adapt_results_cesd20, gcomp_model_list = gcomp_model_list)

# save chart and CSV file
results_report(adapt_results_ci = adapt_results_ci_cesd20, file_path = file_path, label = "cesd20")



#### Appendix 1: Table 10 and Figure 6: Age and cash transfer-only adjustment set ####

## Develop simpler g-computation models with only age and randomized conditional cash transfer (arm) g-computation models
gcomp_model_double_orphan <- depression ~ arm_name + double_orphan + yw_age + yw_age^2 + as.factor(visit) + yw_age*double_orphan + as.factor(visit)*yw_age + as.factor(visit)*double_orphan
gcomp_model_not_care <- depression ~ arm_name + not_care + yw_age + yw_age^2 + as.factor(visit) + yw_age*not_care + as.factor(visit)*yw_age + as.factor(visit)*not_care 
gcomp_model_food_insecure <- depression ~ arm_name + food_insecure + yw_age + yw_age^2 + as.factor(visit) + yw_age*food_insecure + as.factor(visit)*yw_age + as.factor(visit)*food_insecure
gcomp_model_ipv_physical <- depression ~ arm_name + ipv_physical + yw_age + yw_age^2 + as.factor(visit) + yw_age*ipv_physical + as.factor(visit)*yw_age + as.factor(visit)*ipv_physical
gcomp_model_sex_violence <- depression ~ arm_name + sex_violence + yw_age + yw_age^2 + as.factor(visit) + yw_age*sex_violence + as.factor(visit)*yw_age + as.factor(visit)*sex_violence 
gcomp_model_school_violence <- depression ~ arm_name + school_violence + yw_age + yw_age^2 + as.factor(visit) + yw_age*school_violence + as.factor(visit)*yw_age + as.factor(visit)*school_violence
gcomp_model_ipv_sexviolence <- depression ~ arm_name + sex_violence + ipv_physical + ipv_physical*sex_violence + yw_age + yw_age^2 + as.factor(visit) + yw_age*sex_violence + yw_age*ipv_physical + as.factor(visit)*yw_age + as.factor(visit)*sex_violence + as.factor(visit)*ipv_physical

# develop model list
gcomp_simple_model_list <- list(gcomp_model_double_orphan, gcomp_model_not_care, gcomp_model_food_insecure, gcomp_model_ipv_physical, gcomp_model_sex_violence, gcomp_model_school_violence, gcomp_model_ipv_sexviolence)

# reset analytic dataset
adapt_analytic_simple <- dataset_construction(uid_list = adapt_visit201_nona$uid, adapt_visit201_nona, adapt_depression)

# remove baseline visit (201), no outcomes reported
adapt_analytic_simple <- adapt_analytic_simple |>
  filter(visit != 201)

# calculate point estimates
adapt_results_simple <- point_estimates_gcomp(adapt_no201 = adapt_analytic_simple, aces = aces, gcomp_model_list = gcomp_simple_model_list)

# set seed for bootstrap
set.seed(1939)

# bootstrap analysis
adapt_results_ci_simple <- bootstrap_analysis(reps = reps, adapt_visit201 = adapt_visit201_nona, adapt_depression = adapt_depression, aces, analysis = point_estimates_gcomp, point_estimates = adapt_results_simple, gcomp_model_list = gcomp_simple_model_list)

# save chart and CSV file
results_report(adapt_results_ci = adapt_results_ci_simple, file_path = file_path, label = "simple")


#### Appendix 1: Table 11 and Figure 7: Inverse probability weighting analysis ####

## Define IPTW denominator models
iptw_double_orphan <- double_orphan ~ arm_name + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + as.factor(visit)*yw_age 

iptw_not_care <- not_care ~ arm_name + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + as.factor(visit)*yw_age 

iptw_food_insecure <- food_insecure ~ arm_name + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + as.factor(visit)*yw_age 

iptw_ipv_physical <- ipv_physical ~ arm_name + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + as.factor(visit)*yw_age 

iptw_sex_violence <- sex_violence ~ arm_name + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 +  yw_age + yw_age^2 + as.factor(visit) + as.factor(visit)*yw_age 

iptw_school_violence <- school_violence ~ arm_name + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + as.factor(visit)*yw_age

iptw_ipv_sexviolence <- sex_violence ~ arm_name + ipv_physical + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + yw_age*ipv_physical + as.factor(visit)*yw_age + as.factor(visit)*ipv_physical 

## Create list of IPTW denominator models
iptw_model_list <- list(iptw_double_orphan, iptw_not_care, iptw_food_insecure, iptw_ipv_physical, iptw_sex_violence, iptw_school_violence, iptw_ipv_sexviolence)

# Define IPTW numerator models (for stabilized weights)
iptw_double_orphan_st <- double_orphan ~ 1
iptw_not_care_st <- not_care ~ 1
iptw_food_insecure_st <- food_insecure ~ 1
iptw_ipv_physical_st <- ipv_physical ~ 1
iptw_sex_violence_st <- sex_violence ~ 1
iptw_school_violence_st <- school_violence ~ 1

## Create list of IPTW numerator models
iptw_model_list_st <- list(iptw_double_orphan_st, iptw_not_care_st, iptw_food_insecure_st, iptw_ipv_physical_st, iptw_sex_violence_st, iptw_school_violence_st)

# Define IPCW denominator models
ipcw_double_orphan <- depression_nonmissing ~ arm_name + double_orphan + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + yw_age*double_orphan + as.factor(visit)*yw_age + as.factor(visit)*double_orphan

ipcw_not_care <- depression_nonmissing ~ arm_name + not_care + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + yw_age*not_care + as.factor(visit)*yw_age + as.factor(visit)*not_care 

ipcw_food_insecure <- depression_nonmissing ~ arm_name + food_insecure + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + yw_age*food_insecure + as.factor(visit)*yw_age + as.factor(visit)*food_insecure

ipcw_ipv_physical <- depression_nonmissing ~ arm_name + ipv_physical + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + yw_age*ipv_physical + as.factor(visit)*yw_age + as.factor(visit)*ipv_physical
ipcw_sex_violence <- depression_nonmissing ~ arm_name + sex_violence + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + yw_age*sex_violence + as.factor(visit)*yw_age + as.factor(visit)*sex_violence 
ipcw_school_violence <- depression_nonmissing ~ arm_name + school_violence + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 + yw_age + yw_age^2 + as.factor(visit) + yw_age*school_violence + as.factor(visit)*yw_age + as.factor(visit)*school_violence
ipcw_ipv_sexviolence <- depression_nonmissing ~ arm_name + sex_violence + ipv_physical + ipv_physical*sex_violence + someprimaryorless_mother + consumption + consumption_rs1 + consumption_rs2 + consumption_rs3 +  yw_age + yw_age^2 + as.factor(visit) + yw_age*sex_violence + yw_age*ipv_physical + as.factor(visit)*yw_age + as.factor(visit)*sex_violence + as.factor(visit)*ipv_physical

## Create list of IPCW denominator models
ipcw_model_list <- list(ipcw_double_orphan, ipcw_not_care, ipcw_food_insecure, ipcw_ipv_physical, ipcw_sex_violence, ipcw_school_violence, ipcw_ipv_sexviolence)

## Define IPCW numerator model (for stabilized weight)
ipcw_st <- depression_nonmissing ~ 1

# reset analytic dataset
adapt_analytic_ipw <- dataset_construction(uid_list = adapt_visit201_nona$uid, adapt_visit201_nona, adapt_depression)

# remove baseline visit (201), no outcomes reported
adapt_analytic_ipw <- adapt_analytic_ipw |>
  filter(visit != 201)

# calculate point estimates
adapt_results_ipw <- point_estimates_ipw(adapt_no201 = adapt_analytic_main, aces = aces, iptw_model_list = iptw_model_list, iptw_model_list_st = iptw_model_list_st, ipcw_model_list = ipcw_model_list, ipcw_st = ipcw_st)

# set seed for bootstrap
set.seed(1939)

# bootstrap analysis
adapt_results_ci_ipw <- bootstrap_analysis(reps = reps, adapt_visit201 = adapt_visit201_nona, adapt_depression = adapt_depression, aces, analysis = point_estimates_ipw, point_estimates = adapt_results_ipw, iptw_model_list = iptw_model_list, iptw_model_list_st = iptw_model_list_st, ipcw_model_list = ipcw_model_list, ipcw_st = ipcw_st)

# save chart and CSV file
results_report(adapt_results_ci = adapt_results_ci_ipw, file_path = file_path, label = "ipw")

#### Appendix 1: Table 13, E-values for prevalence ratio point estimates ####

# calculate e-value of risk ratio estimate
adapt_report_evalue <- adapt_results_ci_main |>
  mutate(prev_ratio = ifelse(prev_ratio > 1, prev_ratio, 1 / prev_ratio)) |>
  mutate(e_value = prev_ratio + sqrt(prev_ratio * (prev_ratio - 1))) |>
  mutate(e_value_report = round(e_value, digits = 2)) |>
  select(visit_number, ace, e_value_report) |>
  pivot_wider(names_from = visit_number, values_from = e_value_report) |>
  arrange(factor(ace, levels = c("double_orphan", "not_care", "food_insecure", "ipv_physical", "sex_violence", "ipv_sexual_violence", "school_violence")))

write_csv(adapt_report_evalue, paste0(file_path, "Output Data/e-value_rr.csv"))



#### Appendix 1: Table 14, E-values for lower confidence limits of prevalence ratio point estimates ####

# calculate e-value of lower confidence interval limit
adapt_report_evalue_lower <- adapt_results_ci_main |>
  mutate(prev_ratio_lower = ifelse(prev_ratio_lower > 1, prev_ratio_lower, 1 / prev_ratio_lower)) |>
  mutate(e_value_lower = prev_ratio_lower + sqrt(prev_ratio_lower * (prev_ratio_lower - 1))) |>
  mutate(e_value_lower_report = round(e_value_lower, digits = 2)) |>
  select(visit_number, ace, e_value_lower_report) |>
  pivot_wider(names_from = visit_number, values_from = e_value_lower_report) |>
  arrange(factor(ace, levels = c("double_orphan", "not_care", "food_insecure", "ipv_physical", "sex_violence", "ipv_sexual_violence", "school_violence")))


write_csv(adapt_report_evalue_lower, paste0(file_path, "Output Data/e-value_rr_lower.csv"))

#### Appendix 1: Table 15, RR for confounders and ACEs

# Calculate risk ratio between maternal education and each ACE
confounder_exposure_table_some <- adapt_visit201_nona |>
  group_by() |>
  reframe(confounder = "maternal education",
          double_orphan = sum(double_orphan*someprimaryorless_mother)/sum(someprimaryorless_mother) / (sum(double_orphan*(1 - someprimaryorless_mother))/sum(1 - someprimaryorless_mother)),
          not_care = (sum(not_care*someprimaryorless_mother)/sum(someprimaryorless_mother)) / (sum(not_care*(1 - someprimaryorless_mother))/sum(1 - someprimaryorless_mother)),
          food_insecure = (sum(food_insecure*someprimaryorless_mother)/sum(someprimaryorless_mother)) / (sum(food_insecure*(1 - someprimaryorless_mother))/sum(1 - someprimaryorless_mother)),
          ipv_physical = (sum(ipv_physical*someprimaryorless_mother)/sum(someprimaryorless_mother)) / (sum(ipv_physical*(1 - someprimaryorless_mother))/sum(1 - someprimaryorless_mother)),
          sex_violence = (sum(sex_violence*someprimaryorless_mother)/sum(someprimaryorless_mother)) / (sum(sex_violence*(1 - someprimaryorless_mother))/sum(1 - someprimaryorless_mother)),
          school_violence = (sum(school_violence*someprimaryorless_mother)/sum(someprimaryorless_mother)) / (sum(school_violence*(1 - someprimaryorless_mother))/sum(1 - someprimaryorless_mother))
  )

# Calculate risk ratio between age at enrollment (16-17 vs. 13-15) and each ACE
confounder_exposure_table_age <- adapt_visit201_nona |>
  group_by() |>
  reframe(confounder = "enrollment age",
          double_orphan = (sum(double_orphan*(age_cat == "16-17"))/sum((age_cat == "16-17"))) / (sum(double_orphan*(1 - (age_cat == "16-17")))/sum(1 - (age_cat == "16-17"))),
          not_care = (sum(not_care*(age_cat == "16-17"))/sum((age_cat == "16-17"))) / (sum(not_care*(1 - (age_cat == "16-17")))/sum(1 - (age_cat == "16-17"))),
          food_insecure = (sum(food_insecure*(age_cat == "16-17"))/sum((age_cat == "16-17"))) / (sum(food_insecure*(1 - (age_cat == "16-17")))/sum(1 - (age_cat == "16-17"))),
          ipv_physical = (sum(ipv_physical*(age_cat == "16-17"))/sum((age_cat == "16-17"))) / (sum(ipv_physical*(1 - (age_cat == "16-17")))/sum(1 - (age_cat == "16-17"))),
          sex_violence = (sum(sex_violence*(age_cat == "16-17"))/sum((age_cat == "16-17"))) / (sum(sex_violence*(1 - (age_cat == "16-17")))/sum(1 - (age_cat == "16-17"))),
          school_violence = (sum(school_violence*(age_cat == "16-17"))/sum((age_cat == "16-17"))) / (sum(school_violence*(1 - (age_cat == "16-17")))/sum(1 - (age_cat == "16-17")))
  )

# Calculate risk ratio between consumption (below vs. above median) and each ACE
confounder_exposure_table_consumption <- adapt_visit201_nona |>
  group_by() |>
  reframe(confounder = "consumption",
          double_orphan = (sum(double_orphan*(consumption_cat == "Below median"))/sum((consumption_cat == "Below median"))) / (sum(double_orphan*(1 - (consumption_cat == "Below median")))/sum(1 - (consumption_cat == "Below median"))),
          not_care = (sum(not_care*(consumption_cat == "Below median"))/sum((consumption_cat == "Below median"))) / (sum(not_care*(1 - (consumption_cat == "Below median")))/sum(1 - (consumption_cat == "Below median"))),
          food_insecure = (sum(food_insecure*(consumption_cat == "Below median"))/sum((consumption_cat == "Below median"))) / (sum(food_insecure*(1 - (consumption_cat == "Below median")))/sum(1 - (consumption_cat == "Below median"))),
          ipv_physical = (sum(ipv_physical*(consumption_cat == "Below median"))/sum((consumption_cat == "Below median"))) / (sum(ipv_physical*(1 - (consumption_cat == "Below median")))/sum(1 - (consumption_cat == "Below median"))),
          sex_violence = (sum(sex_violence*(consumption_cat == "Below median"))/sum((consumption_cat == "Below median"))) / (sum(sex_violence*(1 - (consumption_cat == "Below median")))/sum(1 - (consumption_cat == "Below median"))),
          school_violence = (sum(school_violence*(consumption_cat == "Below median"))/sum((consumption_cat == "Below median"))) / (sum(school_violence*(1 - (consumption_cat == "Below median")))/sum(1 - (consumption_cat == "Below median")))
  )


# Bind rows together and round
confounder_exposure_table <- rbind(confounder_exposure_table_some, confounder_exposure_table_age, confounder_exposure_table_consumption) |>
  mutate_if(is.numeric, round, digits = 2)

# Save CSV file
write.csv(confounder_exposure_table, paste0(file_path, "Output Data/confounder_exposure.csv"))


#### Appendix 1: Table 16, RR for confounders and outcome

# Calculate risk ratio between each confounder and depression at each follow-up time and round
confounder_outcome_table <- adapt_analytic_main |>
  group_by(visit_number) |>
  reframe(maternal_education = (sum(depression*someprimaryorless_mother, na.rm = TRUE)/sum(someprimaryorless_mother)) / (sum(depression*(1 - someprimaryorless_mother), na.rm = TRUE)/sum(1 - someprimaryorless_mother)),
          enrollment_age = (sum(depression*(age_cat == "16-17"), na.rm = TRUE)/sum((age_cat == "16-17"))) / (sum(depression*(1 - (age_cat == "16-17")), na.rm = TRUE)/sum(1 - (age_cat == "16-17"))),
          consumption = (sum(depression*(consumption_cat == "Below median"), na.rm = TRUE)/sum((consumption_cat == "Below median"))) / (sum(depression*(1 - (consumption_cat == "Below median")), na.rm = TRUE)/sum(1 - (consumption_cat == "Below median")))) |>
  data.table::transpose(keep.names = "confounder", make.names = "visit_number") |>
  mutate_if(is.numeric, round, digits = 2)
  
# Save CSV file
write.csv(confounder_outcome_table, paste0(file_path, "Output Data/confounder_outcome.csv"))

