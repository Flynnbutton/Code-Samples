# This script was developed as part of a group research project evaluating the
# impact of Peru's Internet Para Todos (IPT) initiative on local labour market
# outcomes. Using a district-level panel compiled from Peru's ENAHO household
# survey (constructed separately) this script matches treated and control
# districts via propensity score matching, producing matched datasets for
# subsequent difference-in-differences and event study analysis in Stata.
#
# Though completed as part of a group project, I am the sole author of this
# script.
#
#
#
#
# =============================================================================# =============================================================================
# Diff-in-Diff + Propensity Score Matching
#
# Purpose:  Evaluates the labour market effects of Peru's Internet Para Todos
#           (IPT) rural broadband rollout using a PSM-DiD design. Districts in
#           departments receiving IPT by 2018 (early treated) are matched to
#           comparable control districts on 2016 baseline characteristics.
#           Matching is performed on two samples: the full district panel and a
#           restricted sample excluding second-round IPT departments.
#           Matched district IDs are exported to Stata for DiD estimation and
#           Event Study Analysis
#
# Author:  Flynn Button
#
# Date:     2026
# Data:     ENAHO household survey aggregated to district level, 2016-2022
#
# Output:   psm_matched_ids.dta   — full sample matched districts
#           psm_matched_ids_2.dta — restricted sample matched districts
#           Balance tables and propensity score overlap plots (LaTeX/PNG)
#
# Packages: dplyr, tidyr, fixest, ggplot2, kableExtra, purrr, broom
# =============================================================================

rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(fixest)       
library(Rfast)
library(cobalt)       
library(modelsummary) 
library(ggplot2)
library(kableExtra)   
library(broom)
library(purrr)
library(haven)


# --------------------------------------------------------------------------- #
# PATHS  — adjust as needed
# --------------------------------------------------------------------------- #
root   <- "/Users/Flynn/Desktop/BSE/Semester 2/Peru Project"
data   <- file.path(root, "Data")
code   <- file.path(root, "Code")
output <- file.path(root, "Outputs")



# --------------------------------------------------------------------------- #
# Standardized Mean Difference Function
# --------------------------------------------------------------------------- #
compute_smd <- function(df, vars, treat_var) {
  sapply(vars, function(v) {
    m1 <- mean(df[[v]][df[[treat_var]] == 1], na.rm = TRUE)
    s1 <- sd(df[[v]][df[[treat_var]] == 1],   na.rm = TRUE)
    m0 <- mean(df[[v]][df[[treat_var]] == 0], na.rm = TRUE)
    s0 <- sd(df[[v]][df[[treat_var]] == 0],   na.rm = TRUE)
    (m1 - m0) / sqrt((s1^2 + s0^2) / 2)
  })
}

# --------------------------------------------------------------------------- #
# COVARIATES
# --------------------------------------------------------------------------- #
all_covs <- c("employed", "formal", "hours_week",
              "urban_hh", "hh_electricity_hh", "uses_internet", "internet_home",
              "wage_income", "profit_income", "labor_income",
              "education", "age")

cov_labels <- c(
  employed       = "Employment rate",
  formal         = "Formal employment rate",
  hours_week     = "Weekly hours worked",
  urban_hh       = "Urban household share",
  hh_electricity_hh = "Household electricity access",
  uses_internet  = "Internet usage rate",
  internet_home  = "Home internet access",
  wage_income    = "Mean wage income (soles)",
  profit_income  = "Mean profit income (soles)",
  labor_income   = "Mean labour income (soles)",
  education      = "Mean years of education",
  age            = "Mean age"
)

controls_formula <- ~ profit_income + wage_income + age + education + urban_hh


# =============================================================================
# SECTION 1 — FULL SAMPLE (includes intermediate departments)
# =============================================================================

# --------------------------------------------------------------------------- #
# 1a. Pre-matching balance (full sample, 2016 baseline)
# --------------------------------------------------------------------------- #
df_raw <- read_csv(file.path(data, "district_combined.csv"))

df_raw <- df_raw |>
  mutate(
    dept_code        = floor(district_code / 10000),
    EarlyRolloutDept = as.integer(dept_code %in% c(3, 5, 8, 9, 10))
  )



baseline_full <- df_raw |>
  filter(year == 2016) |>
  distinct(district_code, .keep_all = TRUE)
# Baseline_full is a cross-section of 1,256 districts at 2016 baseline.
# EarlyRolloutDept = 1 for departments 3, 5, 8, 9, 10 (received IPT by 2018).



# T-test balance table (pre-matching)
balance_pre_full <- map_dfr(all_covs, function(v) {
  tt <- t.test(baseline_full[[v]] ~ baseline_full$EarlyRolloutDept,
               var.equal = FALSE)
  tibble(
    variable   = v,
    mean_early = tt$estimate[2],
    mean_late  = tt$estimate[1],
    p_value    = formatC(tt$p.value, format = "f", digits = 3)
  )
}) |>
  mutate(variable = cov_labels[variable])


# Pre-matching balance: significant pre-existing conditions in Early treated vs Control
# Regions. Statistically significantly lower formality rate, income, percent urban etc.
# Imbalance is consistent with IPT targeting underdeveloped areas first



# Export pre-matching balance table (LaTeX)
balance_pre_full |>
  kbl(
    format   = "latex", booktabs = TRUE, digits = 3,
    col.names = c("Variable", "Early Treated", "Control", "p-value"),
    caption  = "Pre-Matching Balance --- Baseline 2016 (Full Sample) \\label{tab:prebalance2}"
  ) |>
  footnote(general = paste(
    "Full sample including intermediate departments.",
    "District-level means at 2016 baseline.",
    "Early treated = departments receiving IPT by 2018.",
    "p-values from two-sample t-tests with unequal variances."
  )) |>
  save_kable(file.path(output, "balance_prematching.tex"))



# SMD pre-matching (full sample)
smd_pre_full <- compute_smd(baseline_full, all_covs, "EarlyRolloutDept")
cat("Pre-matching SMDs (full sample):\n"); print(round(smd_pre_full, 3))

# Pre-matching standard mean errors in full sample show substantial imbalance across most covariates.
# Justifies matching


# --------------------------------------------------------------------------- #
# 1b. PSM — full sample
# --------------------------------------------------------------------------- #

psm_formula_full <- as.formula(paste(
  "EarlyRolloutDept ~",
  paste(all_covs, collapse = " + ")
))


# Drop districts with missing covariates
baseline_full_clean <- baseline_full |>
  filter(if_all(all_of(all_covs), ~ is.finite(.) & !is.na(.)))

cat("Districts dropped due to missing covariates:", 
    nrow(baseline_full) - nrow(baseline_full_clean), "\n")
cat("Districts remaining:", nrow(baseline_full_clean), "\n")

# 87 districts dropped due to missing income covariates
# 1,169 districts retained for matching.




# Estimate propensity scores via logit
ps_model <- glm(psm_formula_full, 
                data   = baseline_full_clean, 
                family = binomial("logit"))

baseline_full_clean$pscore <- predict(ps_model, type = "response")



# Split into treated and control
treated <- baseline_full_clean |> filter(EarlyRolloutDept == 1)
control <- baseline_full_clean |> filter(EarlyRolloutDept == 0)



# Nearest neighbour matching without replacement
set.seed(123)
matched_control_idx <- sapply(treated$pscore, function(ps) {
  which.min(abs(control$pscore - ps))
})



# Apply caliper 0.05 on raw PS scale
ps_diffs        <- abs(treated$pscore - control$pscore[matched_control_idx])
within_caliper  <- ps_diffs <= 0.05

treated_matched <- treated[within_caliper, ]
control_matched <- control[matched_control_idx[within_caliper], ]

cat("Treated matched:", nrow(treated_matched), "\n")
cat("Control matched:", nrow(control_matched), "\n")

matched_full <- bind_rows(treated_matched, control_matched)


# Logit PS model estimated on all 12 covariates. Nearest-neighbour 1:1
# matching without replacement, caliper = 0.05 on raw PS scale.
# 296 treated and 296 control districts matched within caliper.



# Propensity score overlap plot
baseline_full_clean |>
  mutate(Group = ifelse(EarlyRolloutDept == 1, "Early Treated", "Control")) |>
  ggplot(aes(x = pscore, colour = Group, linetype = Group)) +
  geom_density() +
  scale_colour_manual(values = c("Early Treated" = "navy", "Control" = "firebrick")) +
  labs(title = "Propensity Score Overlap (Full Sample)",
       x     = "P(Early Treated)",
       y     = "Density") +
  theme_minimal()

ggsave(file.path(output, "pscore_overlap.png"), width = 7, height = 4)

# Propensity score overlap looks reasonable — common support region ~0.1 to 0.6.
# Control districts concentrated near 0 (low probability of early treatment)
# reflecting their better baseline characteristics. Early treated peak around 0.4.
# Sufficient overlap for PSM matching; caliper will trim poor matches at the tails.


# --------------------------------------------------------------------------- #
# 1c. Post-matching balance (full sample)
# --------------------------------------------------------------------------- #

balance_post_full <- map_dfr(all_covs, function(v) {
  tt <- t.test(matched_full[[v]] ~ matched_full$EarlyRolloutDept,
               var.equal = FALSE)
  tibble(
    variable   = v,
    mean_early = tt$estimate[2],
    mean_late  = tt$estimate[1],
    p_value    = formatC(tt$p.value, format = "f", digits = 3)
  )
}) |>
  mutate(variable = cov_labels[variable])

# Post-matching balance (full sample): PSM successfully eliminates pre-existing differences. 
# All 12 covariates are well-balanced post-matching 
# Most significant remaining gap is hours worked (p = 0.080), 
# but difference of 1 hour per week is marginal




# Post-matching SMD table
smd_post_full <- compute_smd(matched_full, all_covs, "EarlyRolloutDept")
print(smd_post_full)
# Post-matching, all covariates  below 0.2 threshold; most below 0.1. 
# Hours worked (0.144) and formality (0.081) are the largest remaining gaps
# Propensity score matching produced a well-balanced matched sample.






# Combined pre/post SMD table
smd_table_full <- tibble(
  Variable        = cov_labels[all_covs],
  `Pre-Matching`  = round(smd_pre_full,  3),
  `Post-Matching` = round(smd_post_full, 3)
)

smd_table_full |>
  kbl(
    format  = "latex", booktabs = TRUE,
    caption = "Standardized Mean Differences Before and After PSM (Full Sample) \\label{tab:smd2}"
  ) |>
  footnote(general = "SMD = (mean$_{treated}$ - mean$_{control}$) / pooled SD. Values below 0.1 indicate good balance.") |>
  save_kable(file.path(output, "smd_table.tex"))



# Combined balance table
balance_combined_full <- balance_pre_full |>
  rename(early_pre = mean_early, late_pre = mean_late, p_pre = p_value) |>
  left_join(
    balance_post_full |>
      rename(early_post = mean_early, late_post = mean_late, p_post = p_value),
    by = "variable"
  )

balance_combined_full |>
  kbl(
    format    = "latex", booktabs = TRUE, digits = 3,
    col.names = c("Variable",
                  "Early", "Control", "p-val",
                  "Early", "Control", "p-val"),
    caption   = "Covariate Balance Before and After PSM (Full Sample) \\label{tab:balance2}"
  ) |>
  add_header_above(c(" " = 1, "Pre-Matching" = 3, "Post-Matching" = 3)) |>
  footnote(general = paste(
    "Full sample including intermediate departments.",
    "Early treated = departments receiving IPT by 2018.",
    "p-values from two-sample t-tests with unequal variances.",
    "Post-matching uses PSM nearest-neighbour weights."
  )) |>
  save_kable(file.path(output, "balance_full.tex"))



# --------------------------------------------------------------------------- #
# 1d. Save matched district IDs
# --------------------------------------------------------------------------- #

psm_ids_full <- matched_full |>
  select(district_code, EarlyRolloutDept, pscore)

write_dta(psm_ids_full, file.path(data, "psm_matched_ids.dta"))







# =============================================================================
# SECTION 2 — MATCHING ON RESTRICTED SAMPLE
#
# In this round we restrict control group to districts in departments that are
# never treated in our panel by removing second round IPT prioritizes departments
#
# =============================================================================

second_round_depts <- c(2, 4, 12, 13, 14, 19, 20)


df_restricted <- df_raw |>
  filter(!dept_code %in% second_round_depts)


baseline_restr <- df_restricted |>
  filter(year == 2016) |>
  distinct(district_code, .keep_all = TRUE)



# --------------------------------------------------------------------------- #
# 2a. Pre-matching balance (restricted sample)
# --------------------------------------------------------------------------- #

balance_pre_restr <- map_dfr(all_covs, function(v) {
  tt <- t.test(baseline_restr[[v]] ~ baseline_restr$EarlyRolloutDept,
               var.equal = FALSE)
  tibble(
    variable   = v,
    mean_early = tt$estimate[2],
    mean_late  = tt$estimate[1],
    p_value    = formatC(tt$p.value, format = "f", digits = 3)
  )
}) |>
  mutate(variable = cov_labels[variable])

balance_pre_restr |>
  kbl(
    format    = "latex", booktabs = TRUE, digits = 3,
    col.names = c("Variable", "Early Treated", "Control", "p-value"),
    caption   = "Pre-Matching Balance --- Baseline 2016 \\label{tab:prebalance}"
  ) |>
  footnote(general = paste(
    "District-level means at 2016 baseline.",
    "Early treated = departments receiving IPT by 2018.",
    "p-values from two-sample t-tests with unequal variances."
  )) |>
  save_kable(file.path(output, "balance_prematching_2.tex"))

# Pre-matching balance (restricted): similar pattern to full sample but
# slightly smaller gaps. Electricity and age are not significantly different pre-matching.
# All income variables, formality, urban share and internet are strongly imbalanced


smd_pre_restr <- compute_smd(baseline_restr, all_covs, "EarlyRolloutDept")
cat("\nPre-matching SMDs (restricted):\n"); print(round(smd_pre_restr, 3))


# Pre-matching SMDs (restricted sample): similar pattern to full sample.
# Large gaps in profit income, formality, and internet



# --------------------------------------------------------------------------- #
# 2b. PSM — restricted sample
# --------------------------------------------------------------------------- #




baseline_restr_clean <- baseline_restr |>
  filter(if_all(all_of(all_covs), ~ is.finite(.) & !is.na(.)))

cat("Districts dropped due to missing covariates:",
    nrow(baseline_restr) - nrow(baseline_restr_clean), "\n")
cat("Districts remaining:", nrow(baseline_restr_clean), "\n")






ps_model_restr <- glm(psm_formula_full,
                      data   = baseline_restr_clean,
                      family = binomial("logit"))

baseline_restr_clean$pscore <- predict(ps_model_restr, type = "response")

treated_restr <- baseline_restr_clean |> filter(EarlyRolloutDept == 1)
control_restr <- baseline_restr_clean |> filter(EarlyRolloutDept == 0)

set.seed(123)
matched_control_idx_restr <- sapply(treated_restr$pscore, function(ps) {
  which.min(abs(control_restr$pscore - ps))
})

ps_diffs_restr       <- abs(treated_restr$pscore - control_restr$pscore[matched_control_idx_restr])
within_caliper_restr <- ps_diffs_restr <= 0.05

treated_matched_restr <- treated_restr[within_caliper_restr, ]
control_matched_restr <- control_restr[matched_control_idx_restr[within_caliper_restr], ]

cat("Treated matched:", nrow(treated_matched_restr), "\n")
cat("Control matched:", nrow(control_matched_restr), "\n")
# 296 Treated and Control districs


matched_restr <- bind_rows(treated_matched_restr, control_matched_restr)

# Overlap plot
baseline_restr_clean |>
  mutate(Group = ifelse(EarlyRolloutDept == 1, "Early Treated", "Control")) |>
  ggplot(aes(x = pscore, colour = Group, linetype = Group)) +
  geom_density() +
  scale_colour_manual(values = c("Early Treated" = "navy", "Control" = "firebrick")) +
  labs(title = "Propensity Score Overlap (Restricted Sample)",
       x = "P(Treated)", y = "Density") +
  theme_minimal()

ggsave(file.path(output, "pscore_overlap_2.png"), width = 7, height = 4)

# Restricted sample overlap plot: broader common support than full sample,
# Early treated districts concentrated at a higher propensity score than control
# Removing second-round departments shifts the early
# treated distribution rightward. These districts look more distinct from
# controls, which is expected when cleaner treatment/control groups are used.




# --------------------------------------------------------------------------- #
# 2c. Post-matching balance (restricted)
# --------------------------------------------------------------------------- #



balance_post_restr <- map_dfr(all_covs, function(v) {
  tt <- t.test(matched_restr[[v]] ~ matched_restr$EarlyRolloutDept,
               var.equal = FALSE)
  tibble(
    variable   = v,
    mean_early = tt$estimate[2],
    mean_late  = tt$estimate[1],
    p_value    = formatC(tt$p.value, format = "f", digits = 3)
  )
}) |>
  mutate(variable = cov_labels[variable])


# Post-matching balance (restricted): No variable significant at even the 10% level 
# Lowest P value for Employment Rate at 0.279
# Income variables are particularly well balanced
# Matching performance is strong across both samples.



smd_post_restr <- compute_smd(matched_restr, all_covs, "EarlyRolloutDept")

smd_table_restr <- tibble(
  Variable        = cov_labels[all_covs],
  `Pre-Matching`  = round(smd_pre_restr,  3),
  `Post-Matching` = round(smd_post_restr, 3)
)




#    Save Balance Tables
smd_table_restr |>
  kbl(
    format  = "latex", booktabs = TRUE,
    caption = "Standardized Mean Differences Before and After PSM \\label{tab:smd}"
  ) |>
  footnote(general = "SMD = (mean$_{treated}$ - mean$_{control}$) / pooled SD. Values below 0.1 indicate good balance.") |>
  save_kable(file.path(output, "smd_table_2.tex"))



balance_combined_restr <- balance_pre_restr |>
  rename(early_pre = mean_early, late_pre = mean_late, p_pre = p_value) |>
  left_join(
    balance_post_restr |>
      rename(early_post = mean_early, late_post = mean_late, p_post = p_value),
    by = "variable"
  )

balance_combined_restr |>
  kbl(
    format    = "latex", booktabs = TRUE, digits = 3,
    col.names = c("Variable",
                  "Early Treated", "Control", "p-val",
                  "Early Treated", "Control", "p-val"),
    caption   = "Covariate Balance Before and After PSM \\label{tab:balance}"
  ) |>
  add_header_above(c(" " = 1, "Pre-Matching" = 3, "Post-Matching" = 3)) |>
  footnote(general = paste(
    "District-level means at 2016 baseline.",
    "Early treated = departments receiving IPT by 2018.",
    "Control = departments not receiving IPT in sample period.",
    "p-values from two-sample t-tests with unequal variances.",
    "Post-matching uses PSM nearest-neighbour weights."
  )) |>
  save_kable(file.path(output, "balance_full_2.tex"))






# --------------------------------------------------------------------------- #
# 2d. Save matched IDs (restricted)
# --------------------------------------------------------------------------- #

psm_ids_restr <- matched_restr |>
  select(district_code, EarlyRolloutDept, pscore)

write_dta(psm_ids_restr, file.path(data, "psm_matched_ids_2.dta"))

