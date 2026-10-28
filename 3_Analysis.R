# ============================================
# Libraries
# ============================================
### Install ptestr
#install.packages("devtools")
#library(devtools)
#install_github('CoDe-Neuro/ptestR', force = TRUE)

# Load packages
library(haven)
library(readxl)
library(tidyr)
library(dplyr)
library(purrr)
library(stringdist)
library(lmerTest)
library(influence.ME)
library(ptestr)
library(emmeans)
library(r2glmm)
library(broom.mixed)
library(psych)
library(ggfortify)
library(ggplot2)
library(patchwork)

# ============================================
# Load data
# ============================================
df1 <- read_sav("~/Postdoc/CannTeen/Misho_Acute_Baseline_cleaned_v2.4_10.5.22.sav")
df2 <- read_excel("~/Postdoc/CannTeen/Cannteen_Teenagers_Values.xlsx", skip = 1)
df3 <- read_excel("~/Postdoc/CannTeen/Cannteen_Adults_Values.xlsx", skip = 1)
df_drugs <- read_excel("~/Postdoc/CannTeen/CannTeen_acute_visits_for_Uzma_updated.xlsx")
df_other <- read_sav("~/Postdoc/CannTeen/Misho_Invicro123_final_v1.17_Dec22.sav")

# Convert drug df from wide to long format
df_drugs_long <- df_drugs %>%
  pivot_longer(cols = starts_with("Visit"), names_to = "Visit") %>%
  transmute(Sub, Visit = as.integer(gsub("\\D", "", Visit)), condition = value)

# Remove rows from drug df that will be missing from the merged df
# These are:
# 16018: no V3
# 16025: no V2
# 26023: no V1
# 26027: no V3
df_drugs_long <- df_drugs_long[-c(45, 65, 130, 144), ]

# Prep df's for merging
df1 <- df1 %>% mutate(Participant_ID = as.character(Participant_ID))
df2 <- df2 %>% mutate(subid = as.character(subid))
df3 <- df3 %>% mutate(subid = as.character(subid))

# Merge df's 1, 2 and 3
df23 <- bind_rows(df2, df3)
merged_df <- left_join(df1, df23, by = c("Participant_ID" = "subid"))

# Merge df123 with drug df
merged_df_with_drugs <- cbind(df_drugs_long[, c("Sub", "Visit", "condition")], merged_df)

# ============================================
# Clean data
# ============================================
# Convert all labelled or character variables that should be categorical to factors
merged_df_with_drugs <- merged_df_with_drugs %>%
  mutate(across(where(is.labelled), haven::as_factor)) %>%   # SPSS labelled
  mutate(across(where(is.character), as.factor))            # character

# Rename funky Adjusted Glu or Glx variable
names(merged_df_with_drugs)[names(merged_df_with_drugs) == "*Glx = M*(wm + 1.21*gm + 1.55*csf)/(wm+gm)"] <- "Glx_adjusted"

###################################
# Fix typos in the condition column and set the correct levels

# Convert to plain character
vals <- as.character(merged_df_with_drugs$condition)
vals <- trimws(vals)         # remove extra spaces
vals <- tolower(vals)         # make lowercase

# Define the 3 correct target levels
valid_levels <- c("placebo", "thc", "thc+cbd")

# For each entry, find the closest match among valid_levels
corrected <- sapply(vals, function(x) {
  distances <- stringdist::stringdist(x, valid_levels, method = "osa")
  valid_levels[which.min(distances)]
})

# Capitalize nicely
corrected <- ifelse(corrected == "placebo", "Placebo",
                    ifelse(corrected == "thc", "THC", "THC+CBD"))

# Create a review table of old vs corrected values
review_table <- data.frame(
  Original = merged_df_with_drugs$condition,
  Cleaned  = corrected
)

# Show unique mappings for inspection
unique(review_table)

# Finally, assign the cleaned values back to your dataframe
merged_df_with_drugs$condition <- factor(corrected)

# Check the final levels
levels(merged_df_with_drugs$condition)

###################################
# Fix biological sex levels 
merged_df_with_drugs$Biological_Sex <- droplevels(merged_df_with_drugs$Biological_Sex)
# Check the final levels
levels(merged_df_with_drugs$Biological_Sex)

###################################
# Set all variables that will be used in the analysis to the correct types
merged_df_with_drugs <- merged_df_with_drugs %>%
  dplyr::mutate(
    # Modeling predictors
    condition      = factor(condition, levels = c("Placebo", "THC", "THC+CBD")),
    AgeGroup       = if (is.factor(AgeGroup)) AgeGroup else factor(AgeGroup),
    # Random effect grouping
    Participant_ID = as.factor(Participant_ID)
  )

# ============================================
# MRS QC
# ============================================
# Remove session that don't meet QC criteria
merged_df_clean <- merged_df_with_drugs %>%
  filter(`FWHM (ppm)` <= 0.1,
         `Glu+Gln %SD` < 20,
         `S/N` >= 5,
         `Glx_adjusted` >= 1)

# ============================================
# Check for outliers and check residuals 
# ============================================
# Run LME
lme_model <- lmer(Glx_adjusted ~ condition * AgeGroup + (1|Participant_ID), data = merged_df_clean)
summary(lme_model)

# Calculate influence measures by level of the grouping factor (e.g. 'Participant_ID')
infl <- influence(lme_model, group = "Participant_ID")

# Cook’s distance for each subject
cooks <- cooks.distance(infl)

# Plot Cook's distances
plot(cooks, type="h", main="Cook's distance by subject")
abline(h = 4/length(cooks), col = "red", lty = 2)  # rough cutoff

# Annotate max and top 5 influential points with IDs and indices
subj_labels <- names(cooks)
if (is.null(subj_labels) || length(subj_labels) != length(cooks)) {
  subj_labels <- as.character(seq_along(cooks))
}
max_idx <- which.max(cooks)
points(max_idx, cooks[max_idx], col = "red", pch = 19)
text(max_idx, cooks[max_idx], labels = paste0(subj_labels[max_idx], " (", max_idx, ")"),
  pos = 3, cex = 0.8, col = "red")

top_idx <- order(cooks, decreasing = TRUE)[1:min(5, length(cooks))]
points(top_idx, cooks[top_idx], col = "blue", pch = 19)
text(top_idx, cooks[top_idx], labels = paste0(subj_labels[top_idx], " (", top_idx, ")"),
  pos = 3, cex = 0.7, col = "blue")

# REMOVE only the THC condition for ID 16021 and refit model
merged_df_clean_final <- merged_df_clean %>%
  dplyr::filter(!(Sub == 16021 & condition == "THC"))

# Run LME again
lme_model <- lmer(Glx_adjusted ~ condition * AgeGroup + (1|Participant_ID), data = merged_df_clean_final)
summary(lme_model)

# Calculate influence measures by level of the grouping factor (e.g. 'Participant_ID') again
infl <- influence(lme_model, group = "Participant_ID")

# Cook’s distance for each subject again
cooks <- cooks.distance(infl)

# Plot Cook's distances again
plot(cooks, type="h", main="Cook's distance by subject")
abline(h = 4/length(cooks), col = "red", lty = 2)  # rough cutoff

# Annotate max and top 5 influential points with IDs and indices (after refit)
subj_labels <- names(cooks)
if (is.null(subj_labels) || length(subj_labels) != length(cooks)) {
  subj_labels <- as.character(seq_along(cooks))
}
max_idx <- which.max(cooks)
points(max_idx, cooks[max_idx], col = "red", pch = 19)
text(max_idx, cooks[max_idx], labels = paste0(subj_labels[max_idx], " (", max_idx, ")"),
  pos = 3, cex = 0.8, col = "red")

top_idx <- order(cooks, decreasing = TRUE)[1:min(5, length(cooks))]
points(top_idx, cooks[top_idx], col = "blue", pch = 19)
text(top_idx, cooks[top_idx], labels = paste0(subj_labels[top_idx], " (", top_idx, ")"),
  pos = 3, cex = 0.7, col = "blue")

# Residual vs fitted plot
plot(lme_model)

# QQ plot of residuals
qqnorm(residuals(lme_model))
qqline(residuals(lme_model))

# Scale-location plot
plot(fitted(lme_model), sqrt(abs(residuals(lme_model))), main = "Scale-Location", xlab = "Fitted values", ylab = "sqrt(|Residuals|)")

# Save merged_df_clean_final to CSV for future reference
out_dir <- "~/Postdoc/CannTeen/outputs"
utils::write.csv(merged_df_clean_final, file.path(out_dir, "merged_df_clean_final.csv"), row.names = FALSE)

# ============================================
# Analysis
# ============================================
#### 1. Descriptive stats #### 
describe(merged_df_clean_final$Glx_adjusted)
describeBy(merged_df_clean_final$Glx_adjusted, merged_df_clean_final$condition)
describeBy(
  merged_df_clean_final$Glx_adjusted,
  list(
    Condition = merged_df_clean_final$condition,
    AgeGroup  = merged_df_clean_final$AgeGroup
  ),
  mat = TRUE
)
## Neat tables: mean and variance by condition (overall) and by condition within AgeGroup

# By AgeGroup and Condition
tbl_cond_by_age <- merged_df_clean_final %>%
  dplyr::group_by(AgeGroup, condition) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(Glx_adjusted, na.rm = TRUE),
    variance = var(Glx_adjusted, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(AgeGroup, condition) %>%
  dplyr::mutate(dplyr::across(c(mean, variance), ~ round(.x, 2)))

# Overall by Condition (not split by AgeGroup)
tbl_cond <- merged_df_clean_final %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(Glx_adjusted, na.rm = TRUE),
    variance = var(Glx_adjusted, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(condition) %>%
  dplyr::mutate(dplyr::across(c(mean, variance), ~ round(.x, 2)))

# Print tables
tbl_cond_by_age
tbl_cond

# Save tables to CSVs
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
utils::write.csv(tbl_cond_by_age, file.path(out_dir, "Glx_summary_by_age_condition.csv"), row.names = FALSE)
utils::write.csv(tbl_cond,       file.path(out_dir, "Glx_summary_by_condition.csv"),       row.names = FALSE)

## Missing data table
# Sessions missing after QC/removal, by AgeGroup and condition
base_counts <- merged_df_with_drugs %>%
  dplyr::count(AgeGroup, condition, name = "total_sessions")

present_counts <- merged_df_clean_final %>%
  dplyr::count(AgeGroup, condition, name = "present_sessions")

missing_table <- dplyr::full_join(base_counts, present_counts,
                                  by = c("AgeGroup", "condition")) %>%
  dplyr::mutate(
    total_sessions = tidyr::replace_na(total_sessions, 0L),
    present_sessions = tidyr::replace_na(present_sessions, 0L),
    missing_sessions = total_sessions - present_sessions,
    missing_pct = dplyr::if_else(total_sessions > 0,
                                 round(100 * missing_sessions / total_sessions, 1),
                                 NA_real_)
  ) %>%
  dplyr::arrange(AgeGroup, condition)

# Print missing data table
missing_table

# Save to CSV alongside other outputs
utils::write.csv(missing_table,
                 file.path(out_dir, "Glx_missing_sessions_by_age_condition.csv"),
                 row.names = FALSE)

#### 2. (Re-)Fit your models ####
contrasts(merged_df_clean_final$condition)
contrasts(merged_df_clean_final$AgeGroup)

# Use custom Helmert contrasts for condition:
# (i) ActiveAvg(THC, THC+CBD) vs Placebo; (ii) THC+CBD vs THC
merged_df_clean_final$condition <- factor(
  merged_df_clean_final$condition,
  levels = c("Placebo", "THC", "THC+CBD")
)
custom_contr <- cbind(
  ActiveAvg_vs_Placebo = c(-2, 1, 1),
  THC_CBD_vs_THC       = c( 0, -1, 1)
)
rownames(custom_contr) <- c("Placebo", "THC", "THC+CBD")
contrasts(merged_df_clean_final$condition) <- custom_contr

# Use Helmert for AgeGroup
contrasts(merged_df_clean_final$AgeGroup) <- contr.helmert(nlevels(merged_df_clean_final$AgeGroup))

# Show contrasts
contrasts(merged_df_clean_final$condition)
contrasts(merged_df_clean_final$AgeGroup)

lme_model <- lmer(Glx_adjusted ~ condition * AgeGroup + (1|Participant_ID), data = merged_df_clean_final)
summary(lme_model)

# Run LME with permutations
lme_perm_model <- grouped_perm_glmm(merged_df_clean_final, formla = Glx_adjusted ~ condition * AgeGroup + (1|Participant_ID), var_to_perm='Glx_adjusted', permNum = 5000, seed = 42)

#### 3. Fixed effects table (betas, SEs, t, CI) ####
fixef_table <- broom.mixed::tidy(lme_model, effects = "fixed", conf.int = TRUE) # Do not use these p-values

#### 4. Get permutation p-values ####
perm_pvals <- lme_perm_model$p.perm
names(perm_pvals) <- lme_perm_model$term
# add permutation p-values to table
fixef_table$p_perm <- perm_pvals[fixef_table$term]

#### 5. Remove unnecessary columns and re-arrange the rest ####
results_table <- fixef_table %>%
  select(term, estimate, std.error, 
         CI_low = conf.low, CI_high = conf.high, 
         df, statistic, p_perm)

# Save to CSV alongside other outputs
results_table_rounded <- results_table %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))
utils::write.csv(results_table_rounded,
                 file.path(out_dir, "Glx_by_age_condition.csv"),
                 row.names = FALSE)

#### 6. Partial R² per predictor ####
r2_table <- r2beta(lme_model, method = "nsj", partial = TRUE)
r2_table <- as.data.frame(r2_table)
r2_table$term <- rownames(r2_table)
r2_table <- r2_table %>% filter(term != "1")
# Reorder columns nicely
r2_table <- r2_table %>%
  select(Effect, Rsq, lower.CL, upper.CL)
colnames(r2_table) <- c("term", "partialR2", "CI_low", "CI_high")
# !!!!!! R2 is only estimated for the 4 original terms so can't be merged with the other table !!!!!!
# Save to CSV alongside other outputs
r2_table_rounded <- r2_table %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))
utils::write.csv(r2_table_rounded,
                 file.path(out_dir, "Glx_by_age_condition_R2.csv"),
                 row.names = FALSE)

#### 7. Given no significant interactions, collapse across age groups ####
# Check model fit - complex vs simple
m_full <- lmer(
  Glx_adjusted ~ condition * AgeGroup + (1 | Participant_ID),
  data = merged_df_clean_final,
  REML = FALSE
)

m_no_int <- lmer(
  Glx_adjusted ~ condition + AgeGroup + (1 | Participant_ID),
  data = merged_df_clean_final,
  REML = FALSE
)

m_no_age <- lmer(
  Glx_adjusted ~ condition + (1 | Participant_ID),
  data = merged_df_clean_final,
  REML = FALSE
)

anova(m_no_int, m_full)
anova(m_no_age, m_full)

#####################                       #####################
#### 8. Run model only with main effect of condition ####
lme_model_cond <- lmer(Glx_adjusted ~ condition + (1|Participant_ID), data = merged_df_clean_final)
summary(lme_model_cond)

lme_perm_model_cond <- grouped_perm_glmm(merged_df_clean_final, formla = Glx_adjusted ~ condition + (1|Participant_ID), var_to_perm='Glx_adjusted', permNum = 5000, seed = 42)

#### 9. Fixed effects table (betas, SEs, t, CI) ####
fixef_table_cond <- broom.mixed::tidy(lme_model_cond, effects = "fixed", conf.int = TRUE) # Do not use these p-values

#### 10. Get permutation p-values ####
perm_pvals_cond <- lme_perm_model_cond$p.perm
names(perm_pvals_cond) <- rownames(summary(lme_model_cond)$coefficients)
#perm_pvals_clean <- perm_pvals_cond[rownames(summary(lme_model_cond)$coefficients)]
# add permutation p-values to table
fixef_table_cond$p_perm <- perm_pvals_cond[fixef_table_cond$term]

#### 11. Remove unnecessary columns and re-arrange the rest ####
results_table_cond <- fixef_table_cond %>%
  select(term, estimate, std.error, 
         CI_low = conf.low, CI_high = conf.high, 
         df, statistic, p_perm)

# Save to CSV alongside other outputs
results_table_cond_rounded <- results_table_cond %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))
utils::write.csv(results_table_cond_rounded,
                 file.path(out_dir, "Glx_by_condition.csv"),
                 row.names = FALSE)

#### 12. Partial R² per predictor ####
r2_table_cond <- r2beta(lme_model_cond, method = "nsj", partial = TRUE)
r2_table_cond <- as.data.frame(r2_table_cond)
r2_table_cond$term <- rownames(r2_table_cond)
r2_table_cond <- r2_table_cond %>% filter(term != "1")
# Reorder columns nicely
r2_table_cond <- r2_table_cond %>%
  select(Effect, Rsq, lower.CL, upper.CL)
colnames(r2_table_cond) <- c("term", "partialR2", "CI_low", "CI_high")

# Save to CSV alongside other outputs
r2_table_cond_rounded <- r2_table_cond %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))
utils::write.csv(r2_table_cond_rounded,
                 file.path(out_dir, "Glx_by_condition_R2.csv"),
                 row.names = FALSE)

#### 13. Estimate marginal means (variances) ####
emm_cond_cond <- emmeans(lme_model_cond, ~ condition)

# Focused planned contrasts (do not recompute THC+CBD vs THC here):
# (i) THC vs Placebo, (ii) THC+CBD vs Placebo
contrast_defs_placebo <- list(
  "THC vs Placebo"       = c("Placebo" = -1, "THC" = 1, "THC+CBD" = 0),
  "THC+CBD vs Placebo "   = c("Placebo" = -1, "THC" =  0, "THC+CBD" = 1)
)
contrast_placebo <- contrast(emm_cond_cond, method = contrast_defs_placebo, adjust = "holm")

# Tidy and print (add variance = SE^2, rounded like earlier tables)
emm_cond_cond_tbl <- as.data.frame(emm_cond_cond) %>%
  dplyr::mutate(variance = SE^2) %>%
  dplyr::mutate(dplyr::across(c(emmean, variance), ~ round(.x, 2))) %>%
  dplyr::select(condition, emmean, variance, SE, df, lower.CL, upper.CL)
contrast_placebo_tbl <- broom::tidy(contrast_placebo, conf.int = TRUE)
contrast_placebo_tbl_rounded <- contrast_placebo_tbl %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))
emm_cond_cond_tbl
contrast_placebo_tbl

# Save to CSV
utils::write.csv(emm_cond_cond_tbl,           file.path(out_dir, "EMMs_condition_collapsed.csv"), row.names = FALSE)
utils::write.csv(contrast_placebo_tbl_rounded, file.path(out_dir, "Contrasts_placebo_vs_active_collapsed.csv"), row.names = FALSE)

#### 14. Post-hoc pairwise tests ####

# Paired tests per condition pair without global case dropping
# Decision rule: Shapiro-Wilk on within-subject differences; if normal -> paired t-test
# else -> Wilcoxon signed-rank test (paired). p-values adjusted (Holm).

pair_test <- function(df, cond_a, cond_b) {
  wide <- df %>%
    dplyr::filter(condition %in% c(cond_a, cond_b)) %>%
    dplyr::select(Participant_ID, condition, Glx_adjusted) %>%
    tidyr::pivot_wider(names_from = condition, values_from = Glx_adjusted)

  if (!all(c(cond_a, cond_b) %in% names(wide))) {
    return(tibble::tibble(
      contrast = paste(cond_a, "vs", cond_b),
      method = NA_character_,
      n = 0L, stat = NA_real_, df = NA_real_,
      p = NA_real_, shapiro_p = NA_real_
    ))
  }

  wide <- wide %>% tidyr::drop_na(dplyr::all_of(c(cond_a, cond_b)))
  if (nrow(wide) < 2) {
    return(tibble::tibble(
      contrast = paste(cond_a, "vs", cond_b),
      method = NA_character_,
      n = nrow(wide), stat = NA_real_, df = NA_real_,
      p = NA_real_, shapiro_p = NA_real_
    ))
  }

  x <- wide[[cond_a]]
  y <- wide[[cond_b]]
  diff <- x - y
  sh <- tryCatch(stats::shapiro.test(diff), error = function(e) NULL)
  sh_p <- if (!is.null(sh)) sh$p.value else NA_real_

  if (!is.na(sh_p) && sh_p >= 0.05) {
    tst <- stats::t.test(x, y, paired = TRUE)
    tibble::tibble(
      contrast = paste(cond_a, "vs", cond_b),
      method = "paired t-test",
      n = length(diff), stat = unname(tst$statistic),
      df = unname(tst$parameter), p = tst$p.value,
      shapiro_p = sh_p
    )
  } else {
    tst <- stats::wilcox.test(x, y, paired = TRUE, exact = FALSE)
    tibble::tibble(
      contrast = paste(cond_a, "vs", cond_b),
      method = "Wilcoxon signed-rank",
      n = length(diff), stat = unname(tst$statistic),
      df = NA_real_, p = tst$p.value, shapiro_p = sh_p
    )
  }
}

pairs_to_test <- list(
  c("Placebo", "THC"),
  c("Placebo", "THC+CBD"),
  c("THC",     "THC+CBD")
)

pairwise_tests_tbl <- purrr::map_dfr(pairs_to_test, ~ pair_test(merged_df_clean_final, .x[1], .x[2])) %>%
  dplyr::mutate(p_adj_holm = stats::p.adjust(p, method = "holm"))

# Print and save
pairwise_tests_tbl
pairwise_tests_tbl_rounded <- pairwise_tests_tbl %>%
  dplyr::mutate(dplyr::across(where(is.double), ~ round(.x, 3)))
utils::write.csv(pairwise_tests_tbl_rounded, file.path(out_dir, "Pairwise_condition_paired_tests.csv"), row.names = FALSE)

#### 15. Plot Glx across conditions ####

## Violin plot with jitter (uses all data in merged_df_clean_final)
## Map AgeGroup to point color and shape
gg_violin <- ggplot(merged_df_clean_final, aes(x = condition, y = Glx_adjusted, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(aes(color = AgeGroup, shape = AgeGroup), width = 0.15, alpha = 0.85, size = 1.9) +
  labs(title = "Glx Across Conditions", x = "Condition", y = "Glx", color = "Age Group", shape = "Age Group") +
  theme_classic() +
  guides(fill = "none")

print(gg_violin)
ggplot2::ggsave(filename = file.path(out_dir, "Glx_conditions_violin.png"), plot = gg_violin, width = 7, height = 4.5, dpi = 300)

#### 16. Correlation between Glx and THC (and CBD) absorption ####

# Build per-subject Glx shifts (THC - Placebo; THC+CBD - Placebo) using PNum_Simple
glx_shifts <- merged_df_clean_final %>%
  dplyr::select(Participant_ID, PNum_Simple, AgeGroup, condition, Glx_adjusted) %>%
  tidyr::pivot_wider(
    id_cols = c(Participant_ID, PNum_Simple, AgeGroup),
    names_from = condition,
    values_from = Glx_adjusted
  ) %>%
  dplyr::mutate(
    THC_minus_placebo = THC - Placebo,
    THCCBD_minus_placebo = `THC+CBD` - Placebo
  )

# Merge THC absorption deltas from df_other using PNum_Simple
glx_shifts <- glx_shifts %>% dplyr::mutate(Participant_ID = as.character(Participant_ID), PNum_Simple = as.character(PNum_Simple))
df_other_join <- df_other %>%
  dplyr::mutate(PNum_Simple = as.character(PNum_Simple)) %>%
  dplyr::select(
    PNum_Simple,
    T2_THCconc_THC_minusPLA,
    T2_THCconc_THCCBD_minusPLA,
    T2_CBD_conc_THC_CBD
  ) %>%
  dplyr::mutate(
    T2_THCconc_THC_minusPLA = suppressWarnings(as.numeric(T2_THCconc_THC_minusPLA)),
    T2_THCconc_THCCBD_minusPLA = suppressWarnings(as.numeric(T2_THCconc_THCCBD_minusPLA)),
    T2_CBD_conc_THC_CBD = suppressWarnings(as.numeric(T2_CBD_conc_THC_CBD))
  )

glx_absorption <- dplyr::left_join(
  glx_shifts,
  df_other_join,
  by = "PNum_Simple"
)

# Correlation helper: choose Pearson if both variables ~normal, else Spearman
corr_with_assumptions <- function(df, shift_col, absorp_col, label) {
  dat <- df %>% dplyr::select(dplyr::all_of(c(shift_col, absorp_col))) %>% tidyr::drop_na()
  n <- nrow(dat)
  if (n < 3) {
    return(tibble::tibble(
      comparison = label,
      method = NA_character_, n = n, r = NA_real_, p = NA_real_,
      shapiro_p_shift = NA_real_, shapiro_p_absorption = NA_real_
    ))
  }
  sh_x <- tryCatch(stats::shapiro.test(dat[[shift_col]])$p.value, error = function(e) NA_real_)
  sh_y <- tryCatch(stats::shapiro.test(dat[[absorp_col]])$p.value, error = function(e) NA_real_)
  use_pearson <- (!is.na(sh_x) && sh_x >= 0.05) && (!is.na(sh_y) && sh_y >= 0.05)
  method <- if (use_pearson) "pearson" else "spearman"
  ct <- suppressWarnings(stats::cor.test(dat[[shift_col]], dat[[absorp_col]], method = method, exact = FALSE))
  tibble::tibble(
    comparison = label,
    method = method,
    n = n,
    r = unname(ct$estimate),
    p = ct$p.value,
    shapiro_p_shift = sh_x,
    shapiro_p_absorption = sh_y
  )
}

# Permutation test for correlation: fix absorption, permute shift across subjects
perm_test_corr <- function(df, shift_col, absorp_col, method, permNum = 5000, seed = 42) {
  dat <- df %>% dplyr::select(dplyr::all_of(c(shift_col, absorp_col))) %>% tidyr::drop_na()
  n <- nrow(dat)
  if (n < 3) {
    return(list(
      p_perm = NA_real_, r_obs = NA_real_, permNum = permNum,
      perm_df = tibble::tibble(perm_index = integer(0), r_perm = numeric(0))
    ))
  }
  x <- dat[[shift_col]]
  y <- dat[[absorp_col]]
  r_obs <- suppressWarnings(stats::cor(x, y, method = method))
  set.seed(seed)
  r_perm <- replicate(permNum, {
    xp <- sample(x, replace = FALSE)
    suppressWarnings(stats::cor(xp, y, method = method))
  })
  # Two-sided empirical p-value with +1 correction
  p_perm <- (1 + sum(abs(r_perm) >= abs(r_obs))) / (permNum + 1)
  list(
    p_perm = p_perm,
    r_obs = r_obs,
    permNum = permNum,
    perm_df = tibble::tibble(perm_index = seq_along(r_perm), r_perm = as.numeric(r_perm))
  )
}

# Bootstrap percentile CI for correlation: resample subject rows with replacement
boot_ci_corr <- function(df, shift_col, absorp_col, method, B = 2000, seed = 123) {
  dat <- df %>% dplyr::select(dplyr::all_of(c(shift_col, absorp_col))) %>% tidyr::drop_na()
  n <- nrow(dat)
  if (n < 3) {
    return(list(
      ci_low = NA_real_, ci_high = NA_real_, B = B,
      boot_df = tibble::tibble(boot_index = integer(0), r_boot = numeric(0))
    ))
  }
  set.seed(seed)
  r_boot <- replicate(B, {
    idx <- sample.int(n, size = n, replace = TRUE)
    x <- dat[[shift_col]][idx]
    y <- dat[[absorp_col]][idx]
    suppressWarnings(stats::cor(x, y, method = method))
  })
  ci <- stats::quantile(as.numeric(r_boot), probs = c(0.025, 0.975), na.rm = TRUE)
  list(
    ci_low = as.numeric(ci[1]),
    ci_high = as.numeric(ci[2]),
    B = B,
    boot_df = tibble::tibble(boot_index = seq_along(r_boot), r_boot = as.numeric(r_boot))
  )
}

# Run correlations
corr_thc <- corr_with_assumptions(
  glx_absorption,
  "THC_minus_placebo",
  "T2_THCconc_THC_minusPLA",
  "THC shift vs THC absorption"
)

corr_thccbd <- corr_with_assumptions(
  glx_absorption,
  "THCCBD_minus_placebo",
  "T2_THCconc_THCCBD_minusPLA",
  "THC+CBD shift vs THC+CBD absorption"
)

# Additional correlation: THC+CBD Glx shift vs CBD absorption shift
corr_cbd <- corr_with_assumptions(
  glx_absorption,
  "THCCBD_minus_placebo",
  "T2_CBD_conc_THC_CBD",
  "THC+CBD shift vs CBD absorption"
)

# Permutation tests (fix Δabsorption, permute ΔGlx across subjects)
perm_thc <- perm_test_corr(
  glx_absorption,
  "THC_minus_placebo",
  "T2_THCconc_THC_minusPLA",
  method = corr_thc$method[1],
  permNum = 5000,
  seed = 42
)

perm_thccbd <- perm_test_corr(
  glx_absorption,
  "THCCBD_minus_placebo",
  "T2_THCconc_THCCBD_minusPLA",
  method = corr_thccbd$method[1],
  permNum = 5000,
  seed = 42
)

perm_cbd <- perm_test_corr(
  glx_absorption,
  "THCCBD_minus_placebo",
  "T2_CBD_conc_THC_CBD",
  method = corr_cbd$method[1],
  permNum = 5000,
  seed = 42
)

# Bootstrap CIs (percentile) for r by resampling subjects
boot_thc <- boot_ci_corr(
  glx_absorption,
  "THC_minus_placebo",
  "T2_THCconc_THC_minusPLA",
  method = corr_thc$method[1],
  B = 5000,
  seed = 123
)

boot_thccbd <- boot_ci_corr(
  glx_absorption,
  "THCCBD_minus_placebo",
  "T2_THCconc_THCCBD_minusPLA",
  method = corr_thccbd$method[1],
  B = 5000,
  seed = 123
)

boot_cbd <- boot_ci_corr(
  glx_absorption,
  "THCCBD_minus_placebo",
  "T2_CBD_conc_THC_CBD",
  method = corr_cbd$method[1],
  B = 5000,
  seed = 123
)

perm_summary <- dplyr::bind_rows(
  tibble::tibble(
    comparison = "THC shift vs THC absorption",
    p_perm = perm_thc$p_perm,
    permNum = perm_thc$permNum,
    r_obs_perm = perm_thc$r_obs
  ),
  tibble::tibble(
    comparison = "THC+CBD shift vs THC+CBD absorption",
    p_perm = perm_thccbd$p_perm,
    permNum = perm_thccbd$permNum,
    r_obs_perm = perm_thccbd$r_obs
  ),
  tibble::tibble(
    comparison = "THC+CBD shift vs CBD absorption",
    p_perm = perm_cbd$p_perm,
    permNum = perm_cbd$permNum,
    r_obs_perm = perm_cbd$r_obs
  )
)

boot_summary <- dplyr::bind_rows(
  tibble::tibble(
    comparison = "THC shift vs THC absorption",
    boot_CI_low = boot_thc$ci_low,
    boot_CI_high = boot_thc$ci_high,
    boot_B = boot_thc$B
  ),
  tibble::tibble(
    comparison = "THC+CBD shift vs THC+CBD absorption",
    boot_CI_low = boot_thccbd$ci_low,
    boot_CI_high = boot_thccbd$ci_high,
    boot_B = boot_thccbd$B
  ),
  tibble::tibble(
    comparison = "THC+CBD shift vs CBD absorption",
    boot_CI_low = boot_cbd$ci_low,
    boot_CI_high = boot_cbd$ci_high,
    boot_B = boot_cbd$B
  )
)

corr_results <- dplyr::bind_rows(corr_thc, corr_thccbd, corr_cbd) %>%
  dplyr::left_join(perm_summary, by = "comparison") %>%
  dplyr::left_join(boot_summary, by = "comparison") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))

# Save subject-level dataset and correlation summaries
utils::write.csv(glx_absorption, file.path(out_dir, "Glx_absorption_shifts_per_subject.csv"), row.names = FALSE)
utils::write.csv(corr_results, file.path(out_dir, "Correlation_Glx_shifts_vs_THC_absorption.csv"), row.names = FALSE)

# Save permutation distributions
utils::write.csv(perm_thc$perm_df, file.path(out_dir, "Permutations_THCshift_vs_THCabs.csv"), row.names = FALSE)
utils::write.csv(perm_thccbd$perm_df, file.path(out_dir, "Permutations_THCCBDshift_vs_THCCBDabs.csv"), row.names = FALSE)
utils::write.csv(perm_cbd$perm_df, file.path(out_dir, "Permutations_THCCBDshift_vs_CBDabs.csv"), row.names = FALSE)

# Save bootstrap distributions
utils::write.csv(boot_thc$boot_df, file.path(out_dir, "Bootstrap_THCshift_vs_THCabs.csv"), row.names = FALSE)
utils::write.csv(boot_thccbd$boot_df, file.path(out_dir, "Bootstrap_THCCBDshift_vs_THCCBDabs.csv"), row.names = FALSE)
utils::write.csv(boot_cbd$boot_df, file.path(out_dir, "Bootstrap_THCCBDshift_vs_CBDabs.csv"), row.names = FALSE)

# Plot helpers
plot_corr <- function(df, shift_col, absorp_col, method, label, absorption_type, x_label, y_label = "ΔGlx") {
  dfp <- df %>% dplyr::select(Participant_ID, AgeGroup, dplyr::all_of(c(shift_col, absorp_col))) %>% tidyr::drop_na()
  ggplot(dfp, aes_string(x = absorp_col, y = shift_col#, color = "AgeGroup", shape = "AgeGroup"
                         )) +
    geom_point(alpha = 0.85, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
    labs(
      title = paste0("ΔGlx vs ", absorption_type, " Absorption (", label, ")"),
      x = x_label,
      y = y_label#,
      #color = "Age Group",
      #shape = "Age Group"
    ) +
    theme_classic()
}

# Create and save plots
plot_thc <- plot_corr(glx_absorption, "THC_minus_placebo", "T2_THCconc_THC_minusPLA", corr_thc$method[1], "THC - Placebo", absorption_type = "THC", x_label = "THC Concentration")
plot_thccbd <- plot_corr(glx_absorption, "THCCBD_minus_placebo", "T2_THCconc_THCCBD_minusPLA", corr_thccbd$method[1], "THC+CBD - Placebo", absorption_type = "THC", x_label = "THC Concentration")
plot_thccbd_cbd <- plot_corr(glx_absorption, "THCCBD_minus_placebo", "T2_CBD_conc_THC_CBD", corr_cbd$method[1], "THC+CBD - Placebo", absorption_type = "CBD", x_label = "CBD Concentration")

print(plot_thc)
print(plot_thccbd)
print(plot_thccbd_cbd)
ggplot2::ggsave(filename = file.path(out_dir, "Correlation_THCshift_vs_THCabs.png"), plot = plot_thc, width = 7, height = 4.5, dpi = 300)
ggplot2::ggsave(filename = file.path(out_dir, "Correlation_THCCBDshift_vs_THCCBDabs.png"), plot = plot_thccbd, width = 7, height = 4.5, dpi = 300)
ggplot2::ggsave(filename = file.path(out_dir, "Correlation_THCCBDshift_vs_CBDabs.png"), plot = plot_thccbd_cbd, width = 7, height = 4.5, dpi = 300)

# Combined correlation plot: THC plots on top row, CBD plot on bottom row (with a spacer row for extra breathing room)
bottom_row <- patchwork::plot_spacer() + plot_thccbd_cbd + patchwork::plot_spacer() + patchwork::plot_layout(widths = c(1, 2, 1))
plot_corr_combined <- (plot_thc | plot_thccbd) / patchwork::plot_spacer() / bottom_row + patchwork::plot_layout(heights = c(1, 0.1, 1))
print(plot_corr_combined)
ggplot2::ggsave(filename = file.path(out_dir, "Correlation_combined.png"), plot = plot_corr_combined, width = 10, height = 8, dpi = 300)

# Permutation histograms with observed r overlay
hist_thc <- ggplot(perm_thc$perm_df, aes(x = r_perm)) +
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  geom_vline(xintercept = perm_thc$r_obs, color = "firebrick", size = 1) +
  labs(title = "Permutation null: THC shift vs THC absorption",
       x = "r (permuted)", y = "Count",
       subtitle = paste0("Observed r = ", round(perm_thc$r_obs, 3), ", p_perm = ", round(perm_thc$p_perm, 3))) +
  theme_classic()

hist_thccbd <- ggplot(perm_thccbd$perm_df, aes(x = r_perm)) +
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  geom_vline(xintercept = perm_thccbd$r_obs, color = "firebrick", size = 1) +
  labs(title = "Permutation null: THC+CBD shift vs THC+CBD absorption",
       x = "r (permuted)", y = "Count",
       subtitle = paste0("Observed r = ", round(perm_thccbd$r_obs, 3), ", p_perm = ", round(perm_thccbd$p_perm, 3))) +
  theme_classic()

hist_cbd <- ggplot(perm_cbd$perm_df, aes(x = r_perm)) +
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  geom_vline(xintercept = perm_cbd$r_obs, color = "firebrick", size = 1) +
  labs(title = "Permutation null: THC+CBD shift vs CBD absorption",
    x = "r (permuted)", y = "Count",
    subtitle = paste0("Observed r = ", round(perm_cbd$r_obs, 3), ", p_perm = ", round(perm_cbd$p_perm, 3))) +
  theme_classic()

print(hist_thc)
print(hist_thccbd)
print(hist_cbd)
ggplot2::ggsave(filename = file.path(out_dir, "PermHist_THCshift_vs_THCabs.png"), plot = hist_thc, width = 7, height = 4.5, dpi = 300)
ggplot2::ggsave(filename = file.path(out_dir, "PermHist_THCCBDshift_vs_THCCBDabs.png"), plot = hist_thccbd, width = 7, height = 4.5, dpi = 300)
ggplot2::ggsave(filename = file.path(out_dir, "PermHist_THCCBDshift_vs_CBDabs.png"), plot = hist_cbd, width = 7, height = 4.5, dpi = 300)

# Bootstrap histograms with observed r and CI overlay
boot_hist_thc <- ggplot(boot_thc$boot_df, aes(x = r_boot)) +
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  geom_vline(xintercept = corr_thc$r[1], color = "firebrick", size = 1) +
  geom_vline(xintercept = boot_thc$ci_low, linetype = "dashed") +
  geom_vline(xintercept = boot_thc$ci_high, linetype = "dashed") +
  labs(title = "Bootstrap: THC shift vs THC absorption",
       x = "r (bootstrap)", y = "Count") +
  theme_classic()

boot_hist_thccbd <- ggplot(boot_thccbd$boot_df, aes(x = r_boot)) +
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  geom_vline(xintercept = corr_thccbd$r[1], color = "firebrick", size = 1) +
  geom_vline(xintercept = boot_thccbd$ci_low, linetype = "dashed") +
  geom_vline(xintercept = boot_thccbd$ci_high, linetype = "dashed") +
  labs(title = "Bootstrap: THC+CBD shift vs THC+CBD absorption",
       x = "r (bootstrap)", y = "Count") +
  theme_classic()

boot_hist_cbd <- ggplot(boot_cbd$boot_df, aes(x = r_boot)) +
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  geom_vline(xintercept = corr_cbd$r[1], color = "firebrick", size = 1) +
  geom_vline(xintercept = boot_cbd$ci_low, linetype = "dashed") +
  geom_vline(xintercept = boot_cbd$ci_high, linetype = "dashed") +
  labs(title = "Bootstrap: THC+CBD shift vs CBD absorption",
       x = "r (bootstrap)", y = "Count") +
  theme_classic()

print(boot_hist_thc)
print(boot_hist_thccbd)
print(boot_hist_cbd)
ggplot2::ggsave(filename = file.path(out_dir, "BootstrapHist_THCshift_vs_THCabs.png"), plot = boot_hist_thc, width = 7, height = 4.5, dpi = 300)
ggplot2::ggsave(filename = file.path(out_dir, "BootstrapHist_THCCBDshift_vs_THCCBDabs.png"), plot = boot_hist_thccbd, width = 7, height = 4.5, dpi = 300)
ggplot2::ggsave(filename = file.path(out_dir, "BootstrapHist_THCCBDshift_vs_CBDabs.png"), plot = boot_hist_cbd, width = 7, height = 4.5, dpi = 300)

# Spaghetti plots of Glx trajecories per subject (Condition on x-axis)
spag_cond <- merged_df_clean_final %>%
  dplyr::mutate(
    condition = factor(condition, levels = c("Placebo", "THC", "THC+CBD")),
    x_idx = as.integer(condition),
    Participant_ID = as.character(Participant_ID)
  ) %>%
  dplyr::arrange(Participant_ID, x_idx) %>%
  dplyr::group_by(Participant_ID) %>%
  dplyr::mutate(
    x_prev = dplyr::lag(x_idx),
    y_prev = dplyr::lag(Glx_adjusted),
    delta  = Glx_adjusted - dplyr::lag(Glx_adjusted),
    change_dir = dplyr::case_when(
      !is.na(delta) & delta > 0  ~ "Increase",
      !is.na(delta) & delta <= 0 ~ "Decrease",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::ungroup() %>%
  tidyr::drop_na(condition, Glx_adjusted)

# Build per-segment deltas to color increases/decreases
spag_segments <- spag_cond %>%
  dplyr::filter(!is.na(x_prev), !is.na(y_prev), !is.na(change_dir))

# Plot points (for all conditions per subject)
gg_spaghetti_pts <- ggplot(spag_cond, aes(x = x_idx, y = Glx_adjusted,
                                          group = Participant_ID)) +
  geom_point(aes(color = change_dir), size = 1.6, alpha = 0.85, na.rm = TRUE) +
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("Placebo", "THC", "THC+CBD")) +
  scale_color_manual(values = c("Increase" = "red", "Decrease" = "blue"), na.translate = FALSE) +
  labs(title = "Subject trajectories of Glx_adjusted across conditions",
       x = "Condition", y = "Glx_adjusted", color = "Change") +
  theme_classic()

# Overlay colored segments (red = increase, blue = decrease)
gg_spaghetti <- gg_spaghetti_pts +
  geom_segment(data = spag_segments,
               aes(x = x_prev, xend = x_idx,
                   y = y_prev, yend = Glx_adjusted,
                   color = change_dir),
               alpha = 0.6, size = 0.8)

print(gg_spaghetti)
ggplot2::ggsave(filename = file.path(out_dir, "Spaghetti_Glx_by_condition_subject.png"),
                plot = gg_spaghetti, width = 8, height = 5, dpi = 300)

gg_spaghetti_age <- gg_spaghetti + facet_wrap(~ AgeGroup, ncol = 1)
print(gg_spaghetti_age)
ggplot2::ggsave(filename = file.path(out_dir, "Spaghetti_Glx_by_condition_subject_facet_age.png"),
                plot = gg_spaghetti_age, width = 8, height = 6, dpi = 300)
