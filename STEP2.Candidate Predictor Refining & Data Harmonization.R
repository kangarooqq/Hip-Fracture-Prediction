# ==============================================================================
# Part 2: Feature Engineering, Data Harmonization, & Class Standardization
# Multidomain Risk Prediction Model for Hip Fracture (CHARLS Cohort)
# ==============================================================================

library(dplyr)
library(haven)

# ------------------------------------------------------------------------------
# Step 21: Define Initial Dual-Model Candidate Pools
# ------------------------------------------------------------------------------

# Basic Model: Non-laboratory candidate variables
candidate_basic <- setdiff(candidate_basic, lab_vars)

# Extended Model: Basic Model + Laboratory Biomarkers
candidate_extended <- unique(c(candidate_basic, lab_vars))

message(sprintf("Initial Basic Candidate Pool: %d predictors", length(candidate_basic)))
message(sprintf("Initial Extended Candidate Pool: %d predictors", length(candidate_extended)))

# ------------------------------------------------------------------------------
# Step 22 - 28: Anthropometric Variable Cleaning & Plausibility Standardisation
# ------------------------------------------------------------------------------

# Clean physiologically impossible anthropometric measurements:
# 1. Height < 1.2m or > 2.0m -> NA
# 2. Weight < 10kg -> NA
# 3. Waist circumference < 40cm -> NA
# 4. Recalculate BMI based on harmonized metrics to eliminate extreme anomalies

analysis_data <- analysis_data %>%
  mutate(
    mheight_clean = ifelse(mheight < 1.2 | mheight > 2.0, NA_real_, mheight),
    mweight_clean = ifelse(mweight < 10, NA_real_, mweight),
    mwaist_clean  = ifelse(mwaist < 40, NA_real_, mwaist),
    bmi_clean     = mweight_clean / (mheight_clean^2)
  )

# Update Candidate Predictor Pools: Replace raw anthropometrics with cleaned variables
anthro_raw <- c("mheight", "mweight", "mwaist", "bmi")
anthro_clean <- c("mheight_clean", "mwaist_clean", "bmi_clean")

candidate_basic <- setdiff(candidate_basic, anthro_raw)
candidate_basic <- unique(c(candidate_basic, anthro_clean))

candidate_extended <- setdiff(candidate_extended, anthro_raw)
candidate_extended <- unique(c(candidate_extended, anthro_clean))

# ------------------------------------------------------------------------------
# Step 29 - 31: Health Insurance Redundancy Reduction
# ------------------------------------------------------------------------------

# Consolidate 9 granular insurance variables (ea001s1-ea001s11) into primary 'ins'
ea_vars <- c("ea001s1", "ea001s2", "ea001s3", "ea001s4", 
             "ea001s5", "ea001s6", "ea001s7", "ea001s8", "ea001s11")

candidate_basic <- setdiff(candidate_basic, ea_vars)
candidate_extended <- setdiff(candidate_extended, ea_vars)

# Ensure composite 'ins' variable is retained
candidate_basic <- unique(c(candidate_basic, "ins"))
candidate_extended <- unique(c(candidate_extended, "ins"))

# ------------------------------------------------------------------------------
# Step 32 - 33: Social Activity Dimension Reduction
# ------------------------------------------------------------------------------

# Collapse sparse itemized social activities (social1-social11) into a 3-level score
social_vars <- paste0("social", 1:11)

analysis_data$social_count <- apply(
  analysis_data[, social_vars],
  1,
  function(x) ifelse(all(is.na(x)), NA_real_, sum(x == 1, na.rm = TRUE))
)

analysis_data$social_level <- cut(
  analysis_data$social_count,
  breaks = c(-Inf, 0, 1, Inf),
  labels = c("0", "1", "2+"),
  right = TRUE
)

# Update pools with composite social activity level
candidate_basic <- setdiff(candidate_basic, social_vars)
candidate_basic <- unique(c(candidate_basic, "social_level"))

candidate_extended <- setdiff(candidate_extended, social_vars)
candidate_extended <- unique(c(candidate_extended, "social_level"))

# ------------------------------------------------------------------------------
# Step 34 - 35: Lifestyle Behaviors Harmonization (Smoking & Drinking)
# ------------------------------------------------------------------------------

# Construct standardized 3-category behavioral indicators (never, former, current)
analysis_data <- analysis_data %>%
  mutate(
    drinking_status = factor(
      case_when(
        drinkl == 1 ~ "current",
        drinkev == 1 & drinkl == 0 ~ "former",
        drinkev == 0 & drinkl == 0 ~ "never",
        TRUE ~ NA_character_
      ),
      levels = c("never", "former", "current")
    ),
    smoking_status = factor(
      case_when(
        smoken == 1 ~ "current",
        smokev == 1 & smoken == 0 ~ "former",
        smokev == 0 & smoken == 0 ~ "never",
        TRUE ~ NA_character_
      ),
      levels = c("never", "former", "current")
    )
  )

old_behavior_vars <- c("drinkev", "drinkl", "smokev", "smoken")

candidate_basic <- setdiff(candidate_basic, old_behavior_vars)
candidate_basic <- unique(c(candidate_basic, "drinking_status", "smoking_status"))

candidate_extended <- setdiff(candidate_extended, old_behavior_vars)
candidate_extended <- unique(c(candidate_extended, "drinking_status", "smoking_status"))

# ------------------------------------------------------------------------------
# Step 36 - 38: Healthcare Utilization Restructuring
# ------------------------------------------------------------------------------

# Categorize right-skewed healthcare utilization counts into clinical categories
analysis_data <- analysis_data %>%
  mutate(
    hospital_use = factor(
      case_when(
        hospital_time >= 2 ~ "2+",
        hospital_time == 1 ~ "1",
        hospital_time == 0 & hospital == 1 ~ "1", # Resolve minor logic conflict
        hospital_time == 0 ~ "0",
        TRUE ~ NA_character_
      ),
      levels = c("0", "1", "2+")
    ),
    doctor_use = factor(
      case_when(
        doctor_time >= 5 ~ "5+",
        doctor_time >= 2 & doctor_time <= 4 ~ "2-4",
        doctor_time == 1 ~ "1",
        doctor_time == 0 ~ "0",
        TRUE ~ NA_character_
      ),
      levels = c("0", "1", "2-4", "5+")
    )
  )

old_medical_vars <- c("hospital", "hospital_time", "hspnite", "doctor", "doctor_time")

candidate_basic <- setdiff(candidate_basic, old_medical_vars)
candidate_basic <- unique(c(candidate_basic, "hospital_use", "doctor_use"))

candidate_extended <- setdiff(candidate_extended, old_medical_vars)
candidate_extended <- unique(c(candidate_extended, "hospital_use", "doctor_use"))

# ------------------------------------------------------------------------------
# Step 39 - 41: Screen Ultra-Low Prevalence Binary Predictors
# ------------------------------------------------------------------------------

# Remove 'hear_aid' (ultra-sparse positive events and redundant with 'hear' status)
# Retain 'cancre' due to high clinical relevance
candidate_basic <- setdiff(candidate_basic, "hear_aid")
candidate_extended <- setdiff(candidate_extended, "hear_aid")

# ------------------------------------------------------------------------------
# Step 42 - 43: Remove Predictors with High Missingness (>20%) & Low Information Gain
# ------------------------------------------------------------------------------

# Remove income_total, pension, and nation due to high missingness and low predictive gain
exclude_high_missing_low_value <- c("income_total", "pension", "nation")

candidate_basic <- setdiff(candidate_basic, exclude_high_missing_low_value)
candidate_extended <- setdiff(candidate_extended, exclude_high_missing_low_value)

# Note: Retained 'total_cognition' and 'wspeed' for MICE imputation due to critical clinical relevance

# ------------------------------------------------------------------------------
# Step 44 - 46: Communication Cost Removal
# ------------------------------------------------------------------------------

# Remove non-essential 'communication' variable containing non-standard missing codes
candidate_basic <- setdiff(candidate_basic, "communication")
candidate_extended <- setdiff(candidate_extended, "communication")

# ------------------------------------------------------------------------------
# Step 47 - 49: Stata Tagged Missing Values Normalization
# ------------------------------------------------------------------------------

# Standardize Stata tagged NAs (e.g., .a, .b) to standard R NAs
tagged_vars <- c("wspeed", "lgrip", "rgrip")
for (v in tagged_vars) {
  if (v %in% names(analysis_data)) {
    x <- analysis_data[[v]]
    x[haven::is_tagged_na(x)] <- NA
    analysis_data[[v]] <- x
  }
}

# ------------------------------------------------------------------------------
# Step 50 - 52: Factor Type Standardization (Binary, Ordinal, & Unordered)
# ------------------------------------------------------------------------------

# 1. Standardize Binary Predictors (0/1 -> Factor)
binary_vars <- candidate_basic[
  sapply(candidate_basic, function(v) {
    x <- analysis_data[[v]]
    vals <- sort(unique(x[!is.na(x)]))
    length(vals) == 2 && all(vals %in% c(0, 1))
  })
]

for (v in binary_vars) {
  analysis_data[[v]] <- factor(as.numeric(analysis_data[[v]]), levels = c(0, 1))
}

# 2. Standardize Ordinal Predictors (Ordered Factor)
ordinal_vars <- c("hope", "edu", "srh", "satlife", "eyesight_distance", "eyesight_close", "hear")

for (v in intersect(ordinal_vars, names(analysis_data))) {
  # Preserve natural numeric ordering
  vals <- sort(unique(analysis_data[[v]][!is.na(analysis_data[[v]])]))
  analysis_data[[v]] <- factor(analysis_data[[v]], levels = vals, ordered = TRUE)
}

# 3. Standardize Unordered Categorical Predictors
unordered_vars <- c("glass")
for (v in intersect(unordered_vars, names(analysis_data))) {
  vals <- sort(unique(analysis_data[[v]][!is.na(analysis_data[[v]])]))
  analysis_data[[v]] <- factor(analysis_data[[v]], levels = vals, ordered = FALSE)
}

# Clean workspace attributes if necessary
message(sprintf("Refined Candidate Predictors - Basic Model: %d", length(candidate_basic)))
message(sprintf("Refined Candidate Predictors - Extended Model: %d", length(candidate_extended)))