# ==============================================================================
# Part 1: Data Preprocessing, Cohort Selection, and Predictor Pool Setup
# Multidomain Risk Prediction Model for Hip Fracture (CHARLS Cohort)
# ==============================================================================

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
library(haven)
library(dplyr)

# ------------------------------------------------------------------------------
# Step 1: Load Longitudinal Dataset
# Note: Ensure 'charls.dta' is placed in the './data/' directory.
# ------------------------------------------------------------------------------
data_path <- "./data/charls.dta"
if (!file.exists(data_path)) {
  stop("Dataset not found at target path: ", data_path)
}

charls_raw <- read_dta(data_path)

# Verify core identifiers and wave structure
required_vars <- c("ID", "wave", "iwy", "age", "hip", "province")
missing_vars <- setdiff(required_vars, names(charls_raw))
if (length(missing_vars) > 0) {
  stop("Missing critical variables in raw data: ", paste(missing_vars, collapse = ", "))
}

# ------------------------------------------------------------------------------
# Step 2: Define 2011 Baseline Eligible Cohort
# Criteria: Wave 1 (2011), Age >= 60, No baseline hip fracture history
# ------------------------------------------------------------------------------
baseline_2011 <- charls_raw %>%
  filter(wave == 1)

baseline_eligible <- baseline_2011 %>%
  filter(
    !is.na(age),
    age >= 60,
    !is.na(hip),
    hip == 0
  )

message(sprintf("Initial 2011 Baseline N = %d", nrow(baseline_2011)))
message(sprintf("Eligible Baseline Cohort N = %d", nrow(baseline_eligible)))

# ------------------------------------------------------------------------------
# Step 3: Construct Interval-Censored and Follow-Up Outcome Variables
# ------------------------------------------------------------------------------
eligible_ids <- baseline_eligible$ID

followup_long <- charls_raw %>%
  filter(ID %in% eligible_ids) %>%
  arrange(ID, iwy)

# Identify first incident report and overall follow-up boundary
interval_outcome <- followup_long %>%
  group_by(ID) %>%
  summarise(
    baseline_year = min(iwy[wave == 1], na.rm = TRUE),
    n_followup_hip = sum(wave > 1 & !is.na(hip)),
    event = any(wave > 1 & hip == 1, na.rm = TRUE),
    event_right_year = ifelse(
      event,
      min(iwy[wave > 1 & hip == 1], na.rm = TRUE),
      NA_real_
    ),
    last_followup_year = ifelse(
      n_followup_hip > 0,
      max(iwy[wave > 1 & !is.na(hip)], na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  filter(n_followup_hip > 0)

# Determine last event-free interview year prior to event onset
last_negative <- followup_long %>%
  inner_join(
    interval_outcome %>% select(ID, event_right_year),
    by = "ID"
  ) %>%
  filter(
    hip == 0,
    !is.na(iwy),
    is.na(event_right_year) | iwy < event_right_year
  ) %>%
  group_by(ID) %>%
  summarise(
    event_left_year = max(iwy),
    .groups = "drop"
  )

# Calculate interval bounds and censor times relative to baseline year
interval_outcome <- interval_outcome %>%
  left_join(last_negative, by = "ID") %>%
  mutate(
    left_time = event_left_year - baseline_year,
    right_time = ifelse(event, event_right_year - baseline_year, Inf),
    censor_time = last_followup_year - baseline_year
  )

# ------------------------------------------------------------------------------
# Step 4: Generate Primary Analysis (Cox Midpoint) and Sensitivity Outcomes
# ------------------------------------------------------------------------------
outcome_cox <- interval_outcome %>%
  mutate(
    followup_time = ifelse(event, (left_time + right_time) / 2, censor_time),
    event_cox = as.integer(event)
  ) %>%
  select(ID, followup_time, event_cox)

outcome_interval <- interval_outcome %>%
  transmute(
    ID,
    event = as.integer(event),
    interval_left = ifelse(event, left_time, censor_time),
    interval_right = ifelse(event, right_time, Inf)
  )

# Merge back into single master analysis dataset
analysis_data <- baseline_eligible %>%
  inner_join(outcome_cox, by = "ID") %>%
  inner_join(outcome_interval %>% select(ID, interval_left, interval_right), by = "ID")

message(sprintf("Final Analytical Cohort N = %d", nrow(analysis_data)))
message(sprintf("Total Incident Hip Fractures = %d", sum(analysis_data$event_cox)))

# ------------------------------------------------------------------------------
# Step 5: Screen and Refine Candidate Predictor Pool
# ------------------------------------------------------------------------------

# 1. Structural exclusions: IDs, outcome information, geographic stratifiers, raw repeated metrics
exclude_structural <- c(
  "ID", "wave", "householdID", "communityID", "iwy", "iwm", "hip",
  "event_cox", "followup_time", "interval_left", "interval_right",
  "province", "city",
  "wspeed1", "wspeed2", "systo1", "systo2", "systo3",
  "diasto1", "diasto2", "diasto3", "pulse1", "pulse2", "pulse3",
  "lgrip1", "lgrip2", "rgrip1", "rgrip2", "puff1", "puff2", "puff3"
)

candidate_pool <- setdiff(names(analysis_data), exclude_structural)

# 2. Exclude predictors with missingness > 50%
missing_pct <- sapply(analysis_data[candidate_pool], function(x) mean(is.na(x)) * 100)
exclude_high_missing <- names(missing_pct[missing_pct > 50])
candidate_pool <- setdiff(candidate_pool, exclude_high_missing)

# 3. Exclude redundant/derived variables (keep composite or more accurate indicators)
exclude_redundant <- c(
  "hibpe",          # Retain hibpe_1 (includes physical/blood exam data)
  "diabe",          # Retain diabe_1 (includes physical/blood exam data)
  "memeory",        # Retain total_cognition
  "executive",      # Retain total_cognition
  "tyg_bmi",        # Retain individual bmi + tyg
  "chronic",        # Retain specific chronic conditions
  "act_1", "act_2", "act_3", "act_4", 
  "act_5", "act_6", "act_7", "act_8" # Retain comprehensive social1-social11
)
candidate_pool <- setdiff(candidate_pool, exclude_redundant)

# 4. Construct pain burden score (sum of 15 pain sites) and replace individual site items
pain_vars <- paste0("da042s", 1:15)
if (all(pain_vars %in% names(analysis_data))) {
  analysis_data$pain_burden <- apply(
    analysis_data[, pain_vars],
    1,
    function(x) ifelse(all(is.na(x)), NA, sum(x == 1, na.rm = TRUE))
  )
  candidate_pool <- setdiff(candidate_pool, pain_vars)
  candidate_pool <- unique(c(candidate_pool, "pain_burden"))
}

# 5. Correct logical structural zeros in outpatient cost data
analysis_data <- analysis_data %>%
  mutate(
    oopdoc1m = ifelse(doctor == 0 & is.na(oopdoc1m), 0, oopdoc1m),
    totdoc1m = ifelse(doctor == 0 & is.na(totdoc1m), 0, totdoc1m)
  )

# Exclude highly skewed, healthcare-system dependent medical cost variables
exclude_costs <- c("oophos1y", "tothos1y", "oopdoc1m", "totdoc1m")
candidate_pool <- setdiff(candidate_pool, exclude_costs)

# ------------------------------------------------------------------------------
# Step 6: Stratify Candidate Pools (Basic Model vs. Extended Model)
# ------------------------------------------------------------------------------
lab_vars <- grep("^bl_|^tyg$|^bloodweight$", candidate_pool, value = TRUE)

candidate_basic <- setdiff(candidate_pool, lab_vars)
candidate_extended <- candidate_pool

message(sprintf("Candidate Predictors - Basic Model: %d", length(candidate_basic)))
message(sprintf("Candidate Predictors - Extended Model: %d", length(candidate_extended)))

# Save cleaned dataset for subsequent modeling steps
saveRDS(analysis_data, file = "./data/analysis_data_cleaned.rds")