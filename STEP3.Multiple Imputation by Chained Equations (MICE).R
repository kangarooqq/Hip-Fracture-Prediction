# ==============================================================================
# Part 3: Multiple Imputation by Chained Equations (MICE) Workflow
# Multidomain Risk Prediction Model for Hip Fracture (CHARLS Cohort)
# ==============================================================================

library(mice)
library(dplyr)

# ------------------------------------------------------------------------------
# Step 53: Construct Nelson-Aalen Cumulative Hazard Estimator
# Note: Nelson-Aalen hazard estimator and event indicator MUST be included 
# in the imputation model to preserve the target time-to-event outcome structure.
# ------------------------------------------------------------------------------

# Compute Nelson-Aalen estimate for the primary analysis cohort
na_fit <- survival::basehaz(
  survival::coxph(survival::Surv(followup_time, event_cox) ~ 1, data = analysis_data)
)

analysis_data$nelson_aalen <- approx(
  x = na_fit$time, 
  y = na_fit$hazard, 
  xout = analysis_data$followup_time, 
  rule = 2
)$y

# Prepare specific data matrix for MICE workflow
mice_basic_vars <- c(candidate_basic, "event_cox", "nelson_aalen")
mice_basic_data <- analysis_data[, mice_basic_vars]

# ------------------------------------------------------------------------------
# Step 54: Configure Imputation Methods and Predictor Matrix
# ------------------------------------------------------------------------------

# 1. Initialize imputation method vector
meth_basic <- rep("", ncol(mice_basic_data))
names(meth_basic) <- names(mice_basic_data)

# Assign variable-type-specific imputation methods
for (v in candidate_basic) {
  x <- mice_basic_data[[v]]
  
  if (!anyNA(x)) {
    meth_basic[v] <- ""
  } else if (is.ordered(x)) {
    meth_basic[v] <- "polr"      # Proportional odds logistic regression
  } else if (is.factor(x) && nlevels(x) == 2) {
    meth_basic[v] <- "logreg"    # Binary logistic regression
  } else if (is.factor(x) && nlevels(x) > 2) {
    meth_basic[v] <- "polyreg"   # Polytomous logistic regression
  } else if (is.numeric(x) || is.integer(x)) {
    meth_basic[v] <- "pmm"       # Predictive mean matching
  }
}

# Ensure outcome metrics are NOT imputed themselves
meth_basic["event_cox"] <- ""
meth_basic["nelson_aalen"] <- ""

# 2. Setup Predictor Matrix
pred_basic <- make.predictorMatrix(mice_basic_data)

# Self-prediction set to zero
diag(pred_basic) <- 0

# Outcomes must not be imputed by other predictors
pred_basic["event_cox", ] <- 0
pred_basic["nelson_aalen", ] <- 0

# Outcome metrics MUST predict all candidate predictors with missing values
pred_basic[candidate_basic, c("event_cox", "nelson_aalen")] <- 1

# ------------------------------------------------------------------------------
# Step 55 - 56: Formal Execution of Multiple Imputation (Initial Phase)
# ------------------------------------------------------------------------------

message("Starting initial MICE process (m = 10, maxit = 5)...")
set.seed(20260913)

imp_basic <- mice(
  data = mice_basic_data,
  m = 10,
  maxit = 5,
  method = meth_basic,
  predictorMatrix = pred_basic,
  seed = 20260913,
  printFlag = FALSE
)

# Check for logged numerical instability events
if (!is.null(imp_basic$loggedEvents)) {
  warning("MICE logged potential collinearity/prediction warnings:")
  print(imp_basic$loggedEvents)
} else {
  message("Initial MICE run completed smoothly without logged errors.")
}

# ------------------------------------------------------------------------------
# Step 57 - 58: Convergence Assessment & Chain Extension
# ------------------------------------------------------------------------------

message("Extending MICE iterations (additional 5 iterations for convergence)...")
set.seed(20260913)

# Extend iterations from 5 to 10 to ensure stationarity of chain means
imp_basic_final <- mice::mice.mids(
  imp_basic,
  maxit = 5,
  printFlag = FALSE
)

message(sprintf("MICE finalized: %d imputations with %d total iterations per chain.", 
                imp_basic_final$m, imp_basic_final$iteration))

# ------------------------------------------------------------------------------
# Step 59: Export Stacked Long-Format Analysis Dataset
# ------------------------------------------------------------------------------

# Convert multiply imputed mids object into long-format dataframe
basic_long <- complete(
  imp_basic_final,
  action = "long",
  include = FALSE
)

# Synchronize original subject follow-up times into long-format data
basic_long$followup_time <- analysis_data$followup_time[basic_long$.id]

message("--------------------------------------------------")
message(sprintf("Long-Format Imputed Data Created: N = %d rows", nrow(basic_long)))
message(sprintf("Number of Imputations (.imp): %d", length(unique(basic_long$.imp))))
message(sprintf("Total Incident Events Per Imputation: %d", sum(basic_long$event_cox[basic_long$.imp == 1])))
message("--------------------------------------------------")

# Save imputed longitudinal stacked object for downstream modeling (LASSO / Elastic Net)
saveRDS(basic_long, file = "./data/mice_imputed_basic_long.rds")