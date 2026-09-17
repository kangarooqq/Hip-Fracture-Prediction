# ==============================================================================
# Part 4: Design Matrix Construction & Hyperparameter Tuning (Elastic Net)
# Multidomain Risk Prediction Model for Hip Fracture (CHARLS Cohort)
# ==============================================================================

library(dplyr)
library(glmnet)
library(survival)

# ------------------------------------------------------------------------------
# Step 60: Construct Unified Treatment-Coded Design Matrix
# Note: Forced conversion of ordered factors to standard unordered factors to 
# prevent model.matrix() from generating polynomial contrasts (.L, .Q, .C).
# ------------------------------------------------------------------------------

design_data <- basic_long %>%
  select(all_of(candidate_basic)) %>%
  mutate(
    across(
      where(is.ordered),
      ~ factor(as.character(.x), levels = levels(.x), ordered = FALSE)
    )
  )

# Verify zero remaining ordered factors
remaining_ordered <- sum(sapply(design_data, is.ordered))
if (remaining_ordered > 0) {
  stop("Error: Ordered factors detected prior to design matrix generation.")
}

# Generate treatment-coded dummy predictor matrix (excluding intercept)
X_basic <- model.matrix(~ ., data = design_data)[, -1, drop = FALSE]

# Quality Checks on Design Matrix
if (anyNA(X_basic) || any(!is.finite(X_basic))) {
  stop("Error: Design matrix contains NA or infinite values.")
}

zero_var_cols <- colnames(X_basic)[apply(X_basic, 2, function(x) length(unique(x)) <= 1)]
if (length(zero_var_cols) > 0) {
  stop("Error: Zero-variance columns detected: ", paste(zero_var_cols, collapse = ", "))
}

message(sprintf("Design Matrix Generated Successfully: %d observations x %d features", 
                nrow(X_basic), ncol(X_basic)))

# ------------------------------------------------------------------------------
# Step 61: Construct Subject-Level Stratified 10-Fold Cross-Validation Split
# Note: Ensures all 10 imputed versions of the same subject belong to the same fold.
# ------------------------------------------------------------------------------

set.seed(20260913)
K_folds <- 10
n_subject <- nrow(analysis_data)
fold_subject <- integer(n_subject)

# Stratify subject IDs by primary event outcome
event_idx <- which(analysis_data$event_cox == 1)
nonevent_idx <- which(analysis_data$event_cox == 0)

fold_subject[event_idx] <- sample(rep(1:K_folds, length.out = length(event_idx)))
fold_subject[nonevent_idx] <- sample(rep(1:K_folds, length.out = length(nonevent_idx)))

# Map subject-level fold IDs to the stacked longitudinal imputed dataset
foldid_basic <- fold_subject[as.integer(basic_long$.id)]

# Validate subject leakage across folds
fold_leakage <- sum(tapply(foldid_basic, basic_long$.id, function(x) length(unique(x))) > 1)
if (fold_leakage > 0) {
  stop("Error: Data leakage detected across CV folds for identical subjects.")
}

message("Subject-level stratified 10-fold CV partitions established successfully.")

# ------------------------------------------------------------------------------
# Step 62 - 64: Elastic Net Grid Search & 1-SE Model Selection Rule
# ------------------------------------------------------------------------------

# Define Cox outcome
y_basic <- Surv(basic_long$followup_time, basic_long$event_cox)

alpha_grid <- c(0.00, 0.25, 0.50, 0.75, 1.00)
cv_results_list <- vector("list", length(alpha_grid))

alpha_summary <- data.frame(
  alpha = alpha_grid,
  best_C = NA_real_,
  SE_at_best = NA_real_,
  lambda_min = NA_real_,
  n_nonzero_min = NA_integer_
)

message("Executing Elastic Net grid search across alpha parameters...")

for (i in seq_along(alpha_grid)) {
  a <- alpha_grid[i]
  set.seed(20260913)
  
  fit_cv <- cv.glmnet(
    x = X_basic,
    y = y_basic,
    family = "cox",
    alpha = a,
    foldid = foldid_basic,
    type.measure = "C",
    standardize = TRUE,
    grouped = TRUE,
    nlambda = 50
  )
  
  cv_results_list[[i]] <- fit_cv
  
  idx_min <- which.min(abs(fit_cv$lambda - fit_cv$lambda.min))
  coef_min <- as.matrix(coef(fit_cv, s = "lambda.min"))
  
  alpha_summary$best_C[i] <- fit_cv$cvm[idx_min]
  alpha_summary$SE_at_best[i] <- fit_cv$cvsd[idx_min]
  alpha_summary$lambda_min[i] <- fit_cv$lambda.min
  alpha_summary$n_nonzero_min[i] <- sum(coef_min != 0)
}

# Determine optimal alpha via 1-SE parsimony criterion
overall_best_idx <- which.max(alpha_summary$best_C)
max_C <- alpha_summary$best_C[overall_best_idx]
max_SE <- alpha_summary$SE_at_best[overall_best_idx]
C_1se_threshold <- max_C - max_SE

# Filter models performing within 1-SE threshold and choose highest sparsity
eligible_alphas <- alpha_summary %>%
  filter(best_C >= C_1se_threshold)

selected_model <- eligible_alphas %>%
  filter(n_nonzero_min == min(n_nonzero_min)) %>%
  slice(1)

selected_alpha <- selected_model$alpha
selected_lambda <- selected_model$lambda_min

message("--------------------------------------------------")
message(sprintf("Hyperparameter Optimization Results:"))
message(sprintf("Maximum Cross-Validated C-index: %.4f (SE: %.4f)", max_C, max_SE))
message(sprintf("Selected Alpha (1-SE Rule Criteria): %.2f", selected_alpha))
message(sprintf("Selected Lambda (lambda.min): %.6f", selected_lambda))
message(sprintf("Non-zero Feature Count at Selected Model: %d", selected_model$n_nonzero_min))
message("--------------------------------------------------")

# Save hyperparameter tuning artifacts and design matrix
saveRDS(list(
  X_basic = X_basic,
  foldid_basic = foldid_basic,
  selected_alpha = selected_alpha,
  selected_lambda = selected_lambda,
  cv_summary = alpha_summary
), file = "./data/elastic_net_tuning_results.rds")