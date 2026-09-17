# ==============================================================================
# Step 7 Dynamic Predictive Modeling Pipeline: Inner CV Selection & Internal Validation
# ==============================================================================

# Load required libraries
library(glmnet)
library(survival)
library(dplyr)
library(timeROC)

# Set seed for reproducible results across operations
set.seed(20260914)

# Define core variable list for evaluation
final_vars_10 <- c(
  "fall_down", "age", "adlab_c", "arthre", "cesd10",
  "hear", "total_cognition", "wspeed", "puff", "teeth"
)

final_formula <- as.formula(
  paste("Surv(followup_time, event_cox) ~", paste(final_vars_10, collapse = " + "))
)

# ------------------------------------------------------------------------------
# SECTION 1: Outer Fold 1 - Inner CV Tuning (Alpha & Lambda)
# ------------------------------------------------------------------------------

alpha_grid <- c(0, 0.25, 0.50, 0.75, 1)
cv_alpha_models <- vector("list", length(alpha_grid))

alpha_summary <- data.frame(
  alpha = alpha_grid,
  C_min = NA_real_,
  SE_min = NA_real_,
  lambda_min = NA_real_,
  n_nonzero_min = NA_integer_
)

# Cross-validation over alpha grid
for (j in seq_along(alpha_grid)) {
  a <- alpha_grid[j]
  
  # Standardized Cox Elastic-Net model
  cvfit <- cv.glmnet(
    x = X_train1,
    y = y_train1,
    family = "cox",
    alpha = a,
    foldid = inner_foldid1,
    type.measure = "C",
    standardize = TRUE
  )
  
  cv_alpha_models[[j]] <- cvfit
  idx_min <- which.min(abs(cvfit$lambda - cvfit$lambda.min))
  
  alpha_summary$C_min[j] <- cvfit$cvm[idx_min]
  alpha_summary$SE_min[j] <- cvfit$cvsd[idx_min]
  alpha_summary$lambda_min[j] <- cvfit$lambda.min
  alpha_summary$n_nonzero_min[j] <- sum(as.numeric(coef(cvfit, s = "lambda.min")) != 0)
}

# Select best model using 1-SE rule (sparse-first strategy)
best_idx <- which.max(alpha_summary$C_min)
threshold_1se <- alpha_summary$C_min[best_idx] - alpha_summary$SE_min[best_idx]
eligible <- which(alpha_summary$C_min >= threshold_1se)

min_nzero <- min(alpha_summary$n_nonzero_min[eligible])
eligible_sparse <- eligible[alpha_summary$n_nonzero_min[eligible] == min_nzero]
chosen_idx <- eligible_sparse[which.max(alpha_summary$alpha[eligible_sparse])]

chosen_alpha <- alpha_summary$alpha[chosen_idx]

# ------------------------------------------------------------------------------
# SECTION 2: Backward CV Elimination on Training Subset (15 -> 10 Variables)
# ------------------------------------------------------------------------------

# Helper function: Calculate 10-fold x 10-MI concordance index
cv_C_MI <- function(vars) {
  C_mat <- matrix(NA_real_, nrow = 10, ncol = 10)
  
  for (m in 1:10) {
    dat_m <- train1_long[train1_long$.imp == m, ]
    dat_m[vars] <- lapply(dat_m[vars], function(x) {
      if (is.ordered(x)) factor(as.character(x), levels = levels(x), ordered = FALSE) else x
    })
    
    formula_m <- as.formula(paste("Surv(followup_time, event_cox) ~", paste(vars, collapse = " + ")))
    fold_m <- unname(inner_fold_map1[as.character(dat_m$ID)])
    
    for (f in 1:10) {
      train_idx <- fold_m != f
      valid_idx <- fold_m == f
      
      fit_mf <- try(coxph(formula_m, data = dat_m[train_idx, ], ties = "efron"), silent = TRUE)
      if (inherits(fit_mf, "try-error") || any(!is.finite(coef(fit_mf)))) next
      
      lp_valid <- predict(fit_mf, newdata = dat_m[valid_idx, ], type = "lp")
      C_mat[m, f] <- concordance(
        Surv(followup_time, event_cox) ~ lp_valid,
        data = dat_m[valid_idx, ],
        reverse = TRUE
      )$concordance
    }
  }
  
  fold_C <- colMeans(C_mat, na.rm = TRUE)
  c(mean_C = mean(fold_C, na.rm = TRUE), SE_C = sd(fold_C, na.rm = TRUE) / sqrt(10))
}

# ------------------------------------------------------------------------------
# SECTION 3: Bootstrap Internal Validation (Discrimination & Global Slope)
# ------------------------------------------------------------------------------

B <- 500
boot_summary <- data.frame(
  imputation = 1:10,
  apparent_C = NA_real_,
  optimism_C = NA_real_,
  corrected_C = NA_real_,
  apparent_slope = NA_real_,
  optimism_slope = NA_real_,
  corrected_slope = NA_real_
)

for (m in 1:10) {
  dat_m <- basic_long[basic_long$.imp == m, ]
  dat_m$hear <- factor(as.character(dat_m$hear), levels = c("1", "2", "3", "4-5"), ordered = FALSE)
  n_m <- nrow(dat_m)
  
  # Fit base model
  fit_orig <- coxph(final_formula, data = dat_m, ties = "efron")
  C_app <- concordance(fit_orig)$concordance
  
  lp_orig <- predict(fit_orig, newdata = dat_m, type = "lp")
  slope_app <- coef(coxph(Surv(followup_time, event_cox) ~ lp_orig, data = dat_m))[["lp_orig"]]
  
  opt_C_b <- numeric(B)
  opt_slope_b <- numeric(B)
  
  for (b in 1:B) {
    boot_idx <- sample(seq_len(n_m), size = n_m, replace = TRUE)
    dat_boot <- dat_m[boot_idx, ]
    
    fit_boot <- try(coxph(final_formula, data = dat_boot, ties = "efron"), silent = TRUE)
    if (inherits(fit_boot, "try-error") || any(!is.finite(coef(fit_boot)))) next
    
    # Evaluate C-Index Optimism
    C_boot <- concordance(fit_boot)$concordance
    lp_test <- predict(fit_boot, newdata = dat_m, type = "lp")
    C_test <- concordance(
      Surv(followup_time, event_cox) ~ lp_test,
      data = dat_m,
      reverse = TRUE
    )$concordance
    opt_C_b[b] <- C_boot - C_test
    
    # Evaluate Slope Optimism
    lp_boot <- predict(fit_boot, newdata = dat_boot, type = "lp")
    s_boot <- coef(coxph(Surv(followup_time, event_cox) ~ lp_boot, data = dat_boot))[["lp_boot"]]
    s_test <- coef(coxph(Surv(followup_time, event_cox) ~ lp_test, data = dat_m))[["lp_test"]]
    opt_slope_b[b] <- s_boot - s_test
  }
  
  # Store MI metrics
  boot_summary$apparent_C[m] <- C_app
  boot_summary$optimism_C[m] <- mean(opt_C_b, na.rm = TRUE)
  boot_summary$corrected_C[m] <- C_app - mean(opt_C_b, na.rm = TRUE)
  
  boot_summary$apparent_slope[m] <- slope_app
  boot_summary$optimism_slope[m] <- mean(opt_slope_b, na.rm = TRUE)
  boot_summary$corrected_slope[m] <- slope_app - mean(opt_slope_b, na.rm = TRUE)
}

# ------------------------------------------------------------------------------
# SECTION 4: Time-Specific IPCW Calibration & Time-Dependent AUC
# ------------------------------------------------------------------------------

eval_time_calibration <- function(fit, eval_data, times = c(3, 5, 7)) {
  eval_data <- as.data.frame(eval_data)
  eval_data$hear <- factor(as.character(eval_data$hear), levels = c("1", "2", "3", "4-5"), ordered = FALSE)
  
  bh <- basehaz(fit, centered = FALSE)
  lp <- predict(fit, newdata = eval_data, type = "lp", reference = "zero")
  
  # Censoring survival distribution G(t)
  fit_G <- survfit(Surv(followup_time, 1 - event_cox) ~ 1, data = eval_data)
  G_time <- fit_G$time
  G_surv <- fit_G$surv
  
  G_at <- function(t) {
    idx <- findInterval(t, G_time)
    ifelse(idx == 0, 1, G_surv[pmax(idx, 1)])
  }
  
  out <- lapply(times, function(t) {
    H0_t <- approx(bh$time, bh$hazard, xout = t, method = "constant", rule = 2, f = 0)$y
    pred_risk <- pmin(pmax(1 - exp(-H0_t * exp(lp)), 1e-6), 1 - 1e-6)
    x_cloglog <- log(-log(1 - pred_risk))
    
    known_event <- eval_data$event_cox == 1 & eval_data$followup_time <= t
    known_nonevent <- eval_data$followup_time > t
    use <- known_event | known_nonevent
    
    y <- as.integer(known_event[use])
    x_use <- x_cloglog[use]
    
    # IPCW Weighting
    w <- numeric(sum(use))
    event_used <- known_event[use]
    w[event_used] <- 1 / G_at(eval_data$followup_time[use][event_used] - 1e-8)
    w[!event_used] <- 1 / G_at(t)
    
    cal_dat <- data.frame(y = y, x = x_use, w = w)
    
    fit_int <- suppressWarnings(glm(y ~ 1 + offset(x), data = cal_dat, family = binomial(link = "cloglog"), weights = w))
    fit_slope <- suppressWarnings(glm(y ~ x, data = cal_dat, family = binomial(link = "cloglog"), weights = w))
    
    c(intercept = unname(coef(fit_int)[1]), slope = unname(coef(fit_slope)["x"]))
  })
  
  do.call(rbind, out)
}

# Diagnostic verification for Time-Dependent AUC (MI = 1)
dat_auc <- basic_long[basic_long$.imp == 1, ]
dat_auc$hear <- factor(as.character(dat_auc$hear), levels = c("1", "2", "3", "4-5"), ordered = FALSE)

fit_auc <- coxph(final_formula, data = dat_auc, ties = "efron")
lp_auc <- predict(fit_auc, newdata = dat_auc, type = "lp")

auc_res <- timeROC(
  T = dat_auc$followup_time,
  delta = dat_auc$event_cox,
  marker = lp_auc,
  cause = 1,
  weighting = "marginal",
  times = c(3, 5, 7),
  iid = TRUE
)

# Output summary metrics
print(boot_summary)
print(data.frame(time = c(3, 5, 7), AUC = round(as.numeric(auc_res$AUC), 4)))