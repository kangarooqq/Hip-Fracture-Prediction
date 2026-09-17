# ==============================================================================
# Step 8 Dynamic Predictive Modeling Pipeline: Time-Dependent Discrimination, IPCW Brier/IPA Evaluation, and Decision Curve Analysis (DCA)
# ==============================================================================

# Load required libraries
library(survival)
library(timeROC)
library(dplyr)

# Set global seed for reproducible bootstrap operations
set.seed(20260914)

# Define core predictor set and model formula
final_vars_10 <- c(
  "fall_down", "age", "adlab_c", "arthre", "cesd10",
  "hear", "total_cognition", "wspeed", "puff", "teeth"
)

final_formula <- as.formula(
  paste("Surv(followup_time, event_cox) ~", paste(final_vars_10, collapse = " + "))
)

# Configuration parameters
times_eval <- c(3, 5, 7)
n_imputations <- 10
n_bootstraps <- 500

# Helper function: Re-evaluate dynamic IPCW Brier score and IPA on evaluation data
eval_brier_ipa <- function(fit, eval_data, times = c(3, 5, 7)) {
  eval_data <- as.data.frame(eval_data)
  eval_data$hear <- factor(
    as.character(eval_data$hear), 
    levels = c("1", "2", "3", "4-5"), 
    ordered = FALSE
  )
  
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
  
  G_left <- function(t) {
    idx <- findInterval(t - 1e-8, G_time)
    ifelse(idx == 0, 1, G_surv[pmax(idx, 1)])
  }
  
  # Reference null KM distribution
  km <- survfit(Surv(followup_time, event_cox) ~ 1, data = eval_data)
  
  out <- lapply(times, function(t) {
    H0_t <- approx(bh$time, bh$hazard, xout = t, method = "constant", rule = 2, f = 0)$y
    pred <- 1 - exp(-H0_t * exp(lp))
    
    S0_t <- summary(km, times = t, extend = TRUE)$surv
    null_pred <- 1 - S0_t
    
    event_before <- eval_data$event_cox == 1 & eval_data$followup_time <= t
    event_free <- eval_data$followup_time > t
    
    w_event <- rep(0, nrow(eval_data))
    w_free <- rep(0, nrow(eval_data))
    
    w_event[event_before] <- 1 / G_left(eval_data$followup_time[event_before])
    w_free[event_free] <- 1 / G_at(t)
    
    BS_model <- mean(w_event * (1 - pred)^2 + w_free * pred^2)
    BS_null <- mean(w_event * (1 - null_pred)^2 + w_free * null_pred^2)
    
    c(Brier = BS_model, IPA = 1 - BS_model / BS_null)
  })
  
  do.call(rbind, out)
}

# Helper function: Calculate IPCW Net Benefit across risk thresholds
eval_net_benefit <- function(fit, eval_data, t0 = 5, thresholds = seq(0.015, 0.060, by = 0.005)) {
  eval_data <- as.data.frame(eval_data)
  eval_data$hear <- factor(
    as.character(eval_data$hear), 
    levels = c("1", "2", "3", "4-5"), 
    ordered = FALSE
  )
  
  bh <- basehaz(fit, centered = FALSE)
  lp <- predict(fit, newdata = eval_data, type = "lp", reference = "zero")
  
  H0_t <- approx(bh$time, bh$hazard, xout = t0, method = "constant", rule = 2, f = 0)$y
  risk <- 1 - exp(-H0_t * exp(lp))
  
  fit_G <- survfit(Surv(followup_time, 1 - event_cox) ~ 1, data = eval_data)
  G_time <- fit_G$time
  G_surv <- fit_G$surv
  
  G_at <- function(t) {
    idx <- findInterval(t, G_time)
    ifelse(idx == 0, 1, G_surv[pmax(idx, 1)])
  }
  
  G_left <- function(t) {
    idx <- findInterval(t - 1e-8, G_time)
    ifelse(idx == 0, 1, G_surv[pmax(idx, 1)])
  }
  
  event_t <- eval_data$event_cox == 1 & eval_data$followup_time <= t0
  nonevent_t <- eval_data$followup_time > t0
  
  w_event <- rep(0, nrow(eval_data))
  w_nonevent <- rep(0, nrow(eval_data))
  
  w_event[event_t] <- 1 / G_left(eval_data$followup_time[event_t])
  w_nonevent[nonevent_t] <- 1 / G_at(t0)
  
  N <- nrow(eval_data)
  
  sapply(thresholds, function(pt) {
    positive <- risk >= pt
    sum(w_event * positive) / N - sum(w_nonevent * positive) / N * (pt / (1 - pt))
  })
}

# ------------------------------------------------------------------------------
# SECTION 1: Bootstrap Validation for Time-Dependent AUC & Brier / IPA Score
# ------------------------------------------------------------------------------

auc_boot_summary <- data.frame()
brier_boot_summary <- data.frame()

for (m in 1:n_imputations) {
  dat_m <- basic_long[basic_long$.imp == m, ]
  dat_m$hear <- factor(as.character(dat_m$hear), levels = c("1", "2", "3", "4-5"), ordered = FALSE)
  n_m <- nrow(dat_m)
  
  # Apparent base model fitting
  fit_orig <- coxph(final_formula, data = dat_m, ties = "efron", x = TRUE, model = TRUE)
  lp_orig <- predict(fit_orig, newdata = dat_m, type = "lp")
  
  roc_orig <- timeROC(
    T = dat_m$followup_time, delta = dat_m$event_cox, marker = lp_orig,
    cause = 1, weighting = "marginal", times = times_eval, iid = FALSE
  )
  auc_app <- as.numeric(roc_orig$AUC)
  app_perf_brier <- eval_brier_ipa(fit_orig, dat_m, times_eval)
  
  opt_auc <- matrix(NA_real_, nrow = n_bootstraps, ncol = length(times_eval))
  opt_BS <- matrix(NA_real_, nrow = n_bootstraps, ncol = length(times_eval))
  opt_IPA <- matrix(NA_real_, nrow = n_bootstraps, ncol = length(times_eval))
  
  for (b in 1:n_bootstraps) {
    boot_idx <- sample(seq_len(n_m), size = n_m, replace = TRUE)
    dat_boot <- dat_m[boot_idx, , drop = FALSE]
    
    fit_boot <- try(coxph(final_formula, data = dat_boot, ties = "efron", x = TRUE, model = TRUE), silent = TRUE)
    if (inherits(fit_boot, "try-error") || any(!is.finite(coef(fit_boot)))) next
    
    # 1. AUC Bootstrap Optimism Assessment
    lp_boot <- predict(fit_boot, newdata = dat_boot, type = "lp")
    roc_boot <- try(timeROC(T = dat_boot$followup_time, delta = dat_boot$event_cox, marker = lp_boot, cause = 1, weighting = "marginal", times = times_eval, iid = FALSE), silent = TRUE)
    
    lp_test <- predict(fit_boot, newdata = dat_m, type = "lp")
    roc_test <- try(timeROC(T = dat_m$followup_time, delta = dat_m$event_cox, marker = lp_test, cause = 1, weighting = "marginal", times = times_eval, iid = FALSE), silent = TRUE)
    
    if (!inherits(roc_boot, "try-error") && !inherits(roc_test, "try-error")) {
      if (all(is.finite(roc_boot$AUC)) && all(is.finite(roc_test$AUC))) {
        opt_auc[b, ] <- as.numeric(roc_boot$AUC) - as.numeric(roc_test$AUC)
      }
    }
    
    # 2. Brier and IPA Bootstrap Optimism Assessment
    perf_boot <- try(eval_brier_ipa(fit_boot, dat_boot, times_eval), silent = TRUE)
    perf_test <- try(eval_brier_ipa(fit_boot, dat_m, times_eval), silent = TRUE)
    
    if (!inherits(perf_boot, "try-error") && !inherits(perf_test, "try-error")) {
      opt_BS[b, ] <- perf_test[, "Brier"] - perf_boot[, "Brier"]
      opt_IPA[b, ] <- perf_boot[, "IPA"] - perf_test[, "IPA"]
    }
  }
  
  # Summarize AUC metrics for Imputation m
  for (j in seq_along(times_eval)) {
    ok_auc <- is.finite(opt_auc[, j])
    auc_boot_summary <- rbind(
      auc_boot_summary,
      data.frame(
        imputation = m, time = times_eval[j],
        apparent_AUC = auc_app[j], optimism = mean(opt_auc[ok_auc, j]),
        corrected_AUC = auc_app[j] - mean(opt_auc[ok_auc, j])
      )
    )
    
    ok_brier <- is.finite(opt_BS[, j]) & is.finite(opt_IPA[, j])
    brier_boot_summary <- rbind(
      brier_boot_summary,
      data.frame(
        imputation = m, time = times_eval[j],
        apparent_Brier = app_perf_brier[j, "Brier"],
        corrected_Brier = app_perf_brier[j, "Brier"] + mean(opt_BS[ok_brier, j]),
        apparent_IPA = app_perf_brier[j, "IPA"],
        corrected_IPA = app_perf_brier[j, "IPA"] - mean(opt_IPA[ok_brier, j])
      )
    )
  }
}

# ------------------------------------------------------------------------------
# SECTION 2: IPCW Decision Curve Analysis (DCA) Internal Validation (5-Year Horizon)
# ------------------------------------------------------------------------------

dca_thresholds <- c(0.015, 0.020, 0.025, 0.030, 0.040, 0.050, 0.060)
dca_boot_summary <- data.frame()

for (m in 1:n_imputations) {
  dat_m <- basic_long[basic_long$.imp == m, ]
  dat_m$hear <- factor(as.character(dat_m$hear), levels = c("1", "2", "3", "4-5"), ordered = FALSE)
  N_m <- nrow(dat_m)
  
  fit_orig <- coxph(final_formula, data = dat_m, ties = "efron", x = TRUE, model = TRUE)
  NB_app <- eval_net_benefit(fit_orig, dat_m, t0 = 5, thresholds = dca_thresholds)
  
  # Calculate benchmark Treat-All strategy Net Benefit
  fit_G <- survfit(Surv(followup_time, 1 - event_cox) ~ 1, data = dat_m)
  G_time <- fit_G$time; G_surv <- fit_G$surv
  
  G_at <- function(t) { idx <- findInterval(t, G_time); ifelse(idx == 0, 1, G_surv[pmax(idx, 1)]) }
  G_left <- function(t) { idx <- findInterval(t - 1e-8, G_time); ifelse(idx == 0, 1, G_surv[pmax(idx, 1)]) }
  
  event5 <- dat_m$event_cox == 1 & dat_m$followup_time <= 5
  nonevent5 <- dat_m$followup_time > 5
  
  w_event <- rep(0, N_m); w_nonevent <- rep(0, N_m)
  w_event[event5] <- 1 / G_left(dat_m$followup_time[event5])
  w_nonevent[nonevent5] <- 1 / G_at(5)
  
  NB_all <- sapply(dca_thresholds, function(pt) {
    sum(w_event) / N_m - sum(w_nonevent) / N_m * (pt / (1 - pt))
  })
  
  opt_NB <- matrix(NA_real_, nrow = n_bootstraps, ncol = length(dca_thresholds))
  
  for (b in 1:n_bootstraps) {
    boot_idx <- sample(seq_len(N_m), size = N_m, replace = TRUE)
    dat_boot <- dat_m[boot_idx, , drop = FALSE]
    
    fit_boot <- try(coxph(final_formula, data = dat_boot, ties = "efron", x = TRUE, model = TRUE), silent = TRUE)
    if (inherits(fit_boot, "try-error") || any(!is.finite(coef(fit_boot)))) next
    
    NB_boot <- try(eval_net_benefit(fit_boot, dat_boot, t0 = 5, thresholds = dca_thresholds), silent = TRUE)
    NB_test <- try(eval_net_benefit(fit_boot, dat_m, t0 = 5, thresholds = dca_thresholds), silent = TRUE)
    
    if (!inherits(NB_boot, "try-error") && !inherits(NB_test, "try-error")) {
      opt_NB[b, ] <- NB_boot - NB_test
    }
  }
  
  for (j in seq_along(dca_thresholds)) {
    ok <- is.finite(opt_NB[, j])
    dca_boot_summary <- rbind(
      dca_boot_summary,
      data.frame(
        imputation = m, threshold = dca_thresholds[j],
        apparent_NB = NB_app[j], optimism = mean(opt_NB[ok, j]),
        corrected_NB = NB_app[j] - mean(opt_NB[ok, j]),
        NB_all = NB_all[j],
        corrected_advantage_vs_all = (NB_app[j] - mean(opt_NB[ok, j])) - NB_all[j]
      )
    )
  }
}

# ------------------------------------------------------------------------------
# SECTION 3: Aggregate Pipeline Results across Imputed Sets
# ------------------------------------------------------------------------------

final_auc_summary <- auc_boot_summary %>%
  group_by(time) %>%
  summarise(
    mean_apparent_AUC = round(mean(apparent_AUC), 4),
    mean_corrected_AUC = round(mean(corrected_AUC), 4),
    .groups = "drop"
  )

final_brier_summary <- brier_boot_summary %>%
  group_by(time) %>%
  summarise(
    mean_corrected_Brier = round(mean(corrected_Brier), 5),
    mean_corrected_IPA = round(mean(corrected_IPA), 5),
    .groups = "drop"
  )

final_dca_summary <- dca_boot_summary %>%
  group_by(threshold) %>%
  summarise(
    threshold_pct = threshold * 100,
    mean_corrected_NB = round(mean(corrected_NB), 5),
    mean_NB_all = round(mean(NB_all), 5),
    mean_advantage_vs_all = round(mean(corrected_advantage_vs_all), 5),
    .groups = "drop"
  ) %>%
  select(threshold_pct, mean_corrected_NB, mean_NB_all, mean_advantage_vs_all)

# Display Execution Summaries
print(as.data.frame(final_auc_summary))
print(as.data.frame(final_brier_summary))
print(as.data.frame(final_dca_summary))