# ==============================================================================
# Step 9 Dynamic Predictive Modeling Pipeline: Geographic Internal-External, Cross-Validation (IECV) & Interval-Censoring Sensitivity Analysis
# ==============================================================================

# Load required libraries
library(survival)
library(dplyr)
library(mice)

# Set global seed for reproducible bootstrap operations
set.seed(20260914)

# Define core predictor set and formal Cox formula
final_vars_10 <- c(
  "fall_down", "age", "adlab_c", "arthre", "cesd10",
  "hear", "total_cognition", "wspeed", "puff", "teeth"
)

final_formula_iecv <- Surv(followup_time, event_cox) ~ 
  fall_down + age + adlab_c + arthre + cesd10 + 
  hear + total_cognition + wspeed + puff + teeth

# ------------------------------------------------------------------------------
# SECTION 1: Geographic IECV - Independent Evaluation in Test Province (Sichuan)
# ------------------------------------------------------------------------------

# Define test geographic region (Sichuan: N = 657, Events = 23)
test_province <- "四川省"

iecv_train <- analysis_data %>% filter(province != test_province)
iecv_test  <- analysis_data %>% filter(province == test_province)

# 1. Fit Nelson-Aalen hazard estimation exclusively on TRAIN subset
fit_na_iecv <- coxph(Surv(followup_time, event_cox) ~ 1, data = iecv_train)
bh_iecv     <- basehaz(fit_na_iecv, centered = FALSE)

iecv_train$nelson_aalen_train <- approx(
  x = bh_iecv$time, y = bh_iecv$hazard,
  xout = iecv_train$followup_time, method = "constant", rule = 2, f = 0
)$y
iecv_train$nelson_aalen_train[iecv_train$followup_time < min(bh_iecv$time)] <- 0

# 2. Setup and run TRAIN-only MICE (m = 10, maxit = 10)
mice_train_data <- iecv_train %>% 
  select(all_of(candidate_basic), event_cox, nelson_aalen_train)

meth_train <- rep("", ncol(mice_train_data))
names(meth_train) <- names(mice_train_data)

for (v in candidate_basic) {
  x <- mice_train_data[[v]]
  if (!anyNA(x)) {
    meth_train[v] <- ""
  } else if (is.ordered(x)) {
    meth_train[v] <- "polr"
  } else if (is.factor(x) && nlevels(x) == 2) {
    meth_train[v] <- "logreg"
  } else if (is.factor(x) && nlevels(x) > 2) {
    meth_train[v] <- "polyreg"
  } else if (is.numeric(x) || is.integer(x)) {
    meth_train[v] <- "pmm"
  }
}

pred_train <- make.predictorMatrix(mice_train_data)
diag(pred_train) <- 0
pred_train["event_cox", ] <- 0
pred_train["nelson_aalen_train", ] <- 0
pred_train[candidate_basic, c("event_cox", "nelson_aalen_train")] <- 1

imp_iecv_train <- mice(
  mice_train_data, m = 10, maxit = 10, method = meth_train,
  predictorMatrix = pred_train, seed = 20260914, printFlag = FALSE
)

# 3. Leakage-free TEST Imputation (Excluding outcomes and using ignore logic)
iecv_imp_data <- bind_rows(
  iecv_train %>% select(all_of(candidate_basic)),
  iecv_test  %>% select(all_of(candidate_basic))
)

n_train <- nrow(iecv_train)
n_test  <- nrow(iecv_test)
ignore_test <- c(rep(FALSE, n_train), rep(TRUE, n_test))

meth_test <- rep("", length(candidate_basic))
names(meth_test) <- candidate_basic

for (v in candidate_basic) {
  x_all <- iecv_imp_data[[v]]
  x_tr  <- iecv_train[[v]]
  if (!anyNA(x_all)) {
    meth_test[v] <- ""
  } else if (is.ordered(x_tr)) {
    meth_test[v] <- "polr"
  } else if (is.factor(x_tr) && nlevels(x_tr) == 2) {
    meth_test[v] <- "logreg"
  } else if (is.factor(x_tr) && nlevels(x_tr) > 2) {
    meth_test[v] <- "polyreg"
  } else {
    meth_test[v] <- "pmm"
  }
}

pred_test <- make.predictorMatrix(iecv_imp_data)
diag(pred_test) <- 0

imp_iecv_test <- mice(
  iecv_imp_data, m = 10, maxit = 10, method = meth_test,
  predictorMatrix = pred_test, ignore = ignore_test,
  seed = 20260914, printFlag = FALSE
)

# Extract completed test imputation lists
iecv_test_list <- lapply(1:10, function(m) {
  comp_m <- complete(imp_iecv_test, action = m)
  test_m <- comp_m[(n_train + 1):(n_train + n_test), , drop = FALSE]
  test_m$ID <- iecv_test$ID
  test_m$followup_time <- iecv_test$followup_time
  test_m$event_cox <- iecv_test$event_cox
  test_m
})

# 4. Evaluate discrimination (C-Index) across 10 TRAIN models x 10 TEST imputed sets
iecv_train_long <- complete(imp_iecv_train, action = "long", include = FALSE)
iecv_train_long$followup_time <- iecv_train$followup_time[iecv_train_long$.id]
iecv_train_long$hear <- factor(
  ifelse(as.character(iecv_train_long$hear) %in% c("4", "5"), "4-5", as.character(iecv_train_long$hear)),
  levels = c("1", "2", "3", "4-5"), ordered = FALSE
)

iecv_train_fits <- lapply(1:10, function(m) {
  dat_m <- iecv_train_long[iecv_train_long$.imp == m, ]
  coxph(final_formula_iecv, data = dat_m, ties = "efron")
})

lp_sichuan_mat <- matrix(NA_real_, nrow = n_test, ncol = 100)
col_idx <- 1

for (m_tr in 1:10) {
  for (m_te in 1:10) {
    dat_te <- iecv_test_list[[m_te]]
    dat_te$hear <- factor(
      ifelse(as.character(dat_te$hear) %in% c("4", "5"), "4-5", as.character(dat_te$hear)),
      levels = c("1", "2", "3", "4-5"), ordered = FALSE
    )
    lp_sichuan_mat[, col_idx] <- predict(iecv_train_fits[[m_tr]], newdata = dat_te, type = "lp", reference = "zero")
    col_idx <- col_idx + 1
  }
}

lp_sichuan_mean <- rowMeans(lp_sichuan_mat, na.rm = TRUE)
conc_sichuan <- concordance(
  Surv(followup_time, event_cox) ~ lp_sichuan_mean,
  data = iecv_test, reverse = TRUE
)

C_sichuan  <- conc_sichuan$concordance
SE_sichuan <- sqrt(conc_sichuan$var)

# ------------------------------------------------------------------------------
# SECTION 2: Bootstrap Calibration Evaluation for Independent Geographic Test
# ------------------------------------------------------------------------------

# Generate expected 5-year risk on test set
t_eval <- 5
risk5_mat <- matrix(NA_real_, nrow = n_test, ncol = 100)
col_idx <- 1

for (m_tr in 1:10) {
  fit_m <- iecv_train_fits[[m_tr]]
  bh_m  <- basehaz(fit_m, centered = FALSE)
  H0_5  <- approx(bh_m$time, bh_m$hazard, xout = t_eval, method = "constant", rule = 2, f = 0)$y
  
  for (m_te in 1:10) {
    dat_te <- iecv_test_list[[m_te]]
    dat_te$hear <- factor(
      ifelse(as.character(dat_te$hear) %in% c("4", "5"), "4-5", as.character(dat_te$hear)),
      levels = c("1", "2", "3", "4-5"), ordered = FALSE
    )
    lp <- predict(fit_m, newdata = dat_te, type = "lp", reference = "zero")
    risk5_mat[, col_idx] <- 1 - exp(-H0_5 * exp(lp))
    col_idx <- col_idx + 1
  }
}

risk5_mean <- rowMeans(risk5_mat, na.rm = TRUE)

# Bootstrap calibration metrics (B = 1000)
B_cal <- 1000
boot_cal_results <- matrix(NA_real_, nrow = B_cal, ncol = 3, dimnames = list(NULL, c("obs_risk", "intercept", "slope")))

for (b in 1:B_cal) {
  b_idx  <- sample(seq_len(n_test), size = n_test, replace = TRUE)
  dat_b  <- iecv_test[b_idx, , drop = FALSE]
  risk_b <- pmin(pmax(risk5_mean[b_idx], 1e-6), 1 - 1e-6)
  x_b    <- log(-log(1 - risk_b))
  
  fit_G_b <- try(survfit(Surv(followup_time, 1 - event_cox) ~ 1, data = dat_b), silent = TRUE)
  if (inherits(fit_G_b, "try-error")) next
  
  G_time_b <- fit_G_b$time; G_surv_b <- fit_G_b$surv
  G_at_b   <- function(t) { i <- findInterval(t, G_time_b); ifelse(i == 0, 1, G_surv_b[pmax(i, 1)]) }
  G_left_b <- function(t) { i <- findInterval(t - 1e-8, G_time_b); ifelse(i == 0, 1, G_surv_b[pmax(i, 1)]) }
  
  event_b    <- dat_b$event_cox == 1 & dat_b$followup_time <= t_eval
  nonevent_b <- dat_b$followup_time > t_eval
  use_b      <- event_b | nonevent_b
  
  y_b <- as.integer(event_b[use_b])
  w_b <- numeric(sum(use_b))
  
  ev_use_b <- event_b[use_b]
  w_b[ev_use_b]  <- 1 / G_left_b(dat_b$followup_time[use_b][ev_use_b])
  w_b[!ev_use_b] <- 1 / G_at_b(t_eval)
  
  cal_b <- data.frame(y = y_b, x = x_b[use_b], w = w_b)
  
  fit_int_b   <- try(suppressWarnings(glm(y ~ 1 + offset(x), data = cal_b, family = binomial(link = "cloglog"), weights = w)), silent = TRUE)
  fit_slope_b <- try(suppressWarnings(glm(y ~ x, data = cal_b, family = binomial(link = "cloglog"), weights = w)), silent = TRUE)
  km_b        <- try(survfit(Surv(followup_time, event_cox) ~ 1, data = dat_b), silent = TRUE)
  
  if (inherits(fit_int_b, "try-error") || inherits(fit_slope_b, "try-error") || inherits(km_b, "try-error")) next
  
  obs_b <- 1 - summary(km_b, times = t_eval, extend = TRUE)$surv
  vals  <- c(obs_risk = obs_b, intercept = unname(coef(fit_int_b)[1]), slope = unname(coef(fit_slope_b)["x"]))
  
  if (all(is.finite(vals))) boot_cal_results[b, ] <- vals
}

ok_cal <- complete.cases(boot_cal_results)

# ------------------------------------------------------------------------------
# SECTION 3: Sensitivity Analysis - Interval-Censoring Event Timing (25%, 50%, 75%)
# ------------------------------------------------------------------------------

basic_long_sens <- basic_long
basic_long_sens$interval_left  <- analysis_data$interval_left[basic_long_sens$.id]
basic_long_sens$interval_right <- analysis_data$interval_right[basic_long_sens$.id]

basic_long_sens <- basic_long_sens %>%
  mutate(
    time_q25 = ifelse(event_cox == 1, interval_left + 0.25 * (interval_right - interval_left), followup_time),
    time_q50 = followup_time,
    time_q75 = ifelse(event_cox == 1, interval_left + 0.75 * (interval_right - interval_left), followup_time)
  )

# Display IECV & Geographic Results
cat("=== Independent Geographic Validation (Sichuan) ===\n")
cat("C-Index:", round(C_sichuan, 4), "95% CI:", round(max(0, C_sichuan - 1.96 * SE_sichuan), 4), "-", round(min(1, C_sichuan + 1.96 * SE_sichuan), 4), "\n")
cat("Mean Predicted 5-Year Risk:", round(mean(risk5_mean), 4), "\n")
cat("Observed 5-Year Risk 95% CI:", round(quantile(boot_cal_results[ok_cal, "obs_risk"], c(0.025, 0.975)), 4), "\n")
cat("Calibration Intercept 95% CI:", round(quantile(boot_cal_results[ok_cal, "intercept"], c(0.025, 0.975)), 4), "\n")
cat("Calibration Slope 95% CI:", round(quantile(boot_cal_results[ok_cal, "slope"], c(0.025, 0.975)), 4), "\n")