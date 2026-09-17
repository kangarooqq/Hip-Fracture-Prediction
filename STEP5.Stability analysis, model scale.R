# ==============================================================================
# Step 5 Feature selection via LASSO, stability analysis, model scale
# leave-one-variable-out comparison, final 10-predictor model fitting, and 2-year landmark sensitivity analysis.
# ==============================================================================

library(glmnet)
library(survival)
library(mice)
library(dplyr)

# ============================================================
# Step 65: Extract nonzero coefficients from LASSO
# ============================================================
idx_alpha1 <- which(alpha_grid == 1)
fit_selected <- cv_alpha[[idx_alpha1]]
selected_lambda <- fit_selected$lambda.min

cat("Selected alpha =", selected_alpha, "\n")
cat("Selected lambda =", selected_lambda, "\n")

coef_selected <- as.matrix(coef(fit_selected, s = "lambda.min"))
selected_coef <- data.frame(
  term = rownames(coef_selected),
  coefficient = as.numeric(coef_selected[, 1])
)
selected_coef <- selected_coef[selected_coef$coefficient != 0, ]
selected_coef <- selected_coef[order(-abs(selected_coef$coefficient)), ]
rownames(selected_coef) <- NULL

cat("Number of nonzero coefficients =", nrow(selected_coef), "\n\n")
print(selected_coef, row.names = FALSE)

# ============================================================
# Step 66: LASSO selection stability across 10 imputations
# ============================================================
selection_list <- vector("list", 10)

for (m in 1:10) {
  cat("Running imputation", m, "...\n")
  idx <- basic_long$.imp == m
  X_m <- X_basic[idx, , drop = FALSE]
  y_m <- Surv(basic_long$followup_time[idx], basic_long$event_cox[idx])
  fold_m <- foldid_basic[idx]
  
  set.seed(20260913)
  cv_m <- cv.glmnet(
    x = X_m, y = y_m, family = "cox", alpha = 1,
    foldid = fold_m, type.measure = "C",
    standardize = TRUE, grouped = TRUE, nlambda = 50
  )
  coef_m <- as.matrix(coef(cv_m, s = "lambda.min"))
  selection_list[[m]] <- rownames(coef_m)[coef_m[, 1] != 0]
}

all_terms <- colnames(X_basic)
selection_frequency <- data.frame(
  term = all_terms,
  selected_n = sapply(all_terms, function(term) {
    sum(sapply(selection_list, function(x) term %in% x))
  })
)
selection_frequency$selected_pct <- selection_frequency$selected_n / 10 * 100
selection_frequency <- selection_frequency[order(-selection_frequency$selected_n, selection_frequency$term), ]
selection_frequency_nonzero <- selection_frequency[selection_frequency$selected_n > 0, ]

print(selection_frequency_nonzero, row.names = FALSE)
cat("\nTerms selected in >=8/10 imputations =", sum(selection_frequency$selected_n >= 8), "\n")
cat("Terms selected in 10/10 imputations =", sum(selection_frequency$selected_n == 10), "\n")

# ============================================================
# Step 67: Map dummy terms to raw predictors & calculate variable-level selection frequency
# ============================================================
mm_full <- model.matrix(~ ., data = design_data)
assign_vec <- attr(mm_full, "assign")
term_labels <- attr(terms(~ ., data = design_data), "term.labels")

term_map <- data.frame(
  term = colnames(mm_full)[-1],
  variable = term_labels[assign_vec[-1]],
  stringsAsFactors = FALSE
)

variable_selection <- lapply(1:10, function(m) {
  selected_terms <- selection_list[[m]]
  unique(term_map$variable[term_map$term %in% selected_terms])
})

variable_frequency <- data.frame(
  variable = candidate_basic,
  selected_n = sapply(candidate_basic, function(v) {
    sum(sapply(variable_selection, function(x) v %in% x))
  })
)
variable_frequency$selected_pct <- variable_frequency$selected_n * 10
variable_frequency <- variable_frequency[variable_frequency$selected_n > 0, ]
variable_frequency <- variable_frequency[order(-variable_frequency$selected_n, variable_frequency$variable), ]

print(variable_frequency, row.names = FALSE)

# ============================================================
# Step 68: Define core and borderline variables
# ============================================================
core_vars <- variable_frequency$variable[variable_frequency$selected_n >= 8]
borderline_vars <- variable_frequency$variable[variable_frequency$selected_n >= 5 & variable_frequency$selected_n < 8]
unstable_vars <- variable_frequency$variable[variable_frequency$selected_n < 5]

cat("Core variables (>=8/10):\n")
print(core_vars)
cat("\nNumber of core variables =", length(core_vars), "\n")
cat("\nBorderline variables (5-7/10):\n")
print(borderline_vars)

# ============================================================
# Step 69: Verify definition and calculation of wspeed
# ============================================================
for (v in c("wspeed", "wspeed1", "wspeed2")) {
  cat("\n====================\n", v, "\n")
  cat("Label:\n"); print(attr(analysis_data[[v]], "label"))
  cat("Summary:\n"); print(summary(analysis_data[[v]]))
}

wspeed_calc <- rowMeans(analysis_data[, c("wspeed1", "wspeed2")], na.rm = TRUE)
wspeed_calc[is.na(analysis_data$wspeed1) & is.na(analysis_data$wspeed2)] <- NA
cat("\nCorrelation with mean(wspeed1, wspeed2) =", cor(analysis_data$wspeed, wspeed_calc, use = "complete.obs"), "\n")

# ============================================================
# Step 70 & 71: Fix variable addition order for model size curve
# ============================================================
plateau_order_final <- c(
  "fall_down",        # Mandatory inclusion
  "age",              # Age
  "adlab_c",          # ADL
  "arthre",           # Arthritis
  "srh",              # Self-rated health
  "cesd10",           # Depressive symptoms
  "hear",             # Hearing
  "total_cognition",  # Cognition
  "wspeed",           # 2.5m walking speed (time)
  "puff",             # Peak flow test
  "teeth",            # Loss of teeth
  "disability",       # Disability status
  "pain_burden",      # Pain burden
  "hchild",           # Alive children count
  "hospital_use"      # Healthcare utilization
)

# ============================================================
# Step 74 & 75: Fix separation by collapsing hear levels (4 & 5 -> 4-5)
# ============================================================
basic_long$hear <- factor(
  ifelse(as.character(basic_long$hear) %in% c("4", "5"), "4-5", as.character(basic_long$hear)),
  levels = c("1", "2", "3", "4-5")
)

# ============================================================
# Step 76: Re-generate design matrix after collapsing hear
# ============================================================
design_data <- basic_long %>%
  select(all_of(candidate_basic)) %>%
  mutate(across(where(is.ordered), ~ factor(as.character(.x), levels = levels(.x), ordered = FALSE)))

X_basic <- model.matrix(~ ., data = design_data)[, -1, drop = FALSE]
mm_full <- model.matrix(~ ., data = design_data)
assign_vec <- attr(mm_full, "assign")
term_labels <- attr(terms(~ ., data = design_data), "term.labels")

term_map <- data.frame(
  term = colnames(mm_full)[-1],
  variable = term_labels[assign_vec[-1]],
  stringsAsFactors = FALSE
)

# ============================================================
# Step 77: Model size vs. CV C-index curve (5 -> 15 variables)
# ============================================================
model_sizes <- 5:15
plateau_results2 <- data.frame(n_vars = model_sizes, mean_C = NA_real_, SE_C = NA_real_)
fold_C_list2 <- vector("list", length(model_sizes))

for (s in seq_along(model_sizes)) {
  k <- model_sizes[s]
  vars_k <- plateau_order_final[1:k]
  cols_k <- term_map$term[term_map$variable %in% vars_k]
  fold_C <- numeric(10)
  
  for (f in 1:10) {
    test_ids <- which(fold_subject == f)
    pred_MI <- matrix(NA_real_, nrow = length(test_ids), ncol = 10)
    
    for (m in 1:10) {
      rows_m <- which(basic_long$.imp == m)
      ids_m <- basic_long$.id[rows_m]
      train_rows <- rows_m[fold_subject[ids_m] != f]
      test_rows <- rows_m[match(test_ids, ids_m)]
      
      train_df <- data.frame(
        time = basic_long$followup_time[train_rows],
        event = basic_long$event_cox[train_rows],
        X_basic[train_rows, cols_k, drop = FALSE],
        check.names = TRUE
      )
      test_df <- data.frame(X_basic[test_rows, cols_k, drop = FALSE], check.names = TRUE)
      
      fit_k <- coxph(Surv(time, event) ~ ., data = train_df, ties = "efron")
      pred_MI[, m] <- predict(fit_k, newdata = test_df, type = "lp")
    }
    
    lp_mean <- rowMeans(pred_MI)
    fold_C[f] <- concordance(
      Surv(analysis_data$followup_time[test_ids], analysis_data$event_cox[test_ids]) ~ lp_mean,
      reverse = TRUE
    )$concordance
  }
  
  fold_C_list2[[s]] <- fold_C
  plateau_results2$mean_C[s] <- mean(fold_C)
  plateau_results2$SE_C[s] <- sd(fold_C) / sqrt(10)
}

print(plateau_results2, row.names = FALSE)

# ============================================================
# Step 80: Leave-one-variable-out comparison (11 candidates to 10)
# ============================================================
candidate_11 <- c("fall_down", "age", "adlab_c", "arthre", "srh", "cesd10", "hear", "total_cognition", "wspeed", "puff", "teeth")
drop_candidates <- setdiff(candidate_11, "fall_down")

loo_results <- data.frame(dropped_variable = drop_candidates, mean_C = NA_real_, SE_C = NA_real_)
loo_fold_C <- vector("list", length(drop_candidates))

for (j in seq_along(drop_candidates)) {
  drop_v <- drop_candidates[j]
  vars_10 <- setdiff(candidate_11, drop_v)
  cols_10 <- term_map$term[term_map$variable %in% vars_10]
  fold_C <- numeric(10)
  
  for (f in 1:10) {
    test_ids <- which(fold_subject == f)
    pred_MI <- matrix(NA_real_, nrow = length(test_ids), ncol = 10)
    
    for (m in 1:10) {
      rows_m <- which(basic_long$.imp == m)
      ids_m <- basic_long$.id[rows_m]
      train_rows <- rows_m[fold_subject[ids_m] != f]
      test_rows <- rows_m[match(test_ids, ids_m)]
      
      train_df <- data.frame(
        time = basic_long$followup_time[train_rows],
        event = basic_long$event_cox[train_rows],
        X_basic[train_rows, cols_10, drop = FALSE],
        check.names = TRUE
      )
      test_df <- data.frame(X_basic[test_rows, cols_10, drop = FALSE], check.names = TRUE)
      
      fit_10 <- coxph(Surv(time, event) ~ ., data = train_df, ties = "efron")
      pred_MI[, m] <- predict(fit_10, newdata = test_df, type = "lp")
    }
    
    lp_mean <- rowMeans(pred_MI)
    fold_C[f] <- concordance(
      Surv(analysis_data$followup_time[test_ids], analysis_data$event_cox[test_ids]) ~ lp_mean,
      reverse = TRUE
    )$concordance
  }
  
  loo_fold_C[[j]] <- fold_C
  loo_results$mean_C[j] <- mean(fold_C)
  loo_results$SE_C[j] <- sd(fold_C) / sqrt(10)
}

loo_results <- loo_results[order(-loo_results$mean_C), ]
print(loo_results, row.names = FALSE)

# ============================================================
# Step 82: 5 Repeated 10-fold CV for Top-3 10-variable candidates
# ============================================================
top_drop <- c("srh", "total_cognition", "wspeed")
repeat_seeds <- c(20260913, 20260914, 20260915, 20260916, 20260917)
repeat_results <- list()
counter <- 1

for (seed_i in repeat_seeds) {
  set.seed(seed_i)
  fold_rep <- integer(nrow(analysis_data))
  event_idx <- which(analysis_data$event_cox == 1)
  nonevent_idx <- which(analysis_data$event_cox == 0)
  
  fold_rep[event_idx] <- sample(rep(1:10, length.out = length(event_idx)))
  fold_rep[nonevent_idx] <- sample(rep(1:10, length.out = length(nonevent_idx)))
  
  for (drop_v in top_drop) {
    vars_10 <- setdiff(candidate_11, drop_v)
    cols_10 <- term_map$term[term_map$variable %in% vars_10]
    fold_C <- numeric(10)
    
    for (f in 1:10) {
      test_ids <- which(fold_rep == f)
      pred_MI <- matrix(NA_real_, nrow = length(test_ids), ncol = 10)
      
      for (m in 1:10) {
        rows_m <- which(basic_long$.imp == m)
        ids_m <- basic_long$.id[rows_m]
        train_rows <- rows_m[fold_rep[ids_m] != f]
        test_rows <- rows_m[match(test_ids, ids_m)]
        
        train_df <- data.frame(
          time = basic_long$followup_time[train_rows],
          event = basic_long$event_cox[train_rows],
          X_basic[train_rows, cols_10, drop = FALSE],
          check.names = TRUE
        )
        test_df <- data.frame(X_basic[test_rows, cols_10, drop = FALSE], check.names = TRUE)
        
        fit_10 <- coxph(Surv(time, event) ~ ., data = train_df, ties = "efron")
        pred_MI[, m] <- predict(fit_10, newdata = test_df, type = "lp")
      }
      
      lp_mean <- rowMeans(pred_MI)
      fold_C[f] <- concordance(
        Surv(analysis_data$followup_time[test_ids], analysis_data$event_cox[test_ids]) ~ lp_mean,
        reverse = TRUE
      )$concordance
    }
    
    repeat_results[[counter]] <- data.frame(seed = seed_i, dropped_variable = drop_v, mean_C = mean(fold_C))
    counter <- counter + 1
  }
}

repeat_summary <- aggregate(mean_C ~ dropped_variable, data = do.call(rbind, repeat_results), FUN = function(x) c(mean = mean(x), sd = sd(x)))
print(repeat_summary)

# ============================================================
# Step 83-85: Lock Final 10-Variable Model & Fit Cox Models across 10 Imputations
# ============================================================
final_vars_10 <- c("fall_down", "age", "adlab_c", "arthre", "cesd10", "hear", "total_cognition", "wspeed", "puff", "teeth")

final_formula <- Surv(followup_time, event_cox) ~ fall_down + age + adlab_c + arthre + cesd10 + hear + total_cognition + wspeed + puff + teeth

final_cox_fits <- vector("list", 10)
for (m in 1:10) {
  dat_m <- basic_long[basic_long$.imp == m, ]
  final_cox_fits[[m]] <- coxph(final_formula, data = dat_m, ties = "efron", x = TRUE)
}

final_mira <- as.mira(final_cox_fits)
final_pool <- pool(final_mira)
final_pool_summary <- summary(final_pool, conf.int = TRUE, exponentiate = TRUE)

cat("\n--- Final Main Model Pooled Summary ---\n")
print(final_pool_summary, row.names = FALSE)

# ============================================================
# Step 86: Check Proportional Hazards Assumption
# ============================================================
ph_results <- list()
for (m in 1:10) {
  zph_m <- cox.zph(final_cox_fits[[m]], transform = "km", terms = TRUE, singledf = TRUE)
  tmp <- as.data.frame(zph_m$table)
  tmp$term <- rownames(tmp)
  tmp$imputation <- m
  ph_results[[m]] <- tmp
}

ph_terms <- do.call(rbind, ph_results)[do.call(rbind, ph_results)$term != "GLOBAL", ]
ph_summary <- aggregate(p ~ term, data = ph_terms, FUN = function(x) c(median_p = median(x), min_p = min(x)))
print(ph_summary)

# ============================================================
# Step 88B-90: 2-Year Landmark Cohort Construction & Setup
# ============================================================
landmark_check2 <- interval_outcome %>%
  mutate(
    landmark_group = case_when(
      event == TRUE & right_time <= 2 ~ "early_event_0_2y",
      event == TRUE & left_time < 2 & right_time > 2 ~ "ambiguous_crossing_2y",
      event == TRUE & left_time >= 2 ~ "eligible_later_event",
      event == FALSE & censor_time > 2 ~ "eligible_nonevent",
      event == FALSE & censor_time == 2 ~ "censored_at_landmark",
      event == FALSE & censor_time < 2 ~ "censored_before_2y",
      TRUE ~ "other"
    )
  )

landmark_outcome <- landmark_check2 %>%
  filter(landmark_group %in% c("eligible_later_event", "eligible_nonevent")) %>%
  transmute(
    ID,
    event_landmark = as.integer(landmark_group == "eligible_later_event"),
    followup_landmark = ifelse(event_landmark == 1, (left_time + right_time) / 2 - 2, censor_time - 2)
  )

landmark_data <- analysis_data %>% inner_join(landmark_outcome, by = "ID")
landmark_data$hear <- factor(
  ifelse(as.character(landmark_data$hear) %in% c("4", "5"), "4-5", as.character(landmark_data$hear)),
  levels = c("1", "2", "3", "4-5"), ordered = TRUE
)

fit_na_landmark <- coxph(Surv(followup_landmark, event_landmark) ~ 1, data = landmark_data)
bh_landmark <- basehaz(fit_na_landmark, centered = FALSE)

landmark_data$nelson_aalen_landmark <- approx(
  x = bh_landmark$time, y = bh_landmark$hazard,
  xout = landmark_data$followup_landmark, method = "constant", rule = 2, f = 0
)$y

mice_landmark_data <- landmark_data %>% select(all_of(candidate_basic), event_landmark, nelson_aalen_landmark)

meth_landmark <- rep("", ncol(mice_landmark_data))
names(meth_landmark) <- names(mice_landmark_data)

for (v in candidate_basic) {
  x <- mice_landmark_data[[v]]
  if (!anyNA(x)) { meth_landmark[v] <- "" }
  else if (is.ordered(x)) { meth_landmark[v] <- "polr" }
  else if (is.factor(x) && nlevels(x) == 2) { meth_landmark[v] <- "logreg" }
  else if (is.factor(x) && nlevels(x) > 2) { meth_landmark[v] <- "polyreg" }
  else if (is.numeric(x) || is.integer(x)) { meth_landmark[v] <- "pmm" }
}

pred_landmark2 <- make.predictorMatrix(mice_landmark_data)
diag(pred_landmark2) <- 0
pred_landmark2["event_landmark", ] <- 0
pred_landmark2["nelson_aalen_landmark", ] <- 0
pred_landmark2[candidate_basic, c("event_landmark", "nelson_aalen_landmark")] <- 1
pred_landmark2["total_cognition", "nelson_aalen_landmark"] <- 0

# ============================================================
# Step 93-95: Landmark MICE Imputation (10 iterations)
# ============================================================
set.seed(20260913)
imp_landmark2 <- mice(
  mice_landmark_data, m = 10, maxit = 5, method = meth_landmark,
  predictorMatrix = pred_landmark2, seed = 20260913, printFlag = FALSE
)

set.seed(20260913)
imp_landmark_final <- mice::mice.mids(imp_landmark2, maxit = 5, printFlag = FALSE)

# ============================================================
# Step 96-97: Landmark Cox Fitting & Comparison with Main Analysis
# ============================================================
landmark_long <- complete(imp_landmark_final, action = "long", include = FALSE)
landmark_long$followup_landmark <- landmark_data$followup_landmark[landmark_long$.id]
landmark_long$hear <- factor(as.character(landmark_long$hear), levels = c("1", "2", "3", "4-5"), ordered = FALSE)

landmark_formula <- Surv(followup_landmark, event_landmark) ~ fall_down + age + adlab_c + arthre + cesd10 + hear + total_cognition + wspeed + puff + teeth

landmark_cox_fits <- vector("list", 10)
for (m in 1:10) {
  dat_m <- landmark_long[landmark_long$.imp == m, ]
  landmark_cox_fits[[m]] <- coxph(landmark_formula, data = dat_m, ties = "efron", x = TRUE)
}

landmark_mira <- as.mira(landmark_cox_fits)
landmark_pool <- pool(landmark_mira)
landmark_pool_summary <- summary(landmark_pool, conf.int = TRUE, exponentiate = TRUE)

main_compare <- final_pool_summary %>%
  select(term, main_HR = estimate, main_low = conf.low, main_high = conf.high, main_p = p.value)

landmark_compare <- landmark_pool_summary %>%
  select(term, landmark_HR = estimate, landmark_low = conf.low, landmark_high = conf.high, landmark_p = p.value)

comparison_landmark <- merge(main_compare, landmark_compare, by = "term", all = TRUE)
comparison_landmark$HR_change_pct <- (comparison_landmark$landmark_HR / comparison_landmark$main_HR - 1) * 100
comparison_landmark <- comparison_landmark %>% mutate(across(where(is.numeric), ~ round(.x, 4)))

cat("\n--- Main vs. 2-Year Landmark Sensitivity Analysis Comparison ---\n")
print(comparison_landmark, row.names = FALSE)