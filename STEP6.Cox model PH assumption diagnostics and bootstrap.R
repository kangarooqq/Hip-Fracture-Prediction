# ==============================================================================
# Step 6 Landmark sensitivity Cox model PH assumption diagnostics, 
#       time-varying interaction sensitivity, apparent C-index baseline, 
#       fixed-model bootstrap optimism-correction, and leakage-free
#       outer 10-fold cross-validation setup for full-pipeline validation.
# ==============================================================================

library(survival)
library(mice)
library(dplyr)

# ============================================================
# Step 98: 2-Year Landmark Cox model PH assumption test
# ============================================================
ph_landmark <- lapply(
  landmark_cox_fits,
  function(fit) {
    cox.zph(fit, transform = "km", terms = TRUE, singledf = TRUE)
  }
)

# Collect p-values across imputations
ph_p_landmark <- do.call(
  cbind,
  lapply(ph_landmark, function(x) x$table[, "p"])
)
colnames(ph_p_landmark) <- paste0("imp", 1:10)

# Summary
ph_summary_landmark <- data.frame(
  term = rownames(ph_p_landmark),
  median_p = apply(ph_p_landmark, 1, median, na.rm = TRUE),
  min_p = apply(ph_p_landmark, 1, min, na.rm = TRUE),
  n_p_lt_0.05 = apply(ph_p_landmark, 1, function(x) sum(x < 0.05, na.rm = TRUE))
)

ph_summary_landmark$median_p <- round(ph_summary_landmark$median_p, 4)
ph_summary_landmark$min_p <- round(ph_summary_landmark$min_p, 4)

cat("\n--- Landmark Model PH Assumption Summary ---\n")
print(ph_summary_landmark, row.names = FALSE)

cat("\nGLOBAL p-values across imputations:\n")
print(round(ph_p_landmark["GLOBAL", ], 4))

# ============================================================
# Step 99: total_cognition PH signal consistency
# ============================================================
ph_landmark_coef <- lapply(
  landmark_cox_fits,
  function(fit) {
    cox.zph(fit, transform = "km", terms = FALSE, global = TRUE)
  }
)

cat("\n--- Total Cognition Results Across 10 Imputations ---\n")
tc_ph <- do.call(
  rbind,
  lapply(1:10, function(i) {
    tab <- ph_landmark_coef[[i]]$table
    data.frame(
      imputation = i,
      chisq = tab["total_cognition", "chisq"],
      p = tab["total_cognition", "p"]
    )
  })
)

tc_ph$chisq <- round(tc_ph$chisq, 4)
tc_ph$p <- round(tc_ph$p, 4)
print(tc_ph, row.names = FALSE)

cat("\nP < 0.05 =", sum(tc_ph$p < 0.05), "/ 10\n")
cat("Median P =", round(median(tc_ph$p), 4), "\n")

# ============================================================
# Step 100: Check total_cognition Schoenfeld residual trend direction
# ============================================================
tc_trend <- data.frame(
  imputation = 1:10,
  correlation = sapply(ph_landmark_coef, function(z) {
    cor(z$x, z$y[, "total_cognition"], use = "complete.obs")
  }),
  p = tc_ph$p
)
tc_trend$correlation <- round(tc_trend$correlation, 4)

cat("\n--- Total Cognition Residual Trend Direction ---\n")
print(tc_trend, row.names = FALSE)

cat("\nPositive trend =", sum(tc_trend$correlation > 0), "/ 10\n")
cat("Negative trend =", sum(tc_trend$correlation < 0), "/ 10\n")
cat("Median correlation =", round(median(tc_trend$correlation), 4), "\n")

cat("\nSignificant imputations only (P < 0.05):\n")
print(tc_trend[tc_trend$p < 0.05, ], row.names = FALSE)

# ============================================================
# Step 101: Total cognition time-varying effect sensitivity model
# ============================================================
landmark_tv_formula <- Surv(followup_landmark, event_landmark) ~
  fall_down + age + adlab_c + arthre + cesd10 + hear + 
  total_cognition + tt(total_cognition) + wspeed + puff + teeth

landmark_tv_fits <- vector("list", 10)

for (m in 1:10) {
  dat_m <- landmark_long[landmark_long$.imp == m, ]
  landmark_tv_fits[[m]] <- coxph(
    landmark_tv_formula,
    data = dat_m,
    ties = "efron",
    tt = function(x, t, ...) { x * log(t) }
  )
}

bad_tv <- sapply(landmark_tv_fits, function(fit) any(!is.finite(coef(fit))))
cat("\nModels with non-finite coefficients =", sum(bad_tv), "\n")

cat("\nTime-varying cognition coefficients across 10 imputations:\n")
print(round(sapply(landmark_tv_fits, function(fit) coef(fit)["tt(total_cognition)"]), 4))

# ============================================================
# Step 102: Rubin's rules pooling for cognition time-varying effect
# ============================================================
landmark_tv_mira <- as.mira(landmark_tv_fits)
landmark_tv_pool <- pool(landmark_tv_mira)
landmark_tv_summary <- summary(landmark_tv_pool, conf.int = TRUE, exponentiate = FALSE)

tc_tv_result <- landmark_tv_summary[
  landmark_tv_summary$term %in% c("total_cognition", "tt(total_cognition)"),
  c("term", "estimate", "std.error", "df", "p.value", "conf.low", "conf.high")
]

tc_tv_result[sapply(tc_tv_result, is.numeric)] <- lapply(
  tc_tv_result[sapply(tc_tv_result, is.numeric)],
  round,
  4
)

cat("\n--- Pooled Time-Varying Cognition Interaction Results ---\n")
print(tc_tv_result, row.names = FALSE)

# ============================================================
# Step 103: Main model apparent C-index across 10 imputations
# ============================================================
main_cindex <- sapply(final_cox_fits, function(fit) {
  concordance(fit)$concordance
})

cat("\n--- Main Model Apparent C-index ---\n")
cat("C-index in each imputation:\n")
print(round(main_cindex, 4))
cat("\nMean apparent C-index =", round(mean(main_cindex), 4), "\n")
cat("SD =", round(sd(main_cindex), 4), "\n")
cat("Range =", paste(round(range(main_cindex), 4), collapse = " - "), "\n")

# ============================================================
# Step 104: Bootstrap optimism correction smoke test (m = 1, B = 20)
# ============================================================
dat1 <- basic_long[basic_long$.imp == 1, ]
dat1$hear <- factor(as.character(dat1$hear), levels = c("1", "2", "3", "4-5"))

final_formula <- Surv(followup_time, event_cox) ~
  fall_down + age + adlab_c + arthre + cesd10 + hear + 
  total_cognition + wspeed + puff + teeth

fit_original <- coxph(final_formula, data = dat1, ties = "efron", x = TRUE)
C_apparent <- concordance(fit_original)$concordance

set.seed(20260914)
B_test <- 20
boot_result <- matrix(
  NA_real_,
  nrow = B_test,
  ncol = 3,
  dimnames = list(NULL, c("C_boot", "C_test", "optimism"))
)
n <- nrow(dat1)

for (b in 1:B_test) {
  idx <- sample(seq_len(n), size = n, replace = TRUE)
  dat_boot <- dat1[idx, ]
  
  fit_boot <- coxph(final_formula, data = dat_boot, ties = "efron")
  C_boot <- concordance(fit_boot)$concordance
  
  lp_test <- predict(fit_boot, newdata = dat1, type = "lp")
  C_test <- concordance(
    Surv(followup_time, event_cox) ~ lp_test,
    data = dat1,
    reverse = TRUE
  )$concordance
  
  boot_result[b, ] <- c(C_boot, C_test, C_boot - C_test)
}

cat("\n--- Bootstrap Smoke Test Summary (Imputation 1) ---\n")
cat("Apparent C-index =", round(C_apparent, 4), "\n")
cat("Successful bootstraps =", sum(complete.cases(boot_result)), "/", B_test, "\n")
cat("Mean optimism =", round(mean(boot_result[, "optimism"], na.rm = TRUE), 4), "\n")
cat("Optimism-corrected C-index =", round(C_apparent - mean(boot_result[, "optimism"], na.rm = TRUE), 4), "\n")

# ============================================================
# Step 105: Fixed-model Bootstrap (10 imputations x 100 bootstraps)
# ============================================================
set.seed(20260914)
B <- 100

boot_summary_MI <- data.frame(
  imputation = 1:10,
  apparent_C = NA_real_,
  mean_optimism = NA_real_,
  corrected_C = NA_real_,
  successful_B = NA_integer_
)

for (m in 1:10) {
  dat_m <- basic_long[basic_long$.imp == m, ]
  dat_m$hear <- factor(as.character(dat_m$hear), levels = c("1", "2", "3", "4-5"))
  
  fit_orig <- coxph(final_formula, data = dat_m, ties = "efron")
  C_app <- concordance(fit_orig)$concordance
  
  optimism_b <- rep(NA_real_, B)
  n_m <- nrow(dat_m)
  
  for (b in 1:B) {
    idx <- sample(seq_len(n_m), size = n_m, replace = TRUE)
    dat_boot <- dat_m[idx, ]
    
    fit_boot <- try(coxph(final_formula, data = dat_boot, ties = "efron"), silent = TRUE)
    
    if (inherits(fit_boot, "try-error") || any(!is.finite(coef(fit_boot)))) next
    
    C_boot <- concordance(fit_boot)$concordance
    lp_test <- predict(fit_boot, newdata = dat_m, type = "lp")
    C_test <- concordance(
      Surv(followup_time, event_cox) ~ lp_test,
      data = dat_m,
      reverse = TRUE
    )$concordance
    
    optimism_b[b] <- C_boot - C_test
  }
  
  mean_opt <- mean(optimism_b, na.rm = TRUE)
  boot_summary_MI[m, ] <- c(m, C_app, mean_opt, C_app - mean_opt, sum(is.finite(optimism_b)))
}

boot_summary_MI$apparent_C <- round(boot_summary_MI$apparent_C, 4)
boot_summary_MI$mean_optimism <- round(boot_summary_MI$mean_optimism, 4)
boot_summary_MI$corrected_C <- round(boot_summary_MI$corrected_C, 4)

cat("\n--- Fixed-Model Bootstrap Summary (10 Imputations x 100 Bootstraps) ---\n")
print(boot_summary_MI, row.names = FALSE)
cat("\nMean corrected C-index =", round(mean(boot_summary_MI$corrected_C), 4), "\n")
cat("Range corrected C-index =", paste(round(range(boot_summary_MI$corrected_C), 4), collapse = " - "), "\n")

# ============================================================
# Step 106: Full-pipeline validation outer 10-fold setup
# ============================================================
set.seed(20260914)

outer_base <- analysis_data %>%
  select(ID, event_cox, followup_time) %>%
  group_by(event_cox) %>%
  mutate(
    random_order = sample(seq_len(n())),
    outer_fold = rep(1:10, length.out = n())[random_order]
  ) %>%
  ungroup() %>%
  select(-random_order)

cat("\n--- Outer Fold Assignment Check ---\n")
cat("Subjects per outer fold:\n")
print(table(outer_base$outer_fold))
cat("\nEvents per outer fold:\n")
print(table(fold = outer_base$outer_fold, event = outer_base$event_cox))
cat("\nIDs assigned to >1 fold =", sum(table(outer_base$ID) > 1), "\n")

# ============================================================
# Step 107: Outer fold 1 data split and isolation check
# ============================================================
outer_fold_now <- 1

train_ids_1 <- outer_base %>% filter(outer_fold != outer_fold_now) %>% pull(ID)
test_ids_1  <- outer_base %>% filter(outer_fold == outer_fold_now) %>% pull(ID)

outer_train1 <- analysis_data %>% filter(ID %in% train_ids_1)
outer_test1  <- analysis_data %>% filter(ID %in% test_ids_1)

cat("\n--- Outer Fold 1 Split Diagnostics ---\n")
cat("TRAIN N =", nrow(outer_train1), "| Events =", sum(outer_train1$event_cox), "\n")
cat("TEST N  =", nrow(outer_test1),  "| Events =", sum(outer_test1$event_cox), "\n")
cat("ID overlap =", length(intersect(outer_train1$ID, outer_test1$ID)), "\n")

train_missing <- sapply(outer_train1[, candidate_basic], function(x) mean(is.na(x)))
test_missing  <- sapply(outer_test1[, candidate_basic], function(x) mean(is.na(x)))
cat("TRAIN max predictor missing =", round(max(train_missing), 4), "\n")
cat("TEST max predictor missing  =", round(max(test_missing), 4), "\n")

# ============================================================
# Step 108: Outer fold 1 - TRAIN-only Nelson-Aalen & MICE data
# ============================================================
fit_na_train1 <- coxph(Surv(followup_time, event_cox) ~ 1, data = outer_train1)
bh_train1 <- basehaz(fit_na_train1, centered = FALSE)

outer_train1$nelson_aalen_train <- approx(
  x = bh_train1$time,
  y = bh_train1$hazard,
  xout = outer_train1$followup_time,
  method = "constant",
  rule = 2,
  f = 0
)$y

outer_train1$nelson_aalen_train[outer_train1$followup_time < min(bh_train1$time)] <- 0

mice_train1 <- outer_train1 %>%
  select(all_of(candidate_basic), event_cox, nelson_aalen_train)

# ============================================================
# Step 109: Outer fold 1 - TRAIN-only MICE setup
# ============================================================
meth_train1 <- rep("", ncol(mice_train1))
names(meth_train1) <- names(mice_train1)

for (v in candidate_basic) {
  x <- mice_train1[[v]]
  if (!anyNA(x)) { meth_train1[v] <- "" }
  else if (is.ordered(x)) { meth_train1[v] <- "polr" }
  else if (is.factor(x) && nlevels(x) == 2) { meth_train1[v] <- "logreg" }
  else if (is.factor(x) && nlevels(x) > 2) { meth_train1[v] <- "polyreg" }
  else if (is.numeric(x) || is.integer(x)) { meth_train1[v] <- "pmm" }
}

meth_train1["event_cox"] <- ""
meth_train1["nelson_aalen_train"] <- ""

pred_train1 <- make.predictorMatrix(mice_train1)
diag(pred_train1) <- 0
pred_train1["event_cox", ] <- 0
pred_train1["nelson_aalen_train", ] <- 0
pred_train1[candidate_basic, c("event_cox", "nelson_aalen_train")] <- 1

# ============================================================
# Step 110: Outer fold 1 - Leakage-free TEST imputation smoke test
# ============================================================
cv_imp_data1 <- bind_rows(
  outer_train1 %>% select(all_of(candidate_basic)),
  outer_test1  %>% select(all_of(candidate_basic))
)

n_train1 <- nrow(outer_train1)
n_test1  <- nrow(outer_test1)

ignore_test1 <- c(rep(FALSE, n_train1), rep(TRUE, n_test1))

meth_cv1 <- rep("", length(candidate_basic))
names(meth_cv1) <- candidate_basic

for (v in candidate_basic) {
  x <- outer_train1[[v]]
  if (!anyNA(x)) { meth_cv1[v] <- "" }
  else if (is.ordered(x)) { meth_cv1[v] <- "polr" }
  else if (is.factor(x) && nlevels(x) == 2) { meth_cv1[v] <- "logreg" }
  else if (is.factor(x) && nlevels(x) > 2) { meth_cv1[v] <- "polyreg" }
  else if (is.numeric(x) || is.integer(x)) { meth_cv1[v] <- "pmm" }
}

pred_cv1 <- make.predictorMatrix(cv_imp_data1)
diag(pred_cv1) <- 0

set.seed(20260914)
imp_cv1_smoke <- mice(
  cv_imp_data1,
  m = 2, maxit = 2,
  method = meth_cv1, predictorMatrix = pred_cv1,
  ignore = ignore_test1, printFlag = FALSE, seed = 20260914
)

# ============================================================
# Step 111: Outer fold 1 - TRAIN-only official MICE
# ============================================================
set.seed(20260914)
imp_train1 <- mice(
  mice_train1,
  m = 10, maxit = 10,
  method = meth_train1, predictorMatrix = pred_train1,
  n.core = 5, seed = 20260914, printFlag = FALSE
)

# ============================================================
# Step 112: Outer fold 1 - Official TEST multiple imputation
# ============================================================
set.seed(20260914)
imp_test1 <- mice(
  cv_imp_data1,
  m = 10, maxit = 10,
  method = meth_cv1, predictorMatrix = pred_cv1,
  ignore = ignore_test1, n.core = 5, seed = 20260914, printFlag = FALSE
)

test_imp_list1 <- lapply(1:10, function(m) {
  comp_m <- complete(imp_test1, action = m)
  test_m <- comp_m[(n_train1 + 1):(n_train1 + n_test1), , drop = FALSE]
  test_m$ID <- outer_test1$ID
  test_m$followup_time <- outer_test1$followup_time
  test_m$event_cox <- outer_test1$event_cox
  test_m
})

# ============================================================
# Step 113: Outer fold 1 - TRAIN stacked matrix + inner 10-fold
# ============================================================
train1_long <- complete(imp_train1, action = "long", include = FALSE)
train1_long$ID <- outer_train1$ID[train1_long$.id]
train1_long$followup_time <- outer_train1$followup_time[train1_long$.id]

train1_long$hear <- factor(
  ifelse(as.character(train1_long$hear) %in% c("4", "5"), "4-5", as.character(train1_long$hear)),
  levels = c("1", "2", "3", "4-5"), ordered = FALSE
)

design_train1 <- train1_long %>% select(all_of(candidate_basic))
design_train1[] <- lapply(design_train1, function(x) {
  if (is.ordered(x)) { factor(as.character(x), levels = levels(x), ordered = FALSE) } else { x }
})

X_train1 <- model.matrix(~ ., data = design_train1)[, -1, drop = FALSE]
y_train1 <- Surv(train1_long$followup_time, train1_long$event_cox)

set.seed(20260914)
inner_subjects1 <- outer_train1 %>%
  select(ID, event_cox) %>%
  group_by(event_cox) %>%
  mutate(
    random_order = sample(seq_len(n())),
    inner_fold = rep(1:10, length.out = n())[random_order]
  ) %>%
  ungroup() %>%
  select(-random_order)

inner_fold_map1 <- setNames(inner_subjects1$inner_fold, inner_subjects1$ID)
inner_foldid1 <- unname(inner_fold_map1[as.character(train1_long$ID)])

cat("\n--- Step 113 Completed Successfully ---\n")
cat("TRAIN long dimensions =", nrow(train1_long), "x", ncol(train1_long), "\n")
cat("X dimensions =", nrow(X_train1), "x", ncol(X_train1), "\n")
cat("NA in X =", sum(is.na(X_train1)), "| Infinite in X =", sum(!is.finite(X_train1)), "\n")
cat("IDs appearing in >1 inner fold =", sum(tapply(inner_foldid1, train1_long$ID, function(x) length(unique(x))) > 1), "\n")