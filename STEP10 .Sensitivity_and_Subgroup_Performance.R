# ==============================================================================
# Step 10 Dynamic Predictive Modeling Pipeline: Sensitivity Analyses & Subgroup Performance
# Sensitivity_and_Subgroup_Performance
# ==============================================================================

library(survival)
library(mice)
library(dplyr)

# 确保全局随机种子，保证分析完全可复现
set.seed(20260914)

# ------------------------------------------------------------------------------
# SECTION 1: Event-Time Location Sensitivity (Q25 / Q50 / Q75 Interval Sensitivity)
# ------------------------------------------------------------------------------

# 基于已有的 10 个插补数据集拟合锁定 10 变量 Cox 模型
fit_sensitivity_q <- function(time_var) {
  fits <- lapply(1:10, function(m) {
    dat_m <- basic_long %>% filter(.imp == m)
    form_m <- as.formula(paste0("Surv(", time_var, ", event_cox) ~ ", paste(final_vars_10, collapse = " + ")))
    coxph(form_m, data = dat_m, ties = "efron")
  })
  
  pooled <- pool(as.mira(fits))
  summary(pooled, conf.int = TRUE, exponentiate = TRUE)
}

res_q25 <- fit_sensitivity_q("time_q25")
res_q50 <- fit_sensitivity_q("time_q50")
res_q75 <- fit_sensitivity_q("time_q75")

# 汇总与对比 Pooled HR
compare_q <- data.frame(
  term   = res_q50$term,
  HR_Q25 = res_q25$estimate,
  HR_Q50 = res_q50$estimate,
  HR_Q75 = res_q75$estimate
) %>%
  mutate(
    change_Q25_vs_Q50_pct = 100 * (HR_Q25 / HR_Q50 - 1),
    change_Q75_vs_Q50_pct = 100 * (HR_Q75 / HR_Q50 - 1),
    same_direction        = ((HR_Q25 > 1 & HR_Q50 > 1 & HR_Q75 > 1) | 
                               (HR_Q25 < 1 & HR_Q50 < 1 & HR_Q75 < 1))
  )

# ------------------------------------------------------------------------------
# SECTION 2: Semi-Parametric Interval-Censored PH Model (Single Dataset Check)
# ------------------------------------------------------------------------------

if (requireNamespace("icenReg", quietly = TRUE)) {
  dat_ic <- basic_long %>% filter(.imp == 1)
  
  form_ic <- Surv(interval_left, interval_right, type = "interval2") ~
    fall_down + age + adlab_c + arthre + cesd10 + 
    hear + total_cognition + wspeed + puff + teeth
  
  # 拟合区间删失 PH 模型
  fit_ic <- icenReg::ic_sp(formula = form_ic, data = dat_ic, model = "ph", bs_samples = 0)
  
  # 同一数据集的中点 Cox 模型
  fit_mid <- coxph(
    Surv(followup_time, event_cox) ~ fall_down + age + adlab_c + arthre + cesd10 + 
      hear + total_cognition + wspeed + puff + teeth,
    data = dat_ic, ties = "efron"
  )
  
  beta_ic  <- fit_ic$coefficients
  beta_mid <- coef(fit_mid)
  common_terms <- intersect(names(beta_ic), names(beta_mid))
  
  compare_ic <- data.frame(
    term        = common_terms,
    HR_midpoint = exp(beta_mid[common_terms]),
    HR_interval = exp(beta_ic[common_terms])
  ) %>%
    mutate(
      change_pct     = 100 * (HR_interval / HR_midpoint - 1),
      same_direction = sign(beta_ic[common_terms]) == sign(beta_mid[common_terms])
    )
}

# ------------------------------------------------------------------------------
# SECTION 3: Main Text & Supplementary Model Performance Summaries
# ------------------------------------------------------------------------------

# 1. 论文主文精简性能表
main_performance_table <- data.frame(
  Outcome_time          = c("Overall", "3 years", "5 years", "7 years"),
  Discrimination        = c("C-index = 0.6696", "AUC = 0.6886", "AUC = 0.6838", "AUC = 0.6942"),
  Calibration_intercept = c(NA, 0.0112, 0.0482, 0.0619),
  Calibration_slope     = c(0.9238, 0.9268, 1.0207, 1.0688),
  Brier_score           = c(NA, 0.0295, 0.0306, 0.0431),
  IPA                   = c(NA, 0.0094, 0.0131, 0.0203)
)

# 2. 论文主文模型变量标准化标签字典
final_variable_labels <- c(
  "fall_down1"      = "History of falls",
  "age"             = "Age, years",
  "adlab_c"         = "ADL limitations",
  "arthre1"         = "Arthritis",
  "cesd10"          = "Depressive symptoms (CESD-10)",
  "hear2"           = "Hearing: level 2 vs level 1",
  "hear3"           = "Hearing: level 3 vs level 1",
  "hear4-5"         = "Hearing: levels 4–5 vs level 1",
  "total_cognition" = "Global cognition score",
  "wspeed"          = "2.5-m walking time, seconds",
  "puff"            = "Peak expiratory flow",
  "teeth1"          = "Tooth loss / dental status"
)

# ------------------------------------------------------------------------------
# SECTION 4: Subgroup Performance Assessment (Age, Sex, Residence)
# ------------------------------------------------------------------------------

n_obs <- nrow(analysis_data)
lp_mat <- matrix(NA_real_, nrow = n_obs, ncol = 10)
risk5_mat <- matrix(NA_real_, nrow = n_obs, ncol = 10)

# 从 10 个模型中提取个体预测 LP 与 5 年绝对风险
for (m in 1:10) {
  dat_m <- basic_long %>% 
    filter(.imp == m) %>% 
    arrange(as.integer(.id))
  
  dat_m$hear <- factor(
    ifelse(as.character(dat_m$hear) %in% c("4", "5"), "4-5", as.character(dat_m$hear)),
    levels = c("1", "2", "3", "4-5"), ordered = FALSE
  )
  
  id_m <- as.integer(dat_m$.id)
  lp_m <- predict(final_cox_fits[[m]], newdata = dat_m, type = "lp", reference = "zero")
  
  bh_m <- basehaz(final_cox_fits[[m]], centered = FALSE)
  H0_5_m <- approx(x = bh_m$time, y = bh_m$hazard, xout = 5, method = "constant", rule = 2, f = 0)$y
  
  lp_mat[id_m, m]    <- lp_m
  risk5_mat[id_m, m] <- 1 - exp(-H0_5_m * exp(lp_m))
}

analysis_subgroup <- analysis_data %>%
  mutate(
    lp_mean    = rowMeans(lp_mat, na.rm = TRUE),
    risk5_mean = rowMeans(risk5_mat, na.rm = TRUE),
    age_group  = cut(age, breaks = c(60, 70, 80, Inf), right = FALSE, labels = c("60–69", "70–79", "≥80"))
  )

# 计算亚组 C-index 及其解析 95% CI，以及 5 年预测 vs 观察风险
compute_subgroup_metrics <- function(data, variable, dimension_name) {
  groups <- unique(as.character(data[[variable]]))
  groups <- groups[!is.na(groups)]
  
  res <- lapply(groups, function(g) {
    d <- data[as.character(data[[variable]]) == g, , drop = FALSE]
    
    # C-index
    conc <- concordance(Surv(followup_time, event_cox) ~ lp_mean, data = d, reverse = TRUE)
    C    <- conc$concordance
    SE   <- sqrt(conc$var)
    
    # Kaplan-Meier 观察风险 (5年)
    km   <- survfit(Surv(followup_time, event_cox) ~ 1, data = d)
    obs5 <- 1 - summary(km, times = 5, extend = TRUE)$surv
    
    data.frame(
      Dimension    = dimension_name,
      Group        = g,
      N            = nrow(d),
      Events       = sum(d$event_cox),
      C_index      = C,
      CI_low       = max(0, C - 1.96 * SE),
      CI_high      = min(1, C + 1.96 * SE),
      Predicted_5y = mean(d$risk5_mean),
      Observed_5y  = obs5
    )
  })
  
  bind_rows(res)
}

# 提取各临床维度性能
subgroup_results_raw <- bind_rows(
  compute_subgroup_metrics(analysis_subgroup, "age_group", "Age group"),
  compute_subgroup_metrics(analysis_subgroup, "gender",    "Sex"),
  compute_subgroup_metrics(analysis_subgroup, "rural",     "Place of residence")
)

# 格式化亚组表格（使用标准化展示标签）
subgroup_table_final <- subgroup_results_raw %>%
  mutate(
    Group = case_when(
      Dimension == "Sex" & Group == "0" ~ "Female",
      Dimension == "Sex" & Group == "1" ~ "Male",
      Dimension == "Place of residence" & Group == "0" ~ "Urban",
      Dimension == "Place of residence" & Group == "1" ~ "Rural",
      TRUE ~ as.character(Group)
    ),
    `C-index (95% CI)`      = sprintf("%.3f (%.3f–%.3f)", C_index, CI_low, CI_high),
    `Predicted 5-y risk, %` = sprintf("%.2f", 100 * Predicted_5y),
    `Observed 5-y risk, %`  = sprintf("%.2f", 100 * Observed_5y)
  ) %>%
  select(
    Dimension, Group, N, Events, 
    `C-index (95% CI)`, `Predicted 5-y risk, %`, `Observed 5-y risk, %`
  )

# ------------------------------------------------------------------------------
# SECTION 5: Output Formatted Results to Console
# ------------------------------------------------------------------------------

cat("=== SECTION 1: Q25 / Q50 / Q75 Hazard Ratio Sensitivity ===\n")
cat("Directional consistency across all variables:", sum(compare_q$same_direction), "/", nrow(compare_q), "\n")
cat("Max absolute HR shift vs midpoint (%):", round(max(abs(c(compare_q$change_Q25_vs_Q50_pct, compare_q$change_Q75_vs_Q50_pct))), 2), "%\n\n")

if (exists("compare_ic")) {
  cat("=== SECTION 2: Semi-Parametric Interval-Censored PH Check ===\n")
  cat("Directional consistency with Cox:", sum(compare_ic$same_direction), "/", nrow(compare_ic), "\n")
  cat("Max absolute HR change (%):", round(max(abs(compare_ic$change_pct)), 2), "%\n\n")
}

cat("=== SECTION 3: Main Model Performance Table ===\n")
print(main_performance_table, row.names = FALSE)
cat("\n")

cat("=== SECTION 4: Subgroup Performance & Consistency Table ===\n")
print(as.data.frame(subgroup_table_final), row.names = FALSE)