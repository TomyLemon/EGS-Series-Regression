# ==============================================================================
# Empirical Analysis: Nonlinear Market Exposure via EGS Polynomial Series
# Sample Period: July 2016 - June 2026
# ==============================================================================

library(quantmod)
library(sandwich)
library(ggplot2)
library(tikzDevice)

# --- 1. Global Settings and Data Ingestion ---

output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

data_file <- file.path(output_dir, "market_returns_2016_2026.csv")

if (file.exists(data_file)) {
  data_mat <- read.csv(data_file, row.names = 1)
} else {
  getSymbols(c("AAPL", "NVDA", "^GSPC"), src = "yahoo", 
             from = "2016-07-01", to = "2026-07-01")
  
  r_aapl <- 100 * diff(log(Ad(AAPL)))
  r_nvda <- 100 * diff(log(Ad(NVDA)))
  r_sp   <- 100 * diff(log(Ad(GSPC)))
  
  data_mat <- as.data.frame(na.omit(merge(r_aapl, r_nvda, r_sp)))
  colnames(data_mat) <- c("AAPL", "NVDA", "SP500")
  
  write.csv(data_mat, data_file)
}

n <- nrow(data_mat)
Y_aapl <- as.numeric(data_mat$AAPL)
Y_nvda <- as.numeric(data_mat$NVDA)
R_m    <- as.numeric(data_mat$SP500)

R_min   <- min(R_m)
R_max   <- max(R_m)
R_range <- R_max - R_min
X       <- (R_m - R_min) / R_range


# --- 2. Table 1: Descriptive Statistics ---

calc_stats <- function(v) {
  v_mean <- mean(v)
  s2     <- mean((v - v_mean)^2)
  c(
    Mean            = v_mean,
    SD              = sd(v),
    Min             = min(v),
    Max             = max(v),
    Skewness        = mean((v - v_mean)^3) / (s2^(3 / 2)),
    Excess_Kurtosis = mean((v - v_mean)^4) / (s2^2) - 3
  )
}

table1 <- round(rbind(
  AAPL  = calc_stats(Y_aapl),
  NVDA  = calc_stats(Y_nvda),
  SP500 = calc_stats(R_m)
), 3)

write.csv(table1, file.path(output_dir, "Table1_Summary_Statistics.csv"))


# --- 3. Empirical Gram-Schmidt (EGS) Basis Construction ---

build_egs <- function(x_vals, degree) {
  n_obs <- length(x_vals)
  H <- outer(x_vals, 0:degree, "^")
  
  qr_fit <- qr(H / sqrt(n_obs), tol = 0, LAPACK = FALSE)
  U <- qr.R(qr_fit)
  
  signs <- sign(diag(U))
  signs[signs == 0] <- 1
  U <- diag(signs) %*% U
  
  U_inv <- backsolve(U, diag(degree + 1))
  Phi   <- H %*% U_inv
  
  list(Phi = Phi, U_inv = U_inv)
}


# --- 4. Table 3: Degree Selection Criteria ---

degree_selection <- function(Y, max_degree = 9) {
  res_mat <- matrix(0, nrow = max_degree, ncol = 3)
  colnames(res_mat) <- c("CV", "GCV", "Max_Leverage")
  
  for (d in 1:max_degree) {
    egs <- build_egs(X, d)
    fit <- lm(Y ~ egs$Phi[, -1, drop = FALSE])
    
    res <- residuals(fit)
    h   <- hatvalues(fit)
    
    cv_val  <- mean((res / (1 - h))^2)
    gcv_val <- mean(res^2) / (1 - (d + 1) / n)^2
    
    res_mat[d, ] <- c(cv_val, gcv_val, max(h))
  }
  as.data.frame(res_mat)
}

sel_aapl <- degree_selection(Y_aapl)
sel_nvda <- degree_selection(Y_nvda)

table3 <- data.frame(
  Degree       = 1:9,
  AAPL_CV      = sel_aapl$CV,
  AAPL_GCV     = sel_aapl$GCV,
  NVDA_CV      = sel_nvda$CV,
  NVDA_GCV     = sel_nvda$GCV,
  Max_Leverage = sel_aapl$Max_Leverage
)

write.csv(round(table3, 4), file.path(output_dir, "Table3_Degree_Selection.csv"), row.names = FALSE)

degree_aapl <- which.min(sel_aapl$GCV)
degree_nvda <- which.min(sel_nvda$GCV)


# --- 5. Estimation, HAC Inference, and TikZ Export ---

analyze_asset <- function(asset_name, Y, selected_degree) {
  egs <- build_egs(X, selected_degree)
  fit <- lm(Y ~ egs$Phi[, -1, drop = FALSE])
  coefs <- coef(fit)
  
  hac_cov <- NeweyWest(fit, lag = 5, prewhite = FALSE, adjust = FALSE)
  hac_se  <- sqrt(diag(hac_cov))
  z_stat  <- coefs / hac_se
  p_val   <- 2 * pnorm(-abs(z_stat))
  
  tss <- sum((Y - mean(Y))^2)
  delta_r2 <- c(NA, 100 * n * coefs[-1]^2 / tss)
  cum_r2   <- c(NA, cumsum(delta_r2[-1]))
  
  estimates_table <- data.frame(
    Direction             = paste0("phi_", 0:selected_degree),
    Estimate              = coefs,
    HAC_SE                = hac_se,
    z_Statistic           = z_stat,
    p_Value               = p_val,
    Delta_R2_Percent      = delta_r2,
    Cumulative_R2_Percent = cum_r2
  )
  write.csv(estimates_table, file.path(output_dir, paste0("Table_", asset_name, "_EGS_Estimates.csv")), row.names = FALSE)
  
  c_high <- coefs[-c(1, 2)]
  V_high <- hac_cov[-c(1, 2), -c(1, 2), drop = FALSE]
  
  wald_stat <- as.numeric(t(c_high) %*% solve(V_high) %*% c_high)
  wald_df   <- selected_degree - 1
  wald_p    <- pchisq(wald_stat, df = wald_df, lower.tail = FALSE)
  
  wald_result <- data.frame(
    Asset              = asset_name,
    Selected_Degree    = selected_degree,
    Wald_Statistic     = wald_stat,
    Degrees_of_Freedom = wald_df,
    p_Value            = wald_p
  )
  
  x_grid  <- seq(0, 1, length.out = 200)
  rm_grid <- R_min + x_grid * R_range
  
  H_grid   <- outer(x_grid, 0:selected_degree, "^")
  Phi_grid <- H_grid %*% egs$U_inv
  
  D_raw <- matrix(0, nrow = length(x_grid), ncol = selected_degree + 1)
  for (k in 1:selected_degree) {
    D_raw[, k + 1] <- k * x_grid^(k - 1)
  }
  D_phi       <- D_raw %*% egs$U_inv
  beta_design <- D_phi / R_range
  
  fit_mean <- as.numeric(Phi_grid %*% coefs)
  mean_se  <- sqrt(pmax(rowSums((Phi_grid %*% hac_cov) * Phi_grid), 0))
  
  dyn_beta <- as.numeric(beta_design %*% coefs)
  beta_se  <- sqrt(pmax(rowSums((beta_design %*% hac_cov) * beta_design), 0))
  
  func_df <- data.frame(
    Market_Return = rm_grid,
    Fitted_Mean   = fit_mean,
    Mean_CI_Lower = fit_mean - 1.96 * mean_se,
    Mean_CI_Upper = fit_mean + 1.96 * mean_se,
    Dynamic_Beta  = dyn_beta,
    Beta_CI_Lower = dyn_beta - 1.96 * beta_se,
    Beta_CI_Upper = dyn_beta + 1.96 * beta_se
  )
  write.csv(func_df, file.path(output_dir, paste0("Figure_", asset_name, "_Function_Data.csv")), row.names = FALSE)
  
  fit_linear   <- lm(Y ~ R_m)
  linear_alpha <- coef(fit_linear)[1]
  linear_beta  <- coef(fit_linear)[2]
  
  rm_q01 <- quantile(R_m, 0.01)
  rm_q99 <- quantile(R_m, 0.99)
  
  scatter_df <- data.frame(Rm = R_m, Ri = Y)
  p1 <- ggplot() +
    geom_point(data = scatter_df, aes(x = Rm, y = Ri), alpha = 0.12, size = 0.6, color = "gray30") +
    geom_ribbon(data = func_df, aes(x = Market_Return, ymin = Mean_CI_Lower, ymax = Mean_CI_Upper), fill = "#0072B2", alpha = 0.25) +
    geom_line(data = func_df, aes(x = Market_Return, y = Fitted_Mean, color = "EGS Fit"), linewidth = 0.8) +
    geom_abline(intercept = linear_alpha, slope = linear_beta, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
    scale_color_manual(name = "", values = c("EGS Fit" = "#0072B2")) +
    labs(title = paste0("(a) Fitted Conditional Mean $\\widehat{m}_{\\mathrm{", asset_name, "}}(R_m)$"),
         x = "S\\&P 500 Daily Return (\\%)",
         y = paste0(asset_name, " Daily Return (\\%)")) +
    theme_bw(base_size = 9.5) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
  
  p2 <- ggplot(func_df, aes(x = Market_Return, y = Dynamic_Beta)) +
    geom_ribbon(aes(ymin = Beta_CI_Lower, ymax = Beta_CI_Upper), fill = "#009E73", alpha = 0.2) +
    geom_line(color = "#009E73", linewidth = 0.8) +
    geom_hline(yintercept = linear_beta, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
    geom_vline(xintercept = c(rm_q01, rm_q99), linetype = "dotted", color = "#7B3294", linewidth = 1.5) +
    labs(title = paste0("(b) Dynamic Market Beta $\\widehat{\\beta}_{\\mathrm{", asset_name, "}}^{\\mathrm{GS}}(R_m)$"),
         x = "S\\&P 500 Daily Return (\\%)",
         y = "Dynamic Market Beta") +
    theme_bw(base_size = 9.5) +
    theme(panel.grid.minor = element_blank())
  
  f1_path <- file.path(output_dir, paste0("Figure_", asset_name, "_Conditional_Mean.tex"))
  f2_path <- file.path(output_dir, paste0("Figure_", asset_name, "_Dynamic_Beta.tex"))
  
  tikz(f1_path, width = 6, height = 3, standAlone = FALSE)
  print(p1)
  dev.off()
  
  tikz(f2_path, width = 6, height = 3, standAlone = FALSE)
  print(p2)
  dev.off()
  
  
  list(wald = wald_result)
}

# --- 6. Execution ---

res_aapl <- analyze_asset("AAPL", Y_aapl, degree_aapl)
res_nvda <- analyze_asset("NVDA", Y_nvda, degree_nvda)

table_wald <- rbind(res_aapl$wald, res_nvda$wald)
write.csv(table_wald, file.path(output_dir, "Table_Wald_Linearity_Tests.csv"), row.names = FALSE)
