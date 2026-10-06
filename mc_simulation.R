# ==============================================================================
# Monte Carlo Simulation: Empirical Gram-Schmidt (EGS) Polynomial Series
# ==============================================================================

library(ggplot2)
library(tikzDevice)

# --- 1. Global Settings ---

set.seed(2026)

M <- 1000
sample_sizes <- c(100, 2500)

max_d <- 10
reported_degrees <- 1:9
plot_degrees <- c(3, 5, 7)

n_grid <- 200
x_grid <- seq(0, 1, length.out = n_grid)
true_mean_fn <- function(x) sin(2 * pi * x)

output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# --- 2. Basis Functions and Estimators ---

raw_basis <- function(x, d) {
  outer(x, 0:d, "^")
}

theoretical_basis <- function(x, d) {
  n <- length(x)
  Phi_T <- matrix(0, nrow = n, ncol = d + 1)
  z <- 2 * x - 1
  
  P_prev2 <- rep(1, n)
  Phi_T[, 1] <- P_prev2
  
  if (d >= 1) {
    P_prev1 <- z
    Phi_T[, 2] <- sqrt(3) * P_prev1
  }
  if (d >= 2) {
    for (j in 2:d) {
      P_curr <- ((2 * j - 1) * z * P_prev1 - (j - 1) * P_prev2) / j
      Phi_T[, j + 1] <- sqrt(2 * j + 1) * P_curr
      P_prev2 <- P_prev1
      P_prev1 <- P_curr
    }
  }
  return(Phi_T)
}

empirical_basis <- function(H) {
  n <- nrow(H)
  qr_decomp <- qr(H, tol = 0, LAPACK = FALSE)
  if (!all(qr_decomp$pivot == seq_len(ncol(H)))) {
    stop("QR changed the monomial ordering.")
  }
  
  Q <- qr.Q(qr_decomp)
  R <- qr.R(qr_decomp)
  
  s <- sign(diag(R))
  s[s == 0] <- 1
  Q <- sweep(Q, 2, s, "*")
  
  return(sqrt(n) * Q)
}

fit_model <- function(X, y, is_empirical = FALSE) {
  n <- nrow(X)
  if (is_empirical) {
    coefs <- as.vector(crossprod(X, y) / n)
  } else {
    coefs <- as.vector(qr.solve(X, y))
  }
  fitted_vals <- as.vector(X %*% coefs)
  list(coefs = coefs, fitted = fitted_vals)
}

calc_full_stats <- function(v) {
  c(
    Mean   = mean(v),
    SD     = sd(v),
    Median = median(v),
    Min    = min(v),
    Max    = max(v),
    P02.5  = as.numeric(quantile(v, 0.025)),
    P25    = as.numeric(quantile(v, 0.25)),
    P75    = as.numeric(quantile(v, 0.75)),
    P95    = as.numeric(quantile(v, 0.95)),
    P97.5  = as.numeric(quantile(v, 0.975))
  )
}


# --- 3. Simulation Runner ---

run_scenario <- function(n_obs, M, max_d = 10, plot_degrees = c(3, 5, 7)) {
  cat(sprintf("Simulating n = %d, M = %d...\n", n_obs, M))
  
  Phi_T_grid_all <- theoretical_basis(x_grid, max_d)
  
  log_cond_R <- matrix(0, nrow = M, ncol = max_d)
  log_cond_T <- matrix(0, nrow = M, ncol = max_d)
  log_cond_E <- matrix(0, nrow = M, ncol = max_d)
  
  diff_fit_R <- matrix(0, nrow = M, ncol = max_d)
  diff_fit_T <- matrix(0, nrow = M, ncol = max_d)
  
  diff_seq_R <- matrix(0, nrow = M, ncol = max_d - 1)
  diff_seq_T <- matrix(0, nrow = M, ncol = max_d - 1)
  diff_seq_E <- matrix(0, nrow = M, ncol = max_d - 1)
  
  grid_fits <- array(0, dim = c(M, n_grid, length(plot_degrees)))
  
  for (m in 1:M) {
    X <- runif(n_obs, 0, 1)
    Y <- true_mean_fn(X) + rnorm(n_obs, 0, 1)
    
    H_full     <- raw_basis(X, max_d)
    Phi_T_full <- theoretical_basis(X, max_d)
    
    coef_R_by_d <- vector("list", max_d)
    coef_T_by_d <- vector("list", max_d)
    coef_E_by_d <- vector("list", max_d)
    
    for (d in 1:max_d) {
      H_d     <- H_full[, 1:(d + 1), drop = FALSE]
      Phi_T_d <- Phi_T_full[, 1:(d + 1), drop = FALSE]
      Phi_E_d <- empirical_basis(H_d)
      
      log_cond_R[m, d] <- log10(kappa(crossprod(H_d) / n_obs, exact = TRUE))
      log_cond_T[m, d] <- log10(kappa(crossprod(Phi_T_d) / n_obs, exact = TRUE))
      log_cond_E[m, d] <- log10(kappa(crossprod(Phi_E_d) / n_obs, exact = TRUE))
      
      fit_R <- fit_model(H_d, Y, is_empirical = FALSE)
      fit_T <- fit_model(Phi_T_d, Y, is_empirical = FALSE)
      fit_E <- fit_model(Phi_E_d, Y, is_empirical = TRUE)
      
      coef_R_by_d[[d]] <- fit_R$coefs
      coef_T_by_d[[d]] <- fit_T$coefs
      coef_E_by_d[[d]] <- fit_E$coefs
      
      diff_fit_R[m, d] <- max(abs(fit_R$fitted - fit_E$fitted))
      diff_fit_T[m, d] <- max(abs(fit_T$fitted - fit_E$fitted))
      
      if (d %in% plot_degrees) {
        p_idx <- which(plot_degrees == d)
        Phi_T_grid_d <- Phi_T_grid_all[, 1:(d + 1), drop = FALSE]
        grid_fits[m, , p_idx] <- as.vector(Phi_T_grid_d %*% fit_T$coefs)
      }
    }
    
    for (d in 1:(max_d - 1)) {
      diff_seq_R[m, d] <- max(abs(coef_R_by_d[[d + 1]][1:(d + 1)] - coef_R_by_d[[d]]))
      diff_seq_T[m, d] <- max(abs(coef_T_by_d[[d + 1]][1:(d + 1)] - coef_T_by_d[[d]]))
      diff_seq_E[m, d] <- max(abs(coef_E_by_d[[d + 1]][1:(d + 1)] - coef_E_by_d[[d]]))
    }
  }
  
  list(
    n_obs      = n_obs,
    log_cond_R = log_cond_R,
    log_cond_T = log_cond_T,
    log_cond_E = log_cond_E,
    diff_fit_R = diff_fit_R,
    diff_fit_T = diff_fit_T,
    diff_seq_R = diff_seq_R,
    diff_seq_T = diff_seq_T,
    diff_seq_E = diff_seq_E,
    grid_fits  = grid_fits
  )
}

results <- list(
  uniform_100  = run_scenario(100,  M, max_d = max_d, plot_degrees = plot_degrees),
  uniform_2500 = run_scenario(2500, M, max_d = max_d, plot_degrees = plot_degrees)
)

saveRDS(results, file.path(output_dir, "monte_carlo_results.rds"))


# --- 4. Export Tables (CSV) ---

tbl2_compact_list <- list()
tbl2_full_list    <- list()

for (sc_name in names(results)) {
  res <- results[[sc_name]]
  for (d_val in reported_degrees) {
    s_R <- calc_full_stats(res$log_cond_R[, d_val])
    s_T <- calc_full_stats(res$log_cond_T[, d_val])
    s_E <- calc_full_stats(res$log_cond_E[, d_val])
    
    tbl2_compact_list[[length(tbl2_compact_list) + 1]] <- data.frame(
      Sample_Size = paste0("n = ", res$n_obs),
      Degree      = paste0("d = ", d_val),
      Raw_Median  = sprintf("%.2f [%.2f, %.2f]", s_R["Median"], s_R["P25"], s_R["P75"]),
      Theory_Med  = sprintf("%.2f [%.2f, %.2f]", s_T["Median"], s_T["P25"], s_T["P75"]),
      Emp_Median  = sprintf("%.2f [%.2f, %.2f]", s_E["Median"], s_E["P25"], s_E["P75"])
    )
    
    for (b_name in c("Raw", "Theoretical", "Empirical")) {
      st <- switch(b_name, "Raw" = s_R, "Theoretical" = s_T, "Empirical" = s_E)
      tbl2_full_list[[length(tbl2_full_list) + 1]] <- data.frame(
        Sample_Size = paste0("n = ", res$n_obs),
        Degree      = paste0("d = ", d_val),
        Basis       = b_name,
        t(st)
      )
    }
  }
}

write.csv(do.call(rbind, tbl2_compact_list), file.path(output_dir, "table2_condition_numbers.csv"), row.names = FALSE)
write.csv(do.call(rbind, tbl2_full_list),    file.path(output_dir, "table2_condition_numbers_full.csv"), row.names = FALSE)

tbl3_compact_list <- list()
tbl3_full_list    <- list()

for (sc_name in names(results)) {
  res <- results[[sc_name]]
  for (d_val in reported_degrees) {
    st_fit_R <- calc_full_stats(res$diff_fit_R[, d_val])
    st_fit_T <- calc_full_stats(res$diff_fit_T[, d_val])
    st_seq_R <- calc_full_stats(res$diff_seq_R[, d_val])
    st_seq_T <- calc_full_stats(res$diff_seq_T[, d_val])
    st_seq_E <- calc_full_stats(res$diff_seq_E[, d_val])
    
    tbl3_compact_list[[length(tbl3_compact_list) + 1]] <- data.frame(
      Sample_Size   = paste0("n = ", res$n_obs),
      Degree        = paste0("d = ", d_val),
      Fit_Diff_Raw  = sprintf("%.2e (%.2e)", st_fit_R["Median"], st_fit_R["Max"]),
      Fit_Diff_Th   = sprintf("%.2e (%.2e)", st_fit_T["Median"], st_fit_T["Max"]),
      Seq_Shift_Raw = sprintf("%.2e", st_seq_R["Median"]),
      Seq_Shift_Th  = sprintf("%.2e", st_seq_T["Median"]),
      Seq_Shift_Emp = sprintf("%.2e", st_seq_E["Median"])
    )
    
    diag_items <- list(
      "Fit_Diff_Raw"  = st_fit_R,
      "Fit_Diff_Th"   = st_fit_T,
      "Seq_Shift_Raw" = st_seq_R,
      "Seq_Shift_Th"  = st_seq_T,
      "Seq_Shift_Emp" = st_seq_E
    )
    for (item_name in names(diag_items)) {
      tbl3_full_list[[length(tbl3_full_list) + 1]] <- data.frame(
        Sample_Size = paste0("n = ", res$n_obs),
        Degree      = paste0("d = ", d_val),
        Metric      = item_name,
        t(diag_items[[item_name]])
      )
    }
  }
}

write.csv(do.call(rbind, tbl3_compact_list), file.path(output_dir, "table3_numerical_diagnostics.csv"), row.names = FALSE)
write.csv(do.call(rbind, tbl3_full_list),    file.path(output_dir, "table3_numerical_diagnostics_full.csv"), row.names = FALSE)


# --- 5. Export Figures (TikZ) ---

# Figure 1: Fitted Curves
df_fig1_list <- list()
for (res in results) {
  for (p_idx in seq_along(plot_degrees)) {
    d_val <- plot_degrees[p_idx]
    grid_mat <- res$grid_fits[, , p_idx]
    
    df_fig1_list[[length(df_fig1_list) + 1]] <- data.frame(
      x       = x_grid,
      true_y  = true_mean_fn(x_grid),
      mean_y  = colMeans(grid_mat),
      lower_y = apply(grid_mat, 2, quantile, probs = 0.025),
      upper_y = apply(grid_mat, 2, quantile, probs = 0.975),
      degree  = paste0("$d = ", d_val, "$"),
      n_label = paste0("$n = ", res$n_obs, "$")
    )
  }
}
df_fig1 <- do.call(rbind, df_fig1_list)
df_fig1$n_label <- factor(df_fig1$n_label, levels = c("$n = 100$", "$n = 2500$"))

p1_tikz <- ggplot(df_fig1, aes(x = x)) +
  geom_ribbon(aes(ymin = lower_y, ymax = upper_y, fill = "2.5th--97.5th Percentiles"), alpha = 0.25) +
  geom_line(aes(y = mean_y, color = "Monte Carlo Mean"), linewidth = 0.8, linetype = "dashed") +
  geom_line(aes(y = true_y, color = "True $g(x) = \\sin(2\\pi x)$"), linewidth = 0.8) +
  facet_grid(n_label ~ degree) +
  scale_color_manual(name = "", values = c("True $g(x) = \\sin(2\\pi x)$" = "black", "Monte Carlo Mean" = "#0055D4")) +
  scale_fill_manual(name = "", values = c("2.5th--97.5th Percentiles" = "#0055D4")) +
  labs(x = "Regressor $x$", y = "$g(x)$") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text = element_text(size = 9))

tikz(file.path(output_dir, "fig1_uniform_fits.tex"), width = 6.2, height = 3.6, standAlone = FALSE)
print(p1_tikz)
dev.off()
