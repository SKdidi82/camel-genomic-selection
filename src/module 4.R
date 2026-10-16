
## Module 4: Joint demographic-genetic model + Sobol sensitivity analysis


## ---- Packages -----------------------------------------
pkgs <- c("ggplot2", "scales", "dplyr", "tidyr", "openxlsx")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)
library(openxlsx)

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

## ---- Colour palette ------------------------------

col_primary <- "#2C3E50"   # deep charcoal-navy (Total-order ST)
col_accent  <- "#E67E22"   # warm amber (First-order S1)
col_soft    <- "#7F8C8D"

## ---- theme -----------------------------------------------
theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") %+replace%
    theme(
      axis.title   = element_text(size = base_size, face = "plain", colour = col_primary),
      axis.text    = element_text(size = base_size * 0.85, colour = col_primary),
      axis.line    = element_line(linewidth = 0.45, colour = col_primary),
      axis.ticks   = element_line(linewidth = 0.35, colour = col_primary),
      plot.title   = element_text(size = base_size * 1.0, face = "bold", hjust = 0.5,
                                  colour = col_primary, margin = margin(b = 4)),
      plot.margin  = margin(4, 8, 4, 4),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text  = element_text(size = base_size * 0.8, colour = col_primary),
      legend.key.size = unit(0.4, "cm"),
      panel.grid   = element_blank()
    )
}

## ---- 1. Parameters -------------------------------
pn <- c("h2","Ne","G","n_dau","a_b","tau","alpha","ell","delta","s","j","I","p_pt","p_gs","p_dam")
lo <- c(h2=.10, Ne=50,  G=20, n_dau=10, a_b=4, tau=1.04, alpha=4, ell=.8,  delta=.25, s=.88, j=.70, I=1.5, p_pt=.05, p_gs=.02, p_dam=.3)
hi <- c(h2=.40, Ne=500, G=35, n_dau=40, a_b=7, tau=1.12, alpha=7, ell=1.5, delta=1.0,  s=.96, j=.90, I=3.0, p_pt=.30, p_gs=.20, p_dam=.8)
base <- c(h2=.25, Ne=100, G=30, n_dau=20, a_b=5, tau=1.08, alpha=5, ell=1, delta=.5, s=.93, j=.8, I=2, p_pt=.1, p_gs=.1, p_dam=.5)
decisions <- c("n_dau","a_b","p_pt","p_gs","p_dam")

scale_params <- function(U, nm = pn) {
  P <- list()
  for (i in seq_along(nm)) {
    n <- nm[i]; u <- U[, i]
    P[[n]] <- if (n == "Ne") lo[[n]] * (hi[[n]] / lo[[n]])^u else lo[[n]] + (hi[[n]] - lo[[n]]) * u
  }
  P
}

## ---- 2. Joint model --------------------------------
lam_vec <- function(alpha, s, j, m, it = 80) {
  a <- s + 1e-12; b <- rep(3, length(alpha))
  Ff <- function(L) L^alpha - s * L^(alpha - 1) - m * j * s^(alpha - 1)
  for (k in seq_len(it)) {
    mid <- (a + b) / 2; neg <- Ff(mid) < 0
    a <- ifelse(neg, mid, a); b <- ifelse(neg, b, mid)
  }
  (a + b) / 2
}
i_trunc <- function(p) dnorm(qnorm(1 - p)) / p

model <- function(P, Nref = 1000, fixT = NULL) {
  lam <- lam_vec(P$alpha, P$s, P$j, 0.5 / P$I)
  Tg  <- if (is.null(fixT)) P$alpha + P$s / (lam - P$s) else fixT
  Lpt <- P$a_b + 2 * P$tau + P$alpha + P$ell + P$delta
  Lg  <- P$a_b + P$tau
  Me  <- 2 * P$Ne * P$G / log(4 * P$Ne * P$G); H <- P$h2
  k   <- (4 - P$h2) / P$h2; rpt <- sqrt(P$n_dau / (P$n_dau + k)); h <- sqrt(P$h2)
  ipt <- i_trunc(P$p_pt); ig <- i_trunc(P$p_gs); idm <- i_trunc(P$p_dam)
  rraw  <- ((ipt * rpt + idm * h) * (Lg + Tg) / (Lpt + Tg) - idm * h) / ig
  r     <- pmin(pmax(rraw, 0.02), 0.98)
  Nstar <- Me * r^2 / (H * (1 - r^2))
  rg    <- sqrt(Nref * H / (Nref * H + Me))
  dGgs  <- (2 * ig * rg + 2 * idm * h) / (2 * Lg + 2 * Tg)
  dGpt  <- (2 * ipt * rpt + 2 * idm * h) / (2 * Lpt + 2 * Tg)
  list(logN = log10(Nstar), dG = dGgs, ratio = dGgs / dGpt, rraw = rraw, lam = lam, T = Tg)
}

## ---- 3. Local elasticity check  --------------------
Pb <- as.list(base); m0 <- model(Pb)
cat(sprintf("Baseline: lambda = %.4f, T = %.3f y, required accuracy = %.3f, N* = %.0f\n",
            m0$lam, m0$T, m0$rraw, 10^m0$logN))
elast <- function(par, mode) {
  up <- Pb; dn <- Pb; up[[par]] <- Pb[[par]] * 1.01; dn[[par]] <- Pb[[par]] * 0.99
  fx <- if (mode == "fixed") m0$T else NULL
  (model(up, fixT = fx)$logN - model(dn, fixT = fx)$logN) * log(10) / (log(1.01) - log(0.99))
}
pars <- c("alpha", "s", "j", "I", "a_b", "Ne", "h2")
el_tab <- data.frame(parameter = pars,
                     T_fixed = sapply(pars, elast, mode = "fixed"),
                     T_joint = sapply(pars, elast, mode = "joint"))
cat("\nElasticity of N* (T fixed vs T joint):\n")
print(round(el_tab[, -1], 2))

## ---- 4. Sobol indices  --------------------------------
sobol_indices <- function(f, k, N, nboot = 200, seed = 1) {
  set.seed(seed)
  A <- matrix(runif(N * k), N, k); B <- matrix(runif(N * k), N, k)
  fA <- f(A); fB <- f(B)
  mu <- mean(c(fA, fB))          # centre the output
  fA <- fA - mu; fB <- fB - mu
  V  <- var(c(fA, fB))
  fAB <- matrix(0, N, k)
  for (i in seq_len(k)) { AB <- A; AB[, i] <- B[, i]; fAB[, i] <- f(AB) - mu }
  S1 <- sapply(seq_len(k), function(i) mean(fB * (fAB[, i] - fA)) / V)
  ST <- sapply(seq_len(k), function(i) 0.5 * mean((fA - fAB[, i])^2) / V)
  b1 <- bT <- matrix(0, nboot, k)
  for (b in seq_len(nboot)) {
    ix <- sample.int(N, N, replace = TRUE); Vb <- var(c(fA[ix], fB[ix]))
    b1[b, ] <- sapply(seq_len(k), function(i) mean(fB[ix] * (fAB[ix, i] - fA[ix])) / Vb)
    bT[b, ] <- sapply(seq_len(k), function(i) 0.5 * mean((fA[ix] - fAB[ix, i])^2) / Vb)
  }
  list(S1 = S1, ST = ST,
       S1_lo = apply(b1, 2, quantile, 0.025), S1_hi = apply(b1, 2, quantile, 0.975),
       ST_lo = apply(bT, 2, quantile, 0.025), ST_hi = apply(bT, 2, quantile, 0.975))
}
report <- function(res, nm, top = 8) {
  d <- data.frame(parameter = nm, S1 = res$S1, S1_lo = res$S1_lo, S1_hi = res$S1_hi,
                  ST = res$ST, ST_lo = res$ST_lo, ST_hi = res$ST_hi)
  d <- d[order(-d$ST), ]; cat(sprintf("  sum of S1 = %.2f\n", sum(res$S1)))
  dd <- round(d[seq_len(top), -1], 3); rownames(dd) <- d$parameter[seq_len(top)]
  print(dd)
  invisible(d)
}
N <- 2^14

outs <- c(logN = "log10 N* (required reference size)",
          dG   = "Annual gain, genomic sires (N = 1000)",
          ratio = "Gain ratio (genomic / progeny test)")
cat("\n=== A. All 15 parameters varied ===\n")
resA <- list()
for (o in names(outs)) {
  cat(outs[[o]], "\n")
  resA[[o]] <- sobol_indices(function(U) model(scale_params(U))[[o]], length(pn), N)
  report(resA[[o]], pn)
}

cat("\n=== B. Uncertain parameters only ===\n")
unc <- setdiff(pn, decisions)
scale_unc <- function(U) {
  P <- scale_params(U, unc)
  for (n in decisions) P[[n]] <- rep(base[[n]], nrow(U))
  P
}
resB <- list()
for (o in names(outs)) {
  cat(outs[[o]], "\n")
  resB[[o]] <- sobol_indices(function(U) model(scale_unc(U))[[o]], length(unc), N, nboot = 100)
  report(resB[[o]], unc, top = 6)
}
for (o in names(outs)) {
  write.csv(report(resA[[o]], pn,  top = length(pn)),
            sprintf("tables/sobol_%s_all.csv", o), row.names = FALSE)
  write.csv(report(resB[[o]], unc, top = length(unc)),
            sprintf("tables/sobol_%s_uncertain.csv", o), row.names = FALSE)
}
## ---- 5. Descriptive summary  ------------------------------
set.seed(5); Us <- matrix(runif(2^14 * length(pn)), ncol = length(pn)); ms <- model(scale_params(Us))
cat(sprintf("\nShare of samples with required accuracy outside (0.02, 0.98): %.3f\n",
            mean(ms$rraw <= 0.02 | ms$rraw >= 0.98)))
cat(sprintf("N* median %.0f (5–95%%: %.0f – %.0f)\n", 10^median(ms$logN),
            10^quantile(ms$logN, .05), 10^quantile(ms$logN, .95)))
cat(sprintf("Share with genomic gain > progeny-test gain at N = 1000: %.3f\n",
            mean(ms$ratio > 1)))


## FIGURES

save_pub <- function(p, name, w_j = 3.5, h_j = 5.0, w_p = 10, h_p = 6.5) {
  ggsave(sprintf("figures/%s.pdf", name), p,
         width = w_j, height = h_j, device = cairo_pdf, bg = "white")
  ggsave(sprintf("figures/%s.tiff", name), p,
         width = w_j, height = h_j, dpi = 600, compression = "lzw", bg = "white")
  ggsave(sprintf("figures/%s.png", name), p,
         width = w_j, height = h_j, dpi = 600, bg = "white")
}

make_sobol_df <- function(r, nm) {
  data.frame(
    parameter = nm,
    S1 = pmax(r$S1, 0),
    ST = r$ST
  ) %>%
    pivot_longer(cols = c(S1, ST), names_to = "index", values_to = "value") %>%
    mutate(index = factor(index, levels = c("S1", "ST"),
                          labels = c("First-order (S1)", "Total-order (ST)")))
}

for (o in names(outs)) {
  dfA <- make_sobol_df(resA[[o]], pn)
  dfB <- make_sobol_df(resB[[o]], unc)
  
  # Order by total-order index
  ordA <- dfA %>% filter(index == "Total-order (ST)") %>% arrange(value) %>% pull(parameter)
  ordB <- dfB %>% filter(index == "Total-order (ST)") %>% arrange(value) %>% pull(parameter)
  
  dfA$parameter <- factor(dfA$parameter, levels = ordA)
  dfB$parameter <- factor(dfB$parameter, levels = ordB)
  
  pA <- ggplot(dfA, aes(x = value, y = parameter, fill = index)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(values = c("First-order (S1)" = col_accent, "Total-order (ST)" = col_primary)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = "Sobol index", y = NULL, title = "A. All parameters") +
    theme_pub(8.5) +
    theme(legend.position = "bottom")
  
  pB <- ggplot(dfB, aes(x = value, y = parameter, fill = index)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(values = c("First-order (S1)" = col_accent, "Total-order (ST)" = col_primary)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = "Sobol index", y = NULL, title = "B. Uncertain parameters only") +
    theme_pub(8.5) +
    theme(legend.position = "bottom")
  

  library(patchwork)
  p_journal <- pA / pB +
    plot_annotation(title = outs[[o]],
                    theme = theme(plot.title = element_text(size = 10, face = "bold",
                                                            hjust = 0.5, colour = col_primary)))
  save_pub(p_journal, paste0("fig_m4_sobol_", o), w_j = 3.5, h_j = 5.4)
}


## EXCEL TABLES

wb <- createWorkbook()

## Table 1 – Baseline & local elasticities
tab1a <- data.frame(
  Quantity = c("lambda", "Generation time T (y)", "Required accuracy r", "N*"),
  Value = c(round(m0$lam, 4), round(m0$T, 3), round(m0$rraw, 3), round(10^m0$logN))
)
tab1b <- el_tab %>%
  mutate(T_fixed = round(T_fixed, 2), T_joint = round(T_joint, 2)) %>%
  rename(Parameter = parameter, `T fixed` = T_fixed, `T joint` = T_joint)

addWorksheet(wb, "Table1_Baseline_Elasticities")
writeData(wb, "Table1_Baseline_Elasticities", "Baseline results and local elasticities of N* (illustrative)", startRow = 1)
writeData(wb, "Table1_Baseline_Elasticities", tab1a, startRow = 3)
writeData(wb, "Table1_Baseline_Elasticities", "Local elasticities (T fixed vs T from demographic model)", startRow = 9)
writeData(wb, "Table1_Baseline_Elasticities", tab1b, startRow = 10)

style_header <- createStyle(textDecoration = "bold", border = "bottom", borderStyle = "medium",
                            fontSize = 11, fontColour = col_primary, halign = "center")
style_body   <- createStyle(halign = "center", fontSize = 10)
addStyle(wb, "Table1_Baseline_Elasticities", createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
addStyle(wb, "Table1_Baseline_Elasticities", style_header, rows = 3, cols = 1:2)
addStyle(wb, "Table1_Baseline_Elasticities", style_body, rows = 4:7, cols = 2)
addStyle(wb, "Table1_Baseline_Elasticities", style_header, rows = 10, cols = 1:3)
addStyle(wb, "Table1_Baseline_Elasticities", style_body, rows = 11:17, cols = 2:3, gridExpand = TRUE)
setColWidths(wb, "Table1_Baseline_Elasticities", cols = 1:3, widths = c(28, 12, 12))

## Table 2 – Top Sobol indices (example for logN, All parameters)
sobol_to_df <- function(r, nm) {
  data.frame(
    Parameter = nm,
    S1 = round(r$S1, 3),
    ST = round(r$ST, 3)
  ) %>% arrange(desc(ST))
}

addWorksheet(wb, "Table2_Sobol_logN")
writeData(wb, "Table2_Sobol_logN", "Sobol indices for log10 N* – All parameters (illustrative)", startRow = 1)
writeData(wb, "Table2_Sobol_logN", sobol_to_df(resA$logN, pn), startRow = 3)
addStyle(wb, "Table2_Sobol_logN", createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
addStyle(wb, "Table2_Sobol_logN", style_header, rows = 3, cols = 1:3)
addStyle(wb, "Table2_Sobol_logN", style_body, rows = 4:18, cols = 2:3, gridExpand = TRUE)
setColWidths(wb, "Table2_Sobol_logN", cols = 1:3, widths = c(14, 10, 10))

addWorksheet(wb, "Table3_Sobol_dG")
writeData(wb, "Table3_Sobol_dG", "Sobol indices for annual genomic gain – All parameters (illustrative)", startRow = 1)
writeData(wb, "Table3_Sobol_dG", sobol_to_df(resA$dG, pn), startRow = 3)
addStyle(wb, "Table3_Sobol_dG", createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
addStyle(wb, "Table3_Sobol_dG", style_header, rows = 3, cols = 1:3)
addStyle(wb, "Table3_Sobol_dG", style_body, rows = 4:18, cols = 2:3, gridExpand = TRUE)
setColWidths(wb, "Table3_Sobol_dG", cols = 1:3, widths = c(14, 10, 10))

saveWorkbook(wb, "tables/Camel_Module4_Results.xlsx", overwrite = TRUE)

cat("\n=== All outputs written ===\n")
cat("figures:      figures/\n")
cat("Excel tables:         tables/Camel_Module4_Results.xlsx\n")

