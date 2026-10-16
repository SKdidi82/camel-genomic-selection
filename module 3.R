
## Module 3: Optimal contributions under a coancestry constraint


## ---- Packages -------------------
pkgs <- c("ggplot2", "scales", "dplyr", "openxlsx")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(ggplot2)
library(scales)
library(dplyr)
library(openxlsx)

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

have_qp <- requireNamespace("quadprog", quietly = TRUE)
cat(if (have_qp) "Using quadprog\n" else "quadprog not found: using removal heuristic\n")

## ---- Colour palette ------------------------------

col_primary <- "#2C3E50"   # deep charcoal-navy
col_accent  <- "#E67E22"   # warm amber
col_soft    <- "#7F8C8D"   # soft supporting grey
col_bar1    <- "#2C3E50"   # OCS bars
col_bar2    <- "#E67E22"   # truncation bars

## ---- theme ----------------
theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") %+replace%
    theme(
      axis.title   = element_text(size = base_size, face = "plain", colour = col_primary),
      axis.text    = element_text(size = base_size * 0.9, colour = col_primary),
      axis.line    = element_line(linewidth = 0.45, colour = col_primary),
      axis.ticks   = element_line(linewidth = 0.35, colour = col_primary),
      plot.title   = element_text(size = base_size * 1.05, face = "bold", hjust = 0.5,
                                  colour = col_primary, margin = margin(b = 6)),
      plot.margin  = margin(6, 10, 6, 6),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text  = element_text(size = base_size * 0.85, colour = col_primary),
      legend.key.size = unit(0.45, "cm"),
      panel.grid   = element_blank()
    )
}

## ---- 1. Pedigree -> A matrix ------------------
ped  <- read.csv("pedigree_demo.csv")
n    <- nrow(ped)
sire <- ifelse(ped$sire < 0, NA_integer_, ped$sire + 1L)
dam  <- ifelse(ped$dam  < 0, NA_integer_, ped$dam  + 1L)
A <- matrix(0, n, n)
for (i in seq_len(n)) {
  s <- sire[i]; d <- dam[i]
  if (is.na(s) || is.na(d)) {
    A[i, i] <- 1
  } else {
    if (i > 1) {
      v <- 0.5 * (A[1:(i - 1), s] + A[1:(i - 1), d])
      A[i, 1:(i - 1)] <- v
      A[1:(i - 1), i] <- v
    }
    A[i, i] <- 1 + 0.5 * A[s, d]
  }
}
F_all <- diag(A) - 1
cand  <- which(ped$gen == max(ped$gen))
Ac    <- A[cand, cand]
male  <- ped$sex[cand] == "M"
g     <- 0.6 * ped$tbv[cand]
Fbar  <- mean(F_all[cand])
cat(sprintf("Candidates: %d (%d males); mean inbreeding Fbar = %.4f\n",
            length(cand), sum(male), Fbar))

K_of_dF <- function(dF) 2 * (Fbar + dF * (1 - Fbar))
dF_of_K <- function(K)  (K / 2 - Fbar) / (1 - Fbar)

## ---- 2–3. Relaxed + exact OCS  ----------------------------
relaxed_parts <- function(A, g, male) {
  Q  <- cbind(as.numeric(male), as.numeric(!male)); d <- c(0.5, 0.5)
  Ai <- solve(A); M <- solve(t(Q) %*% Ai %*% Q)
  c0 <- as.numeric(Ai %*% Q %*% M %*% d)
  Pg <- as.numeric(g - Q %*% M %*% t(Q) %*% Ai %*% g)
  v  <- as.numeric(Ai %*% Pg)
  list(c0 = c0, v = v, S = sum(g * v),
       K0 = as.numeric(t(d) %*% M %*% d), G0 = sum(c0 * g))
}

ocs_removal <- function(lam, A, g, male) {
  N <- length(g); keep <- rep(TRUE, N)
  repeat {
    idx <- which(keep)
    rp  <- relaxed_parts(A[idx, idx, drop = FALSE], g[idx], male[idx])
    cs  <- rp$c0 + rp$v / lam
    if (all(cs >= -1e-12)) break
    keep[idx[cs < 0]] <- FALSE
  }
  out <- numeric(N); out[idx] <- pmax(cs, 0); out
}
ocs_fixed <- function(lam, A, g, male) {
  if (!have_qp) return(ocs_removal(lam, A, g, male))
  N  <- length(g)
  Q  <- cbind(as.numeric(male), as.numeric(!male))
  D  <- lam * A + 1e-9 * diag(N)
  Am <- cbind(Q, diag(N)); b0 <- c(0.5, 0.5, rep(0, N))
  sol <- quadprog::solve.QP(D, g, Am, b0, meq = 2)$solution
  sol[sol < 1e-9] <- 0
  sol
}
KG <- function(cc) c(K = sum(cc * (Ac %*% cc)), G = sum(cc * g))

## ---- 4. Closed-form check  -------------------
set.seed(1)
gu <- rnorm(120); mu <- rep(c(TRUE, FALSE), each = 60); Iu <- diag(120)
rp <- relaxed_parts(Iu, gu, mu)
lam_crit <- max(2 * 60 * (mean(gu[mu]) - min(gu[mu])),
                2 * 60 * (mean(gu[!mu]) - min(gu[!mu])))
for (f in c(3, 1.2)) {
  lam <- f * lam_crit; cc <- ocs_fixed(lam, Iu, gu, mu)
  stopifnot(abs(sum(cc * gu) - (rp$G0 + rp$S / lam)) < 1e-3)
}
cat(sprintf("A = I check passed (lambda_crit = %.1f)\n", lam_crit))

## ---- 5. Frontier  ------
lams <- exp(seq(log(0.3), log(3000), length.out = 60))
fr <- as.data.frame(t(sapply(lams, function(l) {
  cc <- ocs_fixed(l, Ac, g, male)
  c(lam = l, K = sum(cc * (Ac %*% cc)), G = sum(cc * g), nz = sum(cc > 1e-6))
})))
fr$dF <- dF_of_K(fr$K)
rl <- relaxed_parts(Ac, g, male)
cat(sprintf("Relaxed: S = %.2f, K0 = %.4f, G0 = %.4f\n", rl$S, rl$K0, rl$G0))

ocs_at_dF <- function(target) {
  Kt <- K_of_dF(target)
  f  <- function(ll) {
    cc <- ocs_fixed(exp(ll), Ac, g, male)
    sum(cc * (Ac %*% cc)) - Kt
  }
  ll <- uniroot(f, c(log(0.3), log(3000)), tol = 1e-10)$root
  cc <- ocs_fixed(exp(ll), Ac, g, male)
  list(lambda0 = exp(ll), c = cc, K = KG(cc)[["K"]], G = KG(cc)[["G"]])
}
trunc_sel <- function(nm, nf) {
  stopifnot(nm <= sum(male),nf <= sum(!male)) #NEW
  cc <- numeric(length(g))
  im <- which(male)[order(-g[male])[1:nm]]
  jf <- which(!male)[order(-g[!male])[1:nf]]
  cc[im] <- 0.5 / nm; cc[jf] <- 0.5 / nf; cc
}

cat("\nOptimal contributions at fixed ΔF:\n")
res <- do.call(rbind, lapply(c(0.005, 0.01, 0.02), function(t) {
  o <- ocs_at_dF(t); K <- K_of_dF(t)
  data.frame(dF_pct = 100 * t, gain = o$G,
             relaxed_bound = rl$G0 + sqrt(rl$S * (K - rl$K0)),
             lambda0 = o$lambda0, shadow_price = o$lambda0 * (1 - Fbar),
             males_used = sum(o$c[male] > 1e-6), females_used = sum(o$c[!male] > 1e-6))
}))
print(round(res, 3), row.names = FALSE)

cat("\nTruncation selection comparison:\n")
tr <- do.call(rbind, lapply(list(c(3, 15), c(10, 40), c(30, 40)), function(p) {
  cc <- trunc_sel(p[1], p[2]); kg <- KG(cc)
  Gocs <- approx(fr$K[order(fr$K)], fr$G[order(fr$K)], xout = kg[["K"]])$y
  data.frame(males = p[1], females = p[2], dF_pct = 100 * dF_of_K(kg[["K"]]),
             gain = kg[["G"]], OCS_gain_same_K = Gocs)
}))
print(round(tr, 3), row.names = FALSE)

## ---- 6. Annual view ----------------------------------
L_gen <- c(6, 10, 14)
cat("\nPer-generation cap equivalent to 0.1% per year:\n")
print(data.frame(L_years = L_gen, dF_gen_pct = 0.1 * L_gen))


## FIGURES

save_pub <- function(p, name, w_j = 3.5, h_j = 2.9, w_p = 9, h_p = 6.3) {
  ggsave(sprintf("figures/%s.pdf", name), p,
         width = w_j, height = h_j, device = cairo_pdf, bg = "white")
  ggsave(sprintf("figures/%s.tiff", name), p,
         width = w_j, height = h_j, dpi = 600, compression = "lzw", bg = "white")
  ggsave(sprintf("figures/%s.png", name), p,
         width = w_j, height = h_j, dpi = 600, bg = "white")
}

## ---- Fig 1: Gain–inbreeding frontier -----------------------
dFg <- seq(max(dF_of_K(rl$K0) + 1e-4, -0.006), 0.06, length.out = 200)
Gr  <- rl$G0 + sqrt(rl$S * (K_of_dF(dFg) - rl$K0))

df_front <- data.frame(dF = 100 * fr$dF, G = fr$G)
df_bound <- data.frame(dF = 100 * dFg, G = Gr)
df_trunc <- do.call(rbind, lapply(list(c(3,15), c(5,20), c(10,40), c(30,40)), function(p) {
  kg <- KG(trunc_sel(p[1], p[2]))
  data.frame(dF = 100 * dF_of_K(kg[["K"]]), G = kg[["G"]])
}))

p1 <- ggplot() +
  geom_line(data = df_front, aes(x = dF, y = G), linewidth = 0.95, colour = col_primary) +
  geom_line(data = df_bound, aes(x = dF, y = G), linewidth = 0.7, linetype = "dashed", colour = col_soft) +
  geom_point(data = df_trunc, aes(x = dF, y = G), size = 2.2, colour = col_accent) +
  geom_vline(xintercept = c(0.5, 1, 2), linetype = "dotted", colour = "grey70", linewidth = 0.5) +
  labs(x = "Rate of inbreeding per generation (%)",
       y = expression(paste("Gain (", sigma[A], " × accuracy)")),
       title = "Gain–inbreeding frontier") +
  annotate("text", x = 5.5, y = max(df_front$G) * 0.22, label = "OCS (exact)", 
           size = 2.7, colour = col_primary, hjust = 1) +
  annotate("text", x = 5.5, y = max(df_front$G) * 0.12, label = "Relaxed bound", 
           size = 2.7, colour = col_soft, hjust = 1) +
  annotate("text", x = 5.5, y = max(df_front$G) * 0.02, label = "Truncation", 
           size = 2.7, colour = col_accent, hjust = 1) +
  scale_x_continuous(limits = c(-0.7, 6)) +
  theme_pub(9)
save_pub(p1, "fig_m3_frontier", h_j = 3.0)

## ---- Fig 2: Top contributions -------------------------------
o1 <- ocs_at_dF(0.01)$c
t1 <- trunc_sel(30, 40)
df_c <- data.frame(
  rank = rep(1:45, 2),
  contribution = c(sort(o1, decreasing = TRUE)[1:45], sort(t1, decreasing = TRUE)[1:45]),
  method = rep(c("OCS", "Truncation (30♂, 40♀)"), each = 45)
)

p2 <- ggplot(df_c, aes(x = rank, y = contribution, fill = method)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  scale_fill_manual(values = c("OCS" = col_bar1, "Truncation (30♂, 60♀)" = col_bar2)) +
  labs(x = "Rank (top 45 animals)", y = "Genetic contribution",
       title = "Top 45 contributions (~1% ΔF per generation)") +
  theme_pub(9) +
  theme(legend.position = "bottom")
save_pub(p2, "fig_m3_contributions", h_j = 3.0)

## ---- Fig 3: Shadow price ---------------------------------------------
df_sh <- data.frame(dF = 100 * fr$dF, shadow = fr$lam * (1 - Fbar))

p3 <- ggplot(df_sh, aes(x = dF, y = shadow)) +
  geom_line(linewidth = 0.95, colour = col_primary) +
  scale_y_log10() +
  labs(x = "Rate of inbreeding per generation (%)",
       y = expression(paste("Shadow price  ", dG/d(dF))),
       title = "Marginal gain per unit of inbreeding") +
  theme_pub(9)
save_pub(p3, "fig_m3_shadow_price", h_j = 2.9)


## EXCEL TABLES

wb <- createWorkbook()

## Table 1 – OCS at fixed ΔF
tab1 <- res %>%
  transmute(
    `ΔF (% per generation)` = dF_pct,
    `Genetic gain` = round(gain, 3),
    `Relaxed bound` = round(relaxed_bound, 3),
    `λ₀` = round(lambda0, 1),
    `Shadow price` = round(shadow_price, 1),
    `Males used` = males_used,
    `Females used` = females_used
  )

addWorksheet(wb, "Table1_OCS_fixed_dF")
writeData(wb, "Table1_OCS_fixed_dF", "Optimal contributions at fixed rate of inbreeding (illustrative)", startRow = 1)
writeData(wb, "Table1_OCS_fixed_dF", tab1, startRow = 3)
style_header <- createStyle(textDecoration = "bold", border = "bottom", borderStyle = "medium",
                            fontSize = 11, fontColour = col_primary, halign = "center")
style_body   <- createStyle(halign = "center", fontSize = 10)
addStyle(wb, "Table1_OCS_fixed_dF", createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
addStyle(wb, "Table1_OCS_fixed_dF", style_header, rows = 3, cols = 1:7)
addStyle(wb, "Table1_OCS_fixed_dF", style_body, rows = 4:6, cols = 1:7, gridExpand = TRUE)
setColWidths(wb, "Table1_OCS_fixed_dF", cols = 1:7, widths = c(18, 12, 14, 10, 12, 11, 12))

## Table 2 – Truncation comparison
tab2 <- tr %>%
  transmute(
    Males = males,
    Females = females,
    `ΔF (% per gen.)` = round(dF_pct, 2),
    `Truncation gain` = round(gain, 3),
    `OCS gain (same K)` = round(OCS_gain_same_K, 3)
  )

addWorksheet(wb, "Table2_Truncation")
writeData(wb, "Table2_Truncation", "Truncation selection vs OCS at the same coancestry (illustrative)", startRow = 1)
writeData(wb, "Table2_Truncation", tab2, startRow = 3)
addStyle(wb, "Table2_Truncation", createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
addStyle(wb, "Table2_Truncation", style_header, rows = 3, cols = 1:5)
addStyle(wb, "Table2_Truncation", style_body, rows = 4:6, cols = 1:5, gridExpand = TRUE)
setColWidths(wb, "Table2_Truncation", cols = 1:5, widths = c(10, 10, 14, 14, 16))

saveWorkbook(wb, "tables/Camel_Module3_Results.xlsx", overwrite = TRUE)

cat("\n=== All outputs written ===\n")
cat("figures:      figures/\n")
cat("Excel tables:         tables/Camel_Module3_Results.xlsx\n")

