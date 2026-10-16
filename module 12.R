
## Dairy-camel breeding theory – Modules 1 & 2


## ---- Required packages ------
pkgs <- c("ggplot2", "scales", "dplyr", "openxlsx")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(ggplot2)
library(scales)
library(dplyr)
library(openxlsx)

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

## ---- Colour palette ------------------

col_primary  <- "#2C3E50"   
col_accent   <- "#E67E22"  
col_soft     <- "#7F8C8D"   
col_bar      <- "#34495E"   
col_grid     <- "#ECF0F1"   

## ---- theme ---------------------
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


## MODULE 1 – Age-structured female population 

make_life <- function(alpha, omega, s_juv, s_ad, I, sex_ratio = 0.5) {
  a   <- 0:omega
  l   <- c(1, s_juv * s_ad^(0:(omega - 1)))
  m   <- ifelse(a >= alpha, sex_ratio / I, 0)
  list(a = a, l = l, m = m, alpha = alpha, omega = omega,
       s = s_ad, j = s_juv, mad = sex_ratio / I)
}

euler_lotka <- function(life) {
  a <- life$a[-1]; lm <- (life$l * life$m)[-1]
  f <- function(L) sum(L^(-a) * lm) - 1
  uniroot(f, c(0.85, 3), tol = 1e-12)$root
}

gen_time <- function(life, lambda) {
  a <- life$a[-1]; lm <- (life$l * life$m)[-1]
  sum(a * lambda^(-a) * lm)
}

leslie_lambda <- function(life) {
  w <- life$omega; n <- w + 1
  s <- life$l[-1] / life$l[-n]
  L <- matrix(0, n, n)
  L[1, 1:w] <- s * life$m[-1]
  L[cbind(2:n, 1:w)] <- s
  ev <- eigen(L)$values
  max(Re(ev[abs(Im(ev)) < 1e-9]))
}

fec_elasticity <- function(life, lambda, T) {
  a <- life$a[-1]; lm <- (life$l * life$m)[-1]
  setNames(lambda^(-a) * lm / T, a)
}

lambda_closed <- function(alpha, s, j, mad)
  uniroot(function(L) L^alpha - s * L^(alpha - 1) - mad * j * s^(alpha - 1),
          c(s + 1e-9, 3), tol = 1e-13)$root
T_closed <- function(alpha, s, lambda) alpha + s / (lambda - s)
dlam_dalpha <- function(alpha, s, lambda)
  -lambda * (lambda - s) * log(lambda / s) / (alpha * (lambda - s) + s)

scen <- list(
  camel_slow  = list(alpha = 5, omega = 25, s_juv = 0.80, s_ad = 0.93, I = 2.0),
  camel_fast  = list(alpha = 4, omega = 25, s_juv = 0.85, s_ad = 0.93, I = 1.5),
  cattle_like = list(alpha = 2, omega = 12, s_juv = 0.90, s_ad = 0.82, I = 1.2)
)
res1 <- do.call(rbind, lapply(names(scen), function(nm) {
  p  <- scen[[nm]]
  lf <- do.call(make_life, p)
  lam <- euler_lotka(lf)
  data.frame(scenario = nm, R0 = sum((lf$l * lf$m)[-1]),
             lambda_EL = lam, lambda_leslie = leslie_lambda(lf),
             T_gen = gen_time(lf, lam))
}))
cat("\n== Module 1 results ==\n")
print(res1, digits = 5, row.names = FALSE)
stopifnot(all(abs(res1$lambda_EL - res1$lambda_leslie) < 1e-8))

p  <- scen$camel_slow
lf <- make_life(p$alpha, 400, p$s_juv, p$s_ad, p$I)
lamN <- euler_lotka(lf); lamC <- lambda_closed(p$alpha, p$s_ad, p$s_juv, 0.5 / p$I)
cat(sprintf("closed-form lambda = %.6f, numeric = %.6f\n", lamC, lamN))
cat(sprintf("closed-form T      = %.4f, numeric = %.4f\n",
            T_closed(p$alpha, p$s_ad, lamC), gen_time(lf, lamN)))

fd <- (lambda_closed(p$alpha + 1e-5, p$s_ad, p$s_juv, 0.5 / p$I) -
         lambda_closed(p$alpha - 1e-5, p$s_ad, p$s_juv, 0.5 / p$I)) / 2e-5
cat(sprintf("d lambda/d alpha: formula %.6f, finite difference %.6f\n",
            dlam_dalpha(p$alpha, p$s_ad, lamC), fd))

tab_alpha <- as.data.frame(t(sapply(3:7, function(al) {
  L <- lambda_closed(al, p$s_ad, p$s_juv, 0.5 / p$I)
  c(alpha = al, lambda = L, T = T_closed(al, p$s_ad, L))
})))
cat("\nEffect of age at first calving:\n")
print(round(tab_alpha, 4))


## MODULE 2 – Genomic selection vs progeny testing 

Me_fun   <- function(Ne, G) 2 * Ne * G / log(4 * Ne * G)
r_gen    <- function(N, H, Me) sqrt(N * H / (N * H + Me))
N_star   <- function(r, H, Me) Me * r^2 / (H * (1 - r^2))
k_fun    <- function(h2) (4 - h2) / h2
r_pt_fun <- function(n, h2) sqrt(n / (n + k_fun(h2)))
R2_bull  <- function(n, h2) n / (n + k_fun(h2))
H_fem    <- function(h2, rep, nrec) nrec * h2 / (1 + (nrec - 1) * rep)
i_trunc  <- function(p) dnorm(qnorm(1 - p)) / p

p0 <- list(h2 = 0.25, rep = 0.5, nrec = 1, Ne = 100, G = 30,
           a_b = 5, tau = 1.08, alpha_f = 5, ell = 1, delta = 0.5,
           n_dau = 20, p_sire_pt = 0.10, p_sire_gs = 0.10, p_dam = 0.5,
           T = T_closed(scen$camel_slow$alpha, scen$camel_slow$s_ad, lamC))   # closed form, no age limit

derive <- function(p) with(p, {
  L_pt <- a_b + 2 * tau + alpha_f + ell + delta
  L_g  <- a_b + tau
  Me   <- Me_fun(Ne, G); H <- H_fem(h2, rep, nrec)
  r_pt <- r_pt_fun(n_dau, h2)
  i_pt <- i_trunc(p_sire_pt); i_g <- i_trunc(p_sire_gs); i_d <- i_trunc(p_dam)
  h    <- sqrt(h2)
  r_path <- (i_pt / i_g) * r_pt * L_g / L_pt
  r_four <- ((i_pt * r_pt + i_d * h) * (L_g + T) / (L_pt + T) - i_d * h) / i_g
  list(L_pt = L_pt, L_g = L_g, Me = Me, H = H, r_pt = r_pt, i_pt = i_pt,
       i_g = i_g, i_d = i_d, h = h, r_path = r_path, r_four = r_four, T = T)
})
N_req <- function(r, H, Me) ifelse(r <= 0, 0, ifelse(r >= 1, Inf, N_star(r, H, Me)))

d0 <- derive(p0)
cat("\n== Module 2 breakeven results ==\n")
cat(sprintf("Me = %.0f | L_pt = %.2f y | L_g = %.2f y | r_pt = %.3f\n",
            d0$Me, d0$L_pt, d0$L_g, d0$r_pt))
cat(sprintf("Required accuracy: path-level %.3f | four-path %.3f\n", d0$r_path, d0$r_four))
cat(sprintf("N* path-level %.0f | four-path %.0f\n",
            N_req(d0$r_path, d0$H, d0$Me), N_req(d0$r_four, d0$H, d0$Me)))

dG_pt <- function(d) (2 * d$i_pt * d$r_pt + 2 * d$i_d * d$h) / (2 * d$L_pt + 2 * d$T)
dG_gs <- function(N, d) (2 * d$i_g * r_gen(N, d$H, d$Me) + 2 * d$i_d * d$h) /
  (2 * d$L_g + 2 * d$T)
Nb <- uniroot(function(N) dG_gs(N, d0) - dG_pt(d0), c(1, 1e7), tol = 1e-8)$root
cat(sprintf("Numerical breakeven N = %.1f\n", Nb))
stopifnot(abs(Nb - N_req(d0$r_four, d0$H, d0$Me)) < 1e-2)
p25 <- p0; p25$T <- res1$T_gen[res1$scenario == "camel_slow"]   # life table to 25 y
d25 <- derive(p25)
cat(sprintf("N* with 25-year life table (T = %.2f): %.0f\n", p25$T, N_req(d25$r_four, d25$H, d25$Me)))
cat("\nEffect of repeated records on N*:\n")
tab_rec <- sapply(1:4, function(k) {
  p <- p0; p$nrec <- k; d <- derive(p); N_req(d$r_four, d$H, d$Me)
})
print(setNames(round(tab_rec), paste0("records=", 1:4)))

Nfun <- function(p) { d <- derive(p); N_req(d$r_four, d$H, d$Me) }
elas <- sapply(c("h2", "Ne", "G", "n_dau", "a_b", "tau", "alpha_f", "ell", "delta", "T"),
               function(nm) {
                 up <- p0; dn <- p0; up[[nm]] <- p0[[nm]] * 1.01; dn[[nm]] <- p0[[nm]] * 0.99
                 (log(Nfun(up)) - log(Nfun(dn))) / (log(1.01) - log(0.99))
               })
cat("\nElasticities of N*:\n"); print(round(elas, 2))


## FIGURES

save_pub <- function(p, name, w_j = 3.5, h_j = 2.8, w_p = 9, h_p = 6.2) {
  ggsave(sprintf("figures/%s.pdf", name), p,
         width = w_j, height = h_j, device = cairo_pdf, bg = "white")
  ggsave(sprintf("figures/%s.tiff", name), p,
         width = w_j, height = h_j, dpi = 600, compression = "lzw", bg = "white")
  ggsave(sprintf("figures/%s.png", name), p,
         width = w_j, height = h_j, dpi = 600, bg = "white")
}

## Fig 1 – Fecundity elasticity
lf0 <- do.call(make_life, scen$camel_slow)
lam0 <- euler_lotka(lf0)
el  <- fec_elasticity(lf0, lam0, gen_time(lf0, lam0))
df1 <- data.frame(age = as.numeric(names(el)), elasticity = as.numeric(el))

p1 <- ggplot(df1, aes(x = age, y = elasticity)) +
  geom_col(fill = col_bar, width = 0.72) +
  labs(x = "Age a (years)",
       y = expression(paste("Elasticity of ", lambda, " to ", m[a])),
       title = "Fecundity elasticity (camel_slow, assumed)") +
  scale_x_continuous(breaks = seq(0, 25, 5)) +
  theme_pub(9)
save_pub(p1, "fig1_fecundity_elasticity", h_j = 2.7)

## Fig 2 – Genomic accuracy vs N
Ns <- seq(10, 3000, by = 10)
df2 <- data.frame(N = Ns, r = r_gen(Ns, d0$H, d0$Me))

p2 <- ggplot(df2, aes(x = N, y = r)) +
  geom_line(linewidth = 0.95, colour = col_primary) +
  geom_hline(yintercept = d0$r_path, linetype = "dashed", colour = col_soft, linewidth = 0.65) +
  geom_hline(yintercept = d0$r_four, linetype = "dashed", colour = col_accent, linewidth = 0.65) +
  labs(x = "Reference females (N)", y = "Expected accuracy r",
       title = "Genomic accuracy vs reference size") +
  annotate("text", x = 2900, y = d0$r_path + 0.045, label = "path-level r*",
           size = 2.7, colour = col_soft, hjust = 1) +
  annotate("text", x = 2900, y = d0$r_four + 0.045, label = "four-path r*",
           size = 2.7, colour = col_accent, hjust = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_pub(9)
save_pub(p2, "fig2_accuracy_vs_N")

## Fig 3 – Annual gain vs N
df3 <- data.frame(N = Ns, dG = sapply(Ns, dG_gs, d = d0))

p3 <- ggplot(df3, aes(x = N, y = dG)) +
  geom_line(linewidth = 0.95, colour = col_primary) +
  geom_hline(yintercept = dG_pt(d0), linetype = "dashed", colour = col_accent, linewidth = 0.65) +
  geom_vline(xintercept = Nb, linetype = "dotted", colour = col_soft, linewidth = 0.65) +
  labs(x = "Reference females (N)",
       y = expression(paste(Delta, "G per year (", sigma[A], ")")),
       title = "Annual gain: genomic vs progeny test") +
  annotate("text", x = 2900, y = dG_pt(d0) + 0.0045, label = "progeny test",
           size = 2.7, colour = col_accent, hjust = 1) +
  theme_pub(9)
save_pub(p3, "fig3_gain_vs_N")
cat(sprintf("Peak fecundity elasticity: %.4f at age %d\n",
            max(df1$elasticity), df1$age[which.max(df1$elasticity)]))
cat(sprintf("Elasticity at age 25: %.4f\n", df1$elasticity[df1$age == 25]))
cat(sprintf("At N = 3000: accuracy %.3f, genomic gain %.4f vs progeny test %.4f\n",
            r_gen(3000, d0$H, d0$Me), dG_gs(3000, d0), dG_pt(d0)))
## Fig 4 – N* heatmap
h2s <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.40)
Nes <- c(50, 100, 200, 500, 1000)
grid <- outer(h2s, Nes, Vectorize(function(h, ne) {
  p <- p0; p$h2 <- h; p$Ne <- ne; Nfun(p)
}))
df4 <- expand.grid(h2 = h2s, Ne = Nes)
df4$Nstar <- as.vector(grid)
df4$logN  <- log10(df4$Nstar)

p4 <- ggplot(df4, aes(x = factor(h2), y = factor(Ne), fill = logN)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = round(Nstar)), size = 2.5, colour = "white") +
  scale_fill_gradient(low = "#6C859E", high = col_primary, name = expression(log[10]~N^"*")) +
  labs(x = expression(h^2), y = expression(N[e]),
       title = expression(log[10]~of~required~reference~size~N^"*")) +
  theme_pub(9) +
  theme(legend.position = "right", legend.key.height = unit(0.85, "cm"))
save_pub(p4, "fig4_Nstar_grid", w_j = 3.7, h_j = 3.1)

## Fig 5 – Elasticities
df5 <- data.frame(parameter = names(elas), elasticity = as.numeric(elas))
df5 <- df5[order(df5$elasticity), ]
df5$parameter <- factor(df5$parameter, levels = df5$parameter)

p5 <- ggplot(df5, aes(x = elasticity, y = parameter)) +
  geom_col(fill = col_bar, width = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = col_primary) +
  labs(x = expression(Elasticity~of~N^"*"), y = NULL,
       title = "Drivers of required reference size") +
  theme_pub(9)
save_pub(p5, "fig5_elasticities", h_j = 3.3)


## EXCEL TABLES 

wb <- createWorkbook()

## ---- Table 1: Module 1 scenarios ------------------------
tab1 <- res1 %>%
  transmute(
    Scenario = c("Camel (slow)", "Camel (fast)", "Cattle-like"),
    R0 = round(R0, 3),
    Lambda = round(lambda_EL, 4),
    `Generation time (y)` = round(T_gen, 2)
  )

addWorksheet(wb, "Table1_Module1")
writeData(wb, "Table1_Module1", tab1, startRow = 2)
style_header <- createStyle(textDecoration = "bold", border = "bottom",
                            borderStyle = "medium", fontSize = 11, fontColour = col_primary)
style_body   <- createStyle(halign = "center", fontSize = 10)
addStyle(wb, "Table1_Module1", style_header, rows = 2, cols = 1:4, gridExpand = TRUE)
addStyle(wb, "Table1_Module1", style_body, rows = 3:5, cols = 2:4, gridExpand = TRUE)
setColWidths(wb, "Table1_Module1", cols = 1:4, widths = c(16, 10, 10, 18))
writeData(wb, "Table1_Module1", "Module 1 – Population growth parameters (illustrative assumptions)",
          startRow = 1, startCol = 1)
addStyle(wb, "Table1_Module1", createStyle(textDecoration = "bold", fontSize = 12),
         rows = 1, cols = 1)

## ---- Table 2: Module 2 key results ---------------
tab2 <- data.frame(
  Quantity = c("Me (effective loci)", "L_pt (progeny-test interval, y)",
               "L_g (genomic interval, y)", "r_pt (progeny-test accuracy)",
               "Required r – path-level", "Required r – four-path",
               "N* path-level", "N* four-path", "Breakeven N"),
  Value = c(round(d0$Me), round(d0$L_pt, 2), round(d0$L_g, 2),
            round(d0$r_pt, 3), round(d0$r_path, 3), round(d0$r_four, 3),
            round(N_req(d0$r_path, d0$H, d0$Me)),
            round(N_req(d0$r_four, d0$H, d0$Me)), round(Nb, 1))
)

addWorksheet(wb, "Table2_Breakeven")
writeData(wb, "Table2_Breakeven", tab2, startRow = 2)
addStyle(wb, "Table2_Breakeven", style_header, rows = 2, cols = 1:2)
addStyle(wb, "Table2_Breakeven", style_body, rows = 3:11, cols = 2)
setColWidths(wb, "Table2_Breakeven", cols = 1:2, widths = c(32, 14))
writeData(wb, "Table2_Breakeven", "Module 2 – Genomic selection breakeven quantities (illustrative)",
          startRow = 1)
addStyle(wb, "Table2_Breakeven", createStyle(textDecoration = "bold", fontSize = 12),
         rows = 1, cols = 1)

## ---- Table 3: Elasticities ------------
tab3 <- data.frame(
  Parameter = names(elas),
  Elasticity = round(as.numeric(elas), 2)
) %>% arrange(desc(abs(Elasticity)))

addWorksheet(wb, "Table3_Elasticities")
writeData(wb, "Table3_Elasticities", tab3, startRow = 2)
addStyle(wb, "Table3_Elasticities", style_header, rows = 2, cols = 1:2)
addStyle(wb, "Table3_Elasticities", style_body, rows = 3:12, cols = 2)
setColWidths(wb, "Table3_Elasticities", cols = 1:2, widths = c(14, 12))
writeData(wb, "Table3_Elasticities", "Local elasticities of N* (four-path)", startRow = 1)
addStyle(wb, "Table3_Elasticities", createStyle(textDecoration = "bold", fontSize = 12),
         rows = 1, cols = 1)

saveWorkbook(wb, "tables/Camel_Modules12_Results.xlsx", overwrite = TRUE)

cat("\n=== All outputs written ===\n")
cat("figures:      figures/\n")
cat("Excel tables:         tables/Camel_Modules12_Results.xlsx\n")

