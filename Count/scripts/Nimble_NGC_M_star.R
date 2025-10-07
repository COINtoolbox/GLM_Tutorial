# run_nb_eiv_mstar_nimble.R
# Negative-binomial regression for GC counts with errors-in-variables (NIMBLE)
# Predictor: x = log10(M*/Msun) derived from MV and (B-V) via Bell+03 (Chabrier IMF)
# Output: NGC_Mstar_nimble.pdf

# ------------------------- Setup --------------------------------------------
set.seed(42)
suppressPackageStartupMessages({
  library(nimble)
  library(tidyverse)
  library(Cairo)
  library(latex2exp)
  library(coda)
})

plog <- scales::pseudo_log_trans(base = 10, sigma = 5)

# ------------------------- Data ---------------------------------------------
df <- read.csv(
  "https://raw.githubusercontent.com/COINtoolbox/NB_GCs/refs/heads/master/Dataset/GCs_full.csv",
  header = TRUE, dec = ".", sep = ""
) %>% filter(!is.na(MV_T))

# Bell+03 (Chabrier) M*/L_V:
#   log10(M*/L_V) = -0.778 + 1.305 * (B-V)_0
#   M_V,0 = M_V - A_V ; M_V,Sun = 4.83
M_Vsun <- 4.83

# Typical intrinsic colours by type (used where individual (B-V)_0 is absent)
BV_typical <- c(E = 0.96, S0 = 0.88, S = 0.70, Irr = 0.45)
# Conservative per-type colour uncertainty (tweak if you have better)
BV_sigma   <- c(E = 0.05, S0 = 0.06, S = 0.08, Irr = 0.10)

df$BV  <- BV_typical[df$Type]
df$sBV <- BV_sigma[df$Type]
df$MV0 <- df$MV_T - df$A_V

# Observed x = log10(M*/Msun)
df$logMstar_obs <- -0.4*(df$MV0 - M_Vsun) - 0.778 + 1.305*df$BV

# Propagate uncertainties (first-order):
# d logM*/d M_V = -0.4 ; d logM*/d(B-V) = 1.305
df$err_logMstar <- sqrt( (0.4 * df$err_MV_T)^2 + (1.305 * df$sBV)^2 )
df$err_logMstar[!is.finite(df$err_logMstar)] <- NA
df$err_logMstar <- ifelse(is.na(df$err_logMstar) | df$err_logMstar <= 0, 0.01, df$err_logMstar)

# Count uncertainties for plotting / measurement layer; floor at 0.5 (rounding)
df$N_GC_err <- pmax(df$N_GC_err, 0.5)

# Keep essentials and ensure integer counts
df_fit <- df %>%
  filter(is.finite(logMstar_obs), is.finite(err_logMstar), !is.na(N_GC)) %>%
  mutate(N_GC = pmax(0L, as.integer(round(N_GC))))

# Prediction grid in x = log10(M*/Msun)
Mgrid  <- 2000L
x_pred <- seq(6.5, 12.5, length.out = Mgrid)

# ------------------------- NIMBLE model -------------------------------------
# NB2 parameterization: mean = mu, var = mu + mu^2/phi
# R/nimble's dnbinom(size, prob): mean = size*(1-p)/p
# Convert NB2(mu,phi) -> size = phi, p = phi/(phi + mu)

code <- nimbleCode({
  # Priors
  beta0 ~ dnorm(0, sd = 10)
  beta1 ~ dnorm(0, sd = 10)
  phi   ~ dgamma(shape = 2, rate = 0.1)   # dispersion (k); larger -> closer to Poisson

  # Errors-in-variables on predictor and two-layer count model
  for(i in 1:N){
    x_true[i] ~ dnorm(mean = x_obs[i], sd = x_err[i])   # EIV (symmetric in Normal)

    eta[i] <- beta0 + beta1 * x_true[i]
    mu[i]  <- exp(eta[i])

    p[i] <- phi / (phi + mu[i])
    N_true[i] ~ dnbinom(size = phi, prob = p[i])        # latent integer count

    # Measurement error on observed counts (Normal around integer N_true)
    N_GC[i] ~ dnorm(mean = N_true[i], sd = N_GC_err[i])
  }

  # Predictions on grid
  for (j in 1:M){
    eta_pred[j] <- beta0 + beta1 * x_pred[j]
    mu_pred[j]  <- exp(eta_pred[j])
    p_pred[j]   <- phi / (phi + mu_pred[j])
    N_pred[j]   ~ dnbinom(size = phi, prob = p_pred[j])
  }
})

constants <- list(
  N = nrow(df_fit),
  M = length(x_pred),
  x_pred = x_pred
)

data_list <- list(
  x_obs    = df_fit$logMstar_obs,
  x_err    = df_fit$err_logMstar,
  N_GC     = as.numeric(df_fit$N_GC),   # Normal observation can be real
  N_GC_err = df_fit$N_GC_err
)

inits <- list(
  beta0 = 0,
  beta1 = 1,
  phi   = 5,
  x_true = df_fit$logMstar_obs,
  N_true = pmax(0, as.integer(round(df_fit$N_GC)))
)

model <- nimbleModel(code, data = data_list, constants = constants, inits = inits, check = FALSE)
cModel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c("beta0","beta1","phi","mu_pred","N_pred"))
mcmc  <- buildMCMC(conf)
cMcmc <- compileNimble(mcmc, project = model)

samples <- runMCMC(cMcmc,
                   niter = 30000, nburnin = 10000, thin = 5,
                   nchains = 3, samplesAsCodaMCMC = TRUE, summary = FALSE)

# ------------------------- Summaries on mu_pred ------------------------------
# Combine chains to a matrix
mat_list <- lapply(samples, as.matrix)
S <- do.call(rbind, mat_list)

# Extract columns for mu_pred[j]
mu_cols <- grep("^mu_pred\\[", colnames(S))
mu_mat  <- S[, mu_cols, drop = FALSE]

# Indices j from column names like "mu_pred[123]"
j_idx <- as.integer(gsub("mu_pred\\[|\\]", "", colnames(mu_mat)))
ord   <- order(j_idx)
mu_mat <- mu_mat[, ord, drop = FALSE]
j_idx  <- j_idx[ord]

summ_mu <- tibble(
  x     = x_pred[j_idx],
  lwr95 = apply(mu_mat, 2, quantile, probs = 0.05),
  lwr50 = apply(mu_mat, 2, quantile, probs = 0.25),
  median= apply(mu_mat, 2, quantile, probs = 0.50),
  upr50 = apply(mu_mat, 2, quantile, probs = 0.75),
  upr95 = apply(mu_mat, 2, quantile, probs = 0.95)
)

# ------------------------- Plot ----------------------------------------------
cairo_pdf("NGC_Mstar_nimble.pdf", width = 7, height = 5)
ggplot() +
  geom_point(
    data = df_fit,
    aes(x = logMstar_obs, y = N_GC),
    size = 1, color = "gray70", alpha = 0.8
  ) +
  geom_errorbar(
    data = df_fit,
    aes(x = logMstar_obs,
        ymin = pmax(N_GC - N_GC_err, 0), ymax = N_GC + N_GC_err),
    width = 0.02, color = "gray80", alpha = 0.7
  ) +
  geom_errorbarh(
    data = df_fit,
    aes(y = N_GC,
        xmin = logMstar_obs - err_logMstar,
        xmax = logMstar_obs + err_logMstar),
    height = 0.03, color = "gray80", alpha = 0.7
  ) +
  geom_ribbon(
    data = summ_mu,
    aes(x = x, ymin = lwr95, ymax = upr95),
    fill = "#FFA630", alpha = 0.9
  ) +
  geom_step(
    data = summ_mu,
    aes(x = x, y = median),
    color = "#0474BA", linewidth = 1
  ) +
  scale_y_continuous(
    trans = plog,
    breaks = c(0, 10, 100, 1e3, 1e4, 1e5),
    labels = c("0", expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5))
  ) +
  labs(x = expression(log[10]~(M['\u2605']/M['☉'])), y = expression(N[GC]))  +
  theme_bw(base_size = 22) +
  theme(legend.position = "none")
dev.off()


