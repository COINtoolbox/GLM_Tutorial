# ------------------------------------------------------------
# Required Libraries
# ------------------------------------------------------------
require(R2jags)
library(rjags)
#load.module("extra")
library(ggplot2)
library(ggthemes)
library(Cairo)
library(ggmcmc)
library(plyr)
library(latex2exp)
source("jagsresults.R")
plog <- scales::pseudo_log_trans(base = 10, sigma = 5)

# ------------------------------------------------------------
# Read Data
# ------------------------------------------------------------
df <- read.csv("https://raw.githubusercontent.com/COINtoolbox/NB_GCs/refs/heads/master/Dataset/GCs_full.csv",
               header = TRUE, dec = ".", sep = "")
df <- subset(df, !is.na(MV_T))
M <- 2000

# ---------- Conversão para log M_star ----------
# Cores típicas por morfologia (pode ajustar se tiver algo melhor)
BV_typical <- c(E = 0.96, S0 = 0.88, S = 0.70, Irr = 0.45)
# Incerteza típica da cor por tipo (hiperparâmetro simples)
BV_sigma   <- c(E = 0.05, S0 = 0.06, S = 0.08, Irr = 0.10)

M_Vsun <- 4.83
df$BV     <- BV_typical[df$Type]
df$sBV    <- BV_sigma[df$Type]

# Corrige por extinção no V
df$MV0 <- df$MV_T - df$A_V

# log M*_obs segundo Bell+03 (Chabrier) usado no paper
df$logMstar_obs <- -0.4*(df$MV0 - M_Vsun) - 0.778 + 1.305*df$BV

# Propagação de incertezas (aprox. 1ª ordem)
# d logM*/d M_V = -0.4 ; d logM*/d(B-V) = 1.305
df$err_logMstar <- sqrt( (0.4 * df$err_MV_T)^2 + (1.305 * df$sBV)^2 )

# Grade do preditor em log M*
Msx <- seq(6.5, 12.5, length.out = M)

df <- filter(df,Type=="E")


# ------------------------------------------------------------
# JAGS data (agora com logMstar)
# ------------------------------------------------------------
jags.data <- list(
  N = nrow(df),
  logMstar_obs = df$logMstar_obs,
  err_logMstar = pmax(df$err_logMstar, 1e-3),   # evita variância zero
  N_GC_err = pmax(df$N_GC_err, 0.5),
  Msx = Msx,
  N_GC = df$N_GC,
  M = M
)

# ------------------------------------------------------------
# Inits
# ------------------------------------------------------------
init_fun <- function() {
  list(beta0 = 1,
       beta1 = 1,
       size  = 1,
       N_pred = rep(1, jags.data$M))
}

# ------------------------------------------------------------
# Modelo NB com preditor log M*
# ------------------------------------------------------------
model.NB <- "
model{
  beta0 ~ dnorm(0,1.0E-2)
  beta1 ~ dnorm(0,1.0E-2)
  size  ~ dunif(0.001,10)

  for(i in 1:N){
    # Latente: logMstar_true
    logMstar_true[i] ~ dunif(7.0, 12.5)
    # Observado com erro gaussiano
    logMstar_obs[i] ~ dnorm(logMstar_true[i], 1/pow(err_logMstar[i],2))

    eta[i] <- beta0 + beta1*logMstar_true[i]
    mu[i]  <- exp( max(-20, min(20, eta[i])) )
    p[i]   <- size/(size+mu[i])

    N_true[i] ~ dnegbin(p[i], size)

    # Likelihood para observação ruidosa do count
    N_GC[i] ~ dnorm(N_true[i], 1/pow(N_GC_err[i], 2))
  }

  for(j in 1:M){
    eta_pred[j] <- beta0 + beta1*Msx[j]
    mu_pred[j]  <- exp( max(-20, min(20, eta_pred[j])) )
    p_pred[j]   <- size/(size+mu_pred[j])
    N_pred[j] ~ dnegbin(p_pred[j], size)
  }
}
"

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------
fit <- jags(
  model   = textConnection(model.NB),
  data    = jags.data,
  inits   = list(init_fun(), init_fun(), init_fun()),
  n.chains= 3,
  n.burnin= 20000,
  n.iter  = 50000,
  parameters = c("beta0","beta1","size","mu_pred","N_pred")
)

# ------------------------------------------------------------
# Summaries para mu_pred e N_pred
# ------------------------------------------------------------
mux <- jagsresults(fit, params=c('mu_pred'))

library(tibble)
mu_df <- tibble(
  Msx   = Msx,
  mean  = mux[,"50%"],
  lwr50 = mux[,"25%"],
  upr50 = mux[,"75%"],
  lwr95 = mux[,"2.5%"],
  upr95 = mux[,"97.5%"]
)

# ------------------------------------------------------------
# Plot: relação média (mu_pred)
# ------------------------------------------------------------
cairo_pdf("NGC_Mstar.pdf", width = 7, height = 9)
ggplot() +
  geom_point(data = df, aes(x = logMstar_obs, y = N_GC), size = 1, color = "gray80",
             alpha=0.7) +
  geom_errorbar(data = df, aes(x = logMstar_obs, y = N_GC,
                               ymin = pmax(N_GC - N_GC_err, 0),
                               ymax = N_GC + N_GC_err),
                width = 0.02, color = "gray80", alpha=0.7) +
  geom_errorbarh(data = df, aes(y = N_GC,
                                xmin = logMstar_obs - err_logMstar,
                                xmax = logMstar_obs + err_logMstar),
                 height = 0.03, color = "gray80", alpha=0.7) +
  geom_ribbon(data = mu_df, aes(x = Msx, y = mean, ymin = lwr95, ymax = upr95), fill = "#FFA630") +
#  geom_ribbon(data = mu_df, aes(x = Msx, y = mean, ymin = lwr50, ymax = upr50), fill = "#F17720") +
  geom_step(data = mu_df, aes(x = Msx, y = mean), color = "#0474BA", size = 1) +
  scale_y_continuous(trans = plog,
                     breaks = c(0, 10, 100, 1e3, 1e4, 1e5),
                     labels = c("0", expression(10^1), expression(10^2),
                                expression(10^3), expression(10^4), expression(10^5))) +
  labs(x = expression(log[10]~(M['\u2605']/M['☉'])), y = expression(N[GC])) +
  theme_bw() +
  scale_fill_viridis_d(name="")+
  scale_color_viridis_d(name="")+
  theme(text = element_text(size = 22), legend.position = "none")
dev.off()
