library(tidyverse)
library(parallel)
library(pscl) # compute pseudo R2 for binary outcome
RNGkind("L'Ecuyer-CMRG")
set.seed(2024)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
n_rep <- 1000
n_cores <- 8


# 1 setup -----------------------------------------------

dim <- 5

mag_bias <- 6 # 6 

prop_unbias <- 0.5

ep <- "binary"

model_spec <- tribble(
  ~sm, ~om,
  T, T,
  T, F,
  F, T,
  F, F
)

sample_size <- tribble(
  ~n1, ~n0, ~m,
  50,  25,  150,
  # 50,  25,  75,
  # 50,  25,  300,
  # 175,  175,  350,
  # 175,  175,  750,
  # 175,  175,  1500
  #168, 167, 335
)

prob <- tribble(
  ~p0, ~p1,
  # 0.2, 0.3
  # 0.25,0.35
  0.3, 0.4
  # 0.35,0.45
  # 0.4, 0.5
  
  # 0.24, 0.3
  # 0.28,0.35
  # 0.32, 0.4
  # 0.36,0.45
  # 0.4, 0.5
)

bias_stru <- c("latent partial")

setup <- expand_grid(
  ep,
  sample_size,
  model_spec,
  bias_stru,
  mag_bias,
  prop_unbias,
  prob,
  dim
) %>% 
  mutate(prop_unbias = ifelse(mag_bias == 0, 0, prop_unbias)) %>% 
  distinct() %>% 
  mutate(case = row_number(), .before = everything()) %>% 
  group_by(ep, n1, n0, m, sm, om) %>% 
  mutate(
    main_case = case[1], 
    R2 = 0 # compute later
  ) %>% 
  ungroup

# if (!dir.exists("data")) {dir.create("data")}
if (!dir.exists(str_glue("data_{prob$p0}_{prob$p1}_response_Xdim_{dim}"))) {dir.create(str_glue("data_{prob$p0}_{prob$p1}_response_Xdim_{dim}"))}


# 2 generate simulation data -------------------------------------------------

for (case in 1:nrow(setup)) {
  # 2.1 load setup parameters
  list2env(setup[case,], envir = .GlobalEnv)
  N <- n1 + n0 + m
  sd_ec <- 0.5
  etaj <- 2
  
  # 2.2 generate (X, S, A, Y0, Y1)
  # if mag_bias == 0, generate data; otherwise, skip
  if (mag_bias == 0) {
    simdata_main <- mclapply(1:n_rep, function(i_rep) {
      # X
      p <- dim
      X <- matrix(runif(N * p, -2, 2), N, p)
      Xt <- scale(exp(X) + 10 * sin(X) * cos(X), center = F)
      
      # S
      if (sm) {
        X_pi <- X
      } else {
        X_pi <- Xt
      }
      #eta <- rep(0.1, p)
      eta <- rep(etaj, p) # eta determine the weight of each covariate
      X_pi_eta <- X_pi %*% eta
      eta0 <- uniroot(function(.x) {
        mean(plogis(.x + X_pi_eta)) - (n1 + n0) / N # eta0 to adjust the prob S=1 be more alike (n1 + n0) / N
      }, interval = c(-100, 100)
      )$root
      pi_S <- plogis(eta0 + X_pi_eta)
      S <- rbinom(N, size = 1, prob = pi_S)
      #S <- rbinom(N, 1, (n1 + n0) / (n1 + n0 + m))
      
      # A
      A <- S 
      A[S == 1] <- rbinom(sum(S == 1), 1, n1 / (n1 + n0)) # S=1 (RCT) has A=1 and A=0; S=0 (EC) all A=0
      
      # Y0, Y1
      if (om) {
        X_mu <- X
      } else {
        X_mu <- Xt
      }
      beta0 <- rep(1, p)
      beta1 <- rep(2, p)
      if (ep == "continuous") {
        mu0 <- as.vector(X_mu %*% beta0)
        mu1 <- as.vector(X_mu %*% beta1) + 0.3
        Y0 <- mu0 + rnorm(N)
        Y1 <- mu1 + rnorm(N)
        R2 <- mean(c(var(mu0) / var(Y0), var(mu1) / var(Y1)))
      } else if (ep == "binary") {
        beta0 <- rep(1, p) # weights for each covariates # scale down the weights for covariates
        mu0 <- as.vector(X_mu %*% beta0)
        mu0_intercept <- uniroot(function(.x) {
          mean(plogis(.x + mu0)) - prob$p0 #0.3 # M1: adjust 0.5 to lower value can reduce P(Y1=1)
        }, interval = c(-100, 100)
        )$root
        mu0 <- mu0_intercept + mu0
        Y0 <- rbinom(N, size = 1, prob = plogis(mu0))

        beta1 <- rep(2, p)
        mu1 <- as.vector(X_mu %*% beta1)
        mu1_intercept <- uniroot(function(.x) {
          mean(plogis(.x + mu1)) - prob$p1 #0.4 # original 0.6
        }, interval = c(-100, 100)
        )$root
        mu1 <- mu1_intercept + mu1
        Y1 <- rbinom(N, size = 1, prob = plogis(mu1))

        R2 <- mean(c(
          pR2(glm(Y0 ~ X_mu, family = binomial))[["McFadden"]],
          pR2(glm(Y1 ~ X_mu, family = binomial))[["McFadden"]]
        ))
      }
      # tau <- mean((Y1 - Y0)[S == 1])
      tau <- tribble(
        ~ tau_rd, ~ tau_rr, ~ tau_or,
        mean(Y1[S == 1]) - mean(Y0[S == 1]),
        mean(Y1[S == 1]) / mean(Y0[S == 1]),
        (mean(Y1[S == 1])/(1-mean(Y1[S == 1]))) / (mean(Y0[S == 1])/(1-mean(Y0[S == 1])))
        )
      # tau <- c(
      #   mean(Y1[S == 1]) - mean(Y0[S == 1]),
      #   mean(Y1[S == 1]) / mean(Y0[S == 1]),
      #   (mean(Y1[S == 1])/(1-mean(Y1[S == 1]))) / (mean(Y0[S == 1])/(1-mean(Y0[S == 1])))
      #   )
      
      # output
      lst(X, Xt, S, A, beta0, beta1, mu0, mu1, Y0, Y1, R2, tau)
    }, mc.cores = n_cores)
    # save data
    # save(simdata_main, file = str_glue("data/simdata_main_new_{case}.RData"))
    save(simdata_main, file = str_glue("data_{prob$p0}_{prob$p1}_response_Xdim_{dim}/simdata_main_new_{case}.RData"))
    
  } else {
    # load(str_glue("data/simdata_main_new_{main_case}.RData"))
    load(str_glue("data_{prob$p0}_{prob$p1}_response_Xdim_{dim}/simdata_main_new_{case}.RData"))
  }
  
  setup$R2[case] <- map_dbl(simdata_main, ~{.$R2}) %>% mean()
  
  # 2.3 generate (Y00, Y, Ynull)
  simdata <- mclapply(1:n_rep, function(i_rep) {
    # load required quantities
    mu0 <- simdata_main[[i_rep]]$mu0
    S <- simdata_main[[i_rep]]$S
    A <- simdata_main[[i_rep]]$A
    Y0 <- simdata_main[[i_rep]]$Y0
    Y1 <- simdata_main[[i_rep]]$Y1
    # Y00
    if (ep == "continuous") {
      if (bias_stru == "latent partial") {
        if (mag_bias == 0) {
          id_unbias <- which(S == 0)
          mu00 <- mu0
        } else {
          id_unbias <- sample(which(S == 0), round(sum(S == 0) * prop_unbias))
          mu00 <- mu0 - mag_bias
          mu00[id_unbias] <- mu0[id_unbias]
        }
      }
      Y00 <- mu00 + rnorm(N, sd = sd_ec)
    } else if (ep == "binary") {
      if (bias_stru == "latent partial") {
        if (mag_bias == 0) {
          id_unbias <- which(S == 0)
          mu00 <- mu0
        } else { # some samples have shifted pattern
          id_unbias <- sample(which(S == 0), round(sum(S == 0) * prop_unbias))
          mu00_intercept <- uniroot(function(.x) {
            mean(plogis(.x + mu0)) - prob$p0 - mag_bias / 20
          }, interval = c(-100, 100)
          )$root
          mu00 <- mu00_intercept + mu0
          mu00[id_unbias] <- mu0[id_unbias]
        }
      }
      Y00 <- rbinom(N, size = 1, prob = plogis(mu00)) # Y for EC
    }
    # (Y, Ynull)
    Y <- S * A * Y1 + S * (1 - A) * Y0 + (1 - S) * Y00
    Ynull <- S * Y0 + (1 - S) * Y00
    # output
    lst(mu00, Y00, Y, Ynull, id_unbias)
  }, mc.cores = n_cores)
  save(simdata, file = str_glue("data_{prob$p0}_{prob$p1}_response_Xdim_{dim}/simdata_new_{case}.RData"))
}

save(setup, file = str_glue("data_{prob$p0}_{prob$p1}_response_Xdim_{dim}/setup.RData"))


