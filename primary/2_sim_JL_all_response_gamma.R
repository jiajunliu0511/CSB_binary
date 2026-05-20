install.packages("vctrs")
install.packages("rlang")
install.packages("devtools")
devtools::install_github("ke-zhu/intFRT")
library(intFRT)

library(tidyverse)
library(parallel)
library(tictoc)

RNGkind("L'Ecuyer-CMRG")
set.seed(2025)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("fun_JL_all.R")
print(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))) {
  test <- TRUE
  if (!exists("case")) {
    case <- 1
  }
  n_c <- 8
  n_rep <- 5 # 100
  n_f <- 2
  n_rep_g <- 2
  n_b <- 2
} else {
  test <- FALSE
  case <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  n_c <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
  n_rep <- 400 #400
  n_f <- 2000 #1000
  n_b <- 1000
  n_rep_g <- 200
}

print(str_glue("test = {test}"))
print(str_glue("case = {case}"))
print(str_glue("n_cores = {n_c}"))
print(str_glue("n_rep = {n_rep}"))


# 1 setup -----------------------------------------------
# load("data/setup.RData")
load("data_0.3_0.4/setup.RData")
list2env(setup[case,], envir = .GlobalEnv)
if (ep == "continuous") {
  fam <- "gaussian"
} else if (ep == "binary") {
  fam <- "binomial"
}
load(str_glue("data_0.3_0.4/simdata_main_new_{main_case}.RData"))
load(str_glue("data_0.3_0.4/simdata_new_{case}.RData"))
# load(str_glue("data/simdata_main_new_{main_case}.RData"))
# load(str_glue("data/simdata_new_{case}.RData"))

gamma_grid_fix <- c(0.2, 0.4, 0.6, 0.8)

# 2 simulation -------------------------------------------------

if (!dir.exists("output_gamma")) {dir.create("output_gamma")}
if (!dir.exists("output_gamma/cache")) {dir.create("output_gamma/cache")}
writeLines("", str_glue("output_gamma/cache/prog_{case}.txt"))

run_csb_grid <- function(Y, A, S, X, gamma_grid, cf_score, fam, n_f, o_f) {
  score_lab <- c(
    LC = "(LC-NN)",
    NN = "(NN)",
    AR = "(AR)"
  )[cf_score]
  
  lapply(gamma_grid, function(g) {
    out <- ec_borrow(
      Y, A, S, X,
      method = "Conformal Selective Borrow AIPW",
      family = fam,
      n_fisher = n_f,
      gamma_sel = g,
      cf = "cv+",
      cf_score = cf_score,
      cf_model = "glm",
      cv_fold = 10,
      output_frt = o_f
    ) %>% add_name(score_lab)
    
    out$gamma_sel <- g
    out
  })
}

# i_rep <- 84
tic()
raw_result <- mclapply(1:n_rep, function(i_rep) {
  # save oracle information
  oracle_info <- tibble(
    tau = simdata_main[[i_rep]]$tau, 
    R2 = simdata_main[[i_rep]]$R2, 
    id_unbias = list(simdata[[i_rep]]$id_unbias)
  )
  # load data
  A <- simdata_main[[i_rep]]$A
  S <- simdata_main[[i_rep]]$S
  X <- simdata_main[[i_rep]]$X
  Y <- simdata[[i_rep]]$Y
  Ynull <- simdata[[i_rep]]$Ynull
  
  # save one sample data & output_frt
  if (i_rep == 1) {
    save(
      Y, A, S, X, Ynull, oracle_info,
      file = str_glue("output_gamma/cache/dt_{case}_{i_rep}.RData")
    )
    o_f <- T
  } else {
    o_f <- F
  }

  ## inference
  res_alt <- res_null <- NULL
  
  if (ep == "binary") {
    tryCatch({
      # result under alternative
      res_alt <- c(
        list(
          ec_borrow(Y, A, S, X, "No Borrow AIPW", fam, n_f, output_frt = o_f),
          ec_borrow(Y, A, S, X, "Borrow AIPW", fam, n_f, output_frt = o_f)
        ),
        run_csb_grid(Y, A, S, X, gamma_grid_fix, "LC", fam, n_f, o_f),
        run_csb_grid(Y, A, S, X, gamma_grid_fix, "NN", fam, n_f, o_f),
        run_csb_grid(Y, A, S, X, gamma_grid_fix, "AR", fam, n_f, o_f)
      )
      
      # result under null
      res_null <- c(
        list(
          ec_borrow(Ynull, A, S, X, "No Borrow AIPW", fam, n_f, output_frt = o_f),
          ec_borrow(Ynull, A, S, X, "Borrow AIPW", fam, n_f, output_frt = o_f)
        ),
        run_csb_grid(Ynull, A, S, X, gamma_grid_fix, "LC", fam, n_f, o_f),
        run_csb_grid(Ynull, A, S, X, gamma_grid_fix, "NN", fam, n_f, o_f),
        run_csb_grid(Ynull, A, S, X, gamma_grid_fix, "AR", fam, n_f, o_f)
      )
    }, error = function(e) {
      # output error
      error_msg <- conditionMessage(e)
      save(
        error_msg, Y, A, S, X, Ynull,
        file = str_glue("output_gamma/cache/err_{case}_{i_rep}.RData")
      )
    })
  }

  # output
  cat(i_rep, "\n", file = str_glue("output_gamma/cache/prog_{case}.txt"), 
      append = TRUE)
  lst(res_alt, res_null, oracle_info)
}, mc.cores = n_c)
total_time0 <- toc()
total_time <- unname(total_time0$toc - total_time0$tic)
count_null <- data.frame(res_null = sum(sapply(seq_along(raw_result), function(i) is.null(raw_result[[i]]$res_null))),
                         res_alt = sum(sapply(seq_along(raw_result), function(i) is.null(raw_result[[i]]$res_alt)))
                         )
save(raw_result, count_null, total_time, file = str_glue("output_gamma/cache/raw_{case}.RData"))
# save(raw_result, total_time, file = str_glue("output_gamma/cache/raw_{case}.RData"))
# sum(sapply(seq_along(raw_result), function(i) is.null(raw_result[[i]]$res_alt))) # test how many nulls in raw results

# organize raw result
inf_result <- imap_dfr(raw_result, function(raw_i, i) {
  if (!is.null(raw_i$res_alt) & !is.null(raw_i$res_null)) {
    # each replication
    map2_dfr(raw_i$res_alt, raw_i$res_null, function(res_alt_j, res_null_j) {
      # each method
      id_ec_i <- res_alt_j$dat_info$id_ec[[1]]
      id_unbias_i <- raw_i$oracle_info$id_unbias[[1]]
      id_sel_j <- res_alt_j$out$id_sel[[1]]
      n_ec <- length(id_ec_i)
      n_unbias <- length(id_unbias_i)
      n_sel <- length(id_sel_j)
      n_TU <- intersect(id_sel_j, id_unbias_i) %>% length()
      n_FU <- n_sel - n_TU
      n_FB <- n_unbias - n_TU
      n_TB <- n_ec - n_TU - n_FU - n_FB
      #c(n_TU, n_FU, n_FB, n_TB) %>% matrix(2, 2, byrow = T)
      res_alt_j$res %>% 
        mutate(
          p_value_null = res_null_j$res$p_value,
          p_value_null_upper = res_null_j$res$p_value_upper,
          p_value_null_lower = res_null_j$res$p_value_lower,
          `FU/Sel` = ifelse(n_sel == 0, 0, n_FU / n_sel),
          `TU/Unbias` = n_TU / n_unbias,
          gamma_sel = res_alt_j$gamma_sel
          # cf_sig = res_alt_j$cf_sig
        )
    }) %>%
      mutate(
        rep = i, 
        tau = raw_i$oracle_info$tau,
        R2 = raw_i$oracle_info$R2,
        .before = everything()
      )
  } else {
    NULL
  }
}) %>%
  mutate(case = case, .before = everything()) %>% 
  mutate(method = as_factor(method))

inf_result <- inf_result %>%
  mutate(
    tau = case_when(
      str_detect(method, "RD") ~ tau$tau_rd,
      str_detect(method, "RR") ~ tau$tau_rr,
      str_detect(method, "OR") ~ tau$tau_or,
      TRUE ~ NA_real_
    )
  )

save(inf_result, file = str_glue("output_gamma/cache/inf_{case}.RData"))


# 3 summary ----------------------------------------------------------------
metrics <- inf_result %>% 
  mutate(
    R2 = mean(R2)
  ) %>%
  group_by(
    case, 
    R2, 
    method,
    gamma_sel
  ) %>%
  summarise(
    tau = mean(tau),
    Bias = mean(est - tau),
    SE = mean(se),
    SD = sd(est),
    Var = var(est),
    MSE = mean((est - tau)^2),
    CP = mean(ci_l <= tau & tau <= ci_u),
    Width = mean(ci_u - ci_l),
    `Type I` = mean(p_value_null <= 0.05), 
    `Type I lower` = mean(p_value_null_lower <= 0.05), 
    `Type I upper` = mean(p_value_null_upper <= 0.05), 
    Power = mean(p_value <= 0.05),
    Power_lower = mean(p_value_lower <= 0.05),
    Power_upper = mean(p_value_upper <= 0.05),
    n_sel = round(mean(n_sel)),
    ess_sel = round(mean(ess_sel)),
    `FU/Sel` = ifelse(is.na(n_sel), NA, mean(`FU/Sel`)),
    `TU/Unbias` = ifelse(is.na(n_sel), NA, mean(`TU/Unbias`)),
    runtime = mean(runtime),
    gamma_sel = mean(gamma_sel), # updated
    .groups = "drop"
  ) %>% 
  mutate(
    estimand = case_when(
      grepl("_RD", method) ~ "RD",
      grepl("_RR", method) ~ "RR",
      grepl("_OR", method) ~ "OR",
      TRUE ~ NA_character_
    ),
    `Bias/SD%` = round(abs(Bias / SD) * 100),
    .after = Bias
  )
metrics <- metrics %>%
  mutate(
    estimand = case_when(
      grepl("_RD", method) ~ "RD",
      grepl("_RR", method) ~ "RR",
      grepl("_OR", method) ~ "OR",
      TRUE ~ NA_character_
    )
  )

ref_metrics <- metrics %>%
  filter(method %in% c("No Borrow AIPW_RD", "No Borrow AIPW_RR", "No Borrow AIPW_OR")) %>%
  select(
    estimand,
    Var_ref = Var,
    MSE_ref = MSE,
    Width_ref = Width
  )
power_ref_frt <- metrics %>%
  filter(method == "No Borrow AIPW_RD+FRT") %>%
  pull(Power)

power_ref_asy <- metrics %>%
  filter(method == "No Borrow AIPW_RD") %>%
  pull(Power)

metrics <- metrics %>%
  left_join(ref_metrics, by = "estimand") %>%
  mutate(
    `Var%` = round(Var / Var_ref * 100),
    `MSE%` = round(MSE / MSE_ref * 100),
    `Width%` = round(Width / Width_ref * 100)
  ) %>%
  select(-Var_ref, -MSE_ref, -Width_ref, -estimand)
metrics <- metrics %>%
  mutate(
    `Pow%` = ifelse(!grepl("FRT", method), round(Power / power_ref_asy * 100), NA),
    `FRT Pow%` = ifelse(grepl("FRT", method), round(Power / power_ref_frt * 100), NA)
  )

save(metrics, file = str_glue("output_gamma/metrics_{case}.RData"))

  # mutate(`Bias/SD%` = round(abs(Bias / SD) * 100), .after = `Bias`) %>%
  # select(-SD) %>% 
  # mutate(`Var%` = round((Var / Var[1]) * 100), .after = Var) %>%
  # mutate(`MSE%` = round((MSE / MSE[1]) * 100), .after = MSE) %>%
  # mutate(`Width%` = round((Width / Width[1]) * 100), .after = Width) %>%
  # mutate(`Pow%` = round((Power / Power[1]) * 100), .after = Power)

# metrics <- inf_result %>% 
#   mutate(
#     R2 = mean(R2),
#     tau = mean(tau)
#   ) %>%
#   group_by(
#     case, 
#     R2, 
#     tau, 
#     method
#   ) %>%
#   summarise(
#     Bias = mean(est - tau),
#     SE = mean(se),
#     SD = sd(est),
#     Var = var(est),
#     MSE = mean((est - tau)^2),
#     CP = mean(ci_l <= tau & tau <= ci_u),
#     Width = mean(ci_u - ci_l),
#     `Type I` = mean(p_value_null <= 0.05), 
#     `Type I lower` = mean(p_value_null_lower <= 0.05), 
#     `Type I upper` = mean(p_value_null_upper <= 0.05), 
#     Power = mean(p_value <= 0.05),
#     Power_lower = mean(p_value_lower <= 0.05),
#     Power_upper = mean(p_value_upper <= 0.05),
#     n_sel = round(mean(n_sel)),
#     ess_sel = round(mean(ess_sel)),
#     `FU/Sel` = ifelse(is.na(n_sel), NA, mean(`FU/Sel`)),
#     `TU/Unbias` = ifelse(is.na(n_sel), NA, mean(`TU/Unbias`)),
#     runtime = mean(runtime),
#     gamma_sel = mean(gamma_sel), # updated
#     .groups = "drop"
#   ) %>% 
#   mutate(`Bias/SD%` = round(abs(Bias / SD) * 100), .after = `Bias`) %>%
#   select(-SD) %>% 
#   mutate(`Var%` = round((Var / Var[1]) * 100), .after = Var) %>%
#   mutate(`MSE%` = round((MSE / MSE[1]) * 100), .after = MSE) %>%
#   mutate(`Width%` = round((Width / Width[1]) * 100), .after = Width) %>%
#   mutate(`Pow%` = round((Power / Power[1]) * 100), .after = Power)


# metrics_narm <- inf_result %>% 
#   mutate(
#     R2 = mean(R2, na.rm = T),
#     tau = mean(tau, na.rm = T)
#   ) %>%
#   group_by(
#     case, 
#     R2, 
#     tau, 
#     method
#   ) %>%
#   summarise(
#     Bias = mean(est - tau, na.rm = T),
#     SE = mean(se, na.rm = T),
#     SD = sd(est, na.rm = T),
#     Var = var(est, na.rm = T),
#     MSE = mean((est - tau)^2, na.rm = T),
#     CP = mean(ci_l <= tau & tau <= ci_u, na.rm = T),
#     Width = mean(ci_u - ci_l, na.rm = T),
#     `Type I` = mean(p_value_null <= 0.05, na.rm = T), 
#     Power = mean(p_value <= 0.05, na.rm = T),
#     n_sel = round(mean(n_sel, na.rm = T)),
#     ess_sel = round(mean(ess_sel, na.rm = T)),
#     `FU/Sel` = ifelse(is.na(n_sel), NA, mean(`FU/Sel`, na.rm = T)),
#     `TU/Unbias` = ifelse(is.na(n_sel), NA, mean(`TU/Unbias`, na.rm = T)),
#     runtime = mean(runtime, na.rm = T),
#     .groups = "drop"
#   ) %>% 
#   mutate(`Bias/SD%` = round(abs(Bias / SD) * 100), .after = `Bias`) %>%
#   select(-SD) %>% 
#   mutate(`Var%` = round((Var / Var[1]) * 100), .after = Var) %>%
#   mutate(`MSE%` = round((MSE / MSE[1]) * 100), .after = MSE) %>%
#   mutate(`Width%` = round((Width / Width[1]) * 100), .after = Width) %>%
#   mutate(`Pow%` = round((Power / Power[1]) * 100), .after = Power)
