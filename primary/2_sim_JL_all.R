library(tidyverse)
library(parallel)
library(tictoc)
# library(intFRT)
RNGkind("L'Ecuyer-CMRG")
set.seed(2025)
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# source("fun_JL_github.R")
# source("fun_JL_github_test.R") # new conformal_p with label conditional score
# source("fun_JL_partial.R") # only change csb, nb aipw
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

# 2 simulation -------------------------------------------------

if (!dir.exists("output")) {dir.create("output")}
if (!dir.exists("output/cache")) {dir.create("output/cache")}
writeLines("", str_glue("output/cache/prog_{case}.txt"))

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
      file = str_glue("output/cache/dt_{case}_{i_rep}.RData")
    )
    o_f <- T
  } else {
    o_f <- F
  }
  
  # adaptive selection threshold: rely on the ec_borrow, ec_borrow rely on conformal_p
  ada_g1 <- compute_ada_gamma(
    Y, A, S, X, family = fam, gamma_grid = seq(0,1,0.05),
    n_rep_gamma = n_rep_g, cf_score = "LC"
  )

  ada_g2 <- compute_ada_gamma(
    Y, A, S, X, family = fam, gamma_grid = seq(0,1,0.05),
    n_rep_gamma = n_rep_g, cf_score = "NN"
  )
  
  # ada_g3 <- compute_ada_gamma(
  #   Y, A, S, X, family = fam, gamma_grid = seq(0,1,0.05),
  #   n_rep_gamma = n_rep_g, cf_score = "NNR"
  # )
  
  ada_g3 <- compute_ada_gamma(
    Y, A, S, X, family = fam, gamma_grid = seq(0,1,0.05),
    n_rep_gamma = n_rep_g, cf_score = "Standardized AR"
  )

  ## inference
  res_alt <- res_null <- NULL
  
  if (ep == "binary") {
    tryCatch({
      # result under alternative
      res_alt <- list(
        # ec_borrow(Y, A, S, X, "No Borrow DiM", fam, n_f, n_b, output_frt = o_f), 
        ec_borrow(Y, A, S, X, "No Borrow AIPW", fam, n_f, output_frt = o_f), 
        # ec_borrow(Y, A, S, X, "Borrow Naive", fam, n_f, output_frt = o_f),
        # ec_borrow(Y, A, S, X, "Borrow IPW", fam, n_f, n_b, output_frt = o_f), 
        # ec_borrow(Y, A, S, X, "Borrow CW", fam, n_f, n_b, output_frt = o_f), 
        # ec_borrow(Y, A, S, X, "Borrow OM", fam, n_f, n_b, output_frt = o_f),
        ec_borrow(Y, A, S, X, "Borrow AIPW", fam, n_f, output_frt = o_f),
        # ec_borrow(Y, A, S, X, "Borrow ACW", fam, n_f, output_frt = o_f),
        
        ec_borrow(Y, A, S, X, method = "Conformal Selective Borrow AIPW",family = fam, n_fisher = n_f,
                  gamma_sel = ada_g1,
                  cf = "cv+", # c("split", "full", "jackknife+", "cv+")
                  cf_score = "LC",
                  cf_model = "glm",
                  #n_boot = n_b,
                  cv_fold = 10, output_frt = o_f) %>% intFRT:::add_name("(LC-NN)"), 
        ec_borrow(Y, A, S, X, method = "Conformal Selective Borrow AIPW",family = fam, n_fisher = n_f,
                  gamma_sel = ada_g2,
                  cf = "cv+", # c("split", "full", "jackknife+", "cv+")
                  cf_score = "NN",
                  cf_model = "glm",
                  # n_boot = n_b,
                  cv_fold = 10, output_frt = o_f) %>% intFRT:::add_name("(NN)"),
        # ec_borrow(Y, A, S, X, method = "Conformal Selective Borrow AIPW",family = fam, n_fisher = n_f,
        #           gamma_sel = ada_g3,
        #           cf = "cv+", # c("split", "full", "jackknife+", "cv+")
        #           cf_score = "NNR",
        #           cf_model = "glm",
        #           n_boot = n_b,
        #           cv_fold = 10, output_frt = o_f) %>% intFRT:::add_name("(NNR)")
        ec_borrow(Y, A, S, X, method = "Conformal Selective Borrow AIPW",family = fam, n_fisher = n_f,
                  gamma_sel = ada_g3,
                  cf = "cv+", # c("split", "full", "jackknife+", "cv+")
                  cf_score = "Standardized AR",
                  cf_model = "glm",
                  # n_boot = n_b,
                  cv_fold = 10,
                  sd_est = "model", output_frt = o_f) %>% intFRT:::add_name("(Standardized AR_model)")
      )
      # result under null
      res_null <- list(
        # ec_borrow(Ynull, A, S, X, "No Borrow DiM", fam, n_f, n_b, output_frt = o_f),
        ec_borrow(Ynull, A, S, X, "No Borrow AIPW", fam, n_f, output_frt = o_f), 
        # ec_borrow(Ynull, A, S, X, "Borrow Naive", fam, n_f, output_frt = o_f),
        # ec_borrow(Ynull, A, S, X, "Borrow IPW", fam, n_f, n_b, output_frt = o_f),
        # ec_borrow(Ynull, A, S, X, "Borrow CW", fam, n_f, n_b, output_frt = o_f),
        # ec_borrow(Ynull, A, S, X, "Borrow OM", fam, n_f, n_b, output_frt = o_f),
        ec_borrow(Ynull, A, S, X, "Borrow AIPW", fam, n_f, output_frt = o_f), 
        # ec_borrow(Ynull, A, S, X, "Borrow ACW", fam, n_f, output_frt = o_f),
        
        ec_borrow(Ynull, A, S, X, method = "Conformal Selective Borrow AIPW",family = fam, n_fisher = n_f,
                  gamma_sel = ada_g1,
                  cf = "cv+", # c("split", "full", "jackknife+", "cv+")
                  cf_score = "LC",
                  cf_model = "glm",
                  # n_boot=n_b,
                  cv_fold = 10, output_frt = o_f) %>% intFRT:::add_name("(LC-NN)"),
        ec_borrow(Ynull, A, S, X, method = "Conformal Selective Borrow AIPW",family = fam, n_fisher = n_f,
                  gamma_sel = ada_g2,
                  cf = "cv+", # c("split", "full", "jackknife+", "cv+")
                  cf_score = "NN",
                  cf_model = "glm",
                  # n_boot=n_b,
                  cv_fold = 10, output_frt = o_f) %>% intFRT:::add_name("(NN)"),
        # ec_borrow(Ynull, A, S, X, method = "Conformal Selective Borrow AIPW",family = fam, n_fisher = n_f,
        #           gamma_sel = ada_g3,
        #           cf = "cv+", # c("split", "full", "jackknife+", "cv+")
        #           cf_score = "NNR",
        #           cf_model = "glm",
        #           n_boot=n_b,
        #           cv_fold = 10, output_frt = o_f) %>% intFRT:::add_name("(NNR)")
        ec_borrow(Ynull, A, S, X, method = "Conformal Selective Borrow AIPW",family = fam, n_fisher = n_f,
                  gamma_sel = ada_g3,
                  cf = "cv+", # c("split", "full", "jackknife+", "cv+")
                  cf_score = "Standardized AR",
                  cf_model = "glm",
                  split_train = 0.75,
                  # n_boot=n_b,
                  sd_est = "model", output_frt = o_f) %>% intFRT:::add_name("(Standardized AR_model)")
        
      )
    }, error = function(e) {
      # output error
      error_msg <- conditionMessage(e)
      save(
        error_msg, Y, A, S, X, Ynull,
        file = str_glue("output/cache/err_{case}_{i_rep}.RData")
      )
    })
  }

  # output
  cat(i_rep, "\n", file = str_glue("output/cache/prog_{case}.txt"), 
      append = TRUE)
  lst(res_alt, res_null, oracle_info)
}, mc.cores = n_c)
total_time0 <- toc()
total_time <- unname(total_time0$toc - total_time0$tic)
count_null <- data.frame(res_null = sum(sapply(seq_along(raw_result), function(i) is.null(raw_result[[i]]$res_null))),
                         res_alt = sum(sapply(seq_along(raw_result), function(i) is.null(raw_result[[i]]$res_alt)))
                         )
save(raw_result, count_null, total_time, file = str_glue("output/cache/raw_{case}.RData"))
# save(raw_result, total_time, file = str_glue("output/cache/raw_{case}.RData"))
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

save(inf_result, file = str_glue("output/cache/inf_{case}.RData"))


# 3 summary ----------------------------------------------------------------
metrics <- inf_result %>% 
  mutate(R2 = mean(R2)) %>%
  group_by(case, R2, method) %>% 
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
    gamma_sel = mean(gamma_sel), 
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
ref_metrics <- metrics %>%
  filter(grepl("No Borrow AIPW", method)) %>%
  group_by(estimand) %>%
  summarise(
    Var_ref       = Var[!grepl("FRT", method)],
    MSE_ref       = MSE[!grepl("FRT", method)],
    Width_ref     = Width[!grepl("FRT", method)],
    Power_asy_ref = Power[!grepl("FRT", method)],
    Power_frt_ref = Power[grepl("FRT", method)],
    .groups = "drop"
  )
metrics <- metrics %>%
  left_join(ref_metrics, by = "estimand") %>%
  mutate(
    `Var%`     = round(Var / Var_ref * 100),
    `MSE%`     = round(MSE / MSE_ref * 100),
    `Width%`   = round(Width / Width_ref * 100),
    `Pow%`     = ifelse(!grepl("FRT", method), round(Power / Power_asy_ref * 100), NA),
    `FRT Pow%` = ifelse(grepl("FRT", method),  round(Power / Power_frt_ref * 100), NA)
  ) %>%
  # Step 3: Remove helper columns and estimand to reach exactly 29 columns
  select(-Var_ref, -MSE_ref, -Width_ref, -Power_asy_ref, -Power_frt_ref, -estimand)

save(metrics, file = str_glue("output/metrics_{case}.RData"))

# metrics <- inf_result %>% 
#   mutate(
#     R2 = mean(R2)
#   ) %>%
#   group_by(
#     case, 
#     R2, 
#     method
#   ) %>%
#   summarise(
#     tau = mean(tau),
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
#   mutate(
#     estimand = case_when(
#       grepl("_RD", method) ~ "RD",
#       grepl("_RR", method) ~ "RR",
#       grepl("_OR", method) ~ "OR",
#       TRUE ~ NA_character_
#     ),
#     `Bias/SD%` = round(abs(Bias / SD) * 100),
#     .after = Bias
#   )
# metrics <- metrics %>%
#   mutate(
#     estimand = case_when(
#       grepl("_RD", method) ~ "RD",
#       grepl("_RR", method) ~ "RR",
#       grepl("_OR", method) ~ "OR",
#       TRUE ~ NA_character_
#     )
#   )
# 
# ref_metrics <- metrics %>%
#   filter(method %in% c("No Borrow AIPW_RD", "No Borrow AIPW_RR", "No Borrow AIPW_OR")) %>%
#   select(
#     estimand,
#     Var_ref = Var,
#     MSE_ref = MSE,
#     Width_ref = Width
#   )
# power_ref_frt <- metrics %>%
#   filter(method == "No Borrow AIPW_RD+FRT") %>%
#   pull(Power)
# 
# power_ref_asy <- metrics %>%
#   filter(method == "No Borrow AIPW_RD") %>%
#   pull(Power)
# 
# metrics <- metrics %>%
#   left_join(ref_metrics, by = "estimand") %>%
#   mutate(
#     `Var%` = round(Var / Var_ref * 100),
#     `MSE%` = round(MSE / MSE_ref * 100),
#     `Width%` = round(Width / Width_ref * 100)
#   ) %>%
#   select(-Var_ref, -MSE_ref, -Width_ref, -estimand)
# metrics <- metrics %>%
#   mutate(
#     `Pow%` = ifelse(!grepl("FRT", method), round(Power / power_ref_asy * 100), NA),
#     `FRT Pow%` = ifelse(grepl("FRT", method), round(Power / power_ref_frt * 100), NA)
#   )


