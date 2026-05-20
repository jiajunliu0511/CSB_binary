fit_outcome_model <- function(dat, family, outcome_model) {
  if (outcome_model == "glm") {
    m1 <- glm(Y ~ X, family = family, dat %>% filter(A == 1)) %>%
      predict(dat, "response")
    m0 <- glm(Y ~ X, family = family, dat %>% filter(A == 0)) %>%
      predict(dat, "response")
  } else if (outcome_model == "rf") {
    rlang::check_installed("randomForest", reason = "to use `randomForest()`")
    if (family == "gaussian") {
      m1 <- randomForest::randomForest(dat$X[dat$A == 1, , drop = FALSE], dat$Y[dat$A == 1]) %>%
        predict(dat$X)
      m0 <- randomForest::randomForest(dat$X[dat$A == 0, , drop = FALSE], dat$Y[dat$A == 0]) %>%
        predict(dat$X)
    } else if (family == "binomial") {
      m1 <- randomForest::randomForest(
        dat$X[dat$A == 1, , drop = FALSE], as.factor(dat$Y[dat$A == 1])
      ) %>%
        predict(dat$X, type = "prob") %>%
        {.[, 2]}
      m0 <- randomForest::randomForest(
        dat$X[dat$A == 0, , drop = FALSE], as.factor(dat$Y[dat$A == 0])
      ) %>%
        predict(dat$X, type = "prob") %>%
        {.[, 2]}
    }
  }
  lst(m1, m0)
}

compute_cw <- function(S, X) {
  n_rct <- sum(S == 1)
  n_ec <- sum(S == 0)

  # calibration weights (entropy balancing) for EC sample
  dat_df <- as.data.frame(X)
  dat_df$S <- S
  W.out <- WeightIt::weightit(
    as.formula(paste("S ~", paste(names(dat_df)[-ncol(dat_df)], collapse = " + "))),
    data = dat_df, estimand = "ATT", method = "ebal", include.obj = T
  )
  w <- W.out$weights / n_ec

  # reproduce calibration weights for all sample, since w[S == 1] is constant
  Xscale <- WeightIt:::.make_closer_to_1(X)
  Xtarget <- map(1:ncol(Xscale), function(j) {
    Xscale[,j] - colMeans(Xscale[S == 1, , drop=F])[j]
  }) %>% sapply(function(x) x)
  ww_init <- exp(-Xtarget %*% W.out$obj$`0`$par)
  ww <- ww_init / sum(ww_init[S == 0])
  # # check ww[S == 0] = w[S == 0]
  # max(abs(ww[S == 0] - w[S == 0]))
  # plot(ww[S == 0], w[S==0])
  # tibble(ww[S == 0], w[S==0])

  # calibration weights for RCT sample
  w[S == 1] <- ww[S == 1]
  qhat <- w * n_rct

  if (any(is.na(qhat) | is.nan(qhat) | is.infinite(qhat))) {
    pS <- glm(S ~ X, family = "binomial", tibble(S, X)) %>%
      predict(tibble(S, X), "response")
    qhat <- pS / (1 - pS)
  }

  # check q_cw vs q_ipw
  # pS <- glm(S ~ X, family = "binomial", tibble(S,X)) %>%
  #   predict(tibble(S,X),"response")
  # tibble(q_cw = qhat, q_ipw = pS / (1 - pS)) %>%
  #   ggplot(aes(q_cw, q_ipw)) +
  #   geom_point()+
  #   geom_abline(slope = 1, intercept = 0, color = "red")
  # max(abs(qhat - pS / (1 - pS)))

  qhat
}
dat <- data.frame(
  ID = 1:100,             # Unique ID for each observation
  A = sample(0:1, 100, replace = TRUE),  # Treatment assignment (0 = control, 1 = treated)
  S = sample(0:1, 100, replace = TRUE),  # Selection indicator (1 = in RCT, 0 = not)
  Y = rbinom(100, 1, prob = 0.5), # Outcome variable (continuous)
  X = rnorm(100, mean = 5, sd = 10)
)
rct_aipw <- function(dat, family, outcome_model, small_n_adj) {
  n_rt <- dat %>% filter(A == 1, S == 1) %>% nrow
  n_rct <- dat %>% filter(S == 1) %>% nrow
  dat_rct <- dat %>% filter(S == 1)

  # outcome model
  m10 <- fit_outcome_model(dat_rct, family, outcome_model)
  m1 <- m10$m1
  m0 <- m10$m0

  pA <- n_rt / n_rct
  d1 <- with(
    dat_rct,
    m1 + A / pA * (Y - m1)
  )
  d0 <- with(
    dat_rct,
    m0 + (1 - A) / (1 - pA) * (Y - m0)
  )
  d_rd <- mean(d1) - mean(d0)
  d_rr <- mean(d1) / mean(d0)
  d_or <- (mean(d1)/(1-mean(d1))) / (mean(d0)/(1-mean(d0)))
  est <- c(d_rd, d_rr, d_or)
  
  # d <- with(
  #   dat_rct,
  #   m1 + A / pA * (Y - m1) - m0 - (1 - A) / (1 - pA) * (Y - m0) # d1-d0, take mean then RR,OR,RD
  # )
  
  dd_rd <- d1 - d0 - d_rd
  dd_rr <- (d1 - d0*d_rr) / mean(d0)
  dd_or <- (1-mean(d0))/(mean(d0)) * ((d1 - mean(d1))/(1-mean(d1))^2 - (d0 - mean(d0))/(1-mean(d0))^2 * d_or)

  if (small_n_adj) {
    dof <- max(n_rct - ncol(dat$X) * 2, 1)

    se_rd <- sqrt(sum((dd_rd - mean(dd_rd))^2) / dof^2)
    se_rr <- sqrt(sum((dd_rr - mean(dd_rr))^2) / dof^2)
    se_or <- sqrt(sum((dd_or - mean(dd_or))^2) / dof^2)

  } else {
    se_rd <- sqrt(sum((dd_rd - mean(dd_rd))^2) / n_rct^2)
    se_rr <- sqrt(sum((dd_rr - mean(dd_rr))^2) / n_rct^2)
    se_or <- sqrt(sum((dd_or - mean(dd_or))^2) / n_rct^2)
  }
  se <- c(se_rd, se_rr, se_or)
  
  tibble(
    # est = mean(d),
    est = est,
    se = se
    # d = list(d0,d1),
    # d = list(d) # newly added in utils
  )
}

rct_ec_aipw_acw <- function(dat, family, outcome_model, max_r, small_n_adj, cw = FALSE) {
  n_rt <- dat %>% filter(A == 1, S == 1) %>% nrow
  n_rc <- dat %>% filter(A == 0, S == 1) %>% nrow
  n_rct <- dat %>% filter(S == 1) %>% nrow
  n_all <- dat %>% nrow

  # outcome model
  m10 <- fit_outcome_model(dat, family, outcome_model)
  m1 <- m10$m1
  m0 <- m10$m0

  # treatment group
  # use true propensity score
  pA <- n_rt / n_rct
  w1 <- with(
    dat,
    S * A / pA
  )
  d1 <- with(
    dat,
    (n_all / n_rct) * (S * m1 +  w1 * (Y - m1))
  )

  # control group
  # compute r
  if (family == "gaussian") {
    r1 <- glm(Y ~ X, family = family, dat %>% filter(A == 0, S == 1)) %>%
      resid(type = "response") %>% var
    r0 <- glm(Y ~ X, family = family, dat %>% filter(S == 0)) %>%
      resid(type = "response") %>% var
    r <- min(r1 / r0, max_r)
  } else if (family == "binomial") {
    # for binary outcome, under exchangeablity assumption, r=1 (Li et al., 2023)
    r <- 1
  }
  # compute qhat
  if (cw) {
    qhat <- compute_cw(dat$S, dat$X)
  } else {
    pS <- glm(S ~ X, family = "binomial", dat) %>% predict(dat, "response")
    qhat <- pS / (1 - pS)
  }
  w0init <- with(
    dat,
    qhat * (S * (1 - A) + (1 - S) * r) / (qhat * (1 - pA) + r)
  )
  w0 <- w0init / sum(w0init) * n_rct
  d0 <- with(
    dat,
    (n_all / n_rct) * (S * m0 +  w0 * (Y - m0))
  )

  # compute est
  d_rd <- mean(d1) - mean(d0)
  d_rr <- mean(d1) / mean(d0)
  d_or <- (mean(d1)/(1-mean(d1))) / (mean(d0)/(1-mean(d0)))
  est <- c(d_rd, d_rr, d_or)
  
  # compute se
  dd_rd <- d1 - d0 - dat$S/(n_rct/n_all)*d_rd
  dd_rr <- (d1 - d0*d_rr) / mean(d0)
  dd_or <- (1-mean(d0))/(mean(d0)) * ((d1 - dat$S/(n_rct/n_all)*mean(d1))/(1-mean(d1))^2 - (d0 - dat$S/(n_rct/n_all)*mean(d0))/(1-mean(d0))^2 * d_or)
  
  if (small_n_adj) {
    dof <- max(n_all - ncol(dat$X) * 3, 1)
    
    se_rd <- sqrt(sum((dd_rd - mean(dd_rd))^2) / dof^2)
    se_rr <- sqrt(sum((dd_rr - mean(dd_rr))^2) / dof^2)
    se_or <- sqrt(sum((dd_or - mean(dd_or))^2) / dof^2)
    
  } else {
    se_rd <- sqrt(sum((dd_rd - mean(dd_rd))^2) / n_all^2)
    se_rr <- sqrt(sum((dd_rr - mean(dd_rr))^2) / n_all^2)
    se_or <- sqrt(sum((dd_or - mean(dd_or))^2) / n_all^2)
  }
  se <- c(se_rd, se_rr, se_or)
  # only one estimand
  # if (small_n_adj) {
  #   dof <- max(n_all - ncol(dat$X) * 3, 1)
  #   se <- sqrt(sum((d - mean(d))^2) / dof^2)
  # } else {
  #   se <- sqrt(sum((d - mean(d))^2) / n_all^2)
  # }
  # se <- NA
  
  # output
  # tibble(est, se, ess_sel = max(0, ESS(w0) - n_rc))
  tibble(est, 
         se, 
         ess_sel = max(0, ESS(w0) - n_rc)
         # d = list(d0, d1)
         ) # newly added
}

fit_cf_model_mean <- function(dat_train, dat_pred, family, cf_model) {
  if (cf_model == "glm") {
    fit <- glm(y ~ x, family = family, data = dat_train)
    predict(fit, newdata = dat_pred, type = "response")
  } else if (cf_model == "rf") {
    if (family == "gaussian") {
      fit <- randomForest::randomForest(x = dat_train$x, y = dat_train$y)
      predict(fit, newdata = dat_pred$x)
    } else if (family == "binomial") {
      fit <- randomForest::randomForest(x = dat_train$x, y = as.factor(dat_train$y))
      predict(fit, newdata = dat_pred$x, type = "prob")[,2]
    }
  }
}

# changed to binary
fit_cf_model_mean_sd <- function(dat_train, dat_pred) { 
  # mean model
  fit1 <- glm(y ~ x, family = "binomial", data = dat_train)
  ar_train <- abs(resid(fit1))
  # absolute residual model
  fit2 <- glm(ar_train ~ x, family = gaussian(link = "log"), # ar are skewed, log is more helpful
              data = dat_train %>% mutate(ar_train)) 
  yhat <- predict(fit1, newdata = dat_pred, type = "response")
  sighat <- predict(fit2, newdata = dat_pred, type = "response")
  lst(yhat, sighat)
}

fit_cf_model_quantile <- function(dat_train, dat_pred, a, cf_model) { # family = "gaussian"
  if (cf_model == "glm") {
    fit_low <- quantreg::rq(y ~ x, tau = a/2, data = dat_train)
    fit_up <- quantreg::rq(y ~ x, tau = 1 - a/2, data = dat_train)
    cbind(
      predict(fit_low, newdata = dat_pred),
      predict(fit_up, newdata = dat_pred)
    )
  } else if (cf_model == "rf") {
    rlang::check_installed("grf", reason = "to use `quantile_forest()`")
    fit <- grf::quantile_forest(dat_train$x, dat_train$y, quantiles = c(a/2, 1 - a/2))
    predict(fit, dat_pred$x)$predictions
  }
}

conformal_p <- function(x, y, s, family, cf, cf_score, cf_model,
                        split_train, cv_fold, sig_level,sd_est) {
  id_rc <- which(s == 1) # train & cal
  id_test <- which(s == 0) # test
  
  if (cf == "jackknife+") {
    cv_fold <- length(id_rc)
    cf <- "cv+"
  }
  
  if (cf == "full") {
    p_cf <- map_dbl(id_test, function(j) { # for every EC j
      # data augmentation
      id_train <- id_pred <- c(j, id_rc)
      dat_train <- dat_pred <- tibble(y = y[id_train], x = x[id_train, , drop = FALSE])
      # conformal score
      if (cf_score == "AR") {
        # prediction model
        yhat_pred <- fit_cf_model_mean(dat_train, dat_pred, family, cf_model)
        # compute score
        score_pred <- abs(yhat_pred - y[id_pred])
      } else if (cf_score == "local AR") {
        # prediction model
        hat_pred <- fit_cf_model_mean_sd(dat_train, dat_pred)
        # compute score
        score_pred <- abs(hat_pred$yhat - y[id_pred]) / hat_pred$sighat
      } else if (cf_score == "CQR") {
        # prediction model
        qhat_pred <- fit_cf_model_quantile(dat_train, dat_pred, sig_level, cf_model)
        # compute score
        score_pred <- pmax(qhat_pred[,1] - y[id_pred], y[id_pred] - qhat_pred[,2])
      } else if (cf_score == "NN") {
        rlang::check_installed("RANN", reason = "to use `nn2()`")
        # 1-nearest-neighbor
        d1 <- RANN::nn2(
          dat_train %>% filter(y == 1) %>% pull(x),
          k = 2
        )$nn.dists[,2]
        d0 <- RANN::nn2(
          dat_train %>% filter(y == 0) %>% pull(x),
          k = 2
        )$nn.dists[,2]
        # compute score
        score_pred <- rep(0, nrow(dat_pred))
        score_pred[dat_pred$y == 1] <- d1
        score_pred[dat_pred$y == 0] <- d0
      }
      score_test_j <- score_pred[1]
      score_cal <- score_pred[-1]
      # conformal p-values
      mean(c(score_cal >= score_test_j, 1))
    })
  } else if (cf == "split") {
    # data split
    n_train <- ceiling(split_train * length(id_rc))
    n_cal <- length(id_rc) - n_train
    split_id_rc <- split(
      id_rc,
      c(rep(1, n_cal), rep(2, n_train)) %>% sample
    )
    id_cal <- split_id_rc[[1]]
    id_train <- split_id_rc[[2]]
    dat_train <- tibble(y = y[id_train], x = x[id_train, , drop = FALSE])
    id_pred <- c(id_test, id_cal)
    dat_pred <- tibble(y = y[id_pred], x = x[id_pred, , drop = FALSE])
    # conformal score
    if (cf_score == "AR") {
      # prediction model
      yhat_pred <- fit_cf_model_mean(dat_train, dat_pred, family, cf_model)
      # compute score
      score_pred <- abs(yhat_pred - y[id_pred])
    } else if (cf_score == "local AR") {
      # prediction model
      hat_pred <- fit_cf_model_mean_sd(dat_train, dat_pred)
      # compute score
      score_pred <- abs(hat_pred$yhat - y[id_pred]) / hat_pred$sighat
    } else if (cf_score == "CQR") {
      # prediction model
      qhat_pred <- fit_cf_model_quantile(dat_train, dat_pred, sig_level, cf_model)
      # compute score
      score_pred <- pmax(qhat_pred[,1] - y[id_pred], y[id_pred] - qhat_pred[,2])
    } else if (cf_score == "NN") {
      # 1-nearest-neighbor
      d1 <- RANN::nn2(
        dat_train %>% filter(y == 1) %>% pull(x),
        dat_pred %>% filter(y == 1) %>% pull(x),
        k = 1
      )$nn.dists %>% as.vector()
      d0 <- RANN::nn2(
        dat_train %>% filter(y == 0) %>% pull(x),
        dat_pred %>% filter(y == 0) %>% pull(x),
        k = 1
      )$nn.dists %>% as.vector()
      # compute score
      score_pred <- rep(0, nrow(dat_pred))
      score_pred[dat_pred$y == 1] <- d1
      score_pred[dat_pred$y == 0] <- d0
    }
    score_test <- head(score_pred, length(id_test))
    score_cal <- tail(score_pred, length(id_cal))
    # conformal p-values
    p_cf <- map_dbl(score_test, function(score_test_j) {
      mean(c(score_cal >= score_test_j, 1))
    })
  } else if (cf == "cv+") {
    if (cf_score=="LC"){ # label conditional conformal
      # data split
      id_rc_Y0 <- id_rc[which(y[id_rc]==0)] # they should be larger than 3, ow, we cannot do split, and train set is not enough
      id_rc_Y1 <- id_rc[which(y[id_rc]==1)]
      
      # if (length(id_rc_Y0) < 3 || length(id_rc_Y1) < 3) {
      #   stop("Not enough reference controls for Y=0 or Y=1 to perform label-conditional conformal split.")
      # }
      id_test_Y0 <- id_test[which(y[id_test]==0)]
      id_test_Y1 <- id_test[which(y[id_test]==1)]
      
      split_id_rc_Y0 <- split(
        id_rc_Y0,
        (rep(1:cv_fold, ceiling(length(id_rc_Y0) / cv_fold))[1:length(id_rc_Y0)]) %>% sample
      )
      
      split_id_rc_Y1 <- split(
        id_rc_Y1,
        (rep(1:cv_fold, ceiling(length(id_rc_Y1) / cv_fold))[1:length(id_rc_Y1)]) %>% sample
      )
      
      # conformal score for Y=0
      compare_all_Y0 <- map(split_id_rc_Y0, function(id_cal) {
        id_train <- setdiff(id_rc_Y0, id_cal) # RC-Cali = train
        dat_train <- tibble(y = y[id_train], x = x[id_train, , drop = FALSE])
        id_pred <- c(id_test_Y0, id_cal) # need predicted values
        dat_pred <- tibble(y = y[id_pred], x = x[id_pred, , drop = FALSE])
        # 1-nearest-neighbor
        d0 <- RANN::nn2(
          dat_train %>% filter(y == 0) %>% pull(x),
          dat_pred %>% filter(y == 0) %>% pull(x),
          k = 1
        )$nn.dists %>% as.vector()
        # compute score
        score_pred_Y0 <- rep(0, nrow(dat_pred))
        score_pred_Y0[dat_pred$y == 0] <- d0
        
        score_test_k_Y0 <- head(score_pred_Y0, length(id_test_Y0))
        score_cal_k_Y0 <- tail(score_pred_Y0, length(id_cal))
        # compare
        map_dbl(score_test_k_Y0, function(score_test_l) {
          sum(score_cal_k_Y0 >= score_test_l)
        })
      }) %>% sapply(function(x) x)
      # conformal p-values for Y=0
      p_cf_Y0 <- map_dbl(1:length(id_test_Y0), function(l) {
        (sum(compare_all_Y0[l,]) + 1) / (length(id_rc_Y0) + 1)
      })
      
      # conformal score for Y=1
      compare_all_Y1 <- map(split_id_rc_Y1, function(id_cal) {
        id_train <- setdiff(id_rc_Y1, id_cal) # RC-Cali = train
        dat_train <- tibble(y = y[id_train], x = x[id_train, , drop = FALSE])
        id_pred <- c(id_test_Y1, id_cal) # need predicted values
        dat_pred <- tibble(y = y[id_pred], x = x[id_pred, , drop = FALSE])
        # 1-nearest-neighbor
        d1 <- RANN::nn2(
          dat_train %>% filter(y == 1) %>% pull(x),
          dat_pred %>% filter(y == 1) %>% pull(x),
          k = 1
        )$nn.dists %>% as.vector()
        # compute score
        score_pred_Y1 <- rep(0, nrow(dat_pred))
        score_pred_Y1[dat_pred$y == 1] <- d1
        
        score_test_k_Y1 <- head(score_pred_Y1, length(id_test_Y1))
        score_cal_k_Y1 <- tail(score_pred_Y1, length(id_cal))
        # compare
        map_dbl(score_test_k_Y1, function(score_test_l) {
          sum(score_cal_k_Y1 >= score_test_l)
        })
      }) %>% sapply(function(x) x)
      # conformal p-values for Y=1
      p_cf_Y1 <- map_dbl(1:length(id_test_Y1), function(l) {
        (sum(compare_all_Y1[l,]) + 1) / (length(id_rc_Y1) + 1)
      })
      
      # combine conformal p-values
      p_cf <- rep(NA,length(id_test))
      ind_Y0 <- match(id_test_Y0, id_test)
      ind_Y1 <- match(id_test_Y1, id_test)
      p_cf[ind_Y0] <- p_cf_Y0
      p_cf[ind_Y1] <- p_cf_Y1
    } else { # calculate conformal p together
      split_id_rc <- split(
        id_rc,
        (rep(1:cv_fold, ceiling(length(id_rc) / cv_fold))[1:length(id_rc)]) %>% sample
      )
      # conformal score
      compare_all <- map(split_id_rc, function(id_cal) {
        id_train <- setdiff(id_rc, id_cal) # RC-Cali = train
        dat_train <- tibble(y = y[id_train], x = x[id_train, , drop = FALSE])
        id_pred <- c(id_test, id_cal) # need predicted values
        dat_pred <- tibble(y = y[id_pred], x = x[id_pred, , drop = FALSE])
        if (cf_score == "AR") {
          # prediction model
          yhat_pred <- fit_cf_model_mean(dat_train, dat_pred, family, cf_model)
          # compute score
          score_pred <- abs(yhat_pred - y[id_pred])
        } else if (cf_score == "local AR") {
          # prediction model
          hat_pred <- fit_cf_model_mean_sd(dat_train, dat_pred)
          # compute score
          score_pred <- abs(hat_pred$yhat - y[id_pred]) / hat_pred$sighat
        } else if (cf_score == "CQR") {
          # prediction model
          qhat_pred <- fit_cf_model_quantile(dat_train, dat_pred, sig_level, cf_model)
          # compute score
          score_pred <- pmax(qhat_pred[,1] - y[id_pred], y[id_pred] - qhat_pred[,2])
        } else if (cf_score == "NN") {
          # 1-nearest-neighbor
          d1 <- RANN::nn2(
            dat_train %>% filter(y == 1) %>% pull(x),
            dat_pred %>% filter(y == 1) %>% pull(x),
            k = 1
          )$nn.dists %>% as.vector()
          d0 <- RANN::nn2(
            dat_train %>% filter(y == 0) %>% pull(x),
            dat_pred %>% filter(y == 0) %>% pull(x),
            k = 1
          )$nn.dists %>% as.vector()
          # compute score
          score_pred <- rep(0, nrow(dat_pred))
          score_pred[dat_pred$y == 1] <- d1
          score_pred[dat_pred$y == 0] <- d0
        } else if (cf_score == "NNR"){
          # nearest-neighbor-ratio
          d11 <- RANN::nn2(
            dat_train %>% filter(y == 1) %>% pull(x),
            dat_pred %>% filter(y == 1) %>% pull(x),
            k = 1
          )$nn.dists %>% as.vector()
          d00 <- RANN::nn2(
            dat_train %>% filter(y == 0) %>% pull(x),
            dat_pred %>% filter(y == 0) %>% pull(x),
            k = 1
          )$nn.dists %>% as.vector()
          d10 <- RANN::nn2(
            dat_train %>% filter(y == 0) %>% pull(x),
            dat_pred %>% filter(y == 1) %>% pull(x),
            k = 1
          )$nn.dists %>% as.vector()
          d01 <- RANN::nn2(
            dat_train %>% filter(y == 1) %>% pull(x),
            dat_pred %>% filter(y == 0) %>% pull(x),
            k = 1
          )$nn.dists %>% as.vector()
          # compute score
          score_pred <- rep(0, nrow(dat_pred))
          score_pred[dat_pred$y == 1] <- d11/d10
          score_pred[dat_pred$y == 0] <- d00/d01
        } else if (cf_score == "HPS"){
          # high-probability score
          yhat_pred <- fit_cf_model_mean(dat_train, dat_pred, family, cf_model)
          # compute score
          score_pred <- rep(0, nrow(dat_pred))
          score_pred[dat_pred$y == 1] <- -yhat_pred[dat_pred$y == 1]
          score_pred[dat_pred$y == 0] <- -(1-yhat_pred[dat_pred$y == 0])
        } else if (cf_score == "Standardized AR"){
          # standardized absolute residual
          # prediction model
          hat_pred <- fit_cf_model_mean_sd(dat_train, dat_pred)
          
          # compute sd in different ways
          if (sd_est == "model"){
            score_pred <- abs(hat_pred$yhat - y[id_pred]) / hat_pred$sighat
          } else if (sd_est == "direct"){
            # phat(1-phat)
            score_pred <- abs(hat_pred$yhat - y[id_pred]) / sqrt(hat_pred$yhat*(1-hat_pred$yhat))
          }
        } else if (cf_score == "APS"){
          # adaptive prediction sets
          yhat_pred <- fit_cf_model_mean(dat_train, dat_pred, family, cf_model)
          # compute score
          score_pred <- rep(0, nrow(dat_pred))
          score_pred[dat_pred$y == 1] <- -((1-yhat_pred[dat_pred$y == 1]) * as.integer((1-yhat_pred[dat_pred$y == 1]) <= yhat_pred[dat_pred$y == 1]) + yhat_pred[dat_pred$y == 1])
          score_pred[dat_pred$y == 0] <- -(yhat_pred[dat_pred$y == 0] * as.integer(yhat_pred[dat_pred$y == 0] <= (1-yhat_pred[dat_pred$y == 0])) + (1-yhat_pred[dat_pred$y == 0]))
        } else if (cf_score == "RANK"){
          # rank-based
          yhat_pred <- fit_cf_model_mean(dat_train, dat_pred, family, cf_model)
          # compute score
          score_pred <- rep(0, nrow(dat_pred))
          score_pred[dat_pred$y == 1] <- -as.integer((1-yhat_pred[dat_pred$y == 1]) < yhat_pred[dat_pred$y == 1])
          score_pred[dat_pred$y == 0] <- -as.integer(yhat_pred[dat_pred$y == 0] < (1-yhat_pred[dat_pred$y == 0]))
        } else if (cf_score == "BCE"){
          # binary cross-entropy loss
          yhat_pred <- fit_cf_model_mean(dat_train, dat_pred, family, cf_model)
          # compute score
          score_pred <- -(y[id_pred]*log(yhat_pred)+(1-y[id_pred])*log(1-yhat_pred))
        } else if (cf_score == "HL"){
          # hinge loss
          yhat_pred <- fit_cf_model_mean(dat_train, dat_pred, family, cf_model)
          # compute score
          score_pred <- pmax(0,1-y[id_pred]*log(yhat_pred/(1-yhat_pred))+(1-y[id_pred])*log(yhat_pred/(1-yhat_pred)))
        } else if (cf_score == "SHL"){
          # hinge loss
          yhat_pred <- fit_cf_model_mean(dat_train, dat_pred, family, cf_model)
          # compute score
          score_pred <- (pmax(0,1-y[id_pred]*log(yhat_pred/(1-yhat_pred))+(1-y[id_pred])*log(yhat_pred/(1-yhat_pred))))^2
        } else if (cf_score == "MSE"){
          # Mean Squared Error
          yhat_pred <- fit_cf_model_mean(dat_train, dat_pred, family, cf_model)
          # compute score
          score_pred <- (y[id_pred]-yhat_pred)^2
        } else if (cf_score == "EL"){
          # exponential loss
          yhat_pred <- fit_cf_model_mean(dat_train, dat_pred, family, cf_model)
          # compute score
          score_pred <- exp(-(y[id_pred]-(1-y[id_pred]))*log(yhat_pred/(1-yhat_pred)))
        }
        score_test_k <- head(score_pred, length(id_test))
        score_cal_k <- tail(score_pred, length(id_cal))
        # compare
        map_dbl(score_test_k, function(score_test_l) {
          sum(score_cal_k >= score_test_l)
        })
      }) %>% sapply(function(x) x)
      # conformal p-values
      p_cf <- map_dbl(1:length(id_test), function(l) {
        (sum(compare_all[l,]) + 1) / (length(id_rc) + 1)
      })
    }
  }
  p_cf
}

ec_borrow <- function(
    Y, A, S, X, method, 
    family = "gaussian", 
    n_fisher = NULL, 
    # OM/IPW/CW
    n_boot = NULL,
    # conformal selective borrowing
    gamma_sel = NULL,
    cf = "cv+", # c("split", "full", "jackknife+", "cv+")
    cf_score = "AR",
    cf_model = "glm",
    split_train = 0.75,
    cv_fold = 10,
    sd_est = "model",
    n_rc_min = 3,
    # AIPW
    outcome_model = "glm",
    max_r = Inf,
    # testing
    sig_level = 0.05,
    small_n_adj = TRUE, 
    # computing & output
    parallel = FALSE, 
    n_cores = detectCores(logical = FALSE),
    output_frt = FALSE
) {
  
  dat_origin <- tibble(Y, A, S, X)
  rm(Y, A, S, X)
  
  # No Borrowing
  #   No Borrow DiM
  #   No Borrow AIPW
  
  # Borrowing
  #   Borrow OM
  #   Borrow IPW, Borrow staIPW, Borrow CW
  #   Borrow AIPW, Borrow ACW
  
  # Selective Borrowing
  #   AdaLasso Selective Borrow ACW (Chenyin)
  #   Conformal Selective Borrow AIPW (proposed) 
  #   Conformal Selective Borrow ACW (proposed) 
  
  if (identical(method, "No Borrow DiM")) {
    est_fun <- function(dat) {
      fit <- t.test(
        x = dat %>% filter(A == 1, S == 1) %>% pull(Y),
        y = dat %>% filter(A == 0, S == 1) %>% pull(Y)
      )
      
      mean_d1 <- fit$estimate[1] %>% unname
      mean_d0 <- fit$estimate[2] %>% unname
      d_rd <- mean_d1 - mean_d0
      d_rr <- mean_d1 / mean_d0
      d_or <- (mean_d1/(1-mean_d1)) / (mean_d0/(1-mean_d0))
      est <- c(d_rd, d_rr, d_or)
      
      tibble(
        # est = -(fit$estimate %>% diff %>% unname), 
        est = est,
        # se = fit$stderr, # use bootstrap
        ci_l = fit$conf.int[1],
        ci_u = fit$conf.int[2],
        p_value = fit$p.value,
        # borrow no EC
        ess_sel = 0,
        id_sel = list(NULL)
      )
    }
    gamma_sel <- 1
  }
  
  if (identical(method, "No Borrow AIPW")) {
    est_fun <- function(dat) {
      rct_aipw(dat, family, outcome_model, small_n_adj) %>% 
        # borrow no EC
        mutate(
          ess_sel = 0,
          id_sel = list(NULL)
        )
    }
    gamma_sel <- 1
  }
  
  if (identical(method, "Borrow Naive")) {
    est_fun <- function(dat) {
      n_ec <- dat %>% filter(A == 0, S == 0) %>% nrow
      dat_naive <- dat %>% mutate(S = 1)
      rct_aipw(dat_naive, family, outcome_model, small_n_adj) %>% 
        # borrow no EC
        mutate(
          ess_sel = n_ec,
          id_sel = list(which(dat$S == 0))
        )
    }
    gamma_sel <- 0
  }
  
  if (identical(method, "Borrow OM")) {
    est_fun <- function(dat) {
      n_ec <- dat %>% filter(A == 0, S == 0) %>% nrow
      
      m10 <- fit_outcome_model(dat, family, outcome_model)
      m1 <- m10$m1[dat$S == 1]
      m0 <- m10$m0[dat$S == 1]
      
      # m1 <- m10$m1
      # m0 <- m10$m0
      # d <- dat %>% 
      #   mutate(d_i = m1 - m0) %>% 
      #   filter(S == 1) %>% 
      #   pull(d_i)
      
      # compute est
      d_rd <- mean(m1) - mean(m0)
      d_rr <- mean(m1) / mean(m0)
      d_or <- (mean(m1)/(1-mean(m1))) / (mean(m0)/(1-mean(m0)))
      est <- c(d_rd, d_rr, d_or)
      
      tibble(
        # est = mean(d), 
        est = est,
        # borrow all ECs
        ess_sel = n_ec,
        id_sel = list(which(dat$S == 0))
      )
    }
    gamma_sel <- 0
  }
  
  if (identical(method, "Borrow IPW")) {
    est_fun <- function(dat) {
      n_rt <- dat %>% filter(A == 1, S == 1) %>% nrow
      n_rc <- dat %>% filter(A == 0, S == 1) %>% nrow
      n_rct <- dat %>% filter(S == 1) %>% nrow
      n_all <- dat %>% nrow
      
      # treatment group
      # use true propensity score
      pA <- n_rt / n_rct 
      w1 <- with(
        dat, 
        S * A / pA
      )
      d1 <- with(
        dat,
        (n_all / n_rct) * w1 * Y
      )
      
      # control group
      # compute r
      if (family == "gaussian") {
        r1 <- glm(Y ~ X, family = family, dat %>% filter(A == 0, S == 1)) %>% 
          resid(type = "response") %>% var
        r0 <- glm(Y ~ X, family = family, dat %>% filter(S == 0)) %>% 
          resid(type = "response") %>% var
        r <- min(r1 / r0, max_r)
      } else if (family == "binomial") {
        # for binary outcome, under exchangeablity assumption, r=1 (Li et al., 2023)
        r <- 1
      }
      # compute qhat
      pS <- glm(S ~ X, family = "binomial", dat) %>% predict(dat, "response")
      qhat <- pS / (1 - pS)
      w0 <- with(
        dat, 
        qhat * (S * (1 - A) + (1 - S) * r) / (qhat * (1 - pA) + r)
      )
      d0 <- with(
        dat,
        (n_all / n_rct) * w0 * Y
      )
      
      # compute est
      d_rd <- mean(d1) - mean(d0)
      d_rr <- mean(d1) / mean(d0)
      d_or <- (mean(d1)/(1-mean(d1))) / (mean(d0)/(1-mean(d0)))
      est <- c(d_rd, d_rr, d_or)
      
      # output
      tibble(
        # est = mean(d), 
        est = est,
        # borrow all ECs
        ess_sel = max(0, ESS(w0) - n_rc),
        id_sel = list(which(dat$S == 0))
      )
    }
    gamma_sel <- 0
  }
  
  if (identical(method, "Borrow staIPW")) {
    est_fun <- function(dat) {
      n_rt <- dat %>% filter(A == 1, S == 1) %>% nrow
      n_rc <- dat %>% filter(A == 0, S == 1) %>% nrow
      n_rct <- dat %>% filter(S == 1) %>% nrow
      n_all <- dat %>% nrow
      
      # treatment group
      # use true propensity score
      pA <- n_rt / n_rct 
      w1 <- with(
        dat, 
        S * A / pA
      )
      d1 <- with(
        dat,
        (n_all / n_rct) * w1 * Y
      )
      
      # control group
      # compute r
      if (family == "gaussian") {
        r1 <- glm(Y ~ X, family = family, dat %>% filter(A == 0, S == 1)) %>% 
          resid(type = "response") %>% var
        r0 <- glm(Y ~ X, family = family, dat %>% filter(S == 0)) %>% 
          resid(type = "response") %>% var
        r <- min(r1 / r0, max_r)
      } else if (family == "binomial") {
        # for binary outcome, under exchangeablity assumption, r=1 (Li et al., 2023)
        r <- 1
      }
      # compute qhat
      pS <- glm(S ~ X, family = "binomial", dat) %>% predict(dat, "response")
      qhat <- pS / (1 - pS)
      w0init <- with(
        dat, 
        qhat * (S * (1 - A) + (1 - S) * r) / (qhat * (1 - pA) + r)
      )
      w0 <- w0init / sum(w0init) * n_rct
      d0 <- with(
        dat,
        (n_all / n_rct) * w0 * Y
      )
      
      # compute est
      d_rd <- mean(d1) - mean(d0)
      d_rr <- mean(d1) / mean(d0)
      d_or <- (mean(d1)/(1-mean(d1))) / (mean(d0)/(1-mean(d0)))
      est <- c(d_rd, d_rr, d_or)
      
      # output
      tibble(
        # est = mean(d),
        est = est,
        # borrow all ECs
        ess_sel = max(0, ESS(w0) - n_rc),
        id_sel = list(which(dat$S == 0))
      )
    }
    gamma_sel <- 0
  }
  
  if (identical(method, "Borrow CW")) {
    est_fun <- function(dat) {
      n_rt <- dat %>% filter(A == 1, S == 1) %>% nrow
      n_rc <- dat %>% filter(A == 0, S == 1) %>% nrow
      n_rct <- dat %>% filter(S == 1) %>% nrow
      n_all <- dat %>% nrow
      
      # treatment group
      # use true propensity score
      pA <- n_rt / n_rct 
      w1 <- with(
        dat, 
        S * A / pA
      )
      d1 <- with(
        dat,
        (n_all / n_rct) * w1 * Y
      )
      
      # control group
      # compute r
      if (family == "gaussian") {
        r1 <- glm(Y ~ X, family = family, dat %>% filter(A == 0, S == 1)) %>% 
          resid(type = "response") %>% var
        r0 <- glm(Y ~ X, family = family, dat %>% filter(S == 0)) %>% 
          resid(type = "response") %>% var
        r <- min(r1 / r0, max_r)
      } else if (family == "binomial") {
        # for binary outcome, under exchangeablity assumption, r=1 (Li et al., 2023)
        r <- 1
      }
      # compute qhat
      qhat <- compute_cw(dat$S, dat$X)
      w0init <- with(
        dat, 
        qhat * (S * (1 - A) + (1 - S) * r) / (qhat * (1 - pA) + r)
      )
      w0 <- w0init / sum(w0init) * n_rct
      d0 <- with(
        dat,
        (n_all / n_rct) * w0 * Y
      )
      
      # compute est
      d_rd <- mean(d1) - mean(d0)
      d_rr <- mean(d1) / mean(d0)
      d_or <- (mean(d1)/(1-mean(d1))) / (mean(d0)/(1-mean(d0)))
      est <- c(d_rd, d_rr, d_or)
      
      # output
      tibble(
        # est = mean(d), 
        est = est,
        # borrow all ECs
        ess_sel = max(0, ESS(w0) - n_rc),
        id_sel = list(which(dat$S == 0))
      )
    }
    gamma_sel <- 0
  }
  
  if (identical(method, "Borrow AIPW")) {
    est_fun <- function(dat) {
       rct_ec_aipw_acw(dat, family, outcome_model, max_r, small_n_adj) %>% 
        # borrow all ECs
        mutate(id_sel = list(which(dat$S == 0))) 
    }
    gamma_sel <- 0
  }
  
  if (identical(method, "Borrow ACW")) {
    est_fun <- function(dat) {
      rct_ec_aipw_acw(dat, family, outcome_model, max_r, small_n_adj, cw = TRUE) %>% 
        # borrow all ECs
        mutate(id_sel = list(which(dat$S == 0))) 
    }
    gamma_sel <- 0
  }
  
  if (identical(method, "AdaLasso Selective Borrow ACW")) { # Chenyin's method
    est_fun <- function(dat) {
      fit <- with(dat,
                  srEC(
                    data_rt = list(X = X[S == 1,], Y = Y[S == 1], A = A[S == 1]),
                    data_ec = list(list(X = X[S == 0,], Y = Y[S == 0], A = A[S == 0])),
                    method = "glm"
                  )
      )
      tibble(
        est = as.vector(fit$est$ACW.final), 
        se = as.vector(fit$sd$ACW.final) / sqrt(fit$n_c),
        ess_sel = NA,
        id_sel = list(which(dat$S == 0)[fit$subset.idx])
      )
    }
    gamma_sel <- NA
  }
  
  if (identical(method, c("Conformal Selective Borrow AIPW"))) { # proposed method
    # est_fun
    # version 1
    n_rc_y1 <- sum(dat_origin$A == 0 & dat_origin$S == 1 & dat_origin$Y == 1)
    n_rc_y0 <- sum(dat_origin$A == 0 & dat_origin$S == 1 & dat_origin$Y == 0)
    
    if (n_rc_y0 >= n_rc_min & n_rc_y1 >= n_rc_min){
      est_fun <- function(dat) {
        x <- dat %>% filter(A == 0) %>% pull(X)
        y <- dat %>% filter(A == 0) %>% pull(Y)
        s <- dat %>% filter(A == 0) %>% pull(S)
        p_cf <- conformal_p(
          x, y, s, family,
          cf, cf_score, cf_model, split_train, cv_fold, sig_level,sd_est
        )
        # biased or unbiased
        bias_ec <- ifelse(p_cf >= gamma_sel, 0, 1)
        # estimation
        if (sum(bias_ec == 0) < 5) {
          # if n_sel < 5, do not borrow anyone
          rct_aipw(dat, family, outcome_model, small_n_adj) %>%
            mutate(
              ess_sel = 0,
              id_sel = list(NULL)
            )
        } else {
          # if n_sel >= 5, borrow them
          bias <- rep(0, nrow(dat))
          bias[dat$S == 0] <- bias_ec
          dat_sel <- dat %>% filter(bias == 0)
          rct_ec_aipw_acw(dat_sel, family, outcome_model, max_r, small_n_adj) %>%
            mutate(id_sel = list(which(dat$S == 0 & bias == 0)))
        }
      }
    } else {
      # if not enough, then No Borrow AIPW
      est_fun <- function(dat) {
        rct_aipw(dat, family, outcome_model, small_n_adj) %>%
          # borrow no EC
          mutate(
            ess_sel = 0,
            id_sel = list(NULL)
          )
      }
      gamma_sel <- 1
    }
    
    # version 2
    # est_fun <- function(dat) {
    #   n_rc_y1 <- sum(dat$A == 0 & dat$S == 1 & dat$Y == 1)
    #   n_rc_y0 <- sum(dat$A == 0 & dat$S == 1 & dat$Y == 0)
    #   
    #   if (family == "binomial" & n_rc_y1 < n_rc_min & n_rc_y0 < n_rc_min) {
    #     rct_aipw(dat, family, outcome_model, small_n_adj) %>%
    #       mutate(
    #         ess_sel = 0,
    #         id_sel = list(NULL)
    #       )
    #   } else {
    #     x <- dat %>% filter(A == 0) %>% pull(X)
    #     y <- dat %>% filter(A == 0) %>% pull(Y)
    #     s <- dat %>% filter(A == 0) %>% pull(S)
    #     p_cf <- conformal_p(
    #       x, y, s, family,
    #       cf, cf_score, cf_model, split_train, cv_fold, sig_level
    #     )
    #     # biased or unbiased
    #     bias_ec <- ifelse(p_cf >= gamma_sel, 0, 1)
    #     # estimation
    #     if (sum(bias_ec == 0) < 5) {
    #       # if n_sel < 5, do not borrow anyone
    #       rct_aipw(dat, family, outcome_model, small_n_adj) %>%
    #         mutate(
    #           ess_sel = 0,
    #           id_sel = list(NULL)
    #         )
    #     } else {
    #       # if n_sel >= 5, borrow them
    #       bias <- rep(0, nrow(dat))
    #       bias[dat$S == 0] <- bias_ec
    #       dat_sel <- dat %>% filter(bias == 0)
    #       rct_ec_aipw_acw(dat_sel, family, outcome_model, max_r, small_n_adj) %>%
    #         mutate(id_sel = list(which(dat$S == 0 & bias == 0)))
    #     }
    #   }
    # }
  }
  
  if (identical(method, c("Conformal Selective Borrow ACW"))) { # proposed method
    # est_fun
    est_fun <- function(dat) {
      x <- dat %>% filter(A == 0) %>% pull(X)
      y <- dat %>% filter(A == 0) %>% pull(Y)
      s <- dat %>% filter(A == 0) %>% pull(S)
      p_cf <- conformal_p(
        x, y, s, family,
        cf, cf_score, cf_model, split_train, cv_fold, sig_level
      )
      # biased or unbiased
      bias_ec <- ifelse(p_cf >= gamma_sel, 0, 1)
      # estimation
      if (sum(bias_ec == 0) < 5) {
        # if n_sel < 5, do not borrow anyone
        rct_aipw(dat, family, outcome_model, small_n_adj) %>% 
          mutate(
            ess_sel = 0,
            id_sel = list(NULL)
          )
      } else {
        # if n_sel >= 5, borrow them
        bias <- rep(0, nrow(dat))
        bias[dat$S == 0] <- bias_ec
        dat_sel <- dat %>% filter(bias == 0)
        rct_ec_aipw_acw(dat_sel, family, outcome_model, max_r, small_n_adj, cw = TRUE) %>% 
          mutate(id_sel = list(which(dat$S == 0 & bias == 0)))
      }
    }
  }
  
  # 1 Estimation
  # record run time
  runtime <- system.time(
    # raw output
    out <- est_fun(dat_origin) 
  )[3] %>% unname
  
  # 2 Inference
  if (method %in% c("Borrow OM", "Borrow IPW", "Borrow staIPW", "Borrow CW", "No Borrow DiM")) {
    if (!is.null(n_boot)) {
      # compute bootstrap SE
      runtime_boot <- system.time(
        # bootstrap
        if (parallel) {
          cat(paste0("parallel computing enabled with ", n_cores, 
                     " cores for bootstrap SE of ", method, "\n\n"))
          # five estimands
          out_boot <- mclapply(1:n_boot, function(i) {
            dat_boot <- dat_origin %>%
              group_by(A, S) %>%
              slice_sample(prop = 1, replace = TRUE) %>%
              ungroup()

            tryCatch({
              est_fun(dat_boot)$est  # Returns a vector of length 5
            }, error = function(e) {
              rep(NA, 5)  # Ensure error handling returns a consistent vector of length 5
            })
          }, mc.cores = n_cores) %>%
            do.call(rbind, .)
          
          # only one estimand
          # out_boot <- mclapply(1:n_boot, function(i) {
          #   dat_boot <- dat_origin %>%
          #     group_by(A, S) %>%
          #     slice_sample(prop = 1, replace = TRUE) %>%
          #     ungroup()
          #   tryCatch({
          #     est_fun(dat_boot)$est
          #   }, error = function(e) {
          #     NA
          #   })
          # }, mc.cores = n_cores) %>%
          #   map_dbl(~.)
        } else {
          # five estimands
          out_boot <- map(1:n_boot, ~ {
            dat_boot <- dat_origin %>%
              group_by(A, S) %>%
              slice_sample(prop = 1, replace = TRUE) %>%
              ungroup()

            tryCatch({
              est_fun(dat_boot)$est  # Returns a vector of length 5
            }, error = function(e) {
              rep(NA, 5)  # Ensure error handling returns a vector of length 5
            })
          }) %>%
            do.call(rbind, .)
          
          # only one estimand
          # out_boot <- map_dbl(1:n_boot, ~ {
          #   dat_boot <- dat_origin %>%
          #     group_by(A, S) %>%
          #     slice_sample(prop = 1, replace = TRUE) %>%
          #     ungroup()
          #   tryCatch({
          #     est_fun(dat_boot)$est
          #   }, error = function(e) {
          #     NA
          #   })
          # })
        }
      )[3] %>% unname
      # five estimands
      out$se <- apply(out_boot, 2, function(x) sd(x, na.rm = TRUE))
      ## only one estimand
      # out$se <- sd(out_boot, na.rm = TRUE)
      
      # warning for NA
      if (any(is.na(out_boot))) {
        warning(paste0("There are ", sum(is.na(out_boot)), 
                       " NA in bootstrap for ", method))
      }
      # organize
      res <- tibble(
        method = paste0(method, c("_RD","_RR","_OR")),
        # method = method,
        est = out$est, 
        se = out$se, 
        ci_l = est - qnorm(1 - sig_level / 2) * se, 
        ci_u = est + qnorm(1 - sig_level / 2) * se,
        p_value = (1 - pnorm(abs(est / se))) * 2,
        p_value_upper = 1 - pnorm(est / se),
        p_value_lower = pnorm(est / se),
        n_sel = map_dbl(out$id_sel, length),
        ess_sel = out$ess_sel,
        runtime = runtime + runtime_boot
      )
    } else {
      # organize
      res <- tibble(
        method = paste0(method, c("_RD","_RR","_OR")),
        # method = method,
        est = out$est, 
        se = NA, 
        ci_l = NA, 
        ci_u = NA,
        p_value = NA,
        p_value_upper = NA,
        p_value_lower = NA,
        n_sel = map_dbl(out$id_sel, length),
        ess_sel = out$ess_sel,
        runtime = runtime
      )
    }
  } else {
    # for No Borrow AIPW; No Borrow DiM; Borrow ACW; AdaLasso Selective Borrow ACW; Conformal Selective Borrow AIPW
    # organize
    res <- tibble(
      method = paste0(method, c("_RD","_RR","_OR")),
      # method = method,
      est = out$est, 
      se = out$se, 
      ci_l = est - qnorm(1 - sig_level / 2) * se, 
      ci_u = est + qnorm(1 - sig_level / 2) * se,
      p_value = (1 - pnorm(abs(est / se))) * 2,
      p_value_upper = 1 - pnorm(est / se),
      p_value_lower = pnorm(est / se),
      n_sel = map_dbl(out$id_sel, length),
      ess_sel = out$ess_sel,
      runtime
    )
  }
  
  # 3 Fisher randomization test
  if (!is.null(n_fisher)) {
    if (identical(method, "No Borrow DiM") & family == "binomial") {
      # special case of Fisher's exact test
      x1 <- dat_origin %>% filter(S == 1, A == 0, Y == 1) %>% nrow()
      n1 <- dat_origin %>% filter(S == 1, A == 0) %>% nrow()
      x2 <- dat_origin %>% filter(S == 1, A == 1, Y == 1) %>% nrow()
      n2 <- dat_origin %>% filter(S == 1, A == 1) %>% nrow()
      runtime_frt <- system.time(
        fit <- matrix(
          c(x2, x1, n2 - x2, n1 - x1), 2, 2,
          dimnames = list(c("t", "c"), c("Event", "No Event"))
        ) %>% fisher.test()
      )[3] %>% unname
      
      fit_upper <- matrix(
        c(x2, x1, n2 - x2, n1 - x1), 2, 2,
        dimnames = list(c("t", "c"), c("Event", "No Event"))
      ) %>% fisher.test(alternative = "greater")
      fit_lower <- matrix(
        c(x2, x1, n2 - x2, n1 - x1), 2, 2,
        dimnames = list(c("t", "c"), c("Event", "No Event"))
      ) %>% fisher.test(alternative = "less")
      
      res_frt <- tibble(
        method = paste0(method, "+FRT"),
        est = NA, 
        se = NA, 
        ci_l = NA, 
        ci_u = NA,
        p_value = fit$p.value,
        p_value_upper = fit_upper$p.value, # need revision
        p_value_lower = fit_lower$p.value, # need revision
        n_sel = NA,
        ess_sel = NA,
        runtime = runtime_frt
      )
    } else {
      runtime_frt <- system.time(
        if (parallel) {
          cat(paste0("parallel computing enabled with ", n_cores, 
                     " cores for ", method, "+FRT\n\n"))
          # randomization
          out_frt <- mclapply(1:n_fisher, function(i) {
            dat_rand <- dat_origin %>% 
              mutate(A = {A[S == 1] <- sample(A[S == 1]); A})
            tryCatch({
              est_fun(dat_rand)
            }, error = function(e) {
              NULL
            })
          }, mc.cores = n_cores) %>% 
            map_dfr(~.) %>% 
            mutate(
              cond = floor(map_dbl(id_sel, length) / 10) == floor(res$n_sel / 10)
            )
        } else {
          # randomization
          out_frt <- map_dfr(1:n_fisher, ~{
            dat_rand <- dat_origin %>% 
              mutate(A = {A[S == 1] <- sample(A[S == 1]); A})
            tryCatch({
              est_fun(dat_rand)
            }, error = function(e) {
              NULL
            })
          })
        }
      )[3] %>% unname
      
      # warning for NA
      if (nrow(out_frt) < n_fisher) {
        warning(paste0("nrow(out_frt) is ", nrow(out_frt), " for ",
                       method, "+FRT"))
      }
      if (any(is.na(out_frt$est))) {
        warning(paste0("There are ", sum(is.na(out_frt$est)), " NA in ",
                       method, "+FRT"))
      }
      
      # prepare for p-value (five estimands)
      if (nrow(out) == 3){
        out_rep <- out[rep(seq_len(nrow(out)), times=nrow(out_frt)/3), ]
        out_frt$id <- rep(seq(1,3), times=nrow(out_frt)/3)
        out_frt$count <- ifelse(out_frt$id %in% c(2, 3), 
                                (abs(out_frt$est - 1) >= abs(out_rep$est - 1)), # RR and OR
                                (abs(out_frt$est) >= abs(out_rep$est))
                                ) # RD
        
        out_frt$count_one_upper <- (out_frt$est - out_rep$est >= 0) # H1: tau_RD > 0; tau_RR/OR > 1
        out_frt$count_one_lower <- (out_frt$est - out_rep$est <= 0) # H1: tau_RD < 0; tau_RR/OR < 1
        
        frtp_t <- out_frt %>% group_by(id) %>% summarise(prob = mean(c(count,1),na.rm = TRUE))
        frtp <- frtp_t$prob
        
        frtp_tu <- out_frt %>% group_by(id) %>% summarise(prob = mean(c(count_one_upper,1),na.rm = TRUE))
        frtp_u <- frtp_tu$prob
        
        frtp_tl <- out_frt %>% group_by(id) %>% summarise(prob = mean(c(count_one_lower,1),na.rm = TRUE))
        frtp_l <- frtp_tl$prob
        
      }else{
        frtp <- mean(c(abs(out_frt$est) >= abs(out$est), 1), na.rm = T)
        frtp_u <- mean(c(out_frt$est >= out$est, 1), na.rm = T)
        frtp_l <- mean(c(out_frt$est <= out$est, 1), na.rm = T)
      }
      
      # organize
      if (nrow(out_frt)==0){
        res_frt <- tibble(
          # method = paste0(method, "+FRT"),
          method = paste0(method, c("_RD","_RR","_OR"),"+FRT"),
          est = NA,
          se = NA,
          ci_l = NA,
          ci_u = NA,
          p_value = 1,
          p_value_upper = 1,
          p_value_lower = 1,
          n_sel = NA,
          ess_sel = NA,
          runtime = runtime_frt
        )
      }else{
        res_frt <- tibble(
          # method = paste0(method, "+FRT"),
          method = paste0(method, c("_RD","_RR","_OR"),"+FRT"),
          est = NA,
          se = NA,
          ci_l = NA,
          ci_u = NA,
          p_value = frtp,
          p_value_upper = frtp_u,
          p_value_lower = frtp_l,
          n_sel = NA,
          ess_sel = NA,
          runtime = runtime_frt
        )
      }
    }
    res <- rbind(res, res_frt)
  } else {
    out_frt <- NULL
  }
  dat_info <- tibble(
    n_rt = dat_origin %>% filter(A == 1, S == 1) %>% nrow,
    n_rc = dat_origin %>% filter(A == 0, S == 1) %>% nrow,
    n_rct = n_rt + n_rc,
    n_ec = dat_origin %>% filter(A == 0, S == 0) %>% nrow,
    id_ec = list(which(dat_origin$S == 0))
  )
  if (output_frt) {
    lst(res, out, dat_info, gamma_sel, out_frt)
  } else {
    lst(res, out, dat_info, gamma_sel)
  }
}
ESS <- function (w) {
  sum(w)^2 / sum(w^2)
}

add_name <- function(l, new) {
  l$res <- l$res %>% mutate(method = paste0(method, " ", new))
  l
}
compute_ada_gamma <- function(Y, A, S, X, family, 
                              gamma_grid = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                              measure = "mse_hat", n_rep_gamma = 100,
                              parallel = F, 
                              n_cores = detectCores(logical = FALSE),
                              ...) {
  if (!parallel) {n_cores <- 1}
  res_grid <- mclapply(gamma_grid, function(g) {
    dat_rct <- tibble(Y, A, S, X) %>% filter(S == 1)
    dat_ec <- tibble(Y, A, S, X) %>% filter(S == 0)
    n_rt1 <- floor(0.5 * sum(dat_rct$A == 1))
    n_rc1 <- floor(0.5 * sum(dat_rct$A == 0))
    est_rep <- map(1:n_rep_gamma, ~ {
      fold1_id <- c(
        sample(which(dat_rct$A == 1), size = n_rt1),
        sample(which(dat_rct$A == 0), size = n_rc1)
      )
      dat1 <- bind_rows(dat_rct[fold1_id,], dat_ec)
      dat2 <- bind_rows(dat_rct[-fold1_id,], dat_ec)
      est_csb <- ec_borrow(
        dat1$Y, dat1$A, dat1$S, dat1$X, 
        "Conformal Selective Borrow AIPW", family, n_fisher = NULL,
        gamma_sel = g, ...
      )$res$est[1]
      est_nb <- ec_borrow(
        dat2$Y, dat2$A, dat2$S, dat2$X,
        "No Borrow AIPW", family, n_fisher = NULL
      )$res$est[1]
      lst(est_csb, est_nb)
    })
    var_hat <- map_dbl(est_rep, ~ {.$est_csb}) %>% var
    bias2_hat <- mean(map_dbl(est_rep, ~ {.$est_csb - .$est_nb}), na.rm = T)^2
    mse_hat <- var_hat + bias2_hat
    # output
    cat(paste0("For gamma_sel = ", g, 
               ", MSE = ", mse_hat, "\n\n"))
    lst(mse_hat, var_hat, bias2_hat)
  }, mc.cores = n_cores)
  gamma_grid[which.min(map_dbl(res_grid, measure))]
}

