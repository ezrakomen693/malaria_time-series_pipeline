# TIME SERIES PIPELINE: Monthly Malaria Cases (Kenya, 2013-2022)
# Data: Simulated, KHIS/DHIS2-inspired (NOT real surveillance data)
# Pipeline: simulate -> diagnose -> model/validate -> check -> forecast
# Estimand: immediate level change in monthly cases after bed-net
#           distribution (interrupted time series, ITS), true value known
#           because we simulate it.

suppressPackageStartupMessages({
  library(forecast)   # auto.arima, ets, snaive, accuracy, ggAcf, mstl
  library(tseries)    # adf.test, kpss.test
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(lubridate)
})

# CONFIG: every assumption lives here, nothing is hard-coded below
CONFIG <- list(
  seed                = 123,
  n                   = 120,                  # Jan 2013 - Dec 2022
  start               = as.Date("2013-01-01"),
  base                = 100,                  # baseline monthly cases
  slope               = 0.5,                  # cases/month linear trend
  season_amp          = 15,                   # annual seasonal amplitude
  season_peak_month   = 6,                    # assumed peak (post long rains)
  ar_phi              = 0.5,                  # AR(1) in noise (0 = iid)
  noise_sd            = 6,                    # innovation SD
  intervention_t      = 61,                   # first treated month = Jan 2018
  intervention_effect = -8,                   # true level shift
  holdout             = 12,                   # months held out for validation
  horizon             = 12,                   # months forecast ahead
  fig_dir             = "figures"
)

#  SIMULATE
simulate_malaria <- function(cfg) {
  set.seed(cfg$seed)
  t     <- seq_len(cfg$n)
  dates <- seq(cfg$start, by = "month", length.out = cfg$n)
  mon   <- lubridate::month(dates)

  trend  <- cfg$slope * t
  season <- cfg$season_amp * cos(2 * pi * (mon - cfg$season_peak_month) / 12)
  step   <- as.integer(t >= cfg$intervention_t)
  noise  <- as.numeric(arima.sim(model = list(ar = cfg$ar_phi),
                                 n = cfg$n, sd = cfg$noise_sd))

  # Gaussian can go negative; counts cannot. Clip and round.
  cases <- pmax(0, round(cfg$base + trend + season +
                           cfg$intervention_effect * step + noise))

  tibble(date = dates, t = t, month = mon, step = step,
         trend = trend, season = season, noise = noise, cases = cases)
}

to_ts <- function(df) {
  ts(df$cases,
     start     = c(lubridate::year(df$date[1]), lubridate::month(df$date[1])),
     frequency = 12)
}

# xreg must be a named matrix so train/forecast columns match
make_xreg <- function(step) {
  matrix(step, ncol = 1, dimnames = list(NULL, "step"))
}

# DIAGNOSE (before modelling)
diagnose_series <- function(y, cfg) {
  p_series <- autoplot(y) +
    geom_vline(xintercept = time(y)[cfg$intervention_t],
               linetype = "dashed", colour = "firebrick") +
    labs(title = "Monthly malaria cases (simulated)",
         subtitle = "Dashed line = bed-net distribution",
         x = NULL, y = "Cases")

  p_stl  <- autoplot(mstl(y)) + labs(title = "STL decomposition")
  p_acf  <- ggAcf(y,  lag.max = 36) + labs(title = "ACF")
  p_pacf <- ggPacf(y, lag.max = 36) + labs(title = "PACF")

  adf  <- tseries::adf.test(y)
  kpss <- tseries::kpss.test(y, null = "Trend")

  tests <- tibble(
    test      = c("ADF (H0: unit root)", "KPSS (H0: trend-stationary)"),
    statistic = c(unname(adf$statistic), unname(kpss$statistic)),
    p_value   = c(adf$p.value, kpss$p.value)
  )

  list(
    plots = (p_series / p_stl) | (p_acf / p_pacf),
    tests = tests,
    nsdiffs = nsdiffs(y),                      # seasonal differences suggested
    ndiffs  = ndiffs(y)                        # regular differences suggested
  )
}

# MODEL + VALIDATE (holdout on the last `holdout` months)
compare_models <- function(y, step, cfg) {
  n       <- length(y)
  n_train <- n - cfg$holdout

  y_train <- subset(y, end = n_train)
  y_test  <- subset(y, start = n_train + 1)
  x_train <- make_xreg(step[seq_len(n_train)])
  x_test  <- make_xreg(step[(n_train + 1):n])
  h       <- cfg$holdout

  fits <- list(
    snaive     = snaive(y_train, h = h),                     # benchmark
    ets        = forecast(ets(y_train), h = h),
    arima      = forecast(auto.arima(y_train), h = h),
    arima_its  = forecast(auto.arima(y_train, xreg = x_train),
                          h = h, xreg = x_test)
  )

  scores <- bind_rows(lapply(names(fits), function(m) {
    a <- accuracy(fits[[m]], y_test)["Test set", c("RMSE", "MAE", "MASE")]
    tibble(model = m, RMSE = a[["RMSE"]], MAE = a[["MAE"]], MASE = a[["MASE"]])
  })) %>% arrange(RMSE)

  list(scores = scores, forecasts = fits)
}

# RESIDUAL CHECKS
check_residuals <- function(fit, lag = 24) {
  res    <- residuals(fit)
  fitdf  <- sum(fit$arma[1:4])                 # p + q + P + Q
  lb     <- Box.test(res,   lag = lag, type = "Ljung-Box", fitdf = fitdf)
  lb_sq  <- Box.test(res^2, lag = lag, type = "Ljung-Box")  # ARCH-type signal
  sw     <- shapiro.test(as.numeric(res))

  tibble(
    check   = c("Ljung-Box (residual autocorrelation)",
                "Ljung-Box on squared residuals (volatility clustering)",
                "Shapiro-Wilk (normality)"),
    p_value = c(lb$p.value, lb_sq$p.value, sw$p.value)
  )
}

# FINAL FIT, EFFECT ESTIMATE, FORECAST
final_fit <- function(y, step, cfg) {
  fit <- auto.arima(y, xreg = make_xreg(step))

  est <- coef(fit)[["step"]]
  se  <- sqrt(diag(vcov(fit)))[["step"]]
  effect <- tibble(
    estimate = est,
    ci_low   = est - 1.96 * se,
    ci_high  = est + 1.96 * se,
    truth    = cfg$intervention_effect
  )

  # Assumes the intervention effect persists over the forecast horizon
  future_x <- make_xreg(rep(1L, cfg$horizon))
  fc <- forecast(fit, h = cfg$horizon, xreg = future_x)

  p_fc <- autoplot(fc) +
    geom_vline(xintercept = time(y)[cfg$intervention_t],
               linetype = "dashed", colour = "firebrick") +
    labs(title = "12-month forecast (ARIMA with intervention regressor)",
         x = NULL, y = "Cases")

  list(fit = fit, effect = effect, forecast = fc, plot = p_fc)
}

# RUN
run_pipeline <- function(cfg = CONFIG) {
  dir.create(cfg$fig_dir, showWarnings = FALSE)

  df <- simulate_malaria(cfg)
  y  <- to_ts(df)

  diag <- diagnose_series(y, cfg)
  cmp  <- compare_models(y, df$step, cfg)
  fin  <- final_fit(y, df$step, cfg)
  resid_checks <- check_residuals(fin$fit)

  ggsave(file.path(cfg$fig_dir, "01_diagnostics.png"), diag$plots,
         width = 12, height = 8, dpi = 300)
  ggsave(file.path(cfg$fig_dir, "02_forecast.png"), fin$plot,
         width = 8, height = 5, dpi = 300)

  cat("\n--- Stationarity tests ---\n");      print(diag$tests)
  cat("\nSuggested differencing: d =", diag$ndiffs,
      "| D =", diag$nsdiffs, "\n")
  cat("\n--- Holdout accuracy (last", cfg$holdout, "months) ---\n")
  print(cmp$scores)
  cat("\n--- Final model ---\n");             print(fin$fit)
  cat("\n--- Intervention effect (estimate vs truth) ---\n"); print(fin$effect)
  cat("\n--- Residual checks ---\n");         print(resid_checks)

  invisible(list(data = df, ts = y, diagnostics = diag, comparison = cmp,
                 final = fin, residual_checks = resid_checks))
}

results <- run_pipeline()
results$data              # simulated data frame (includes true components)
results$diagnostics$tests # ADF and KPSS
results$comparison$scores # holdout accuracy table
results$final$effect      # intervention estimate vs truth
results$residual_checks   # diagnostic p-values
