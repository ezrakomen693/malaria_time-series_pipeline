# Malaria Time Series Pipeline: Simulated Monthly Cases, Kenya 2013-2022

An R teaching project that walks through a complete time series workflow on **simulated** malaria surveillance data, including an **interrupted time series (ITS)** analysis of a bed-net distribution.

> **Data note:** the data are simulated to resemble the structure of KHIS/DHIS2 monthly facility reporting. They are **not** real surveillance data, and no result here is a statement about malaria in Kenya. The point of simulating is that we know the true answer, so we can check whether our methods recover it.

**Author:** Ezra Komen Kipyegon, BSc Biostatistics, JKUAT

---

## 1. Who this is for

Third-year biostatistics or epidemiology students who have met regression, basic probability, and the idea of a time series (trend, seasonality, autocorrelation), and want to see them used together in one reproducible script.

**Helpful background:** R basics (functions, data frames, `ggplot2`), linear regression, the meaning of ACF/PACF. You do not need prior ARIMA experience.

## 2. Learning objectives

By the end you should be able to:

1. Simulate a time series from explicit components (trend + seasonality + intervention + autocorrelated noise).
2. Decompose a series with STL and read ACF/PACF plots.
3. Explain what ADF and KPSS test, and why a structural break weakens them.
4. Fit and compare seasonal naive, ETS, ARIMA, and ARIMA-with-regressor models using a **holdout**, not just AIC.
5. Estimate an intervention effect as a regression coefficient with ARIMA errors, and compare it with the known truth.
6. Run residual diagnostics and say what a failed check would mean.
7. State clearly what an ITS design can and cannot say about causation.

## 3. The research question

**Objective:** estimate the immediate change in the level of monthly malaria cases after a bed-net distribution.

**Estimand:** the step change in mean monthly cases at the intervention month, holding the trend, seasonality, and error structure fixed. In the simulation the true value is **-8 cases per month**.

**Design:** single-group interrupted time series, with no control series.

## 4. Quick start

### Install packages (once)

```r
install.packages(c("forecast", "tseries", "dplyr", "tibble",
                   "ggplot2", "patchwork", "lubridate"))
```

### Run

```r
source("malaria_ts_pipeline.R")
```

The script runs `run_pipeline()` at the bottom. It prints results to the console, saves two figures to `figures/`, and returns everything invisibly as `results`:

```r
results$data              # simulated data frame (includes true components)
results$diagnostics$tests # ADF and KPSS
results$comparison$scores # holdout accuracy table
results$final$effect      # intervention estimate vs truth
results$residual_checks   # diagnostic p-values
```

## 5. How the pipeline is organised

| Stage | Function | What it does |
|---|---|---|
| 0 | `CONFIG` | Every assumption in one list: no magic numbers in the code |
| 1 | `simulate_malaria()` | Builds the series from named components |
| 2 | `diagnose_series()` | Plot, STL decomposition, ACF/PACF, ADF, KPSS, suggested differencing |
| 3 | `compare_models()` | Fits 4 models on the first 108 months, scores them on the last 12 |
| 4 | `check_residuals()` | Ljung-Box, Ljung-Box on squared residuals, Shapiro-Wilk |
| 5 | `final_fit()` | Refits on all data, extracts the effect with a 95% CI, forecasts 12 months |
| 6 | `run_pipeline()` | Runs everything in order and saves figures |

### The data-generating model

$$Y_t = \beta_0 + \beta_1 t + A\cos\!\left(\tfrac{2\pi(m_t - m^*)}{12}\right) + \delta\,\mathbb{1}(t \ge t_0) + \varepsilon_t,\qquad \varepsilon_t = \phi\,\varepsilon_{t-1} + u_t$$

with $u_t \sim N(0,\sigma^2)$. Then $Y_t$ is rounded and floored at 0 so it behaves like a count.

### Key settings in `CONFIG`

| Setting | Default | Meaning |
|---|---|---|
| `n` | 120 | Months (Jan 2013 - Dec 2022) |
| `base`, `slope` | 100, 0.5 | Baseline level and linear trend |
| `season_amp`, `season_peak_month` | 15, 6 | Annual seasonal amplitude and assumed peak month |
| `ar_phi` | 0.5 | Autocorrelation of the noise (0 = independent errors) |
| `noise_sd` | 6 | Innovation standard deviation |
| `intervention_t` | 61 | First treated month (Jan 2018) |
| `intervention_effect` | -8 | True level shift, our estimand |
| `holdout` | 12 | Months reserved for validation |

## 6. Methods in plain language

**STL decomposition** splits the series into trend, seasonal, and remainder parts. It is a look-before-you-model step, not a test.

**ADF and KPSS** have opposite nulls. ADF's null is a unit root (non-stationary); KPSS's null is stationarity around a trend. Reading them together is more informative than either alone. Both lose power when the series has a structural break (as ours does at the intervention), so treat them as supporting evidence, not as the final word on differencing.

**ETS vs ARIMA.** ETS models smooth level, trend, and season. ARIMA models the autocorrelation structure. Their AICc values are **not comparable**, which is why we compare on holdout error (RMSE, MAE, MASE).

**ARIMA with a regressor** (`arima_its`) is regression with ARIMA errors. The `step` coefficient is the level shift, and the ARIMA part absorbs autocorrelation so the standard error is not falsely small.

**MASE** scales error by a naive benchmark. Below 1 means better than that benchmark, so it is comparable across series.

## 7. Reading the output

Report it the way you would in a methods section:

- **Objective:** what you estimated, in one sentence.
- **Methods:** simulation design, models compared, validation scheme, diagnostics.
- **Results:** the effect estimate with its 95% CI and the truth, the best model on holdout, which residual checks passed.
- **Conclusion:** what the estimate supports, and what it does not.

**Residual checks, and what failure means:**

| Check | Small p-value suggests |
|---|---|
| Ljung-Box (residuals) | Leftover autocorrelation: the model is under-specified, so SEs are unreliable |
| Ljung-Box (squared residuals) | Volatility clustering: variance is not constant |
| Shapiro-Wilk | Non-normal residuals: prediction intervals may be miscalibrated |

With 120 points, a significant p-value can reflect a trivial departure. Always look at the plots as well.

## 8. Assumptions and limitations

- **Association, not causation.** With no control series, anything else that changed at the same time (rainfall, other campaigns, health-worker strikes, DHIS2 reporting-rate changes) would look like a bed-net effect. Real analyses need a comparison series or a strong argument that these are ruled out.
- **Known intervention date.** ITS assumes the change point is specified in advance. Searching for the break in the data is a different, riskier analysis.
- **Persistent effect.** The forecast assumes the shift lasts.
- **Gaussian counts.** Real case counts are overdispersed and non-negative; negative binomial models are the usual choice. The simulation is deliberately simple.
- **No missing data, outliers, or reporting drift.** Real KHIS data have all three.
- **Single-run result.** One seed gives one estimate. Whether the method is unbiased needs many replications (see exercise 3).

## 9. Exercises

1. **Break the assumption.** Set `ar_phi = 0` and then `ar_phi = 0.8`. How does the width of the effect's CI change, and why?
2. **Shift the date.** Set `intervention_t = 60`. Does anything change in the estimate? Which month is now the first treated month?
3. **Simulation study (research-style).** Wrap `simulate_malaria()` and `final_fit()` in a loop over 500 seeds. Compute the mean estimate (bias), the empirical SD, and the proportion of 95% CIs covering -8. Then repeat with a weaker effect (`-3`) and report how power falls.
4. **Misspecification.** Fit a plain linear regression of cases on `step` and `t` (no ARIMA errors). Compare its standard error with the ARIMA one. This is the commonest ITS mistake.
5. **Extension for a real task.** An analyst at a county health office has real DHIS2 monthly counts and a net campaign date. List what you would add to this pipeline before trusting its result (hint: overdispersion, reporting completeness, a control series, seasonal-by-year variation).

## 10. Troubleshooting

| Problem | Likely cause |
|---|---|
| `there is no package called ...` | Re-run the `install.packages()` line |
| `could not find function "ggAcf"` | `forecast` did not load; reinstall it |
| `figures/` is empty | Check your working directory with `getwd()` |
| Warnings from `adf.test` about p-value range | Expected when the statistic is outside the lookup table; read as "p < 0.01" or "p > 0.1" |

## 11. Where to go next

- Replace Gaussian noise with a **negative binomial** model and add **missing months, outliers, and reporting-rate drift**.
- Add a **control series** (a comparison district) for a controlled ITS.
- Model **rainfall or temperature** as external regressors.
- Try **segmented regression with Newey-West errors** and compare with ARIMA.

## 12. Reproducibility

Set seeds via `CONFIG$seed`, record `sessionInfo()` with any results you report, and keep the script and figures under version control.
