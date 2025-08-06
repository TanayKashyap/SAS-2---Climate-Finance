library(StMoMo)
library(stats)
library(ggplot2)
library(forecast)  # for time series forecasting

# --- Helper function to load mortality and exposure CSVs ---
load_data <- function(m_file, e_file) {
  deaths <- read.csv(m_file, header = TRUE, row.names = 1)
  exposure <- read.csv(e_file, header = TRUE, row.names = 1)
  ages <- as.numeric(rownames(deaths))
  years <- as.numeric(sub("X", "", colnames(deaths)))
  Dxt <- as.matrix(sapply(deaths, as.numeric))
  Ext <- as.matrix(sapply(exposure, as.numeric))
  return(list(Dxt = Dxt, Ext = Ext, ages = ages, years = years))
}

# --- Custom LC + Temperature Model Fitting Function ---
fit_LC_with_temp <- function(log_mx, ages, years, temp_vec) {
  n_ages <- length(ages)
  n_years <- length(years)
  
  ax_start <- rowMeans(log_mx)
  kt_start <- colMeans(log_mx)
  bx_start <- rep(0.1, n_ages)
  cx_start <- rep(0.01, n_ages)
  
  par_start <- c(ax_start, bx_start, kt_start, cx_start)
  
  idx_ax <- 1:n_ages
  idx_bx <- (n_ages + 1):(2 * n_ages)
  idx_kt <- (2 * n_ages + 1):(2 * n_ages + n_years)
  idx_cx <- (2 * n_ages + n_years + 1):(3 * n_ages + n_years)
  
  loss_fn <- function(par) {
    ax <- par[idx_ax]
    bx <- par[idx_bx]
    kt <- par[idx_kt]
    cx <- par[idx_cx]
    fitted <- outer(ax, rep(1, n_years)) +
      outer(bx, kt) +
      outer(cx, temp_vec)
    sum((log_mx - fitted)^2, na.rm = TRUE)
  }
  
  fit <- optim(par = par_start, fn = loss_fn, method = "BFGS", control = list(maxit = 1000))
  par_hat <- fit$par
  ax_hat <- par_hat[idx_ax]
  bx_hat <- par_hat[idx_bx]
  kt_hat <- par_hat[idx_kt]
  cx_hat <- par_hat[idx_cx]
  
  fitted_log_mx <- outer(ax_hat, rep(1, n_years)) +
    outer(bx_hat, kt_hat) +
    outer(cx_hat, temp_vec)
  
  return(list(ax = ax_hat, bx = bx_hat, kt = kt_hat, cx = cx_hat,
              ages = ages, years = years, temp = temp_vec,
              fitted_log_mx = fitted_log_mx))
}

# --- Forecast function for extended LC + Temperature model ---
forecast_LC_with_temp <- function(lc_temp_fit, h = 10) {
  ax <- lc_temp_fit$ax
  bx <- lc_temp_fit$bx
  kt <- lc_temp_fit$kt
  cx <- lc_temp_fit$cx
  ages <- lc_temp_fit$ages
  years <- lc_temp_fit$years
  temp <- lc_temp_fit$temp
  
  n_ages <- length(ages)
  n_years <- length(years)
  
  # Forecast kt using ARIMA(0,1,0) with drift (random walk with drift)
  kt_ts <- ts(kt, start = min(years), frequency = 1)
  fit_kt <- forecast::Arima(kt_ts, order = c(0,1,0), include.drift = TRUE)
  kt_forecast <- forecast::forecast(fit_kt, h = h)$mean
  
  # Forecast temperature using auto.arima (automatic selection)
  temp_ts <- ts(temp, start = min(years), frequency = 1)
  fit_temp <- forecast::auto.arima(temp_ts)
  temp_forecast <- forecast::forecast(fit_temp, h = h)$mean
  
  # Combine historical and forecasted kt and temperature
  kt_all <- c(kt, kt_forecast)
  temp_all <- c(temp, temp_forecast)
  
  years_forecast <- seq(max(years) + 1, max(years) + h)
  years_all <- c(years, years_forecast)
  
  # Compute fitted + forecast log mortality rates
  fitted_forecast_log_mx <- outer(ax, rep(1, length(years_all))) +
    outer(bx, kt_all) +
    outer(cx, temp_all)
  
  return(list(
    ages = ages,
    years = years_all,
    ax = ax,
    bx = bx,
    kt = kt_all,
    cx = cx,
    temp = temp_all,
    fitted_log_mx = fitted_forecast_log_mx
  ))
}

# --- Main wrapper function ---
fit_both_models <- function(male_mortality_file, male_exposure_file,
                            female_mortality_file, female_exposure_file,
                            temp_file, h = 10) {
  temp_df <- read.csv(temp_file)
  colnames(temp_df) <- c("Year", "AvgTemp")
  temp_lookup <- setNames(temp_df$AvgTemp, temp_df$Year)
  
  process_sex <- function(m_file, e_file, sex_label) {
    data <- load_data(m_file, e_file)
    log_mx <- log(data$Dxt / data$Ext)
    T_t <- temp_lookup[as.character(data$years)]
    
    # --- Baseline LC using StMoMo ---
    data_obj <- structure(list(
      Dxt = data$Dxt,
      Ext = data$Ext,
      ages = data$ages,
      years = data$years,
      type = "central",
      series = "total",
      label = paste0(sex_label, " Mortality")
    ), class = "StMoMoData")
    
    LC <- lc(link = "log", const = "sum")
    ages_fit <- data$ages[data$ages >= 20 & data$ages <= 89]
    wxt <- genWeightMat(ages = ages_fit, years = data$years, clip = 3)
    LC_fit <- fit(LC, data = data_obj, ages.fit = ages_fit, wxt = wxt)
    LC_forecast <- forecast(LC_fit, h = h)
    
    # --- Extended LC with Temperature ---
    lc_temp_fit <- fit_LC_with_temp(log_mx, data$ages, data$years, T_t)
    lc_temp_forecast <- forecast_LC_with_temp(lc_temp_fit, h = h)
    
    return(list(
      baseline = list(model = LC_fit, forecast = LC_forecast),
      extended = list(fit = lc_temp_fit, forecast = lc_temp_forecast)
    ))
  }
  
  male <- process_sex(male_mortality_file, male_exposure_file, "Male")
  female <- process_sex(female_mortality_file, female_exposure_file, "Female")
  
  return(list(male = male, female = female))
}

# File paths (example)
male_mortality_file <- "C:/Users/tanay/OneDrive/Desktop/Projects/DRP Summer 2025/1991-2023 - Male Mortality Count age 20-89.csv"
male_exposure_file <- "C:/Users/tanay/OneDrive/Desktop/Projects/DRP Summer 2025/1991-2023 - Male Exposure Count 20-89.csv"

female_mortality_file <- "C:/Users/tanay/OneDrive/Desktop/Projects/DRP Summer 2025/1991-2023 - Female Mortality Count age 20-89.csv"
female_exposure_file <- "C:/Users/tanay/OneDrive/Desktop/Projects/DRP Summer 2025/1991-2023 - Female Exposure Count 20-89.csv"

temp_file <- "C:/Users/tanay/OneDrive/Desktop/Projects/DRP Summer 2025/mean_annual_temperature.csv"

# Run the full pipeline with 10-year forecast horizon
results <- fit_both_models(
  male_mortality_file,
  male_exposure_file,
  female_mortality_file,
  female_exposure_file,
  temp_file,
  h = 10
)

plot_params <- function(results, sex = "male", use_forecast = FALSE) {
  baseline_fit <- results[[sex]]$baseline$model
  baseline_forecast <- results[[sex]]$baseline$forecast
  if (use_forecast) {
    extended_fit <- results[[sex]]$extended$forecast
  } else {
    extended_fit <- results[[sex]]$extended$fit
  }
  
  ages <- extended_fit$ages
  
  # 1. ax - Baseline vs Extended (fit only)
  plot(ages, baseline_fit$ax, type = "l", col = "blue", lwd = 2,
       main = paste0("ax (Age Effect) - ", sex),
       xlab = "Age", ylab = "ax")
  lines(ages, extended_fit$ax, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = c("Baseline LC", "Extended LC+Temp"),
         col = c("blue", "red"), lty = c(1,2), bty = "n")
  
  # 2. bx - Baseline vs Extended (fit only)
  plot(ages, baseline_fit$bx, type = "l", col = "blue", lwd = 2,
       main = paste0("bx (Age Sensitivity) - ", sex),
       xlab = "Age", ylab = "bx")
  lines(ages, extended_fit$bx, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = c("Baseline LC", "Extended LC+Temp"),
         col = c("blue", "red"), lty = c(1,2), bty = "n")
  
  # 3. cx - Extended only
  plot(ages, extended_fit$cx, type = "l", col = "darkgreen", lwd = 2,
       main = paste0("cx (Temp Sensitivity) - ", sex),
       xlab = "Age", ylab = "cx")
  
  # 4. kt - Baseline vs Extended over time
  years_orig <- baseline_fit$years
  h <- length(baseline_forecast$kt)
  years_forecast <- seq(max(years_orig) + 1, by = 1, length.out = h)
  years_baseline <- c(years_orig, years_forecast)
  kt_baseline <- c(baseline_fit$kt, baseline_forecast$kt)
  
  years_extended <- extended_fit$years
  kt_extended <- extended_fit$kt
  
  plot(years_baseline, kt_baseline, type = "l", col = "blue", lwd = 2,
       main = paste0("kt (Time Index) - ", sex),
       xlab = "Year", ylab = "kt")
  lines(years_extended, kt_extended, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = c("Baseline LC", "Extended LC+Temp"),
         col = c("blue", "red"), lty = c(1,2), bty = "n")
}


# Plot with forecasts included for kt
plot_params(results, sex = "male", use_forecast = TRUE)
plot_params(results, sex = "female", use_forecast = TRUE)


calculate_fit_metrics <- function(log_mx_true, log_mx_fitted) {
  mse <- mean((log_mx_true - log_mx_fitted)^2, na.rm = TRUE)
  rmse <- sqrt(mse)
  mape <- mean(abs((exp(log_mx_true) - exp(log_mx_fitted)) / exp(log_mx_true)), na.rm = TRUE)
  return(list(RMSE = rmse, MAPE = mape))
}

# Male metrics
male_data <- load_data(male_mortality_file, male_exposure_file)
log_mx_male <- log(male_data$Dxt / male_data$Ext)

baseline_male_fitted_logmx <- fitted(results$male$baseline$model, type = "rates")
extended_male_fitted_logmx <- results$male$extended$fit$fitted_log_mx

male_baseline_metrics <- calculate_fit_metrics(log_mx_male, log(baseline_male_fitted_logmx))
male_extended_metrics <- calculate_fit_metrics(log_mx_male, extended_male_fitted_logmx)

# Female metrics
female_data <- load_data(female_mortality_file, female_exposure_file)
log_mx_female <- log(female_data$Dxt / female_data$Ext)

baseline_female_fitted_logmx <- fitted(results$female$baseline$model, type = "rates")
extended_female_fitted_logmx <- results$female$extended$fit$fitted_log_mx

female_baseline_metrics <- calculate_fit_metrics(log_mx_female, log(baseline_female_fitted_logmx))
female_extended_metrics <- calculate_fit_metrics(log_mx_female, extended_female_fitted_logmx)

# View results
male_baseline_metrics
male_extended_metrics

female_baseline_metrics
female_extended_metrics

# Residuals Diagnostics
plot_residuals <- function(log_mx_true, log_mx_fitted, sex = "male", model = "Extended LC+Temp") {
  residuals <- log_mx_true - log_mx_fitted
  residuals_vec <- as.vector(residuals)
  residuals_vec <- residuals_vec[!is.na(residuals_vec)]
  
  # Histogram
  hist(residuals_vec, breaks = 50, main = paste("Residual Histogram -", model, "-", sex),
       xlab = "Residual", col = "lightblue", border = "white")
  
  # Q-Q Plot
  qqnorm(residuals_vec, main = paste("Q-Q Plot -", model, "-", sex))
  qqline(residuals_vec)
}


# Male - Baseline
plot_residuals(log_mx_male, log(baseline_male_fitted_logmx), sex = "male", model = "Baseline LC")

# Male - Extended
plot_residuals(log_mx_male, extended_male_fitted_logmx, sex = "male", model = "LC+Temp")

# Female - Baseline
plot_residuals(log_mx_female, log(baseline_female_fitted_logmx), sex = "female", model = "Baseline LC")

# Female - Extended
plot_residuals(log_mx_female, extended_female_fitted_logmx, sex = "female", model = "LC+Temp")

# Load temperature data
temp_df <- read.csv("C:/Users/tanay/OneDrive/Desktop/Projects/DRP Summer 2025/mean_annual_temperature.csv")
colnames(temp_df) <- c("Year", "AvgTemp")

# Plot
ggplot(temp_df, aes(x = Year, y = AvgTemp)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(color = "darkred", size = 2) +
  labs(
    title = "Mean Annual Temperature in Toronto (1991–2023)",
    x = "Year",
    y = "Mean Temperature (°C)"
  ) +
  theme_minimal(base_size = 14) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed")




