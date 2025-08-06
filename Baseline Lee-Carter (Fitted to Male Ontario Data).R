# Load necessary libraries
library(StMoMo)
library(ggplot2)

# Function to create Lee-Carter model from CSV data
create_lee_carter_from_csv <- function(mortality_file, exposure_file) {

  # =====================================================
  # STEP 1: READ AND PREPARE DATA FROM CSV
  # =====================================================

  # Read CSV file
  # Expected CSV format:
  # - First column: Age
  # - Remaining columns: Years (e.g., 1950, 1951, ..., 2020)
  # - Each cell mortality rates per 1000 population
  print("Reading CSV file...")
  mortality_data <- read.csv(mortality_file, header = TRUE, row.names = 1)
  exposure_data <- read.csv(exposure_file, header = TRUE, row.names = 1)

  print(paste("Data dimensions:", nrow(mortality_data), "ages x", ncol(mortality_data), "years"))

  # =====================================================
  # STEP 2: CONVERT TO STMOMO FORMAT
  # =====================================================

  # Extract ages and years from the data structure
  ages <- as.numeric(rownames(mortality_data))
  years <- as.numeric(sub("X", "", colnames(mortality_data)))

  mortality_data <- as.data.frame(lapply(mortality_data, function(x) as.numeric(as.character(x))))
  exposure_data <- as.data.frame(lapply(exposure_data, function(x) as.numeric(as.character(x))))

  
  # Convert to matrix format (ages x years)
  Dxt <- as.matrix(mortality_data)  # Deaths matrix
  Ext <- as.matrix(exposure_data) # Exposure matrix

  View(Dxt)  # View the deaths matrix
  View(Ext)  # View the exposure matrix

  print("Data matrices created:")
  print(paste("Ages range:", min(ages), "to", max(ages)))
  print(paste("Years range:", min(years), "to", max(years)))

  # =====================================================  
  # STEP 3: CREATE STMOMODATA OBJECT
  # =====================================================

  # Create StMoMoData object
  mortality_data <- structure(list(
    Dxt = Dxt,
    Ext = Ext, 
    ages = ages,
    years = years,
    type = "central",
    series = "total",
    label = "CSV Mortality Data"
  ), class = "StMoMoData")

  print("StMoMoData object created successfully")

  # =====================================================
  # STEP 4: DEFINE LEE-CARTER MODEL
  # =====================================================

  # Create Lee-Carter model using built-in function
  LC <- lc(link = "log", const = "sum")
  print("Lee-Carter model defined")
  print(LC)

  # =====================================================
  # STEP 5: FIT THE MODEL
  # =====================================================

  print("Fitting Lee-Carter model...")

  # Define age range for fitting (adjust as needed)
  ages_fit <- ages[ages >= 20 & ages <= 90]  # Fit for ages 20-90

  # Create weight matrix to exclude edge cohorts if needed
  wxt <- genWeightMat(ages = ages_fit, years = years, clip = 3)

  # Fit the model
  LC_fit <- fit(LC, 
                data = mortality_data,
                ages.fit = ages_fit,
                wxt = wxt,
                verbose = TRUE)

  print("Model fitting completed!")

  # =====================================================
  # STEP 6: DISPLAY RESULTS
  # =====================================================

  # Display model summary
  print("=== MODEL SUMMARY ===")
  print(LC_fit)

  # Display parameter estimates as tables (data frames)
  cat("\nTable: Age-specific parameters (αx) & Age-specific sensitivity (βx)\n")

  table <- data.frame(Age = LC_fit$ages, Alpha = LC_fit$ax, Beta = LC_fit$bx)
  View(table)  # Show all rows

  cat("\nTable: Time index (κt)\n")

  kappa_table <- data.frame(Year = LC_fit$years, Kappa = LC_fit$kt)
  View((head(kappa_table,1)))  # Show all rows

  # =====================================================
  # STEP 7: PLOT RESULTS
  # =====================================================

  print("Creating plots...")

  # Plot the fitted model
  plot(LC_fit, nCol = 3)

  

  # =====================================================
  # STEP 8: FORECASTING (OPTIONAL)
  # =====================================================

  print("Generating forecasts...")

  # Forecast mortality for next 10 years
  LC_forecast <- forecast(LC_fit, h = 10)
  print("Forecast completed")

  # Plot forecasts
  plot(LC_forecast)


  # =====================================================
  # STEP 9: RETURN RESULTS
  # =====================================================

  return(list(
    model = LC,
    fitted_model = LC_fit,
    forecast = LC_forecast,
    data = mortality_data
  ))
}


# Instructions for use with your own CSV file:
print("INSTRUCTIONS FOR USING WITH YOUR OWN CSV FILE:")
print("1. Prepare your CSV file with:")
print("   - First column: Age values (e.g., 20, 21, 22, ...)")
print("   - Header row: Year values (e.g., 2000, 2001, 2002, ...)")
print("   - Each cell: Death counts or mortality rates for that age-year combination")
print("2. Call: results <- create_lee_carter_from_csv('your_file.csv')")
print("3. The function will automatically fit the Lee-Carter model")

#Input CSV
mortality_file <- file.path("C:","Users","tanay","OneDrive",
"Desktop","Projects","DRP Summer 2025",
"1991-2023 - Male Mortality Count age 20-89.csv")

exposure_file <- file.path("C:","Users","tanay","OneDrive",
"Desktop","Projects","DRP Summer 2025",
"1991-2023 - Male Exposure Count 20-89.csv")

#Fit Model
results <- create_lee_carter_from_csv(mortality_file,exposure_file)
print("\n=== RESULTS SUMMARY ===")
print("Available results:")
print(names(results))
print("Model fitting successful!")
return(results)