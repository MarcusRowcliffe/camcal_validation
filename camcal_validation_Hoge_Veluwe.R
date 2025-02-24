
### Packages
install_and_load <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    library(package_name, character.only = TRUE)
  }
}

packages <- c("tidyverse", "MASS", "MASS", "lattice",
              "pbapply", "lme4", "parallel", "cluster", "tidyverse", "ggplot2", "tidyr")

invisible(lapply(packages, install_and_load))

devtools::source_url("https://raw.githubusercontent.com/MarcusRowcliffe/distanceDF/master/distancedf.r")
devtools::source_url("https://raw.githubusercontent.com/nilanjanchatterjee/camcal_validation/main/CTtracking_err.r")
source("https://raw.githubusercontent.com/MarcusRowcliffe/sbd/main/R/hmean.r")



##############################################################################
######################### Load in the required data  #########################
predval <-read.csv("https://raw.githubusercontent.com/nilanjanchatterjee/camcal_validation/main/camera_deployment_validation_.csv")
Hoge <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HogeVeluwe/HogeVeluweDigidat.csv")
posdat_mov <-read.csv("https://raw.githubusercontent.com/nilanjanchatterjee/camcal_validation/main/Speed_seq_data.csv")

#### Option 2: Read in data from "Realistic" scneario for comparison at a later date. 
flat50 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/flat50.csv")
flat100 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/flat100.csv")
flat200 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/flat200.csv")
slop50 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/slop50.csv")
slop100 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/slop100.csv")
slop200 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/slop200.csv")




# # ##############################################################################################################################################################
# # ############################################################### Hoge Veluwe Validation ####################################################################

# Goal: To estimate how accurately camera traps can measure distance from the camera and to assess how errors 
# in these measurements affect detection distances.
# 
#    Approach:

# Calibration: Uses known distances (from poles in the camera's field of view) to build a model linking image pixels to real-world distances.
# Simulation: Generates random animal positions and calculates "true" distances, then adds errors based on observed biases.
# Validation: Compares the performance of different detection models (half-normal and hazard rate) to understand which handles errors better.



# In this dataset, all poles were digitised with known distances and there was no predefined separation of calibration and test data.
# This means the data needs to be split randomly within the script (as done by oneRepCoefficients),
# simulating the separation of calibration and test sets for model validation.



# Rename 
dat <- Hoge

# Function: extract coefficients
oneRepCoefficients <- function(dat, dep) {
  # Randomly split the data into calibration and test sets
  sq <- 1:nrow(dat)
  i <- sq %in% sample(sq, round(0.5 * nrow(dat)))
  
  # Fit the calibration model
  dmod <- cal.dep(dat[i, ], flex = FALSE)
  names(dmod) <- dep
  
  # Extract calibration coefficients
  coefs <- coef(dmod[[1]]$model)
  
  # Return coefficients as a df
  return(data.frame(deployment = dep, b1 = coefs[1], b2 = coefs[2], b3 = coefs[3]))
}

# Function: Train and extract coefficients for one deployment
oneDepCoefficients <- function(dep, reps) {
  # Filter for the specific deployment
  deployment_data <- filter(dat, folder == dep)
  
  # Perform replicates and combine results
  coefs_list <- replicate(reps, oneRepCoefficients(deployment_data, dep), simplify = FALSE)
  
  # Summarise coefficients by averaging across replicates
  bind_rows(coefs_list) %>%
    group_by(deployment) %>%
    summarise(across(starts_with("b"), mean))
}

# Extract calibration coefficients for all deployments
deps <- unique(dat$folder)
modcoef <- lapply(deps, oneDepCoefficients, reps = 20) %>% bind_rows()

# failsafe to ensure cttracking loaded last
devtools::source_url("https://raw.githubusercontent.com/MarcusRowcliffe/CTtracking/master/CTtracking.r")

# Function to erform one replicate test for deployment to generate predictions
oneRepPredictions <- function(dat, dep) {
  # split the data into calibration and test sets
  sq <- 1:nrow(dat)
  i <- sq %in% sample(sq, round(0.5 * nrow(dat)))
  
  # Fit the calibration model
  dmod <- cal.dep(dat[i, ], flex = FALSE)
  names(dmod) <- dep
  
  # Predict positions for the test set
  dat[!i, ] %>%
    dplyr::rename(x = xg, y = yg) %>%
    predict.pos(dmod, "folder")
}

# Function: Train and test predictions for one deployment across multiple replicates
oneDep <- function(dep, reps) {
  # Filter data for the specific deployment
  deployment_data <- filter(dat, folder == dep)
  
  # Perform replicates and combine results
  bind_rows(replicate(reps, oneRepPredictions(deployment_data, dep), simplify = FALSE))
}

# Generate predictions for all deployments
predictions <- lapply(deps, oneDep, 20)

# Plot predictions for all deployments
par(mfrow = c(5, 4), mar = c(2, 2, 1, 1)) # Set plot layout and margins
lim <- c(0, 20) # Define plot limits
for (i in seq_along(predictions)) {
  with(predictions[[i]], plot(distance, radius))
  legend("topleft", deps[i], bty = "n")
  lines(lim, lim, col = 2) # Add reference line
}




##############################################################################
############### Extracting the distance prediction error #####################

# Create a dataframe with all available columns from predictions
distance_predictions <- data.frame()

# Loop through each deployment and extract data
for (i in seq_along(predictions)) {
  pred <- predictions[[i]]
  pred <- pred %>%
    dplyr::mutate(
      error = radius - distance, 
      folder = deps[i]           # Add deployment identifier
    )
  distance_predictions <- dplyr::bind_rows(distance_predictions, pred)
}

##############################################################################
############### Filtering extreme outliers from error data ###################

# Remove infinite values
filtered_distance_predictions <- distance_predictions[!is.infinite(distance_predictions$radius), ] 

# Calculate IQR for error and define thresholds
iqr_error <- IQR(filtered_distance_predictions$error, na.rm = TRUE)
q1_error <- quantile(filtered_distance_predictions$error, 0.25, na.rm = TRUE)
q3_error <- quantile(filtered_distance_predictions$error, 0.75, na.rm = TRUE)
lower_threshold <- q1_error - 1.5 * iqr_error
upper_threshold <- q3_error + 1.5 * iqr_error

# Filter the dataset to remove outliers
filtered_distance_predictions <- filtered_distance_predictions %>%
  dplyr::filter(error >= lower_threshold & error <= upper_threshold)




# Extract distance estimation statistics for Hoge Veluwe
hv_distance_stats <- list(
  mean_error = mean(filtered_distance_predictions$error, na.rm = TRUE),
  sd_error = sd(filtered_distance_predictions$error, na.rm = TRUE),
  error_range = range(filtered_distance_predictions$error, na.rm = TRUE),
  r_squared = summary(lm(radius ~ distance, data = filtered_distance_predictions))$r.squared,
  p_value = cor.test(filtered_distance_predictions$distance, filtered_distance_predictions$radius)$p.value
)

hv_distance_stats


# plot residuals vs distance
ggplot(filtered_distance_predictions, aes(x = distance, y = error)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  labs(title = "Residuals vs Distance", x = "Distance (m)", y = "Residual Error (m)") +
  theme_minimal()

# linear model to check trend
lm_model <- lm(error ~ distance, data = filtered_distance_predictions)
summary(lm_model)




##############################################################################
############ Compare distance estimation across methods (HH vs HV)  ##########

predval <-read.csv("https://raw.githubusercontent.com/nilanjanchatterjee/camcal_validation/main/camera_deployment_validation_.csv")

head(predval)
head(filtered_distance_predictions)

# distance in predval in cm, so divide by 100
predval$distance <- predval$distance/100
predval$diff <- predval$diff/100



# Add labels
predval$Dataset <- "Hampstead Heath"
filtered_distance_predictions$Dataset <- "Hoge Veluwe"

# Select columns and combine 
hh_data <- predval %>% dplyr::select(distance, radius, Dataset)
hv_data <- filtered_distance_predictions %>% dplyr::select(distance, radius, Dataset)
combined_data <- dplyr::bind_rows(hh_data, hv_data)


##########################################################################################
#### Exploring the distance data 

head(combined_data)
combined_data$error <- combined_data$radius - combined_data$distance
summary(combined_data)

ggplot(combined_data, aes(x = error, fill = Dataset)) +
  geom_density(alpha = 0.5) +
  labs(title = "Error Distribution in Distance Estimation",
       x = "Error (Radius - Distance)",
       y = "Density") +
  theme_minimal()



combined_data %>%
  dplyr::group_by(Dataset) %>%
  dplyr::summarise(
    Mean_Error = mean(error, na.rm = TRUE),
    SD_Error = sd(error, na.rm = TRUE),
    Min_Error = min(error, na.rm = TRUE),
    Max_Error = max(error, na.rm = TRUE),
    Median_Error = median(error, na.rm = TRUE),
    Q1_Error = quantile(error, 0.25, na.rm = TRUE),
    Q3_Error = quantile(error, 0.75, na.rm = TRUE)
  )


##########################################################################################





# Add R-squared values
r2_hh <- summary(lm(radius ~ distance, data = hh_data))$r.squared
r2_hv <- summary(lm(radius ~ distance, data = hv_data))$r.squared

print(r2_hh)
print(r2_hv)

### Figure 3: scatter plots comparing object distances predicted by different models 

ggplot(combined_data, aes(x = distance, y = radius, colour = Dataset)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Dataset, labeller = labeller(Dataset = c("Hoge Veluwe" = "Hoge Veluwe", 
                                                       "Hampstead Heath" = "Hampstead Heath"))) +
  labs(
    x = "Actual Distance (m)",
    y = "Estimated Radius (m)"
  ) +
  scale_colour_manual(values = c("Hoge Veluwe" = "skyblue", "Hampstead Heath" = "orange")) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 28, face = "bold"),
    axis.text.x = element_text(size = 18, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 30),
    legend.position = "none",  # Removes legend
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_equal() +  # Equal axis scales
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "black") +  # x=y line
  geom_text(data = data.frame(Dataset = c("Hampstead Heath", "Hoge Veluwe"),
                              x = c(5, 5),
                              y = c(15, 15),
                              label = c(paste0("R² = ", round(r2_hh, 2)), 
                                        paste0("R² = ", round(r2_hv, 2)))),
            aes(x = x, y = y, label = label),
            size = 6, color = "#B22222", fontface = "bold", inherit.aes = FALSE)


# 
# # Add labels
# predval$Dataset <- "Realistic Scenario"
# filtered_distance_predictions$Dataset <- "Controlled Scenario"
# 
# # Selectcolumns and combine 
# hh_data <- predval %>% dplyr::select(distance, radius, Dataset)
# hv_data <- filtered_distance_predictions %>% dplyr::select(distance, radius, Dataset)
# combined_data <- dplyr::bind_rows(hh_data, hv_data)
# 
# 
# # Add R-squared values
# r2_hh <- summary(lm(radius ~ distance, data = hh_data))$r.squared
# r2_hv <- summary(lm(radius ~ distance, data = hv_data))$r.squared
# 
# print(r2_hh)
# print(r2_hv)
# ### Figure 3: scatter plots comparing object distances predicted by diff models 
# 
# ggplot(combined_data, aes(x = distance, y = radius, colour = Dataset)) +
#   geom_point(alpha = 0.8) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~Dataset, labeller = labeller(Dataset = c("Controlled Scenario" = "Controlled Scenario", 
#                                                        "Realistic Scenario" = "Realistic Scenario"))) +
#   labs(
#     x = "Actual Distance (m)",
#     y = "Estimated Radius (m)"
#   ) +
#   scale_colour_manual(values = c("Controlled Scenario" = "skyblue", "Realistic Scenario" = "orange")) +
#   theme_minimal() +
#   theme(
#     strip.text = element_text(size = 28, face = "bold"),
#     axis.text.x = element_text(size = 18, angle = 0, hjust = 0.5),
#     axis.text.y = element_text(size = 18),
#     axis.title = element_text(size = 30),
#     legend.position = "none",  # Removes legend
#     plot.margin = margin(10, 10, 10, 10)
#   ) +
#   coord_equal() +  # Equal axis scales
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "black") +  # x=y line
#   geom_text(data = data.frame(Dataset = c("Realistic Scenario", "Controlled Scenario"),
#                               x = c(5, 5),
#                               y = c(15, 15),
#                               label = c(paste0("R² = ", round(r2_hh, 2)), 
#                                         paste0("R² = ", round(r2_hv, 2)))),
#             aes(x = x, y = y, label = label),
#             size = 6, color = "#B22222", fontface = "bold", inherit.aes = FALSE)


# # ##############################################################################################################################################################
# # ############################################################### Simulating distance error  ####################################################################



##############################################################################
################### Assign parameters for simulation  #########################

# Set parameters
reps <- 500
b_avg <- apply(modcoef[, 2:4], 2, mean, na.rm = TRUE)  # Average coefficients
maxr <- 25

# Calculate means and standard deviations for x, y, and error
mnx <- mean(Hoge$xg, na.rm = TRUE)  # Mean x-coordinate
mny <- mean(Hoge$yg, na.rm = TRUE)  # Mean y-coordinate
sdx <- sd(Hoge$xg, na.rm = TRUE)    # SD x-coordinate
sdy <- sd(Hoge$yg, na.rm = TRUE)    # SD y-coordinate
mnerr <- mean(filtered_distance_predictions$error, na.rm = TRUE)  # Mean error
sderr <- sd(filtered_distance_predictions$error, na.rm = TRUE)   # SD error





##############################################################################
######################### Define the function ################################

# Function to simulate effective radius estimation
sim_rep <- function(points, b, maxr, mnx, mny, mnerr, sdx, sdy, sderr) {
  
  safe_try <- function(expr) {
    repeat {
      result <- tryCatch(expr, error = function(e) NA)
      if (!is.na(result)) return(result)
    }
  }
  
  # Step 1: Generate simulated x and y coordinates
  xsim <- (rnorm(points, mnx, sdx) / 2000) - 0.5
  ysim <- rnorm(points, mny, sdy) / 1500
  
  # Step 2: Calculate true radius using the calibration model
  radius <- b[1] / (ysim - (b[2] + b[3] * xsim))
  radius <- radius[radius > 0 & radius < maxr]  # Filter valid radius values
  
  if (length(radius) == 0) return(rep(NA, 10))  # Return NA vector if radius is empty
  
  # Step 3: Add observed error to the true radius
  err <- rnorm(length(radius), mnerr, sderr)
  pred <- radius + err
  dat <- data.frame(true = radius, err = pred)
  
  # Step 4: Fit detection functions and retrieve estimates 
  r_hn_true <- safe_try(fitdf(true ~ 1, key = "hn", transect = "point", order = 0, data = dat)$edd$estimate[1])
  r_hr_true <- safe_try(fitdf(true ~ 1, key = "hr", transect = "point", order = 0, data = dat)$edd$estimate[1])
  r_hn_err <- safe_try(fitdf(err ~ 1, key = "hn", transect = "point", order = 0, data = dat)$edd$estimate[1])
  r_hr_err <- safe_try(fitdf(err ~ 1, key = "hr", transect = "point", order = 0, data = dat)$edd$estimate[1])
  
  # Step 5: Calculate differences and AICs 
  diff_hn <- if (!is.na(r_hn_true) && !is.na(r_hn_err)) r_hn_err - r_hn_true else NA
  diff_hr <- if (!is.na(r_hr_true) && !is.na(r_hr_err)) r_hr_err - r_hr_true else NA
  aic_hn_true <- safe_try(fitdf(true ~ 1, key = "hn", transect = "point", order = 0, data = dat)$ddf$criterion[1])
  aic_hr_true <- safe_try(fitdf(true ~ 1, key = "hr", transect = "point", order = 0, data = dat)$ddf$criterion[1])
  aic_hn_err <- safe_try(fitdf(err ~ 1, key = "hn", transect = "point", order = 0, data = dat)$ddf$criterion[1])
  aic_hr_err <- safe_try(fitdf(err ~ 1, key = "hr", transect = "point", order = 0, data = dat)$ddf$criterion[1])
  
  # Return a fixed-length vector of results
  return(c(r_hn_true, r_hr_true, r_hn_err, r_hr_err, diff_hn, diff_hr, aic_hn_true, aic_hr_true, aic_hn_err, aic_hr_err))
}

# Define the column names for the output
column_names <- c("r_hn_true", "r_hr_true", "r_hn_err", "r_hr_err",
                  "diff_hn", "diff_hr", "aic_hn_true", "aic_hr_true",
                  "aic_hn_err", "aic_hr_err")


#####################################################################
#################### Run the function - 50 points ###################

#### Option 1: Run the simulation from scratch
# hoge_50 <- data.frame(t(pbreplicate(
#   reps,
#   suppressMessages(sim_rep(points = 50, b = b_avg, maxr = 25, mnx, mny, mnerr, sdx, sdy, sderr))
# )))
# 
# Assign column names
# colnames(hoge_50) <- column_names

#### Option 2: Read in pre-simulated data to save time
hoge_50 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HogeVeluwe/Simulated_Data/hoge_50.csv")


######################################################################
################ Run the function - 100 points #######################

# ### Option 1: Run the simulation from scratch
# hoge_100 <- data.frame(t(pbreplicate(
#   reps,
#   suppressMessages(sim_rep(points = 100, b = b_avg, maxr = 25, mnx, mny, mnerr, sdx, sdy, sderr))
# )))
# 
# # Assign column names
# colnames(hoge_100) <- column_names


#### Option 2: Read in pre-simulated data to save time
hoge_100 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HogeVeluwe/Simulated_Data/hoge_100.csv")

#####################################################################
#################### Run the function - 200 points ##################

# ### Option 1: Run the simulation from scratch
# hoge_200 <- data.frame(t(pbreplicate(
#   reps,
#   suppressMessages(sim_rep(points = 200, b = b_avg, maxr = 25, mnx, mny, mnerr, sdx, sdy, sderr))
# )))
# 
# # Assign column names
# colnames(hoge_200) <- column_names

#### Option 2: Read in pre-simulated data to save time
hoge_200 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HogeVeluwe/Simulated_Data/hoge_200.csv")


##############################################################################
######################### HN vs HR ##########################################

# Evaluating different detection functions - hazard rate and half normal. 
# Calculate the proportion of simulations where HR has a lower AIC than HN:
mean(hoge_50$aic_hr_true < hoge_50$aic_hn_true)
mean(hoge_50$aic_hr_err < hoge_50$aic_hn_err)




### Figure 4

# Reduce space between plots by adjusting margins
par(mfrow = c(1, 3), mar = c(5, 6, 4, 1), oma = c(1, 1, 1, 1), mgp = c(2.5, 1, 0), cex = 1.2)

y_limits <- c(-1.5, 1.5)

# Plot 1: Hoge Veluwe
boxplot(hoge_50$diff_hn, hoge_100$diff_hn, hoge_200$diff_hn,
        hoge_50$diff_hr, hoge_100$diff_hr, hoge_200$diff_hr,
        names = rep(c(50, 100, 200), 2),
        col = rep(c("#E9967A", "#8FBC8F"), each = 3),
        xlab = "Points", ylab = "Error (m)", main = "Hoge Veluwe",
        cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.2, ylim = y_limits)

# Add reference line at zero
abline(h = 0, col = "magenta", lwd = 2)

# Plot 2: Hampstead Heath (Flat ground)
boxplot(flat50$diff_hn, flat100$diff_hn, flat200$diff_hn,
        flat50$diff_hr, flat100$diff_hr, flat200$diff_hr,
        names = rep(c(50, 100, 200), 2),
        col = rep(c("#E9967A", "#8FBC8F"), each = 3),
        xlab = "Points", ylab = "Error (m)", main = "Hampstead Heath: Flat ground",
        cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.2, ylim = y_limits)

legend(x = -0.75, y = -0.5, c("Half normal", "Hazard rate"),
       fill = c("#E9967A", "#8FBC8F"),
       cex = 1.5,  # Slightly larger legend text
       box.lty = 0,
       x.intersp = 0.2,
       text.font = 2)

abline(h = 0, col = "magenta", lwd = 2)

# Plot 3: Hampstead Heath (Sloping ground)
boxplot(slop50$diff_hn, slop100$diff_hn, slop200$diff_hn,
        slop50$diff_hr, slop100$diff_hr, slop200$diff_hr,
        names = rep(c(50, 100, 200), 2),
        col = rep(c("#E9967A", "#8FBC8F"), each = 3),
        xlab = "Points", ylab = "Error (m)", main = "Hampstead Heath: Sloping ground",
        cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.2, ylim = y_limits)

# Add reference line at zero
abline(h = 0, col = "magenta", lwd = 2)

# Reset layout
par(mfrow = c(1, 1))


##############################################################################
######################### Compare the results across methods #################

# Comparing results between Hampstead Heath data and Hoge Veluwe data
# Here, we'll compare the simulation with 200 points between Hampstead Heath (HH) and Hoge Veluwe (HV). As HH data was split based on flat or 
# sloping, I'll combine the two to compare with HV. 

flat200 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/flat200.csv")
slop200 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/slop200.csv")


# Combine flat and sloping datasets for HH
hh_200 <- rbind(flat200, slop200)

# Dataframe sizes will vary as HH combined flat and sloping simulations. Therefore, randomly sample from HH to match size of HV.
set.seed(123)  
hh_200_sampled <- hh_200[sample(nrow(hh_200), nrow(hoge_200)), ]

# combine dfs for comparison
combined_results <- data.frame(
  Radius = c(hh_200_sampled$r_hn_true, hoge_200$r_hn_true,
             hh_200_sampled$r_hn_err, hoge_200$r_hn_err,
             hh_200_sampled$r_hr_true, hoge_200$r_hr_true,
             hh_200_sampled$r_hr_err, hoge_200$r_hr_err),
  Method = rep(c("Half Normal", "Half Normal", "Hazard Rate", "Hazard Rate"), 
               each = nrow(hoge_200) + nrow(hh_200_sampled)),
  Dataset = rep(c("HH", "NLD"), each = nrow(hh_200_sampled), times = 4)
)

# filter df 
focused_results <- data.frame(
  Radius = c(hh_200_sampled$r_hn_true, hoge_200$r_hn_true,
             hh_200_sampled$r_hn_err, hoge_200$r_hn_err),
  Type = rep(c("True", "Error"), each = nrow(hh_200_sampled) + nrow(hoge_200)),
  Dataset = rep(c("HH", "NLD"), each = nrow(hh_200_sampled), times = 2)
)

# Calculate average differences
avg_diff_HH <- mean(c(hh_200_sampled$r_hn_err - hh_200_sampled$r_hn_true,
                      hh_200_sampled$r_hr_err - hh_200_sampled$r_hr_true), na.rm = TRUE)
avg_diff_NLD <- mean(c(hoge_200$r_hn_err - hoge_200$r_hn_true,
                       hoge_200$r_hr_err - hoge_200$r_hr_true), na.rm = TRUE)

# add labels
text_labels <- data.frame(
  Type = "True",
  Dataset = c("Controlled Scenario", "Realistic Scenario"),
  Label = c(paste("Avg diff:", round(avg_diff_NLD, 2), "m"),
            paste("Avg diff:", round(avg_diff_HH, 2), "m"))
)

# Adjust x positions for labels
text_labels$X_Pos <- c(1.5, 3.5)

# Rename datasets
focused_results$Dataset <- recode(focused_results$Dataset, 
                                  "HH" = "Realistic Scenario", 
                                  "NLD" = "Controlled Scenario")

# Reorder for plot
focused_results$Type_Dataset <- factor(paste(focused_results$Type, focused_results$Dataset),
                                       levels = c("True Controlled Scenario", "Error Controlled Scenario",
                                                  "True Realistic Scenario", "Error Realistic Scenario"))


#### Figure 5


# Update the Dataset column in focused_results
focused_results$Dataset <- ifelse(focused_results$Dataset == "Controlled Scenario", "Hoge Veluwe", "Hampstead Heath")



ggplot(focused_results, aes(x = Type_Dataset, y = Radius, fill = Dataset)) +
  geom_boxplot(outlier.size = 1.5, alpha = 0.8, position = position_dodge(width = 0.8)) +
  scale_x_discrete(labels = c("True", "Error", "True", "Error")) +
  scale_fill_manual(values = c("Hoge Veluwe" = "skyblue", "Hampstead Heath" = "orange")) +
  labs(
    title = "Impact of Error on Estimated Effective Radius",
    x = "Estimate Type",
    y = "Effective Radius (m)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 30),
    plot.title = element_text(size = 30, hjust = 0.5),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20),
    legend.position = c(0.9, 0.78)
  ) +
  guides(fill = guide_legend(title = "Dataset")) +
  geom_text(data = text_labels, 
            aes(x = X_Pos, y = 14.5, label = Label), 
            size = 7, color = "#B22222", fontface = "bold", parse = FALSE)


# get summary of the two methods
summary(hh_200_sampled)
summary(hoge_200)


### Results:

# Results:
#   Dutch Calibration Shows Higher Effective Radii Estimates:
#   Mean True Radii (HN): Dutch 10.80 m, HH 6.06 m; Mean True Radii (HR): Dutch 10.08 m, HH 5.53 m.
# 
# Dutch Calibration Has Lower Bias and Deviation:
#   Average Bias: Dutch +0.24 m, HH +0.42 m.
# 
# 
# Dutch Calibration Has Better Model Fit (Lower AIC):
#   Mean AIC (HN True): Dutch 575.0, HH 846.9; Mean AIC (HR True): Dutch 563.6, HH 795.5.
# 
# Hazard Rate Outperforms Half Normal Across Both Methods:
#   Average Deviation: Dutch +0.24 m, HH +0.42 m; HR has lower AIC and smaller bias across both datasets.



#########################################  Read in other simulated dataframes to make comparisons #########################
###########################################################################################################################

# flat50 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/flat50.csv")
# flat100 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/flat100.csv")
# flat200 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/flat200.csv")
# slop50 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/slop50.csv")
# slop100 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/slop100.csv")
# slop200 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/slop200.csv")
# 
# # summaries for HH (realistic scenario)
# summary(flat50)
# summary(flat100)
# summary(flat200)
# summary(slop50)
# summary(slop100)
# summary(slop200)
# 
# 
# # summaries for HV (controlled scenario)
# summary(hoge_50)
# summary(hoge_100)
# summary(hoge_200)



##############################################################################################################################################################
######################################################## Speed simulation function ##################################################################

# Here, we introduce sequence-level errors (seq_error_updated): a shared error term for all observations within the same sequence. 
# This reflects real-world scenarios where measurements within the same sequence are likely to be influenced by the same systematic biases or 
# environmental factors, creating correlation.
# 
# In the speed simulation below, we have used the Dutch-derived error (error in filtered_distance_predictions), propagating this through the 
# deployment-level mixed-effects model and into the radius and speed calculations.


# ensure required functions for speed_err are loaded
devtools::source_url("https://raw.githubusercontent.com/nilanjanchatterjee/camcal_validation/main/CTtracking_err.r")


# PART 1: Data Preparation

lmrsum_updated <- summary(lmer(as.numeric(error) ~ 1 + (1 | folder), data = filtered_distance_predictions))
depmn_updated <- lmrsum_updated$coefficients[1]
depsd_updated <- sqrt(lmrsum_updated$varcor$folder[1])

deps_updated <- unique(predval$deployment)
probdep_updated <- as.numeric(xtabs(~deployment, predval) / nrow(predval))

posdat_mov_updated <- posdat_mov
posdat_mov_updated$uid <- paste(posdat_mov_updated$siteid, posdat_mov_updated$sequence_id, sep = "-")
names(posdat_mov_updated)[7] <- "seq_id"
names(posdat_mov_updated)[20] <- "sequence_id"

s1_updated <- as.numeric(xtabs(~sequence_id, posdat_mov_updated))
posdat_mov_updated$new_seq <- rep(sample(deps_updated, 675, replace = TRUE, prob = probdep_updated), times = s1_updated)
posdat_mov_updated$new_seq1 <- paste(posdat_mov_updated$sequence_id, posdat_mov_updated$new_seq, sep = "-")
names(posdat_mov_updated)[20] <- "new_seq_id"
names(posdat_mov_updated)[22] <- "sequence_id"

seqdat_updated <- seq.summary(posdat_mov_updated)
seqdat1_updated <- subset(seqdat_updated, seqdat_updated$timediff < 2000)

avgpixdif_updated <- 400
pixdiff_updated <- seq.data(posdat_mov_updated)$pixdiff



### CHecking correlation in error 
summary(filtered_distance_predictions$error)
print(depsd_updated)

cor.test(filtered_distance_predictions$error[-1], filtered_distance_predictions$error[-length(filtered_distance_predictions$error)])
cor.test(predval$diff[-1], predval$diff[-length(predval$diff)])



# # PART 2: Update the error propagation and speed generation with correlated errors

# ### Option 1: Run the simulation from scratch
# speeds_err_500_updated <- pbreplicate(500, {
#   dep_error_updated <- rnorm(length(deps_updated), depmn_updated, depsd_updated) * 0.01
#   seq_error_updated <- rnorm(length(unique(posdat_mov_updated$sequence_id)), mean = 0, sd = 0.1)
# 
#   posdat_try_updated <<- posdat_mov_updated
#   posdat_try_updated$radius_err <- posdat_try_updated$radius +
#     dep_error_updated[match(posdat_try_updated$new_seq, deps_updated)] +
#     seq_error_updated[match(posdat_try_updated$sequence_id, unique(posdat_mov_updated$sequence_id))]
# 
#   # Save posdat_try_updated globally
#   assign("posdat_try_updated", posdat_try_updated, envir = .GlobalEnv)
# 
# 
#   seqdat_try_updated <<- seq.summary(posdat_try_updated)
# 
#   spd_updated <- seqdat_try_updated$speed_err[is.finite(seqdat_try_updated$speed_err) &
#                                                 seqdat_try_updated$pixdiff > 10 &
#                                                 seqdat_try_updated$timediff < 2000]
# 
#   sampled_spd_updated <- sample(1 / spd_updated, size = 500, replace = FALSE)
#   1 / mean(sampled_spd_updated[is.finite(sampled_spd_updated)])
# })
# 
# 
# # Convert the matrix to a df for saving
# speeds_err_500_df <- as.data.frame(t(speeds_err_500_updated))

# Option 2: Read in pre-simulated data to save time
speeds_err_500_df <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HogeVeluwe/Simulated_Data/speeds_err_dutch.csv")
seqdat_try_updated <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HogeVeluwe/Simulated_Data/seqdat_try_dutch.csv")

# Change format 
speeds_err_500_matrix <- t(as.matrix(speeds_err_500_df))


# PART 3: Visualisation and True Speed Calculation
trspd_updated <- 1 / seqdat1_updated$speed[is.finite(seqdat1_updated$speed) & seqdat1_updated$pixdiff > 10]
truespeed_updated <- 1 / mean(trspd_updated[is.finite(trspd_updated)])
hmean_trspd <- hmean(seqdat1_updated$speed[is.finite(seqdat1_updated$speed) & 
                                             seqdat1_updated$speed > 0 & 
                                             seqdat1_updated$pixdiff > 10])




par(mfrow = c(1, 2))
plot(seqdat_try_updated$pixdiff, seqdat_try_updated$speed, log = "xy", xlab = "Pixel difference", ylab = "Speed")
plot(seqdat_try_updated$pixdiff, seqdat_try_updated$speed_err, log = "xy", xlab = "Pixel difference", ylab = "Speed")

boxplot(speeds_err_500_matrix, names = c("500_rep_updated"), ylab = "Estimated speed (m/s)")
abline(h = truespeed_updated, lwd = 2, col = "red")

par(mfrow = c(1, 1))




#######################################################################################################################################################
################################################# Examining differences between Dutch & HH data  ######################################################

# read in HH data
hh_speeds_err_500_df <- read.csv(
  "https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/speeds_err_hh.csv")


# rename dutch data
dutch_speeds_err_500_df <- speeds_err_500_df


### Figure 6 

# Prep data for plot
dutch_speeds_err_500_matrix <- speeds_err_500_matrix
hh_speeds_err_500_matrix <- t(as.matrix(hh_speeds_err_500_df))

combined_data <- list(
  "Dutch" = as.vector(dutch_speeds_err_500_matrix),
  "HH" = as.vector(hh_speeds_err_500_matrix)
)





df <- data.frame(
  Scenario = rep(c("Controlled Scenario", "Realistic Scenario"), each = length(combined_data[[1]])),
  Speed = c(combined_data[[1]], combined_data[[2]])
)

# Rename scenarios
df <- df %>%
  mutate(Scenario = recode(Scenario, 
                           "Controlled Scenario" = "Hoge Veluwe", 
                           "Realistic Scenario" = "Hampstead Heath"))

##########################################################################################
#### Exploring the speed data
summary(df)

library(dplyr)

# Compute summary statistics for speed
speed_summary <- df %>%
  group_by(Scenario) %>%
  summarise(
    Mean_Speed = mean(Speed, na.rm = TRUE),
    SD_Speed = sd(Speed, na.rm = TRUE),
    Min_Speed = min(Speed, na.rm = TRUE),
    Q1_Speed = quantile(Speed, 0.25, na.rm = TRUE),
    Median_Speed = median(Speed, na.rm = TRUE),
    Q3_Speed = quantile(Speed, 0.75, na.rm = TRUE),
    Max_Speed = max(Speed, na.rm = TRUE)
  )

# Print results
print(speed_summary)




##########################################################################################



# Comparison plot speeds
ggplot(df, aes(x = Scenario, y = Speed, fill = Scenario)) +
  geom_boxplot(outlier.size = 1.5, alpha = 0.8, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Hoge Veluwe" = "skyblue", "Hampstead Heath" = "orange")) +
  labs(
    title = "Comparison of Speed Estimation Errors",
    x = "",
    y = "Estimated Speed (m/s)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 22, angle = 0, hjust = 0.5),  # Increased size
    axis.text.y = element_text(size = 22),  # Increased y-axis label size
    axis.title = element_text(size = 30),
    plot.title = element_text(size = 30, hjust = 0.5),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20),
    legend.position = c(0.85, 0.7)
  ) +
  guides(fill = guide_legend(title = "Scenario")) +
  geom_hline(yintercept = truespeed_updated, linetype = "dashed", color = "red", size = 1.5)




# Convert the data frames to matrices and vectors
dutch_speeds <- as.vector(t(as.matrix(dutch_speeds_err_500_df)))
hh_speeds <- as.vector(t(as.matrix(hh_speeds_err_500_df)))

# Calculate summary stats
dutch_summary <- summary(dutch_speeds)
hh_summary <- summary(hh_speeds)

# and standard deviation
dutch_sd <- sd(dutch_speeds, na.rm = TRUE)
hh_sd <- sd(hh_speeds, na.rm = TRUE)

# Combine results for comparison
summary_stats <- data.frame(
  Metric = c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max", "Std. Dev"),
  Dutch = c(dutch_summary, dutch_sd),
  Hampstead_Heath = c(hh_summary, hh_sd)
)

print(summary_stats)

# Results very similar between Dutch and HH datasets. The mean simulated speeds are close to the true speed for both approaches & 
# the differences seem very marginal. 






  

