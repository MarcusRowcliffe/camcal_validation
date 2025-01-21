
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


##############################################################################
######################### Load in the required data  #########################
predval <-read.csv("https://raw.githubusercontent.com/nilanjanchatterjee/camcal_validation/main/camera_deployment_validation_.csv")
Hoge <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HogeVeluwe/HogeVeluweDigidat.csv")
posdat_mov <-read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/posdat_mov.csv")




# # ##############################################################################################################################################################
# # ############################################################### Hoge Veluwe Validation ####################################################################

# Goal: To estimate how accurately camera traps can measure distance from the camera and to assess how errors 
# in these measurements affect detection distances.
# 
#    Approach:
     
      # Calibration: Uses known distances (from poles in the camera's field of view) to build a model linking image pixels to real-world distances.
      # Simulation: Generates random animal positions and calculates "true" distances, then adds errors based on observed biases.
      # Validation: Compares the performance of different detection models (e.g., half-normal and hazard rate) to understand which handles errors better.
      # Analysis: Tests how sample size and terrain type (flat vs sloping) influence the reliability of distance estimates.



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
    rename(x = xg, y = yg) %>%
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
simulated_vals <- data.frame()

# Loop through each deployment and extract data
for (i in seq_along(predictions)) {
  pred <- predictions[[i]]
  pred <- pred %>%
    dplyr::mutate(
      error = radius - distance, 
      folder = deps[i]           # Add deployment identifier
    )
  # keep all columns
  simulated_vals <- dplyr::bind_rows(simulated_vals, pred)
}

##############################################################################
############### Filtering extreme outliers from error data ###################

# Remove infinite values
filtered_simulated_vals <- simulated_vals[!is.infinite(simulated_vals$radius), ] 

# Calculate IQR for error and define thresholds
iqr_error <- IQR(filtered_simulated_vals$error, na.rm = TRUE)
q1_error <- quantile(filtered_simulated_vals$error, 0.25, na.rm = TRUE)
q3_error <- quantile(filtered_simulated_vals$error, 0.75, na.rm = TRUE)
lower_threshold <- q1_error - 1.5 * iqr_error
upper_threshold <- q3_error + 1.5 * iqr_error

# Filter the dataset to remove outliers
filtered_simulated_vals <- filtered_simulated_vals %>%
  dplyr::filter(error >= lower_threshold & error <= upper_threshold)





##############################################################################
################### Assign parameters for simulation #########################

# Set simulation parameters
reps <- 500
b_avg <- apply(modcoef[, 2:4], 2, mean, na.rm = TRUE)  # Average coefficients
maxr <- 25

# Calculate means and standard deviations for x, y, and error
mnx <- mean(Hoge$xg, na.rm = TRUE)  # Mean x-coordinate
mny <- mean(Hoge$yg, na.rm = TRUE)  # Mean y-coordinate
sdx <- sd(Hoge$xg, na.rm = TRUE)    # SD x-coordinate
sdy <- sd(Hoge$yg, na.rm = TRUE)    # SD y-coordinate
mnerr <- mean(filtered_simulated_vals$error, na.rm = TRUE)  # Mean error
sderr <- sd(filtered_simulated_vals$error, na.rm = TRUE)   # SD error


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


##############################################################################
######################### Review the results #################################

# Comparing results between Hampstead Heath data and Hoge Veluwe data
flat50 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/flat50.csv")
hoge_50

# View summary of results
summary(flat50)
summary(hoge_50)


# Combine datasets for plotting
combined_results <- data.frame(
  Radius = c(flat50$r_hn_true, hoge_50$r_hn_true,
             flat50$r_hn_err, hoge_50$r_hn_err,
             flat50$r_hr_true, hoge_50$r_hr_true,
             flat50$r_hr_err, hoge_50$r_hr_err),
  Method = rep(c("Half Normal", "Half Normal", "Hazard Rate", "Hazard Rate"), each = nrow(hoge_50) + nrow(flat50)),
  Dataset = rep(c("HH", "NLD"), each = nrow(flat50), times = 4)
)

# Create the boxplot
ggplot(combined_results, aes(x = interaction(Method, Dataset), y = Radius, fill = Dataset)) +
  geom_boxplot() +
  scale_fill_manual(values = c("orange", "skyblue")) +
  labs(
    title = "Comparison of Effective Radius Estimates",
    x = "Method and Dataset",
    y = "Effective Radius (m)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(title = "Dataset"))











##############################################################################################################################################################
######################################################## Speed simulation function ##################################################################

# Here, we introduce sequence-level errors (seq_error_updated): a shared error term for all observations within the same sequence. 
# This reflects real-world scenarios where measurements within the same sequence are likely to be influenced by the same systematic biases or 
# environmental factors, creating correlation.
# 
# In the speed simulation below, we have used the Dutch-derived error (error in filtered_simulated_vals), propagating this through the 
# deployment-level mixed-effects model and into the radius and speed calculations.


# ensure required functions for speed_err are loaded
devtools::source_url("https://raw.githubusercontent.com/nilanjanchatterjee/camcal_validation/main/CTtracking_err.r")


# PART 1: Data Preparation

lmrsum_updated <- summary(lmer(as.numeric(error) ~ 1 + (1 | folder), data = filtered_simulated_vals))
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

par(mfrow = c(1, 2))
plot(seqdat_try_updated$pixdiff, seqdat_try_updated$speed, log = "xy", xlab = "Pixel difference", ylab = "Speed")
plot(seqdat_try_updated$pixdiff, seqdat_try_updated$speed_err, log = "xy", xlab = "Pixel difference", ylab = "Speed")

boxplot(speeds_err_500_matrix, names = c("500_rep_updated"), ylab = "Estimated speed (m/s)")
abline(h = truespeed_updated, lwd = 2, col = "red")





#######################################################################################################################################################
################################################# Examining differences between Dutch & HH data  ######################################################

# read in HH data
hh_speeds_err_500_df <- read.csv(
  "https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/speeds_err_hh.csv")

# rename dutch data
dutch_speeds_err_500_df <- speeds_err_500_df

# Convert the data frames to matrices for compatibility
dutch_speeds_err_500_matrix <- speeds_err_500_matrix
hh_speeds_err_500_matrix <- t(as.matrix(hh_speeds_err_500_df))

# Combine the data into a single list for the plot
combined_data <- list(
  "Dutch" = as.vector(dutch_speeds_err_500_matrix),
  "HH" = as.vector(hh_speeds_err_500_matrix)
)


# Create the boxplot
boxplot(
  combined_data,
  names = c("Dutch", "HH"),
  ylab = "Estimated speed (m/s)",
  main = "Comparison of Speed Estimation Errors"
)

# Add reference line for true speed
abline(h = truespeed_updated, lwd = 2, col = "red")


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



  

