
### Packages
install_and_load <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    library(package_name, character.only = TRUE)
  }
}

packages <- c("tidyverse", "MASS", "MASS", "lattice",
              "pbapply", "lme4", "parallel")

invisible(lapply(packages, install_and_load))

devtools::source_url("https://raw.githubusercontent.com/MarcusRowcliffe/distanceDF/master/distancedf.r")
devtools::source_url("https://raw.githubusercontent.com/nilanjanchatterjee/camcal_validation/main/CTtracking_err.r")


##############################################################################
################## Data files for calibration validation #####################
predval <-read.csv("https://raw.githubusercontent.com/nilanjanchatterjee/camcal_validation/main/camera_deployment_validation_.csv")
sitecal <-read.csv("https://raw.githubusercontent.com/nilanjanchatterjee/camcal_validation/main/Site_digitisation_data.csv")
modcoef <- read.csv("https://raw.githubusercontent.com/nilanjanchatterjee/camcal_validation/main/pole_11_mod_param.csv")
posdat_mov <-read.csv("https://raw.githubusercontent.com/nilanjanchatterjee/camcal_validation/main/Speed_seq_data.csv")



#########################################################################################
################### Simulation for Effective radius estimation
#Function:
# 1. samples from x,y co-ordinate distributions
# 2. calculates "true" radius from these using the calibration model
# 3. adds error to true radius, sampled from observed error distribution
# 4. fits detection functions to true and error radii using half normal and 
#    hazard rate keys
# 5. returns a vector of radius estimates and AICs for each of the 4 models, 
#    plus differences between true and error estimates for half normal and 
#    hazard rate models

#INPUTS
# points: number of points to sample
# b: calibration coeficients
# maxr: maximum radius to retain in samples
# mnx, mny, mnerr: means of x,y co-ordinate and error distributions
# sdx, sdy, sderr: sds of x,y co-ordinate and error distributions


sim_rep <- function(points, b, maxr, mnx, mny, mnerr, sdx, sdy, sderr) {
  xsim <- (rnorm(points, mnx, sdx) / 2000) - 0.5
  ysim <- rnorm(points, mny, sdy) / 1500
  radius <- b[1] / (ysim - (b[2] + b[3] * xsim))
  radius <- radius[radius > 0 & radius < maxr]
  
  if (length(radius) == 0) return(rep(NA, 10))
  
  err <- rnorm(length(radius), mnerr, sderr)
  pred <- radius + err
  dat <- data.frame(true = radius, err = pred)
  
  # Fit models and retrieve results
  r_hn_true <- tryCatch(fitdf(true ~ 1, key = "hn", transect = "point", order = 0, data = dat)$edd$estimate[1], error = function(e) NA)
  r_hr_true <- tryCatch(fitdf(true ~ 1, key = "hr", transect = "point", order = 0, data = dat)$edd$estimate[1], error = function(e) NA)
  r_hn_err <- tryCatch(fitdf(err ~ 1, key = "hn", transect = "point", order = 0, data = dat)$edd$estimate[1], error = function(e) NA)
  r_hr_err <- tryCatch(fitdf(err ~ 1, key = "hr", transect = "point", order = 0, data = dat)$edd$estimate[1], error = function(e) NA)
  
  # Calculate differences and AICs
  diff_hn <- if (!is.na(r_hn_true) && !is.na(r_hn_err)) r_hn_err - r_hn_true else NA
  diff_hr <- if (!is.na(r_hr_true) && !is.na(r_hr_err)) r_hr_err - r_hr_true else NA
  aic_hn_true <- tryCatch(fitdf(true ~ 1, key = "hn", transect = "point", order = 0, data = dat)$ddf$criterion[1], error = function(e) NA)
  aic_hr_true <- tryCatch(fitdf(true ~ 1, key = "hr", transect = "point", order = 0, data = dat)$ddf$criterion[1], error = function(e) NA)
  aic_hn_err <- tryCatch(fitdf(err ~ 1, key = "hn", transect = "point", order = 0, data = dat)$ddf$criterion[1], error = function(e) NA)
  aic_hr_err <- tryCatch(fitdf(err ~ 1, key = "hr", transect = "point", order = 0, data = dat)$ddf$criterion[1], error = function(e) NA)
  
  # Return a fixed-length result vector
  return(c(r_hn_true, r_hr_true, r_hn_err, r_hr_err, diff_hn, diff_hr, aic_hn_true, aic_hr_true, aic_hn_err, aic_hr_err))
}


# Extract necessary vals from site digitisation data: 
mnx <- mean(sitecal$x) 
mny <- mean(sitecal$y) 
sdx <- sd(sitecal$x) 
sdy <- sd(sitecal$y)

# Calculate error, mean error and sd of error. 
err <- predval$radius - predval$distance/100
mnerr <- mean(err, na.rm=T)
sderr <- sd(err, na.rm=T)


# Model parameters from the non-linear least square and surface details of the 10 camera locations
bflat <- subset(modcoef, location_type=="flat")[, 2:4] %>%
  apply(2, mean)



# Define the column names 
column_names <- c("r_hn_true", "r_hr_true", "r_hn_err", "r_hr_err", 
                  "diff_hn", "diff_hr", "aic_hn_true", "aic_hr_true", 
                  "aic_hn_err", "aic_hr_err")


# Set the number of reps
reps <- 500



#######################################################################################################
################################# Flat Surfaces #######################################################

### Option 1: Run the simulation from scratch

# Flat50
# flat50 <- data.frame(t(pbreplicate(reps, suppressMessages(
#   sim_rep(points=50, b=bflat, maxr=25, mnx, mny, mnerr, sdx, sdy, sderr)
# ))))
# colnames(flat50) <- column_names  

# Flat100
# flat100 <- data.frame(t(pbreplicate(reps, suppressMessages( 
#   sim_rep(points=100, b=bflat, maxr=25, mnx, mny, mnerr, sdx, sdy, sderr)
# ))))
# colnames(flat100) <- column_names  

# Flat200
# new_flat200 <- data.frame(t(pbreplicate(reps, suppressMessages(
#   sim_rep(points=200, b=bflat, maxr=25, mnx, mny, mnerr, sdx, sdy, sderr)
# ))))
# colnames(new_flat200) <- column_names

#### Option 2: Read in pre-simulated data to save time
flat50 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/flat50.csv")
flat100 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/flat100.csv")
flat200 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/flat200.csv")



# Plots
boxplot(flat50$r_hn_true, flat50$r_hn_err, flat50$r_hr_true, flat50$r_hr_err,
        names=rep(c("true", "error"), 2),
        col=rep(c("orange", "skyblue"), each=2),
        ylab="Effective radius (m)", main="Flat ground, 50 points")
legend("topright", c("Half normal", "Hazard rate"), fill=c("orange", "skyblue"))

boxplot(flat50$diff_hn, flat100$diff_hn, flat200$diff_hn, 
        flat50$diff_hr, flat100$diff_hr, flat200$diff_hr,
        names=rep(c(50,100,200), 2),
        col=rep(c("orange", "skyblue"), each=3),
        xlab="Points", ylab="Error (m)", main="Flat ground")
legend("topleft", c("Half normal", "Hazard rate"), fill=c("orange", "skyblue"))
lines(c(0,7), rep(mnerr,2), col="magenta")


# Calculate the proportion of simulations where HR has a lower AIC than HN:
mean(flat50$aic_hr_true < flat50$aic_hn_true)
mean(flat50$aic_hr_err < flat50$aic_hn_err)
mean(flat100$aic_hr_true < flat100$aic_hn_true)
mean(flat100$aic_hr_err < flat100$aic_hn_err)
mean(flat200$aic_hr_true < flat200$aic_hn_true)
mean(flat200$aic_hr_err < flat200$aic_hn_err)
#Hazard rate almost always preferred



#######################################################################################################
################################# Sloping surfaces. ##################################################

bslop <- subset(modcoef, location_type=="sloping")[, 2:4] %>%
  apply(2, mean)


### Option 1: Run the simulation from scratch
# Slop 50
# slop50 <- data.frame(t(pbreplicate(reps, suppressMessages(
#   sim_rep(points=50, b=bslop, maxr=15, mnx, mny, mnerr, sdx, sdy, sderr)
# ))))
# 
# colnames(slop50) <- column_names  

# Slop 100
# slop100 <- data.frame(t(pbreplicate(reps, suppressMessages( 
#   sim_rep(points=100, b=bslop, maxr=15, mnx, mny, mnerr, sdx, sdy, sderr)
# ))))
# 
# colnames(slop100) <- column_names  

# Slop 200
# slop200 <- data.frame(t(pbreplicate(reps, suppressMessages(
#   sim_rep(points=200, b=bslop, maxr=15, mnx, mny, mnerr, sdx, sdy, sderr)
# ))))
# 
# colnames(slop200) <- column_names  

#### Option 2: Read in pre-simulated data to save time
slop50 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/slop50.csv")
slop100 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/slop100.csv")
slop200 <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/slop200.csv")


# Plots
boxplot(slop50$r_hn_true, slop50$r_hn_err, slop50$r_hr_true, slop50$r_hr_err,
        slop100$r_hn_true, slop100$r_hn_err, slop100$r_hr_true, slop100$r_hr_err,
        slop200$r_hn_true, slop200$r_hn_err, slop200$r_hr_true, slop200$r_hr_err,
        names=c("50_true", "50_error","50_true", "50_error",
                "100_true", "100_error","100_true", "100_error",
                "200_true", "200_error","200_true", "200_error"), 
        col=rep(c("orange", "skyblue"), each=2),
        ylab="Effective radius (m)", main="Sloping ground")
legend("topright", c("Half normal", "Hazard rate"), fill=c("orange", "skyblue"))

boxplot(slop50$diff_hn, slop100$diff_hn, slop200$diff_hn, 
        slop50$diff_hr, slop100$diff_hr, slop200$diff_hr,
        names=rep(c(50,100,200), 2),
        col=rep(c("orange", "skyblue"), each=3),
        xlab="Points", ylab="Error (m)", main="Sloping ground")
legend("topleft", c("Half normal", "Hazard rate"), fill=c("orange", "skyblue"))
lines(c(0,7), rep(mnerr,2), col="magenta")



# Calculate the proportion of simulations where HR has a lower AIC than HN:
mean(slop50$aic_hr_true < slop50$aic_hn_true)
mean(slop50$aic_hr_err < slop50$aic_hn_err)
mean(slop100$aic_hr_true < slop100$aic_hn_true)
mean(slop100$aic_hr_err < slop100$aic_hn_err)
mean(slop200$aic_hr_true < slop200$aic_hn_true)
mean(slop200$aic_hr_err < slop200$aic_hn_err)
#Hazard rate almost always preferred



# ##############################################################################################################################################################
# ######################################################## UPDATED: Speed simulation function ##################################################################

# The original script introduced errors at both deployment-level and observational-level. Here, we amend the original to replace observational-level errors with
# sequence-level errors (seq_error_updated): a shared error term for all observations within the same sequence. This reflects real-world scenarios where 
# measurements within the same sequence are likely to be influenced by the same systematic biases or environmental factors, creating correlation.


# diff (error in radius predction) currently in centimetres so convert to metres
print(head(predval$diff))
predval$diff <- predval$diff/100


# PART 1: Data Preparation
lmrsum_updated <- summary(lmer(as.numeric(diff) ~ 1 + (1 | deployment), data = predval))
depmn_updated <- lmrsum_updated$coefficients[1]
depsd_updated <- sqrt(lmrsum_updated$varcor$deployment[1])

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

# PART 2: Update the error propagation and generate speeds

### Option 1: Run the simulation from scratch
# speeds_err_500_updated <- pbreplicate(500, {
#   dep_error_updated <- rnorm(length(deps_updated), depmn_updated, depsd_updated) * 0.01
#   seq_error_updated <- rnorm(length(unique(posdat_mov_updated$sequence_id)), mean = 0, sd = 0.1)
# 
#   posdat_try_updated <<- posdat_mov_updated
#   posdat_try_updated$radius_err <- posdat_try_updated$radius +
#     dep_error_updated[match(posdat_try_updated$new_seq, deps_updated)] +
#     seq_error_updated[match(posdat_try_updated$sequence_id, unique(posdat_mov_updated$sequence_id))]
# 
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


#### Option 2: Read in pre-simulated data to save time
# read in speed err
hh_speeds_err_500_df_alt <- read.csv(
  "https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/speeds_err_hh.csv")

speeds_err_500_updated <- t(as.matrix(hh_speeds_err_500_df))

# read in speed vals
seqdat_try_updated <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/seqdat_try_hh.csv")



# PART 3: Visualisation and True Speed Calculation
trspd_updated <- 1 / seqdat1_updated$speed[is.finite(seqdat1_updated$speed) & seqdat1_updated$pixdiff > 10]
truespeed_updated <- 1 / mean(trspd_updated[is.finite(trspd_updated)])

par(mfrow = c(1, 2))
plot(seqdat_try_updated$pixdiff, seqdat_try_updated$speed, log = "xy", xlab = "Pixel difference", ylab = "Speed")
plot(seqdat_try_updated$pixdiff, seqdat_try_updated$speed_err, log = "xy", xlab = "Pixel difference", ylab = "Speed")

boxplot(speeds_err_500_updated, names = c("500_rep_updated"), ylab = "Estimated speed (m/s)")
abline(h = truespeed_updated, lwd = 2, col = "red")



# ##############################################################################################################################################################
# ############################################################### ORIGINAL - Speed simulation function ####################################################################



# Here, the previous (original) speed simulation function has been included (but commented out) for reference. 
# As a reminder, here, the variance is calculated at two levels: 
# deployment (deployment-level errors) and observational (observation-level errors).

# Outputs key parameters:
# depmn (mean deployment-level error),
# depsd (standard deviation of deployment-level error),
# obssd (observation-level error variance).

# PART 1 - Data preparation
### Mixed effects models for the deployment and observational level error segregation
lmrsum <-summary(lmer(as.numeric(diff) ~ 1+ (1|deployment), data= predval))
lmrsum

### Extract the coefficient from the mixed model to be used in the simulation error
depmn <-lmrsum$coefficients[1]
depsd <-sqrt(lmrsum$varcor$deployment[1])
obssd <-lmrsum$sigma

### Unique deployment level for the deployment level error
deps <- unique(predval$deployment)
probdep <-as.numeric(xtabs(~deployment, predval)/nrow(predval))

### Load the speed sequences data and make a unique column for the sequence with the camera id
head(posdat_mov)
posdat_mov$uid <- paste(posdat_mov$siteid, posdat_mov$sequence_id, sep = "-")
names(posdat_mov)[7] <-"seq_id"
names(posdat_mov)[20] <-"sequence_id"

### Adding a column for matching the folder name and prepare a new column
s1 <-as.numeric(xtabs(~sequence_id, posdat_mov))
posdat_mov$new_seq <-rep(sample(deps, 675,replace = T, prob = probdep), times= s1)
posdat_mov$new_seq1 <- paste(posdat_mov$sequence_id,posdat_mov$new_seq, sep = "-")
names(posdat_mov)[20] <-"new_seq_id"
names(posdat_mov)[22] <-"sequence_id"

seqdat <- seq.summary(posdat_mov)
seqdat1 <-subset(seqdat, seqdat$timediff<2000)

avgpixdif <- 400
pixdiff <- seq.data(posdat_mov)$pixdiff

# # PART 2 - ERROR PROPAGATION AND GENERATING SPEEDS
#
# ### Option 1: Run the simulation from scratch
# speeds_err_500 <- pbreplicate(500, {
#   dep_error <- rnorm(length(deps), depmn, depsd) * 0.01
#   obsd <- (posdat_mov$number_of_points) * obssd * pixdiff / avgpixdif
#   obsd[is.na(obsd)] <- 0
#   obsd[obsd > 1] <- 1
# 
#   posdat_try <<- posdat_mov
#   posdat_try$radius_err <- posdat_try$radius +
#     dep_error[match(posdat_try$new_seq, deps)] +
#     rnorm(nrow(posdat_try), sd = obsd)
# 
#   # Save posdat_try globally
#   assign("posdat_try", posdat_try, envir = .GlobalEnv)
# 
# 
#   seqdat_try <<- seq.summary(posdat_try)
# 
#   spd <- sample(1 / (seqdat_try$speed_err[is.finite(seqdat_try$speed_err) &
#                                             seqdat_try$pixdiff > 10 & seqdat_try$timediff < 2000]), size = 500, replace = FALSE)
#   1 / mean(spd[is.finite(spd)])
# })
# 
# 
# 
# 
# # save the outputs
# speeds_err_500_df_original <- as.data.frame(t(speeds_err_500))
# 
# write.csv(
#   speeds_err_500_df_original,
#   file = "/Users/harryjobanputra/Documents/ZSL - Research assistant/WP3 - Position & Speed in Camera Traps/calculated_data_to_speed_scripts/speeds_err_hh_ORIGINAL.csv",
#   row.names = FALSE)
# 
# 
# write.csv(
#   seqdat_try,
#   file = "/Users/harryjobanputra/Documents/ZSL - Research assistant/WP3 - Position & Speed in Camera Traps/calculated_data_to_speed_scripts/seqdat_try_hh_ORIGINAL.csv",
#   row.names = FALSE
# )


#### Option 2: Read in pre-simulated data to save time
seqdat_try <- read.csv("https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/seqdat_try_hh_ORIGINAL.csv")


speeds_err_500_df_original <- read.csv(
  "https://raw.githubusercontent.com/harryjobann/camcal_validation_2025/refs/heads/main/HampsteadHeath/Simulated_data/speeds_err_hh_ORIGINAL.csv")

speeds_err_500 <- t(as.matrix(speeds_err_500_df_original))


# PART 3 - VISUALISATION & CALCULATING TRUE SPEED
trspd <- 1/seqdat1$speed[is.finite(seqdat1$speed) & seqdat1$pixdiff>10]
truespeed <- 1/(mean(trspd[is.finite(trspd)]))

par(mfrow =c(1,2))
plot(seqdat_try$pixdiff, seqdat_try$speed, log="xy", xlab= "Pixel difference", ylab= "Speed")
plot(seqdat_try$pixdiff, seqdat_try$speed_err, log="xy", xlab= "Pixel difference", ylab= "Speed")

boxplot(speeds_err_500, names= c("500_rep"), ylab= "Estimated speed (m/s)")
abline(h= truespeed, lwd=2, col="red")







# ###############################################################################################################################################
# ###############################################################################################################################################
# ########################################################### VERIFICATION OF UPDATES ###########################################################

# Here, we compare the error propagation methods between the original (deployment level & observation level error) 
# and updated (deployment level & sequence-level error)


# Check speed_err in both dataframes
head(seqdat_try$speed_err)
head(seqdat_try_updated$speed_err)

# Summary of speed errors
summary_original <- summary(seqdat_try$speed_err)
summary_updated <- summary(seqdat_try_updated$speed_err)

# Print summaries to console for quick inspection
cat("Original Method - Speed Error Summary:\n")
print(summary_original)

cat("\nUpdated Method - Speed Error Summary:\n")
print(summary_updated)


par(mfrow =c(1,2))

# Plot 1 - Speed
plot(seqdat_try$pixdiff, seqdat_try$speed, log="xy",
     xlab= "Pixel difference",
     main = "Speed",
     ylab= "Speed")


# Plot 2: Compare Speed Errors 
plot(seqdat_try$pixdiff, seqdat_try$speed_err,
     log = "xy",
     col = "blue",
     pch = 16,
     cex = 0.6,
     xlab = "Pixel Difference",
     ylab = "Speed Error",
     main = "Speed Error Comparison")
points(seqdat_try_updated$pixdiff, seqdat_try_updated$speed_err,
       col = "red",
       pch = 17,
       cex = 0.6)
legend("topright",
       legend = c("Original", "Updated"),
       col = c("blue", "red"),
       pch = c(16, 17),
       bty = "n")



# Filter data to remove NA vals
filtered_original <- subset(seqdat_try, is.finite(speed_err) & is.finite(pixdiff))
filtered_updated <- subset(seqdat_try_updated, is.finite(speed_err) & is.finite(pixdiff))


# summary stats for speed_err in original method
summary_stats_original <- summary(filtered_original$speed_err)
sd_original <- sd(filtered_original$speed_err, na.rm = TRUE)
cat("\nSummary statistics for Original Method:\n")
print(summary_stats_original)
cat("Standard Deviation:", sd_original, "\n")

# summary stats for speed_err in updated method
summary_stats_updated <- summary(filtered_updated$speed_err)
sd_updated <- sd(filtered_updated$speed_err, na.rm = TRUE)
cat("\nSummary statistics for Updated Method:\n")
print(summary_stats_updated)
cat("Standard Deviation:", sd_updated, "\n")


# linear models
lm_original <- lm(speed_err ~ pixdiff, data = filtered_original)
lm_updated <- lm(speed_err ~ pixdiff, data = filtered_updated)

print(summary_original <- summary(lm_original))
print(summary_updated <- summary(lm_updated))



#### Results of amending the error propogation method:
# --> The mean error is lower in the updated method (1.73 vs. 2.44), and the median error shows a substantial reduction (0.11 vs. 0.86).
# --> IQR is significantly smaller in the updated method. 





