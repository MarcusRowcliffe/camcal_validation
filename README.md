# camcal_validation
Scripts and two datasets to validate distance measured from camera traps using calibration tools. Original analysis conducted by Nilan Chaterjee, found here: https://github.com/nilanjanchatterjee/camcal_validation

## Details

# Dataset 1: Hampstead Heath (HH)

1. HH_validation_code.R - R script to simulate the effective detection distance and speed from the camera-trap detections in Hampstead Heath. 
2. Site_digitisation_data.csv - Dataset contians the digitised coordinate from images for a reliastic representation of digitised object location in camera traps
3. pole_11_mod_param.csv - Model parameters from the non-linear least square and surface details of the 10 camera locations for validation purpose
4. camera_deloyment_validation.csv - Dataset with the 10 cameras deployed for the validation of the calibration model


# Dataset 2: Hoge Veluwe (NLD)

1. NLD_validation_code.R - R code to simulate the effective detection distance and  speed from the camera-trap detections in Hoge Veluwe.
2. 

# Shared data and code for both datasets

1. Speed_seq_data.csv - Movement digitisation data to generate speed estimates, used in process of validating both datasets.
2. CTtracking_err.r - Package containing functions to analyse the camera calibration dataset and estimate radius and speed with radius_err and speed_err.
