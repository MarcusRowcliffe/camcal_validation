# camcal_validation
Scripts and two datasets to validate distance measured from camera traps using calibration tools. Original analysis conducted by Nilan Chatterjee, found here: https://github.com/nilanjanchatterjee/camcal_validation

## Details

### Main scripts

1. camcal_validation_HH.R - Amended R script to simulate the effective detection distance and speed from the camera-trap detections in Hampstead Heath.
2. camcal_validation_Hoge_Veluwe.R - R code to simulate the effective detection distance and speed from the camera-trap detections in Hoge Veluwe.

### Dataset 1: Hampstead Heath

1. Site_digitisation_data.csv - Dataset contians the digitised coordinate from images for a reliastic representation of digitised object location in camera traps
2. pole_11_mod_param.csv - Model parameters from the non-linear least square and surface details of the 10 camera locations for validation purpose

### Dataset 2: Hoge Veluwe

1. HogeVeluweDigidat.csv - dataset from HogeVeluwe, containing distances, predicted radius vals, along with x & y coords of digitised poles.

### Shared data required for both datasets

1. Speed_seq_data.csv - Movement digitisation data to generate speed estimates, used in process of validating both datasets.
2. camera_deployment_validation_.csv - distance data used in both sequences


