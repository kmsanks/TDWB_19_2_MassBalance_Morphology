# TDWB_19_2_MassBalance_Morphology
#Introduction
This README file describes kmsanks/TDWB_19_2_MassBalance_Morphology repository of scripts used to calculate changes in morphology and mass balance due to marsh deposition in an experimental setting. The raw experimental data can be found in the Tulane_Sediment_Dynamics_Stratigraphy_TSDS project space at: https://sead2.ncsa.illinois.edu/spaces/5825f529e4b0f3dd19c8d93a. The data used here is TDB-18-1 and TDWB-19-2-Surface-Processes. The processed data used in the code is contained herein. 

The experimental delta data come from two experiments run at the Tulane Sediment Dynamics Laboratory. Both experiments were setup identically, except for TDWB-19-2 (treatment) had marsh proxy deposition, while TDB-18-1 (control) did not. 
For more information on experimental conditions, please see the data repositories hosted at: https://sead2.ncsa.illinois.edu/spaces/5825f529e4b0f3dd19c8d93a.

#What does this repository provide?
This repository contains all of the code relevant to produce the results and figures contained in the manuscript titled "Marsh sedimentation controls delta slope and morphology in an experimental setting (Sanks et al., in review)". The code contained herein was written using MATLAB R2019b. Note: all the "Core Script Files" should run out the box, as relative paths were used in creation of this script.

#Contents
The repository contains three folders:
1. code
2. data
3. figures
Please clone the repository in full in order to use the repo, as the code relies on the .mat files contained in the "data" folder. All data needed to run scripts is contained herein.

#Data: TDWB_19_2_MassBalance_Morphology\data
1. ZD_18.mat - A matrix of size 796x522x560, where 560 is time. Each timestep contains a elevation data collected from the control experiment via LiDAR and post-processed into a 5mmx5mm grid. The data has 796 rows and 522 columns corresponding to basin location.
  *Note that the matrix has 560 timesteps, which corresponds to one LiDAR scan per hour. The first run hour is t = 1, total run time = 560 hours.
2. ZD_19.mat - A matrix of size 750x747x281, where 281 is time. Each timestep contains a elevation data collected from the treatment experiment via LiDAR and post-processed into a 5mmx5mm grid. The data has 750 rows and 747 columns corresponding to basin location. 
  *Note that the matrix only has 281 timesteps because the LiDAR data was collected every other hour. The first run hour is t = 0, total run time = 560 hours.
3. chanMaps_18.mat - A matrix of size 796x522x560. The data contained herein is the same reference frame as described above for ZD_18 (control). This data is a binary matrix, where channel pixels = 1 and non-channel pixels = 0.
4. chanMaps_19.mat - 3. chanMaps_18.mat - A matrix of size 796x522x560. The data contained herein is the same reference frame as described above for ZD_19 (treatment). This data is a binary matrix, where channel pixels = 1 and non-channel pixels = 0.
5. interpolated_marsh_strat_frac.mat - A matrix of 750x747 containing the interpolated fraction of marsh preserved in the stratigraphy (i.e., thickness of marsh deposit at pixel/total stratigraphic thickness of pixel). For interpolation methods see the supporting information from (Sanks et al, in review). 
6. MRD.mat (61x4), GBMD.mat (61x4), mekong.mat (22x4), riogrande.mat (31x4) - Global delta data. These data contain 4 columns. Column 1 is the median elevation of 2 m bins, column 2 is the number of pixels at that elevation, column 3 is the normalized elevation (column 1/channel depth), and column 4 is the probability of a pixel being in that elevation bin (column 2/sum(column 2)). 
  *Note: the data herein comes from Google Earthh Engine ETOPO. The raw data can be found in the corresponding *_raw.csv
  *MRD = Mississippi River Delta, GBMD = Ganges Brahamaputra Meghna River Delta, mekong = Mekong River Delta, and riogrande = Rio Grande River Delta.
  *The elevation data ranges from -1 to 3 channel depths relative to sea level. 

#Core Script Files: TDWB_19_2_MassBalance_Morphology\data
1. area.m - This script calculates area of both the control and treatment experiments for different zones of the delta. 
  -This sciript produces Figure 2a. 
  -Data needed: ZD_18.mat and ZD_19.mat.
2. elevation_histogram.m - This script calculates the distribution of elevations relative to sea level for both the control and treatment experiments.
  -This script produces Figure 2b.
  -Data needed: ZD_18.mat and ZD_19.mat.
3. slope.m - This script calculates the delta top slope for different zones in both the control and treatment experiments using a gradient.
  -This script produces Figure 2c.
  -Data needed: ZD_18.mat and ZD_19.mat
4. radial_elevation.m - This script calculates the mean and standard deviation about the mean for various radial distances from the entrance channel for both the control and treatment experiments.
  -This script produces Figure 2d.
  -Data needed: ZD_18.mat, ZD_19.mat, chanMaps_18.mat, and chanMaps_19.mat.
5. volume.m - This script calculates the 50% delta top, 10% marsh window, and 90% above marsh volume for the control and treatment experiments. See manuscript for further details.
  -This script produces some results in Table 1.
  -Data needed: ZD_18.mat, ZD_19.mat, and interpolated_marsh_strat_frac.mat.
6. trapping_efficiency.m - This script calculates the trapping efficiency for the regions describe above (volume.m), as well as the volume accumulated in the off-shore for both the control and treatment experiments.
  -This script produces some results in Table 1 and Figure 3.
  -Data needed: ZD_18.mat, ZD_19.mat, and interpolated_marsh_strat_frac.mat.
7. delta_hypsometry.m - This script calculates the global delta hypsometry and the experimental delta hypsometry (normalized by channel depth for comparison across scales).
  -This script produces Figure 4.
  -Data needed: MRD.mat, GBMD.mat, mekong.mat, riogrande.mat, ZD_18.mat, and ZD_19.mat.
  
#Supplementary Script Files
1. GEE_delta_data.m - This script creates the matrices of the global delta data (MRD.mat, GBMD.mat, mekong.mat, and riogrande.mat).
  *Note: This script does not needed to be run because the data produced here is already included in the repository.
  *Note: This script will not run out of the box, and you must read through and follow the step by step directions contained in the code in order to replicate the matrices provided. 
    -User must copy/paste data from the delta .csv files included (e.g., MRD_raw.csv, GBMD_raw.csv, mekong_raw.csv, and riogrande_raw.csv).
    -The data here is clipped to -1 to 3 channel depths relative to sea level.
    -Only paste the data within this range. 
    
  #Using this repository
  Clone the repository to your desktop. The scripts can be run in any order, but please note that the "Supplementary Script Files" will not run without some copy and paste from the .csv files provided in the "data" folder. More information can be found in script.  
