# Segment-slender-regions
This repository contains an algorithm to segment slender regions in 3D microstructures

Two sample material microstructures are included namely: FOAM_MICROSTRUCTURE.mat and IWP_MICROSTRUCTURE.mat
Please run the MATLAB file named Main_file.m
The required matlab functions to the 'Main_file' are included.

Tips: 1. Studying FOAM_MICROSTRUCTURE is resource intensive. If the user wishes to study the algorithm and has limited computational resources, please use IWP_MICROSTRUCTURE.
      2. While studying FOAM_MICROSTRUCTURE, the most resource-intensive sections are calculation of variables named cross_area, strut_segments and strut_pixels. If you don't           want to calculate these variables, de-activate the corresponding sections and load the MATLAB data files namely cross_area.mat, strut_segments.mat and strut_pixels.mat
      3. While studying IMP_MICROSTRUCTURE, the above-mentioned variables can be quickly calculated.
      4. Please keep all the files in the same folder.

INPUTS: 1. Load the  MATLAB data file with the material microstructure.
        2. Define the threshold values

OUTPUT: 1. If you want to see the segmented struts in FOAM_MICROSTRUCTURE, run the corresponding section in the Main_file that will plot a chosen strut in the microstructure.
        2. If you want to see the segmented struts in IWP_MICROSTRUCTURE, use 'Volume Viewer' app in MATLAB- Import volume               'ori' which is the microstructure and then import labelled volume 'strut_region' which has the segmented struts.
