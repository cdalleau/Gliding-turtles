# Gliding-turtles
Abstract: This repository contains several scripts about the horizontal and vertical analysis of turtle trajectories. This scripts create several products with the CSV and NetCDF format.

I did a three months internship at ICTS SOCIB
The main goal of our project is to analyze the turtles’ trajectory. For this, we wrote four scripts to define different types of product:
•	the L0 product contains the raw data
•	the L1_1 product contains an analysis of the horizontal trajectory
•	the L1_2 product contains an analysis of the  vertical trajectory
•	the L2 product contains interpolated positions based on a state-space model.
We also wrote two scripts to view the variables of products and plot the track.
Thus a first script allows the creation of a file containing the raw data. The data are imported from the Wildlife Computers data portal or a folder containing the CSV files. The interesting data are selected, merged together and exported in a L0 product. This product contains the raw data about the latitude (Figure 1 and Figure 2), longitude (Figure 1 and Figure 3), sea temperature (Figure 4), depth (Figure 5) and the battery voltage (Figure 6). The locations are sorted by the location class of Argos and recorded in the variable: location quality. 
Then a second script uses this L0 product to analyze the horizontal trajectory. The data are imported and the interesting data are selected. The data are filtered by several functions:
•	a function filtering the positions on land,
•	the function Argosfilter filtering the positions according to a maximum velocity, the angle and the distance between two positions (Freitas et al., 2008),
•	a function filtering the temperature and the depth according to a maximum and a minimum chosen beforehand.
This process interpolates also (1) the latitude (Figure 7 and Figure 8), longitude (Figure 7 and Figure 9) and voltage battery (Figure 10); (2) calculates the velocity (Figure 11) and the cumulative distance (Figure 12). After these steps the data are exported in a L1_1 product.
The analysis of the vertical trajectory is processed in a third script using the L1_1 product. After the importation and selection of the data, the function diveMove (P. Luque, 2007) are used to:
•	detect the individual dives with their different phases: descent, bottom, ascent and surface (Figure 13 and Figure 14)
•	calculate dive statistics for each dive: distance and duration for each phase (Figure 15, 16, 17 and 18)
This data and the data from L1_1 product are exported in a L1_2 product.
	Finally the data from L1_1 product are imported to interpolate positions (Figure 19, Figure 20 and Figure 21) based on a state-space model (SSM) and the behavior mode is estimated (Figure 22) using the algorithm developed by Jonsen et al. (Jonsen et al, 2005). This script allows the exportation of several plots supplied by the SSM and the data in a L2 product. 

