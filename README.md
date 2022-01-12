# PyRES
Python code for repeating earthquakes

# Abstract
The software package PyRES is an open source seismological software. It consists of some utilities for data preparation (utils directory) and the main program (xxxx.py). PyRES is designed to discriminate repeating earthquakes starting from a family of candidate repeating earthquakes, based on the cross-correlation values and S-P time difference between pairs of earthquakes (estimated using cross spectrum). 

# Motivation and significance
The code PyRES is inspired to previously published methods that combine both seismic waveform similarity, using cross-correlation function, and differential S-P travel time measured at each seismic station (Chen et al., 2008 and Shakibay Senobari and Funning, 2019).
The code is versatile and works with and without P and S-wave phase pickings. At the moment the reading of phases in archive Hypoellipse and Hypoinverse formats are implemented. The code has been tested using synthetic and real data, providing accurate results. It contributes to the implementations of open-source Python packages in seismology aiming to support the activities of researchers and the reproducibility of scientific results.

# Requirements
Python installed to run the program (version 3.0 or more) is required [http://python.org]. Dependancies include Obspy, Numpy, Scipy, Scikit-learn, mtspec (©2009-2016, Lion Krischer, Moritz Beyreuther) and Matplotlib libraries. to do

# How to use it
A complete step-by-step working flow is reported in the jupiter notebook named
PyRES.ipynb

# References

Sugan, M., Campanella, S., Vuan, A., Shakibay Senobari, N., (2022). A Python Code for Detecting Repeating Earthquakes from Self-similar Waveforms (PyRES). Submitted to BSSA ?

Shakibay Senobari, N., Funning G. J., (2019), Widespread Fault Creep in the Northern San Francisco Bay Area Revealed by Multistation Cluster Detection of Repeating Earthquakes, Geophysical Research Letters, 46(12), 6425-6434, https://doi.org/10.1029/2019GL082766.

Chen, K. H., Nadeau, R. M. and Rau, R.-J. (2008). Characteristic repeating microearthquakes on an arc-continent collision boundary – the Chihshang fault of eastern Taiwan, Earth Planet. Sci. Lett., 276, 262–272, https://doi.org/10.1016/j.epsl.2008.09.021.

