# WaggleDance
Waggle Dance Analysis scripts used in Sreelakshmi Suresh's MS thesis

######## This script has been adapted for a new dataset by Akul Gulrajani and Sreelakshmi Suresh (ssuresh@illinois.edu) from edits made by Douglas Sponsler and Reed Johnson. 
######## Unused portions of the original script have been commented out but not deleted. Additions are followed by the tag # Sponsler:

# Electronic Supplementary Material 4 - R script to simulate dances from known waggle
# dance durations and headings.
# ------------------------------------------------------------------------------------

# Article title: Incorporating variability in honey bee waggle dance decoding improves
# the mapping of communicated resource locations

# Journal: Journal of Comparative Physiology A

# Authors: Roger Schurch, Margaret J. Couvillon, Dominic D. R. Burns, Kiah
# Tasman, David Waxman and Francis L. W. Ratnieks

# Corresponding author: Roger Schurch, Evolution, Behaviour and
# Environment, School of Life Sciences, University of Sussex, Brighton,
# BN1 9QG, United Kingdom, R.Schuerch@sussex.ac.uk

# script files that go with this script:
# ESM_3.jag

# data files that go with this script:
# ESM_5.csv

# the file will output a comma separated value file and an ASC raster file
# you must create a "data" folder within the folder from which you are running the scripts
# or adapt the paths for the instructions below
