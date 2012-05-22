#!/bin/bash
# Get test data - WMAP data, Simulated Planck IQU covariance matrices.

mkdir wmap
cd wmap

#WMAP beams
wget http://lambda.gsfc.nasa.gov/data/map/dr4/ancillary/beams/wmap_ampl_bl_K1_7yr_v4.txt
wget http://lambda.gsfc.nasa.gov/data/map/dr4/ancillary/beams/wmap_ampl_bl_Ka1_7yr_v4.txt
wget http://lambda.gsfc.nasa.gov/data/map/dr4/ancillary/beams/wmap_ampl_bl_Q1_7yr_v4.txt
wget http://lambda.gsfc.nasa.gov/data/map/dr4/ancillary/beams/wmap_ampl_bl_Q2_7yr_v4.txt
wget http://lambda.gsfc.nasa.gov/data/map/dr4/ancillary/beams/wmap_ampl_bl_V1_7yr_v4.txt
wget http://lambda.gsfc.nasa.gov/data/map/dr4/ancillary/beams/wmap_ampl_bl_V2_7yr_v4.txt
wget http://lambda.gsfc.nasa.gov/data/map/dr4/ancillary/beams/wmap_ampl_bl_W1_7yr_v4.txt
wget http://lambda.gsfc.nasa.gov/data/map/dr4/ancillary/beams/wmap_ampl_bl_W2_7yr_v4.txt
wget http://lambda.gsfc.nasa.gov/data/map/dr4/ancillary/beams/wmap_ampl_bl_W3_7yr_v4.txt
wget http://lambda.gsfc.nasa.gov/data/map/dr4/ancillary/beams/wmap_ampl_bl_W4_7yr_v4.txt

#WMAP mask
wget http://lambda.gsfc.nasa.gov/data/map/dr4/ancillary/masks/wmap_temperature_analysis_mask_r9_7yr_v4.fits

#WMAP IQU MAPS
wget http://lambda.gsfc.nasa.gov/data/map/dr4/skymaps/7yr/raw/wmap_da_iqumap_r9_7yr_K1_v4.fits
wget http://lambda.gsfc.nasa.gov/data/map/dr4/skymaps/7yr/raw/wmap_da_iqumap_r9_7yr_Ka1_v4.fits
#wget http://lambda.gsfc.nasa.gov/data/map/dr4/skymaps/7yr/raw/wmap_da_iqumap_r9_7yr_Q1_v4.fits
#wget http://lambda.gsfc.nasa.gov/data/map/dr4/skymaps/7yr/raw/wmap_da_iqumap_r9_7yr_Q2_v4.fits
#wget http://lambda.gsfc.nasa.gov/data/map/dr4/skymaps/7yr/raw/wmap_da_iqumap_r9_7yr_V1_v4.fits
#wget http://lambda.gsfc.nasa.gov/data/map/dr4/skymaps/7yr/raw/wmap_da_iqumap_r9_7yr_V2_v4.fits
#wget http://lambda.gsfc.nasa.gov/data/map/dr4/skymaps/7yr/raw/wmap_da_iqumap_r9_7yr_W1_v4.fits
#wget http://lambda.gsfc.nasa.gov/data/map/dr4/skymaps/7yr/raw/wmap_da_iqumap_r9_7yr_W2_v4.fits
#wget http://lambda.gsfc.nasa.gov/data/map/dr4/skymaps/7yr/raw/wmap_da_iqumap_r9_7yr_W3_v4.fits
#wget http://lambda.gsfc.nasa.gov/data/map/dr4/skymaps/7yr/raw/wmap_da_iqumap_r9_7yr_W4_v4.fits

cd ..

