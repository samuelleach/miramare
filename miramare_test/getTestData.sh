#!/bin/bash

# Script to download foreground template for miramare simulation mode.

#Unzip other data
cd mask
gunzip *gz
cd ../weight
gunzip *gz
cd ..

cd data

#Dust template in Galactic coords
wget http://people.sissa.it/~leach/miramare/data/dust_carlo_150ghz_512.fits

#Same dust template in RA/Dec
wget http://people.sissa.it/~leach/miramare/data/dust_carlo_150ghz_512_radec.fits

