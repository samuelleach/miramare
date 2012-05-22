#!/bin/bash

# Test make_window with Gaussian beams
./bin/make_window parfiles/make_window.par

./bin/cl2pf ./bin/make_window --input_fwhm_arcmin 30. --input_beam_file '' --nside 512 --target_fwhm_arcmin 60. --signal_beam_file data/sigbeam.txt --variance_beam_file data/varbeam.txt --inversevariance_beam_file data/invvarbeam.txt --verbose T


# Test make_window on the WMAP beams.
./bin/make_window parfiles/make_window_K1.par
./bin/make_window parfiles/make_window_Ka1.par
./bin/make_window parfiles/make_window_Q1.par
./bin/make_window parfiles/make_window_Q2.par
./bin/make_window parfiles/make_window_V1.par
./bin/make_window parfiles/make_window_V2.par
./bin/make_window parfiles/make_window_W1.par
./bin/make_window parfiles/make_window_W2.par
./bin/make_window parfiles/make_window_W3.par
./bin/make_window parfiles/make_window_W4.par

# Test ud_grade_mask on WMAP mask.
./bin/ud_grade_mask parfiles/ud_grade_mask.par

# Test lincom_map on WMAP K and Ka bands.
./bin/lincom_map parfiles/lincom_map.par

# Test combination of masks
./bin/combine_mask parfiles/combine_mask.par
#cl2pf ./bin/combine_mask --nmaps 3 --infile\(1\) miramare_test/weight/hits_150GHz.fits --infile\(2\) miramare_test/weight/hits_250GHz.fits --infile\(3\) miramare_test/weight/hits_410GHz.fits --threshold 15000 --outfile \!miramare_test/mask/composite_2048.fits --verbose  T
