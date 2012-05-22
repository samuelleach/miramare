PRO make_rms

; make rms files for I,Q,U using provided noise matrices
; For rms mode, we are just taking the diagonal elements
; of the I,Q,U matrix i.e. no cross-correlations between I,Q,U
; For I, assume this has been smoothed and white noise added,
; so that this is known beforehand
;
; This procedure also outputs the I,Q,U N_inv 6 column matrix
; 
; So the steps are:-
; 1. Sum up N^-1 subpixels for the required Nside
; 2. Inver each pixel matrix
; 3. Take diagonal parts to get Q,U rms noises
;
; Notes:-
; 1. T and pol matrices are supplied separately by Sara Ricciardi
; 2. Maps must be in NESTED format
; 3. For cov matrix output, TT is set to the noise_level^2
;    and TQ=TU=0
;
; 24-Oct-2008  C. Dickinson  1st go
; 27-Nov-2008  C. Dickinson  Added output of cov matrix
;---------------------------------------------------------------
; set data directories, inputs etc.
indir  = '/attic/cdickins/WG2/FFP2/data/'                     ; input directory
outdir = '/attic/cdickins/WG2/FFP2/proc2/'                     ; output directory
T_files = file_search(indir, 'T_cov_ffp2_*.fits',count=nmaps_T)  ; list of files
pol_files = file_search(indir, 'pol_cov_ffp2_*.fits',count=nmaps_pol)  ; list of files
nside_out = 64L                                               ; output Nside
noise_level = 2.                                              ; noise level for I maps

; make sure we have the same number of T and pol cov files!
IF (nmaps_T NE nmaps_pol) THEN BEGIN
print,''
print, '***NOT THE SAME NUMBER OF T AND POL MATRIX FILES***
print, '***Aborting....'
print, ''
ENDIF



; do each frequency separately

FOR i=0L, nmaps_T-1 DO BEGIN
T_infile = T_files[i]
pol_infile = pol_files[i]
read_fits_map, T_infile, Tcov
read_fits_map, pol_infile, polcov
npix1 = n_elements(Tcov[*,0])   ; number of pixels on input
print, 'Doing map ', i+1, ' out of ', nmaps_T, '.....' 

; construct N for each pixel at a time
invcov1 = dblarr(npix1,6)
FOR j=0L, npix1-1 DO BEGIN
cov1 = [[Tcov[j,0], Tcov[j,1], Tcov[j,2]], [Tcov[j,1], polcov[j,0], polcov[j,1]], [Tcov[j,2],polcov[j,1],polcov[j,2]]] 

; invert this matrix
temp = invert(cov1,status, /double)
invcov1[j,*] = [temp[0],temp[1],temp[2],temp[4],temp[5],temp[8]]
IF (status EQ 1) THEN print, 'Inversion failed....'

ENDFOR

; Now sum up sub-pixels of N^-1 to the output nside (this assumes maps
; are in NESTED format!)
npix2 = 12L * nside_out^2
invcov2 = dblarr(npix2,6)
rat = npix1/npix2
FOR j=0L, npix2-1 DO BEGIN
index = j*rat
invcov2[j,0] = total(invcov1[index:index+rat-1,0],1,/double)
invcov2[j,1] = total(invcov1[index:index+rat-1,1],1,/double)
invcov2[j,2] = total(invcov1[index:index+rat-1,2],1,/double)
invcov2[j,3] = total(invcov1[index:index+rat-1,3],1,/double)
invcov2[j,4] = total(invcov1[index:index+rat-1,4],1,/double)
invcov2[j,5] = total(invcov1[index:index+rat-1,5],1,/double)
ENDFOR

; invert this matrix to get N again at the lower resolution and then
; get elements for the rms map
cov3= dblarr(npix2,6)
rmsmap = dblarr(npix2,3)
FOR j=0L, npix2-1 DO BEGIN
invcov3 = [[invcov2[j,0], invcov2[j,1], invcov2[j,2]], [invcov2[j,1], invcov2[j,3], invcov2[j,4]], [invcov2[j,2],invcov2[j,4],invcov2[j,5]]] 

; invert this matrix
temp = invert(invcov3,status, /double)
cov3[j,*] = [temp[0],temp[1],temp[2],temp[4],temp[5],temp[8]]
IF (status EQ 1) THEN print, 'Inversion failed....'





; put each entry into the rms map, converting to uK_CMB (rms)
temp2 = strpos(T_infile,'GHz')
freq_str = strmid(T_infile,temp2-3,3)
freq = float(freq_str)


rmsmap[j,0] =  noise_level
rmsmap[j,1] =  sqrt(temp[4])  * 1.0e6 * planckcorr(freq)
rmsmap[j,2] =  sqrt(temp[8])  * 1.0e6 * planckcorr(freq)

ENDFOR


;STOP

; output rms map
outfile = outdir + 'rms_TQU_' + freq_str + 'GHz_h' + strcompress(string(nside_out,format='(i3)'),/remove_all) + '.fits'
write_tqu, outfile, rmsmap, /nested, units='uK_CMB'


; output covariance columns with 
; i) I set to the noise_level^2
; ii) TQ and TU are set to 0
cov4 = cov3 * (1.0e6 * planckcorr(freq))^2
cov4[*,0] = noise_level^2
cov4[*,1] = 0.
cov4[*,2] = 0.

mapstruct = create_struct('HDR', [''], 'TT_COV', cov4[*,0], 'TQ_COV', cov4[*,1], 'TU_COV', cov4[*,2], 'QQ_COV', cov4[*,3], 'QU_COV', cov4[*,4], 'UU_COV', cov4[*,5])

outfile = outdir + 'cov_TQU_' + freq_str + 'GHz_h' + strcompress(string(nside_out,format='(i3)'),/remove_all) + '.fits'
write_fits_sb, outfile, 0, mapstruct, coord='G', /nest




ENDFOR








;---------------------------------------------------------------
; Finish up
STOP
END
