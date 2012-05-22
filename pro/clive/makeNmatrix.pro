PRO makeNmatrix

; procedure to run monte-carlo simulations based on a noise matrix
; to smooth a map and create the appropriate matrix
;
; 14-Jan-2009  C. Dickinson   1st go
;------------------------------------------------------------------
; set data directories, inputs etc.
indir  = '/attic/cdickins/WG2/FFP2/data/'                     ; input directory
outdir = '/attic/cdickins/WG2/FFP2/mc/'                     ; output directory
mapfiles = file_search(indir, 'TQU_ffp2_*.fits',count=nmaps)  ; list of files
T_files = file_search(indir, 'T_cov_ffp2_*.fits',count=nmaps_T)  ; list of files
pol_files = file_search(indir, 'pol_cov_ffp2_*.fits',count=nmaps_pol)  ; list of files
nside_out = 64L                                               ; output Nside
inbeams = [52.8,33,39.6,30.6,24,21,14,13.2,10,7.1,5,5,5,5]    ; beam fwhm on input (arcmin)
;inbeams = [5,5]                                               ; for 545/857GHz!!!!
resout = 180.                                                 ; output resolution (for T map)
noise_level = 2.                                              ; white noise level to add (uK_CMB)
niter=500                                                       ; number of MC iterations to get good noise statistics
seed = ''

;--------------------------------------------------------------------- 
; 1. Make a noise realization for each map and smooth it

FOR i=1L, nmaps-1 DO BEGIN
print, ''
print, '***Running map ', i+1, ' out of ', nmaps

    FOR j=0L, niter-1 DO BEGIN
print, '***Running iteration ', j+1, ' out of ', niter
infile = mapfiles[i]
temp = strsplit(infile, '.', /extract)
ntemp = n_elements(temp)
nside_str = temp[ntemp-2]
freq_str = strmid(infile,strpos(infile,'GHz')-3,3) 
outfile = outdir + 'noise' + strcompress(string(j,format='(i6)'),/remove_all) + '_' + freq_str + 'GHz'  + '.' + strcompress(string(nside_str,format='(i5)'),/remove_all) + '.' + strcompress(string(resout,format='(i6)'),/remove_all) + 'arcmin' + '.fits'

;read_fits_map, infile, inmap    ; inmap here is the TQU map
T_infile = T_files[i]
pol_infile = pol_files[i]
read_fits_map, T_infile, Tcov
read_fits_map, pol_infile, polcov
npix_in = n_elements(Tcov[*,0])   ; number of pixels on input

; make a noise realization
nside_in = long(nside_str)
npix_in = 12L * nside_in^2
npol = 3
noise_realization = randomn(seed,npix_in,npol)

; calculate a realization with correct normalization and correlations
noise_map = fltarr(npix_in,npol)
nmatrix = fltarr(npol,npol)
FOR k=0L, npix_in-1 DO BEGIN

; make noise matrix for this pixel
nmatrix[0,0] = Tcov[k,0]
nmatrix[0,1] = Tcov[k,1]
nmatrix[1,0] = Tcov[k,1]
nmatrix[0,2] = Tcov[k,2]
nmatrix[2,0] = Tcov[k,2]
nmatrix[1,1] = polcov[k,0]
nmatrix[1,2] = polcov[k,1]
nmatrix[2,1] = polcov[k,1]
nmatrix[2,2] = polcov[k,2]

; Cholesky decompose and multiply to get noise realization
choldc,nmatrix,rms
nmatrix[0,0] = rms[0]
nmatrix[1,1] = rms[1]
nmatrix[2,2] = rms[2]
noise_map[k,*] = nmatrix[*,*] ## noise_realization[k,*]
ENDFOR

; convert to uK_CMB
noise_map = noise_map * 1.0e6 / planckcorr(float(freq_str))

; smooth the map
write_tqu, outfile, noise_map, /nest, units='uK_CMB'
fwhm = sqrt(resout^2 - inbeams[i]^2)
lmax = (3*nside_out)-1L
npix_out = 12L*nside_out^2
ismoothing, outfile, outfile, fwhm=fwhm, /nest, simul=2, lmax=lmax, tmpdir='/attic/cdickins/tmp/'



; degrade the map
infile = outfile
outfile = outdir + 'noise' + strcompress(string(j,format='(i6)'),/remove_all) + '_' + freq_str + 'GHz'  + '.' + strcompress(string(nside_out,format='(i5)'),/remove_all) + '.' + strcompress(string(resout,format='(i6)'),/remove_all) + 'arcmin' + '.fits'
ud_grade, infile, outfile, nside_out=nside_out,order_in='nest',order_out='nest'


; store the noise realization
read_fits_map, outfile, noise_map
IF (j EQ 0) THEN noise_maps = fltarr(npix_out,npol,niter)
noise_maps[*,*,j] = noise_map[*,*]

ENDFOR

; now calculate the noise matrix based on these simulations
newmatrix = fltarr(npix_out,6)
print, '***Calculating new noise matrix...'
FOR k=0L, npix_out-1 DO BEGIN
newmatrix[k,0] = variance(noise_maps[k,0,*])
newmatrix[k,3] = variance(noise_maps[k,1,*])
newmatrix[k,5] = variance(noise_maps[k,2,*])
newmatrix[k,1] = mean( (noise_maps[k,0,*]-mean(noise_maps[k,0,*])) * (noise_maps[k,1,*]*mean(noise_maps[k,0,*])) )
newmatrix[k,2] = mean( (noise_maps[k,0,*]-mean(noise_maps[k,0,*])) * (noise_maps[k,2,*]*mean(noise_maps[k,2,*])) )
newmatrix[k,4] = mean( (noise_maps[k,1,*]-mean(noise_maps[k,1,*])) * (noise_maps[k,2,*]*mean(noise_maps[k,2,*])) )



ENDFOR


; output the new noise matrix
mapstruct = create_struct('HDR', [''], 'TT_COV', newmatrix[*,0], 'TQ_COV', newmatrix[*,1], 'TU_COV', newmatrix[*,2], 'QQ_COV', newmatrix[*,3], 'QU_COV', newmatrix[*,4], 'UU_COV', newmatrix[*,5])

outfile = outdir + 'new_cov_TQU_' + freq_str + 'GHz_' + strcompress(string(resout,format='(i6)'),/remove_all) + 'arcmin' + '_h' + strcompress(string(nside_out,format='(i3)'),/remove_all) + '.fits'
write_fits_sb, outfile, 0, mapstruct, coord='G', /nest

; remove temporary noise files
str = 'rm ' + outdir + 'noise*.fits'
spawn, str


ENDFOR





















;------------------------------------------------------------------
STOP
END
