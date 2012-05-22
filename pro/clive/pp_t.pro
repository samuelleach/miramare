PRO pp_t

; procedure to do the pre-processing for FFP 2 maps for running with
; Commander; This one is just for total-intensity
;
; 1. Downgrade maps to Nside=64
; 2. Smooth T maps only to 3deg resolution
; 3. Convert to uK_CMB and Add 2uK rms white noise to T 
; 4. Just take diagonal elements to make rms file for Commander
; 5. Make beam files for each frequency
;
; Notes:-
;
; -Requires Healpix v2.10
; -Uses split up maps made by Sara Ricciardi
; -Renamed "corrected" files to default filename
; -Same noise level is added to all maps (re-run if you want different
; ones)
; Be careful with inbeams - this has to be correct for the files being
;                           used (do by hand!)
;
; 21-Oct-2008  C. Dickinson.  1st go
;--------------------------------------------------------------------- 
; set data directories, inputs etc.
;indir  = '/attic/cdickins/WG2/FFP2/data/'                     ; input directory
;outdir = '/attic/cdickins/WG2/FFP2/proc/'                     ; output directory
;mapfiles = file_search(indir, 'TQU_ffp2_*.fits',count=nmaps)  ; list of files

indir  = '/attic/cdickins/PSM_e2e/FFP_v2_all_components_outputs/frequency_maps/'                     ; input directory
outdir = '/attic/cdickins/WG2/FFP2/inputs_3deg/'                     ; output directory
mapfiles = file_search(indir, 'coaddedmap_*.fits',count=nmaps)  ; list of files


nside_out = 64L                                               ; output Nside
inbeams = [52.8,33,39.6,30.6,24,21,14,13.2,10,7.1,5,5,5,5]    ; beam fwhm on input (arcmin)
;inbeams = [5,5]                                               ; for 545/857GHz!!!!
resout = 180.                                                 ; output resolution
noise_level = 0.                                              ; white noise level to add (uK_CMB)


; 1. Degrade maps
FOR i=0L, nmaps-1 DO BEGIN

;IF (i NE 11) THEN CONTINUE     ; ****for picking out 1 frequency only!!!!*****
infile = mapfiles[i]
temp = strsplit(infile, '.', /extract)
ntemp = n_elements(temp)
nside_str = temp[ntemp-2]
freq_str_pos = strpos(infile,'GHz')
IF (freq_str_pos LT 0) THEN freq_str_pos = strpos(infile,'Ghz')
freq_str = strmid(infile,freq_str_pos-3,3) 

; make sure freq_str is correct 
FOR j=0L, 5L DO BEGIN
freq_str_a = strmid(freq_str,0,1) 
test = strmatch(freq_str_a, '0') + strmatch(freq_str_a, '1') + strmatch(freq_str_a, '2') + strmatch(freq_str_a, '3') + strmatch(freq_str_a, '4') + strmatch(freq_str_a, '5') + strmatch(freq_str_a, '6') + strmatch(freq_str_a, '7') + strmatch(freq_str_a, '8') + strmatch(freq_str_a, '9') 
IF (test EQ 0) THEN freq_str = strmid(freq_str,1,strlen(freq_str)-1) ELSE CONTINUE
ENDFOR

outfile = temp[0] + '.' + strcompress(string(nside_out,format='(i5)'),/remove_all) + '.fits'
temp = strsplit(outfile, '/', /extract)
ntemp = n_elements(temp)
outfile = outdir + temp[ntemp-1]

print, 'Degrading map ' + infile
print, 'Nside in  = ', nside_str
print, 'Nside out = ', nside_out 
print, 'Freq  = ', freq_str
print, 'Output file = ', outfile

ud_grade, infile, outfile, nside_out=nside_out

; add Nside keyword! (and check ordering)
read_fits_map, outfile, inmap, ordering=ordering,nside=nside
IF (ordering EQ 'RING') THEN ud_grade,outfile,outfile,nside_out=nside_out,order_in='ring', order_out='nest'
read_fits_map, outfile, inmap, ordering=ordering,nside=nside
IF (n_elements(inmap[0,*]) LE 2) THEN write_fits_map, outfile, inmap, /nest ELSE write_tqu, outfile, inmap, /nest



;--------------------------------------------------------------------- 
; 2. Smooth T maps



infile = outfile
temp = strsplit(infile, '.', /extract)
ntemp = n_elements(temp)
outfile = temp[0] + '.' + strcompress(string(nside_out,format='(i5)'),/remove_all) + '.' + strcompress(string(resout,format='(i6)'),/remove_all) + 'arcmin' + '.fits'
fwhm = sqrt(resout^2 - inbeams[i]^2)
lmax = (3*nside_out)-1L
ismoothing, infile, outfile, fwhm=fwhm, /nest, simul=2, lmax=lmax, tmpdir='/attic/cdickins/tmp/'




;--------------------------------------------------------------------- 
; 3. Convert to microK_CMB, add 2uK white noise and output just the
; temperature maps

infile = outfile
temp = strsplit(infile, '/', /extract)
ntemp = n_elements(temp)
outfile = outdir + 'T_ffp2_' + freq_str + 'GHz'  + '.' + strcompress(string(nside_out,format='(i5)'),/remove_all) + '.' + strcompress(string(resout,format='(i6)'),/remove_all) + 'arcmin' + '.' + strcompress(string(noise_level,format='(i9)'), /remove_all) + 'uK_CMB'+ '.fits'

read_fits_map, infile, inmap,ordering=ordering

outmap = inmap[*,0] * 1.0e6  ; convert T to microK RJ
freq = float(freq_str)
outmap = outmap * planckcorr(freq)  ; ***convert to thermodynamic 
noise = randomn(long(freq*i),12L*nside_out^2)* noise_level  ; add white noise
outmap = outmap + noise

write_fits_map, outfile, outmap, /nest, units='uK_CMB'






;--------------------------------------------------------------------- 
; 4. Make rms maps
outfile = outdir + 'rms_' + strcompress(string(noise_level, format='(i9)'), /remove_all) + 'muK_T.fits'
rms = (outmap * 0.0)+ noise_level
write_fits_map, outfile, rms, /nest, units='uK_CMB'



;--------------------------------------------------------------------- 
; 5. Make beam file for each frequency
outfile = outdir + 'beam_3deg_T.fits'
beam = gaussbeam(resout,lmax)
cl2fits, beam, outfile




ENDFOR
;--------------------------------------------------------------------- 
; Finish up
STOP
END
