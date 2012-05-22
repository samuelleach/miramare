PRO pp_pol

; procedure to do the pre-processing for FFP 2 maps for running with
; Commander; This one is for intensity and polarization
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
; -rms files should be made using make_rms.pro (requires matrices)
;
; 21-Oct-2008  C. Dickinson.  1st go
;--------------------------------------------------------------------- 
; set data directories, inputs etc.
indir  = '/attic/cdickins/WG2/FFP2/data/'                     ; input directory
outdir = '/attic/cdickins/WG2/FFP2/proc2/'                     ; output directory
mapfiles = file_search(indir, 'TQU_ffp2_*.fits',count=nmaps)  ; list of files
nside_out = 32L                                               ; output Nside
inbeams = [52.8,33,39.6,30.6,24,21,14,13.2,10,7.1,5,5,5,5]    ; beam fwhm on input (arcmin)
;inbeams = [5,5]                                               ; for 545/857GHz!!!!
resout = 180.                                                 ; output resolution (for T map)
noise_level = 2.                                              ; white noise level to add (uK_CMB)


; 1. Degrade maps
FOR i=0L, nmaps-1 DO BEGIN
infile = mapfiles[i]
temp = strsplit(infile, '.', /extract)
ntemp = n_elements(temp)
nside_str = temp[ntemp-2]
freq_str = strmid(infile,strpos(infile,'GHz')-3,3) 
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

; add Nside keyword!
read_fits_map, outfile, inmap_TQU    ; inmap here is the TQU map
write_tqu, outfile, inmap_TQU, /nest



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
; 3. Convert to microK_CMB, add 2uK white noise and output TQU maps
; with just T smoothed
 
infile = outfile
temp = strsplit(infile, '/', /extract)
ntemp = n_elements(temp)
outfile = outdir + 'TQU_ffp2_' + freq_str + 'GHz'  + '.' + strcompress(string(nside_out,format='(i5)'),/remove_all) + '.' + 'pp' + '.' + strcompress(string(resout,format='(i6)'),/remove_all) + 'arcmin' + '.' + strcompress(string(noise_level,format='(i9)'), /remove_all) + 'uK_CMB'+ '.fits'

read_fits_map, infile, inmap,ordering=ordering      ; this is a smoothed version of TQU

inmap[*,1:2] = inmap_TQU[*,1:2]       ; replace Q and U with the original unsmoothed maps


outmap = inmap * 1.0e6  ; convert T to microK RJ
freq = float(freq_str)
outmap = outmap * planckcorr(freq)  ; convert to thermodynamic 
noise = randomn(long(freq*i),12L*nside_out^2)* noise_level  ; make white noise map
outmap[*,0] = outmap[*,0] + noise[*]          ;add it only to the T map

write_tqu, outfile, outmap, /nest, units='uK_CMB'






;--------------------------------------------------------------------- 
; 4. Make rms maps - do this using make_rms.pro (it needs the matrices)
;outfile = outdir + 'rms_' + strcompress(string(noise_level, format='(i9)'), /remove_all) + 'muK_T.fits'
;rms = (outmap * 0.0)+ noise_level
;write_fits_map, outfile, rms, /nest, units='uK_CMB'



;--------------------------------------------------------------------- 
; 5. Make beam file for each frequency
outfile = outdir + 'beam_' + strcompress(string(inbeams[i], format='(f9.1)'),/remove_all) + 'arcmin.fits'
beam_T = gaussbeam(resout,lmax)
beam_Q = gaussbeam(inbeams[i],lmax)
beam_U = beam_Q
beam_V = beam_Q * 0.0
beam = [[beam_T],[beam_Q], [beam_U],[beam_V]]
cl2fits, beam, outfile




ENDFOR
;--------------------------------------------------------------------- 
; Finish up
STOP
END
