PRO pp_smooth

; procedure to do the pre-processing for FFP 2 maps for running with
; Commander; This one is for intensity and polarization with smoothing
; A separate routine (calcnoisematrix.pro) calculates the resulting
; noise matrix based on MCs.
;
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
; 14-Jan-2009  C. Dickinson.  1st go
;--------------------------------------------------------------------- 
; set data directories, inputs etc.
indir  = '/attic/cdickins/WG2/FFP2/data/'                     ; input directory
outdir = '/attic/cdickins/WG2/FFP2/proc3/'                     ; output directory
mapfiles = file_search(indir, 'TQU_ffp2_*.fits',count=nmaps)  ; list of files
nside_out = 64L                                               ; output Nside
inbeams = [52.8,33,39.6,30.6,24,21,14,13.2,10,7.1,5,5,5,5]    ; beam fwhm on input (arcmin)
;inbeams = [5,5]                                               ; for 545/857GHz!!!!
resout = 180.                                                 ; output resolution (for T map)
noise_level = 2.                                              ; white noise level to add (uK_CMB) to I maps

;--------------------------------------------------------------------- 
; 1. Smooth maps by required amount

FOR i=0L, nmaps-1 DO BEGIN
infile = mapfiles[i]
temp = strsplit(infile, '.', /extract)
ntemp = n_elements(temp)
nside_str = temp[ntemp-2]
freq_str = strmid(infile,strpos(infile,'GHz')-3,3) 
outfile = outdir + 'TQU_ffp2_' + freq_str + 'GHz'  + '.' + strcompress(string(nside_str,format='(i5)'),/remove_all) + '.' + strcompress(string(resout,format='(i6)'),/remove_all) + 'arcmin' + '.fits'

; add Nside keyword!
read_fits_map, infile, inmap    ; inmap here is the TQU map
write_tqu, outfile, inmap, /nest

fwhm = sqrt(resout^2 - inbeams[i]^2)
lmax = (3*nside_out)-1L
ismoothing, outfile, outfile, fwhm=fwhm, /nest, simul=2, lmax=lmax, tmpdir='/attic/cdickins/tmp/'





;--------------------------------------------------------------------- 
; 2. Convert to microK_CMB and add white noise to I maps only
 
infile = outfile
temp = strsplit(infile, '/', /extract)
ntemp = n_elements(temp)
outfile = outdir + 'TQU_ffp2_' + freq_str + 'GHz'  + '.' + strcompress(string(nside_str,format='(i5)'),/remove_all) + '.' + 'pp' + '.' + strcompress(string(resout,format='(i6)'),/remove_all) + 'arcmin' + '_uK_CMB'+ '.fits'

outmap = inmap * 1.0e6  ; convert T to microK RJ
freq = float(freq_str)
outmap = outmap * planckcorr(freq)  ; convert to thermodynamic 



; add white noise to the T map only
outmap[*,0] = outmap[*,0] + randomn(-i*5L,12L*nside_out^2)*noise_level


write_tqu, outfile, outmap, /nest, units='uK_CMB'


;--------------------------------------------------------------------- 
; 3. Make beam file for each frequency
outfile = outdir + 'beam_' + strcompress(string(inbeams[i], format='(f9.1)'),/remove_all) + 'arcmin.fits'
beam_T = gaussbeam(resout,lmax)
beam_Q = beam_T
beam_U = beam_T
beam_V = beam_Q * 0.0
beam = [[beam_T],[beam_Q], [beam_U],[beam_V]]
cl2fits, beam, outfile




ENDFOR
;--------------------------------------------------------------------- 
; Finish up
STOP
END
