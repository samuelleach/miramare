PRO write_cvbl, fwhm

;  Finds beam window function amplitude CVBL to convolve the variance
;  maps. The variance beam is the square of the temperature beam, so
;  first we transform the temperature beam window function CONVBL to
;  get the beam profile, then square it and transform back to get CVBL.

; Input:
;        fwhm: Required fwhm beamwidth (degrees)

code   = ['Ka1']
;code   = ['K1','Ka1','Q1','Q2','V1','V2','W1','W2','W3','W4']

FOR i =0,n_elements(code)-1 DO BEGIN
; Open beam window function file

    fname = '../data/wmap/wmap_'+code[i]+'_ampl_bl_5yr_v3.txt'

    PRINT, 'Opening file',fname

    OPENR, 1, fname

; Read through header
    aline = ''
    REPEAT READF, 1, aline UNTIL STRMID(aline,0,1) NE '#'

; Translate first line
    READS, aline, multipole, amplitude
    ll = fltarr(2001)
    bl = fltarr(2001)
    ll[0] = multipole
    bl[0] = amplitude

; read rest of file
    count = 0
    WHILE ~ EOF(1) DO BEGIN
        count += 1
        IF count gt 2000 THEN BEGIN
            PRINT, 'Error: too many lines in Bl file'
            GOTO, QUIT
        ENDIF
        READF, 1, multipole, amplitude
        ll[count] = multipole
        bl[count] = amplitude
    ENDWHILE

    CLOSE, 1

; Trim arrays
    nbl = count
    ll = ll[0:nbl]
    bl = bl[0:nbl]

; Generate normalised radial profile of convolving beam:
    convbl = bl[0] * gaussbeam(fwhm*60.,nbl) / bl

; Fudge to kill divide-by zeroes from inadequate-precision bl values:
    mask = WHERE(bl EQ 0.0)
    IF mask[0] NE -1 THEN convbl[mask] = 0.0

; Make Lagrange-polynomial transform of amp to get beam window
; function

; Choose scale to match size of beam. First find half-power point
    junk = WHERE(convbl GT 0.5, lcount)
    lhalf = junk[lcount-1]
    print, 'lhalf ',lhalf
	   
; Calculate beam out to about 20 * half-power radius, roughly (note that
; this gives 50 points out to the half-power point).

    rad = FINDGEN(1001)*10.0*!pi/(lhalf*1000.)

    x = COS(rad)
    sinrad = SIN(rad)
    lgndr = FLTARR(1001,nbl+1)

    FOR l= 0,nbl DO lgndr[*,l] = LEGENDRE(x,l)
    
; Generate radial profile of convolving beam:
    conva = FLTARR(1001)
    FOR j=0,1000 DO conva[j] = TOTAL((ll+0.5)*convbl*lgndr[j,*],/DOUBLE)

    conva = conva / (2.0*!pi)

    PRINT, 'Peak of convolving beam is ', conva[0], MAX(conva)

; Define variance beam amplitude array of same size as convbl
    cvbl = convbl

; Square convolving beam and convert back to window function
    mult = sinrad*conva^2
    FOR l = 0,nbl DO cvbl[l] = INT_TABULATED(rad,mult*lgndr[*,l])

; Put in 2pi normalization factor:
    cvbl = 2.0*!pi*cvbl

    PRINT, 'CVBL[0] =', cvbl[0]

; Write out beam window function
;    filout = code[i]+STRING(fwhm,FORMAT="('_c',F3.1,'_var_window.txt')")
    filout = code[i]+'_c'+strtrim(fwhm,2)+'_var_window.txt'
    PRINT, 'Writing ', filout
    OPENW, 2, filout
    PRINTF, 2, '# Beam window function for convolving variance maps'
;    PRINTF, 2, '# Output beam: ',STRING(fwhm*60.,FORMAT="(F5.1)"),' arcmin'
    PRINTF, 2, '# Output beam: ',STRtrim(fwhm*60.,2)+' arcmin'
    PRINTF, 2, '# Created by get_cvbl.pro at', SYSTIME()
    FOR l=0,nbl DO PRINTF, 2, l, cvbl[l]

    CLOSE, 2
    
ENDFOR


QUIT:

END
