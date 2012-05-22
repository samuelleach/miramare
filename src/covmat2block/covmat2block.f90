! This code reads a pixel-pixel covariance matrix and saves it on disc 
! in expanded (block) form, where correlation of different maps
! is presented in separate blocks

PROGRAM covmat2block
  USE healpix_types
  USE extension, ONLY   : getArgument, nArguments
  USE covmat_util, ONLY : tic, toc

  IMPLICIT NONE

  CHARACTER(len=80)      :: argument, covmat_file, outfilename
  INTEGER(dp)            :: npix, nstokes, nside, imap, istokes, isignal
  INTEGER, PARAMETER     :: covmat_unit=55
  REAL(dp), POINTER      :: buffer(:)
  LOGICAL                :: there
  INTEGER                :: ierr, iarg, nelements, icol

  IF (nArguments() < 3 .OR. nArguments() > 4)  THEN
     STOP 'Usage: covmat2block <covmat> <nside> <nstokes> [<covmat_out>]'
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     READ (argument, *) nside
     WRITE (*,'(a,i4)') 'Nside   == ', nside

     CALL getArgument(3, argument)
     READ (argument, *) nstokes
     WRITE (*,'(a,i4)') 'Nstokes == ', nstokes

     outfilename = 'block_'//TRIM(ADJUSTL(covmat_file))
     IF (nArguments() > 3) THEN
        CALL getArgument(4, argument)
        outfilename = TRIM(ADJUSTL(argument))
     END IF
     WRITE (*,*) ' Writing to ' // TRIM(outfilename)
  END IF

  npix = 12*nside**2
  ALLOCATE(buffer(0:npix*nstokes-1), stat=ierr)
  IF (ierr/=0) STOP 'no room for buffer'

  ! Read the matrix full column at a time
  OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old', &
       form='unformatted', access='direct', recl=npix*nstokes*8)

  ! Write the converted matrix single component map at a time
  OPEN(unit=covmat_unit+1, file=TRIM(ADJUSTL(outfilename)), &
       status='replace', form='unformatted', access='direct', recl=npix*8)

  !write (*,*) '  WARNING converting from CMB uK to antenna K'

  CALL tic
  imap = 1
  DO istokes = 1, nstokes
     DO icol = istokes, nstokes*npix, nstokes
        READ (covmat_unit, rec=icol) buffer
        !buffer = buffer/(1E6*1.1331762)**2 ! from CMB uK to antenna K
        DO isignal = 0, nstokes-1
           WRITE (covmat_unit+1, rec=imap) buffer(isignal:npix*nstokes-1:nstokes)
           imap = imap + 1
        END DO
     END DO
  END DO
  CALL toc('convert to block')

  CLOSE(unit=covmat_unit)
  CLOSE(unit=covmat_unit+1)

END PROGRAM covmat2block