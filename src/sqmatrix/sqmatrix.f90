PROGRAM sqmatrix

   use fits_routines
   use matrix

   implicit none
   integer, parameter :: sp = kind(1.0)
   integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: isp = selected_int_kind(9)
   integer, parameter :: idp = selected_int_kind(18)

   integer      :: iargc, nopar, i, k, m, icol, istep, nobad, nozero, degrade
   integer      :: nbuff, nosteps, ncc, nside_in, nside_out, nrot, rp
   integer(idp) :: offset
   real(dp)     :: a(3,3), d(3), v(3,3)

   character(len=200) :: file_in  = '', file_out = ''
   character(len=20)  :: cnside, ordering

   type(handle_hdu)   :: inhandle, outhandle

   real(dp),allocatable :: cc_in(:,:), cc_out(:,:)
   real(dp),allocatable :: cc(:)
   real(dp),allocatable :: dpbuffer(:)
   real(sp),allocatable :: spbuffer(:)

   nside_out = -1
   nopar = iargc()

   if (nopar==2) then
      call getarg(1,file_in)
      call getarg(2,file_out)
   elseif (nopar==3) then
      call getarg(1,file_in)
      call getarg(2,file_out)
      call getarg(3,cnside)

      if (cnside(1:6).ne.'nside=') then
         write(*,*) 'Error in input parameters'
         write(*,*) 'Usage: sqmatrix file_in file_out [nside=n]'
         stop
      else
         read(cnside(7:20),*) nside_out
      endif
   else
      write(*,*) 'Error in input parameters'
      write(*,*) 'Usage: sqmatrix file_in file_out [nside=n]'
      stop
   endif


! Open the input matrix file
   nside_in = 0
   ordering = ''
   call fits_open(inhandle,file_in)
   call fits_get_key(inhandle,'NSIDE',nside_in)
   call fits_get_key(inhandle,'ORDERING',ordering)

   if (nside_in==0) then
      write(*,*) 'ERROR: keyword NSIDE not found in header.'
      stop
   elseif (ordering=='') then
      write(*,*) 'ERROR: keyword ORDERING not found in header.'
      stop
   elseif (ordering.ne.'NESTED'.and.ordering.ne.'nested') then
      write(*,*) 'ERROR: expect the matrix to be in NESTED pixeling scheme.'
      stop
   endif

   ncc = inhandle%ncols
   if (ncc.ne.1.and.ncc.ne.6) then
      write(*,*) 'ERROR: unexpected number of columns in input file'
      stop
   endif

   if (nside_out.lt.0) then
       nside_out = nside_in
       degrade = 1
       write(*,*) 'nside =',nside_in
   elseif (nside_out.gt.nside_in) then
      nside_out = nside_in
      degrade = 1
      write(*,*) 'Warning: nside_out>nside_in.'
      write(*,*) 'Setting nside_out =',nside_in
   elseif (nside_out.le.nside_in) then
      degrade = nside_in**2/nside_out**2
      write(*,*) 'nside_in  =',nside_in
      write(*,*) 'nside_out =',nside_out
      if (nside_out**2*degrade.ne.nside_in**2) then
         write(*,*) 'Error: nside must be a power of two.'
         stop
      endif
   endif


! Create and initialize the output file
   call fits_create(outhandle,file_out)

   rp = min(1024,nside_out*nside_out)
   if (ncc==6) then
      call fits_initialize_bintab(outhandle,'eeeeee',  &
                   (/'sc(1,1)','sc(1,2)','sc(1,3)',    &
                     'sc(2,2)','sc(2,3)','sc(3,3)'/),  &
                   (/ rp,rp,rp,rp,rp,rp /) )
   else
      call fits_initialize_bintab(outhandle,'e','sc(1,1)',rp)
   endif

   call fits_put_key(outhandle,'CREATOR','sqmatrix',   &
                               'Software creating FITS file')
   call fits_put_comment(outhandle,'----------------------------------')
   call fits_put_comment(outhandle,'   Healpix Map Specific Keywords  ')
   call fits_put_comment(outhandle,'----------------------------------')
   call fits_put_key(outhandle,'PIXTYPE', 'HEALPIX','HEALPIX Pixelisation')
   call fits_put_key(outhandle,'ORDERING','NESTED', 'Pixel ordering scheme')
   call fits_put_key(outhandle,'NSIDE',nside_out,     &
                               'Resolution parameter for HEALPIX')
   call fits_put_key(outhandle,'FIRSTPIX',0,'First pixel # (0 based)')
   call fits_put_key(outhandle,'LASTPIX',12*nside_out*nside_out-1,      &
                               'Last pixel # (0 based)')

   nbuff = nside_out*nside_out
   nosteps = 12*nside_out*nside_out/nbuff

   allocate(cc_in(ncc,nbuff))
   allocate(cc_out(ncc,nbuff))
   allocate(cc(ncc))
   allocate(dpbuffer(nbuff*degrade))
   allocate(spbuffer(nbuff))


! Start computation loop
   offset = 0
   nobad  = 0
   nozero = 0

   do istep = 1,nosteps

      write(*,'(" i =",i4," /",i4)') istep,nosteps

      cc_in = 0.0
! input the input
      write(*,*) 'Reading data'
      do icol = 1,ncc

         call fits_read_bincolumn(inhandle,icol,dpbuffer,nbuff*degrade,  &
                                  offset*degrade)

! Degrade the resolution if needed
         m = 0
         do i = 1,nbuff
            do k = 1,degrade
               m = m+1
               cc_in(icol,i) = cc_in(icol,i) +dpbuffer(m)

            enddo

         enddo

!
      enddo

      write(*,*) 'Eigenvalue decomposition'
! Compute the square root through eigenvalue decomposition
      do i = 1,nbuff

         if (ncc==6) then
            cc = cc_in(:,i)

            a(1,1) = cc(1)
            a(2,1) = cc(2)
            a(3,1) = cc(3)
            a(1,2) = cc(2)
            a(2,2) = cc(4)
            a(3,2) = cc(5)
            a(1,3) = cc(3)
            a(2,3) = cc(5)
            a(3,3) = cc(6)

! Eigenvalue decomposition
            call Jacobi(a,3,d,v,nrot)

! Construct the square root
            if (any(d.lt.0)) then
               cc = -1.e30
               nobad = 0.0
            elseif (any(d.lt.1.e-40)) then
               cc = 1.e30
               nozero = nozero+1
            else
               d = 1.d0/sqrt(d)
               cc(1) = d(1)*v(1,1)*v(1,1)+d(2)*v(1,2)*v(1,2)+d(3)*v(1,3)*v(1,3)
               cc(2) = d(1)*v(1,1)*v(2,1)+d(2)*v(1,2)*v(2,2)+d(3)*v(1,3)*v(2,3)
               cc(3) = d(1)*v(1,1)*v(3,1)+d(2)*v(1,2)*v(3,2)+d(3)*v(1,3)*v(3,3)
               cc(4) = d(1)*v(2,1)*v(2,1)+d(2)*v(2,2)*v(2,2)+d(3)*v(2,3)*v(2,3)
               cc(5) = d(1)*v(2,1)*v(3,1)+d(2)*v(2,2)*v(3,2)+d(3)*v(2,3)*v(3,3)
               cc(6) = d(1)*v(3,1)*v(3,1)+d(2)*v(3,2)*v(3,2)+d(3)*v(3,3)*v(3,3)
            endif

            cc_out(:,i) = cc
         elseif (ncc==1) then

            cc = cc_in(1,i)

            if (cc(1).lt.0) then
               cc(1) = -1.e30
               nobad = nobad+1
            elseif (cc(1).lt.1.e-40) then
               cc(1) = 1.e30
               nozero = nozero+1
            else
               cc(1) = 1.d0/sqrt(cc(1))
            endif
            cc_out(:,i) = cc

         endif
      enddo

! output the output
      write(*,*) 'Writing output'
      do icol = 1,ncc
         spbuffer = cc_out(icol,:)
         call fits_write_bincolumn(outhandle,icol,spbuffer,nbuff,offset)
      enddo

      offset = offset+nbuff

   enddo

   call fits_close(inhandle)
   call fits_close(outhandle)

   write(*,*) 'Output written in file ',trim(file_out)
   write(*,*) nozero,'singular matrices'
   if (nobad.gt.0) then
      write(*,*) 'Warning!'
      write(*,*) nobad, 'pixels with negative eigenvalues'
   endif

   deallocate(cc_in,cc_out,cc)
   deallocate(dpbuffer,spbuffer)


END PROGRAM
