MODULE fits_routines

! Routines for handling of fitsio files
! Partly copied from LevelS

   implicit none
   integer, parameter, private :: sp = kind(1.0)
   integer, parameter, private :: dp = kind(1.0d0)
   integer, parameter, private :: isp = selected_int_kind(9)
   integer, parameter, private :: idp = selected_int_kind(18)

   integer,private :: ID = -1

   character, parameter :: datatype_short  = 'I'
   character, parameter :: datatype_int    = 'J'
   character, parameter :: datatype_long   = 'K'
   character, parameter :: datatype_real   = 'E'
   character, parameter :: datatype_double = 'D'
   character, parameter :: datatype_char   = 'A'

! FITS file access types
   integer, parameter, public :: FITS_READONLY  = 0
   integer, parameter, public :: FITS_READWRITE = 1

! FITS HDU types
   integer, parameter, public :: FITS_IMAGE  = 0
   integer, parameter, public :: FITS_ASCTAB = 1
   integer, parameter, public :: FITS_BINTAB = 2

! constant for yet undefined values
   integer, parameter, public :: FITS_INVALID = -123456

!! data type containing information about the currently connected HDU
   type handle_hdu
      integer                     :: unit    = -1
      integer                     :: nrows   = -1
      integer                     :: ncols   = -1
      integer                     :: hdutype = -1
      type(handle_column),pointer :: columns(:) =>NULL()
   end type

!! data type containing information about a single column in a FITS table
   type handle_column
      character(len=80) :: name = ''
      character(len=80) :: unit = ''
      integer           :: repcount = 1
      character         :: datatype = ''
      integer(idp)      :: length = 0
   end type

   public set_fits_id,            &
          fits_open,              &
          fits_close,             &
          fits_create,            &
          fits_goto_hdu,          &
          fits_initialize_bintab, &
          fits_find_column,       &
          fits_read_bincolumn,    &
          fits_write_bincolumn,   &
          fits_put_key,           &
          fits_get_key,           &
          fits_put_comment

   private check_errors,         &
           write_errors,         &
           check_column_read,    &
           check_column_write,   &
           read_column_dbin,     &
           read_column_sbin,     &
           read_column_ibin,     &
           write_column_dbin,    &
           write_column_sbin,    &
           write_column_ibin,    &
           get_key_int,          &
           get_key_logical,      &
           get_key_string,       &
           get_key_real,         &
           get_key_double,       &
           put_key_int,          &
           put_key_logical,      &
           put_key_string,       &
           put_key_real,         &
           put_key_double

   interface fits_read_bincolumn
      module procedure read_column_dbin,  &
                       read_column_sbin,  &
                       read_column_ibin
   end interface

   interface fits_write_bincolumn
      module procedure write_column_dbin,  &
                       write_column_sbin,  &
                       write_column_ibin
   end interface

   interface fits_initialize_bintab
      module procedure initialize_bintab_simple,     &
                       initialize_bintab_onecolumn,  &
                       initialize_bintab_vec
   end interface

   interface fits_get_key
      module procedure get_key_int,     &
                       get_key_logical, &
                       get_key_string,  &
                       get_key_real,    &
                       get_key_double
   end interface

   interface fits_put_key
      module procedure put_key_int,     &
                       put_key_logical, &
                       put_key_string,  &
                       put_key_real,    &
                       put_key_double
   end interface

CONTAINS

!------------------------------------------------------------------------------


   SUBROUTINE set_fits_id(id_in)
! Set process ID to be printed out with error messages
      integer :: id_in
      ID = id_in

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE fits_open(handle,file,rwmode,hdu)
! Open an existing file

      type(handle_hdu),intent(out)         :: handle
      character(len=*),intent(in)          :: file
      integer,         intent(in),optional :: rwmode, hdu
      integer                              :: rwm, nhdu, stat, block
      logical                              :: there

      inquire(file=trim(file),exist=there)

      if (.not.there) then
         if (ID.ge.0) then
            write(*,'(i3,a)') ID,': ERROR in fits_open: File not found'
            write(*,'(i3,a)') ID, trim(file)
         else
            write(*,*) 'ERROR in fits_open: File not found'
            write(*,*) trim(file)
         endif
         stop
      endif

      rwm = fits_readonly
      if (present(rwmode)) rwm=rwmode

      nhdu = 2
      if (present(hdu)) nhdu=hdu

      stat = 0
      call FTGIOU(handle%unit,stat)
      call FTOPEN(handle%unit,trim(file),rwm,block,stat)
      call check_errors(stat)

      call fits_goto_hdu(handle,nhdu)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE fits_goto_hdu(handle,nhdu)

      type(handle_hdu),intent(inout) :: handle
      integer,         intent(in)    :: nhdu
      integer                        :: stat, i, repcount, tnull
      real(dp)                       :: tscal,tzero
      character(len=80)              :: tdisp, ttype, tunit
      character                      :: datatype

      stat = 0
      call FTMAHD(handle%unit,nhdu,handle%hdutype,stat)
      call check_errors(stat)

      call FTGNCL(handle%unit,handle%ncols,stat)
      call FTGNRW(handle%unit,handle%nrows,stat)
      call check_errors(stat)

      if (associated(handle%columns)) deallocate(handle%columns)

      allocate(handle%columns(handle%ncols))

      stat = 0
      do i = 1,handle%ncols

         call FTGBCL(handle%unit,i,ttype,tunit,datatype,repcount,   &
                     tscal,tzero,tnull,tdisp,stat)

         handle%columns(i)%datatype = uppercase(datatype)
         handle%columns(i)%repcount = repcount
         handle%columns(i)%length   = repcount*handle%nrows
         handle%columns(i)%name     = trim(ttype)
         handle%columns(i)%unit     = trim(tunit)

         call check_errors(stat)
      enddo

  CONTAINS

      FUNCTION uppercase(c) RESULT(cu)
         character :: c, cu
         integer   :: ic

         ic = ichar(c)
         if (ic.ge.ichar('a').and.ic.le.ichar('z')) &
            ic = ic +ichar('A')-ichar('a')
         cu = char(ic)

      END FUNCTION

   END SUBROUTINE


!------------------------------------------------------------------------------


   FUNCTION fits_find_column(handle,colname) RESULT(colnum)
! Find column with given name
! Return -1 if file or column not found.

     integer          :: colnum
     type(handle_hdu) :: handle
     character(len=*) :: colname
     integer          :: i

     colnum = -1
     do i = 1,handle%ncols
        if (handle%columns(i)%name=='colname') then
           colnum = i
           exit
        endif
     enddo

   END FUNCTION


!------------------------------------------------------------------------------


   SUBROUTINE check_column_read(handle,icol,frow,felem,n,offset)

      type(handle_hdu),intent(in)  :: handle
      integer,         intent(out) :: frow, felem
      integer,         intent(in)  :: n, icol
      integer(idp),    intent(in)  :: offset
      integer(idp)                 :: idrep, ilast

      if (icol.gt.handle%ncols) then
         if (ID.ge.0) then
            write(*,'(i3,a)') ID,': ERROR in fits_read_column:' // &
                                 ' Column index too large'
            write(*,'(i3,a,2i5)') ID,': icol, ncols = ',icol,handle%ncols
         else
            write(*,*) 'ERROR in fits_read_column: Column index too large'
            write(*,'(a,2i5)') ': icol, ncols = ',icol,handle%ncols
         endif
         stop
      endif

      idrep = handle%columns(icol)%repcount
      frow  = offset/idrep+1
      felem = mod(offset,idrep)+1

      ilast = offset+n

      if (ilast.gt.handle%columns(icol)%length) then
         if (ID.ge.0) then
            write(*,'(i3,a)') ID,': ERROR in fits_check_column_read:' // &
                                 ' Reading past end of file'
            write(*,'(i3,a,2i9)') ID,': offset, n =',offset,n
            write(*,'(i3,a,i9)') ID,': length    =',handle%columns(icol)%length
         else
            write(*,*) 'ERROR in fits_check_column_read:' // &
                       ' Reading past end of file'
            write(*,'(a,2i9)') 'offset, n =',offset,n
            write(*,'(a,i9)') 'length    =',handle%columns(icol)%length
         endif
         stop
      endif

   END SUBROUTINE



!------------------------------------------------------------------------------


   SUBROUTINE read_column_dbin(handle,icol,dbuffer,n,offset)
! Read n samples from column icol, skipping the first offset samples

      type(handle_hdu),intent(in)          :: handle
      integer,         intent(in)          :: n, icol
      real(dp),        intent(out)         :: dbuffer(n)
      integer(idp),    intent(in),optional :: offset
      integer                              :: stat, frow, felem
      integer(idp)                         :: ioff
      logical                              :: anynul

      if (present(offset)) then
         ioff = offset
      else
         ioff = 0
      endif

      call check_column_read(handle,icol,frow,felem,n,ioff)

      stat = 0
      call FTGCVD(handle%unit,icol,frow,felem,n,0._dp,dbuffer,anynul,stat)

      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE read_column_sbin(handle,icol,sbuffer,n,offset)
! Read n samples from column icol, skipping the first offset samples

      type(handle_hdu),intent(in)          :: handle
      integer,         intent(in)          :: n, icol
      real(sp),        intent(out)         :: sbuffer(n)
      integer(idp),    intent(in),optional :: offset
      integer                              :: stat, frow, felem
      integer(idp)                         :: ioff
      logical                              :: anynul

      if (present(offset)) then
         ioff = offset
      else
         ioff = 0
      endif

      call check_column_read(handle,icol,frow,felem,n,ioff)

      stat = 0
      call FTGCVE(handle%unit,icol,frow,felem,n,0.,sbuffer,anynul,stat)

      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE read_column_ibin(handle,icol,ibuffer,n,offset)
! Read n samples from column icol, skipping the first offset samples

      type(handle_hdu),intent(in)          :: handle
      integer,         intent(in)          :: n, icol
      integer,         intent(out)         :: ibuffer(n)
      integer(idp),    intent(in),optional :: offset
      integer                              :: stat, frow, felem
      integer(idp)                         :: ioff
      logical                              :: anynul

      if (present(offset)) then
         ioff = offset
      else
         ioff = 0
      endif

      call check_column_read(handle,icol,frow,felem,n,ioff)

      stat = 0
      call FTGCVJ(handle%unit,icol,frow,felem,n,0,ibuffer,anynul,stat)

      call check_errors(stat)

   END SUBROUTINE


!-------------------------------------------------------------------------------


   SUBROUTINE fits_close(handle)

     type(handle_hdu),intent(inout) :: handle
     integer                        :: stat

     stat = 0
     call FTCLOS(handle%unit,stat)
     call FTFIOU(handle%unit,stat)
     call check_errors(stat)

     handle%unit = -1
     if (associated(handle%columns)) deallocate(handle%columns)

     handle%columns =>NULL()
     handle%hdutype = -1
     handle%nrows = -1
     handle%ncols = -1

   END SUBROUTINE



!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


   SUBROUTINE fits_create(handle,file)

      type(handle_hdu), intent(out) :: handle
      character(len=*), intent(in)  :: file
      integer :: stat

      stat = 0
      call ftgiou(handle%unit,stat)
      call ftinit(handle%unit,trim(file),1,stat)
      call ftphps(handle%unit,8,0,(/0,0/),stat)
      call ftpdat(handle%unit,stat)

      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE initialize_bintab_simple(handle,datatype)

      type(handle_hdu), intent(inout) :: handle
      character(len=*), intent(in)    :: datatype
      character(len=80),allocatable   :: tform(:), tunit(:), ttype(:)
      integer                         :: stat, i, n

      n = len(datatype)

      handle%ncols = n
      handle%nrows = 0

      allocate(tform(n),tunit(n),ttype(n))
      tunit = ''
      ttype = ''

      if (associated(handle%columns)) deallocate(handle%columns)
      allocate (handle%columns(n))

      handle%columns%name     = ''
      handle%columns%unit     = ''
      handle%columns%repcount = 1

      do i = 1,n
         tform(i) = '1'//datatype(i:i)
         handle%columns(i)%datatype = datatype(i:i)
      enddo

      stat = 0
      call FTIBIN(handle%unit,0,n,ttype,tform,tunit,'xtension',0,stat)
      call check_errors(stat)

      handle%hdutype = FITS_BINTAB

      deallocate(tform,tunit,ttype)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE initialize_bintab_onecolumn(handle,datatype,name,repcount,unit)

      type(handle_hdu), intent(inout)       :: handle
      character,        intent(in)          :: datatype
      character(len=*), intent(in)          :: name
      integer,          intent(in),optional :: repcount
      character(len=*), intent(in),optional :: unit
      character(len=80)                     :: tform(1), tunit(1), ttype(1)
      integer                               :: stat

      handle%ncols = 1
      handle%nrows = 0

      if (associated(handle%columns)) deallocate(handle%columns)
      allocate (handle%columns(1))

      ttype = name
      handle%columns%name = ttype

      if (present(repcount)) then
         write (tform(1),*) repcount
         tform(1) = trim(adjustl(tform(1)))//datatype
         handle%columns%repcount = repcount
      else
         tform(1) = '1' // datatype
         handle%columns%repcount = 1
      endif

      if (present(unit)) then
         tunit = unit
      else
         tunit = ''
      endif
      handle%columns%unit = tunit

      handle%columns%datatype = datatype

      stat = 0
      call FTIBIN(handle%unit,0,1,ttype,tform,tunit,'xtension',0,stat)
      call check_errors(stat)

      handle%hdutype = FITS_BINTAB

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE initialize_bintab_vec(handle,datatype,name,repcount,unit)

      type(handle_hdu), intent(inout)       :: handle
      character(len=*), intent(in)          :: datatype
      character(len=*), intent(in)          :: name(:)
      integer,          intent(in),optional :: repcount(:)
      character(len=*), intent(in),optional :: unit(:)
      character(len=80),allocatable         :: tform(:), tunit(:), ttype(:)
      integer                               :: stat, i, n

      n = len(datatype)

      handle%ncols = n
      handle%nrows = 0

      allocate(tform(n),tunit(n),ttype(n))

      if (associated(handle%columns)) deallocate(handle%columns)
      allocate (handle%columns(n))

      ttype = name
      handle%columns%name = ttype

      if (present(repcount)) then
         do i = 1,n
            write (tform(i),*) repcount(i)
            tform(i) = trim(adjustl(tform(i)))//datatype(i:i)
         enddo
         handle%columns%repcount = repcount
      else
         do i = 1,n
            tform(i) = '1' // datatype(i:i)
         enddo
         handle%columns%repcount = 1
      endif

      if (present(unit)) then
         tunit = unit
      else
         tunit = ''
      endif
      handle%columns%unit = tunit

      do i = 1,n
         handle%columns(i)%datatype = datatype(i:i)
      enddo

      stat = 0
      call FTIBIN(handle%unit,0,n,ttype,tform,tunit,'xtension',0,stat)
      call check_errors(stat)

      handle%hdutype = FITS_BINTAB

      deallocate(tform,tunit,ttype)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE check_column_write(handle,icol,frow,felem,offset)

      type(handle_hdu),intent(in)  :: handle
      integer,         intent(out) :: frow, felem
      integer,         intent(in)  :: icol
      integer(idp),    intent(in)  :: offset
      integer(idp)                 :: idrep

      if (icol.gt.handle%ncols) then
         if (ID.ge.0) then
            write(*,'(i3,a)') ID,': ERROR in fits_write_column:' // &
                                 ' Column index too large'
            write(*,'(i3,a,2i5)') ID,': icol, ncols = ',icol,handle%ncols
         else
            write(*,*) 'ERROR in fits_write_column: Column index too large'
            write(*,'(a,2i5)') ': icol, ncols = ',icol,handle%ncols
         endif
         stop
      endif

      idrep = handle%columns(icol)%repcount
      frow  = offset/idrep+1
      felem = mod(offset,idrep)+1

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE write_column_dbin(handle,colnum,dbuffer,n,offset)

      type(handle_hdu),intent(inout)   :: handle
      integer,         intent(in)      :: colnum, n
      real(dp),        intent(in)      :: dbuffer(n)
      integer(idp),intent(in),optional :: offset
      integer                          :: stat, frow, felem, i

      if (present(offset)) then
         call check_column_write(handle,colnum,frow,felem,offset)
      else
         frow = 1
         felem = 1
      endif

      stat = 0
      call FTPCLD(handle%unit,colnum,frow,felem,n,dbuffer,stat)
      call FTGNRW(handle%unit,handle%nrows,stat)

      call check_errors(stat)

      do i = 1,handle%ncols
         handle%columns(i)%length = handle%columns(i)%repcount*handle%nrows
      enddo

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE write_column_sbin(handle,colnum,sbuffer,n,offset)

      type(handle_hdu),intent(inout)   :: handle
      integer,         intent(in)      :: colnum, n
      real(sp),        intent(in)      :: sbuffer(n)
      integer(idp),intent(in),optional :: offset
      integer                          :: stat, frow, felem, i

      if (present(offset)) then
         call check_column_write(handle,colnum,frow,felem,offset)
      else
         frow = 1
         felem = 1
      endif

      stat = 0
      call FTPCLE(handle%unit,colnum,frow,felem,n,sbuffer,stat)
      call FTGNRW(handle%unit,handle%nrows,stat)

      call check_errors(stat)

      do i = 1,handle%ncols
         handle%columns(i)%length = handle%columns(i)%repcount*handle%nrows
      enddo

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE write_column_ibin(handle,colnum,ibuffer,n,offset)

      type(handle_hdu),intent(inout)   :: handle
      integer,         intent(in)      :: colnum, n
      integer,         intent(in)      :: ibuffer(n)
      integer(idp),intent(in),optional :: offset
      integer                          :: stat, frow, felem, i

      if (present(offset)) then
         call check_column_write(handle,colnum,frow,felem,offset)
      else
         frow = 1
         felem = 1
      endif

      stat = 0
      call FTPCLJ(handle%unit,colnum,frow,felem,n,ibuffer(1:n),stat)
      call FTGNRW(handle%unit,handle%nrows,stat)

      call check_errors(stat)

      do i = 1,handle%ncols
         handle%columns(i)%length = handle%columns(i)%repcount*handle%nrows
      enddo

   END SUBROUTINE


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


   SUBROUTINE get_key_int(handle,keyname,value)

      integer           :: value
      type(handle_hdu)  :: handle
      character(len=*)  :: keyname
      character(len=80) :: comment
      integer           :: stat

      stat = 0
      call FTGKYJ(handle%unit,trim(keyname),value,comment,stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE get_key_logical(handle,keyname,value)

      logical           :: value
      type(handle_hdu)  :: handle
      character(len=*)  :: keyname
      character(len=80) :: comment
      integer           :: stat

      stat = 0
      call FTGKYL(handle%unit,trim(keyname),value,comment,stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE get_key_string(handle,keyname,value)

      character(len=*),intent(out) :: value
      type(handle_hdu)  :: handle
      character(len=*)  :: keyname
      character(len=80) :: comment
      integer           :: stat

      stat = 0
      call FTGKYS(handle%unit,trim(keyname),value,comment,stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE get_key_real(handle,keyname,value)

      real              :: value
      type(handle_hdu)  :: handle
      character(len=*)  :: keyname
      character(len=80) :: comment
      integer           :: stat

      stat = 0
      call FTGKYE(handle%unit,trim(keyname),value,comment,stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE get_key_double(handle,keyname,value)

      double precision  :: value
      type(handle_hdu)  :: handle
      character(len=*)  :: keyname
      character(len=80) :: comment
      integer           :: stat

      stat = 0
      call FTGKYD(handle%unit,trim(keyname),value,comment,stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


   SUBROUTINE put_key_int(handle,keyname,value,comment)

      type(handle_hdu) :: handle
      integer          :: value
      character(len=*) :: keyname,comment
      integer          :: stat

      if (len_trim(keyname).gt.8) keyname=keyname(1:8)
      if (len_trim(comment).gt.47) comment=comment(1:47)

      stat = 0
      call FTUKYJ(handle%unit,trim(keyname),value,trim(comment),stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE put_key_logical(handle,keyname,value,comment)

      type(handle_hdu) :: handle
      logical          :: value
      character(len=*) :: keyname,comment
      integer          :: stat

      if (len_trim(keyname).gt.8)  keyname=keyname(1:8)
      if (len_trim(comment).gt.47) comment=comment(1:47)

      stat = 0
      call FTUKYL(handle%unit,trim(keyname),value,trim(comment),stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE put_key_string(handle,keyname,value,comment)

      type(handle_hdu) :: handle
      character(len=*) :: value
      character(len=*) :: keyname,comment
      integer          :: stat

      if (len_trim(keyname).gt.8)  keyname=keyname(1:8)
      if (len_trim(comment).gt.47) comment=comment(1:47)
      if (len_trim(value).gt.20)   value=value(1:20)

      stat = 0
      call FTUKYS(handle%unit,trim(keyname),value,trim(comment),stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE put_key_real(handle,keyname,value,comment)

      type(handle_hdu) :: handle
      real             :: value
      character(len=*) :: keyname,comment
      integer          :: stat

      if (len_trim(keyname).gt.8)  keyname=keyname(1:8)
      if (len_trim(comment).gt.47) comment=comment(1:47)

      stat = 0
      call FTUKYE(handle%unit,trim(keyname),value,-15,trim(comment),stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE put_key_double(handle,keyname,value,comment)

      type(handle_hdu) :: handle
      double precision :: value
      character(len=*) :: keyname,comment
      integer          :: stat

      if (len_trim(keyname).gt.8) keyname=keyname(1:8)
      if (len_trim(comment).gt.47) comment=comment(1:47)

      stat = 0
      call FTUKYD(handle%unit,trim(keyname),value,-20,trim(comment),stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE fits_put_comment(handle,comment)

      type(handle_hdu) :: handle
      character(len=*) :: comment
      integer          :: stat

      if (len_trim(comment).gt.47) comment=comment(1:47)

      stat = 0
      call FTPCOM(handle%unit,trim(comment),stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE check_errors(stat)
! Print error messages and stop

      integer,intent(inout) :: stat
!      character(len=30)     :: errtext
!      character(len=80)     :: errmessage

      if (stat<=0) return

      call write_errors(stat)
      stop

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE write_errors(stat)
! Print error messages and proceed

      integer,intent(inout) :: stat
      character(len=30)     :: errtext
      character(len=80)     :: errmessage

      if (stat<=0) return

      if (ID.ge.0) then
         call ftgerr(stat,errtext)
         write(*,'(i3,": ",a,i5)') ID,'FITSIO Error Status =',stat
         write(*,'(i3,": ",a)')    ID, errtext

         call ftgmsg(errmessage)
         do while (errmessage /= ' ')
            write(*,'(i3,": ",a)') ID,errmessage
            call ftgmsg(errmessage)
         end do

      else
         call ftgerr(stat,errtext)
         write(*,*) 'FITSIO Error Status =',stat
         write(*,*) errtext

         call ftgmsg(errmessage)
         do while (errmessage /= ' ')
            write(*,*) errmessage
            call ftgmsg(errmessage)
         end do

      endif

   END SUBROUTINE


!------------------------------------------------------------------------------

END MODULE
