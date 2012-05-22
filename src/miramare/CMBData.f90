module CMBData
  use HealpixObj
  
  integer, parameter :: nmaxchannel  = 20
  integer, parameter :: nmaxtemplate = 10

  type HealpixMask
     integer*8 npix,nbad
     integer nside,ordering
     integer*8, dimension(:), allocatable :: good
     integer*8, dimension(:), allocatable :: bad 
  end type HealpixMask

  type RegionalModel
     !Healpix related:
     integer nside
     integer*8 npix
     integer*8, dimension(:), allocatable :: pix
     !Getting ready for reduced memory mode.
     integer*8 firstpix_index
     !Component related:
     integer ncomponent
     character(32), dimension(:), allocatable :: component_list
     integer nparam_component   ! For convenience (but right now information is redundant with component list)
                                !                  Could use just firstparam_index and nparamcomponent
     !Region labels:
     integer number
     character(32) label
     !Likelihood
     real*8 chi2
     integer firstparam_index,lastparam_index
  end type RegionalModel

contains
  
  subroutine ReadMask(mask,maskfile,getbadpix,nside_out)
    type(HealpixMask), intent(out) :: mask
    character(*), intent(in)       :: maskfile
    logical, intent(in), optional  :: getbadpix
    integer, intent(in), optional  :: nside_out
    
    logical getthosebadpix
    type(HealpixMap) mask_in, mask_in_degraded

    integer*8 ii,pp,qq,npix_in

    if (present(getbadpix)) then
       getthosebadpix = getbadpix
    else
       getthosebadpix = .false.
    endif
    
    if (Feedback>0) print *,' Reading mask ', trim(maskfile)
    call HealpixMap_read(mask_in,trim(maskfile))

    if(present(nside_out)) then
      call HealpixMap_udgrade(mask_in, mask_in_degraded, nside_out)
      do pp = 0, mask_in_degraded%npix-1
         if (mask_in_degraded%tqu(pp,1) .lt. 1.) mask_in_degraded%tqu(pp,1) = 0.
      enddo
      mask_in = mask_in_degraded
    endif
    
    npix_in = 0
    do ii = 0, nside2npix(mask_in%nside)-1
      if (mask_in%TQU(ii,1) .gt. 0.5) npix_in = npix_in+1
    end do
    if (Feedback>0) print *,' npix_observed = ',npix_in
    mask%nside    = mask_in%nside
    mask%npix     = npix_in
    mask%ordering = mask_in%ordering

    allocate(mask%good(mask%npix))
 
    if(getthosebadpix) then
       mask%nbad = nside2npix(mask_in%nside)-npix_in
       allocate(mask%bad(1:mask%nbad))
    endif
    
    !--------------------------
    !Get observed pixel numbers
    !--------------------------
    pp = 0
    qq = 0
    do ii=0, nside2npix(mask_in%nside)-1
      if (mask_in%TQU(ii,1) .gt. 0.5d0) then
         pp            = pp+1
         mask%good(pp) = ii
      else if (mask_in%TQU(ii,1) .le. 0.5d0 .and. getthosebadpix) then
         qq           = qq+1
         mask%bad(qq) = ii
      endif
    end do
    
    if (Feedback>0) print *,' Finished reading mask '

 end subroutine ReadMask


  subroutine MakeDummyMask(mask,nside_out,getbadpix)
    type(HealpixMask), intent(out) :: mask
    integer, intent(in)            :: nside_out
    logical, intent(in), optional  :: getbadpix

    logical getthosebadpix
    integer*8 ii,pp,npix

    if (present(getbadpix)) then
       getthosebadpix=getbadpix
    else
       getthosebadpix = .false.
    endif
    

    npix=nside2npix(nside_out)

    mask%nside=nside_out
    mask%npix=npix

    allocate(mask%good(npix))
 
    if(getthosebadpix) then
       mask%nbad = 0
    endif
    
    !--------------------------
    !Get observed pixel numbers
    !--------------------------
    pp = 0
    do ii = 0, npix-1
       pp = pp+1
       mask%good(pp)=ii
    end do

 end subroutine MakeDummyMask


 subroutine ReadRegionalMask(Region,Regionmaskfile,nregion,npixtot,mask)
    type(RegionalModel), intent(out), allocatable :: Region(:)
    character(*), intent(in)                      :: Regionmaskfile
    integer, intent(out)                          :: nregion
    integer, intent(out)                          :: npixtot
    type(HealpixMask), intent(in), optional       :: mask

    type(HealpixMap) regionmask

    integer*8 pp,rr,npix,index
    integer*8, dimension(:), allocatable :: temp_counter

    npixtot = 0

    if(regionmaskfile .ne. '') then
       if (Feedback>0) print *,' Reading regional mask ', trim(regionmaskfile)
       call HealpixMap_read(regionmask,trim(regionmaskfile))
    else
       if (Feedback>0) print *,' All-sky single region.'
       call HealpixMap_init(regionmask,nside2npix(mask%nside),1)
       regionmask%tqu(:,1) = 1.
    end if


    !----------------------------------------------
    !Modify the regional mask using the binary mask
    !----------------------------------------------
    if(present(mask)) then
       if (mask%nbad .gt. 0) then
          regionmask%TQU(mask%bad(:),1) = 0
       end if
    end if

    !-----------------------------------------------------------------------
    !Get number of regions, assumed to be a sequence of consecutive integers
    !Regions labelled 0 are ignored.
    !-----------------------------------------------------------------------
    nregion = maxval(regionmask%TQU(:,1))
    allocate(region(nregion))

    !--------------------------------------------------------------------
    !Get number of pixels per region and allocate memory for pixel values
    !--------------------------------------------------------------------
    index = 0
    do rr = 1 ,nregion
       npix = count(regionmask%TQU(:,1) .eq. rr)
       region(rr)%number = rr
       region(rr)%npix   = npix
       region(rr)%nside  = regionmask%nside
       allocate(region(rr)%pix(npix))
    end do

    !----------------------------
    !Get pixel numbers per region
    !----------------------------
    allocate(temp_counter(nregion))
    temp_counter(:) = 0
    do pp = 0, nside2npix(regionmask%nside)-1
       rr = regionmask%TQU(pp,1)
       if (rr .gt. 0) then
          temp_counter(rr)                 = temp_counter(rr) + 1
          Region(rr)%pix(temp_counter(rr)) = pp
          npixtot                          = npixtot + 1
       endif
    end do

    !----------------------------
    !Make checks on pixel numbers
    !----------------------------
    npixtot = 0
    do rr = 1,nregion
      Region(rr)%firstpix_index = npixtot
      npixtot                = npixtot + Region(rr)%npix
       if(temp_counter(rr) .ne. Region(rr)%npix) then
          print *,'Problem reading in regional mask.'
          stop
       endif
    end do
    deallocate(temp_counter)

    call HealpixMap_free(regionmask)
    if (Feedback>0) print *,' Finished reading regional mask.'

 end subroutine ReadRegionalMask

 subroutine HealpixMap_getpatchfromfile(mdata,filename,nmap,maplist,npix, &
      pixlist,offset,factor)
   integer*8 npix
   integer nmap
   real, intent(out)  :: mdata(nmap,npix)
   integer*8 :: pixlist(npix)
   integer :: maplist(nmap)
   character(LEN=512) filename
   real  :: offset,factor

   integer*8 ss,pp
   type(HealpixMap) Hmap

   if (feedback .gt. 0) print *, ' Reading ',trim(filename)
   call HealpixMap_read(Hmap,trim(filename))
   
   if (factor == 1.0 .and. offset == 0.0) then
      do ss = 1, nmap
         do pp=1, npix
            mData(ss,pp) = Hmap%TQU(pixlist(pp),maplist(ss))
         end do
      end do
   else
      do ss = 1, nmap
         do pp = 1, npix
            mData(ss,pp) = Hmap%TQU(pixlist(pp),maplist(ss))*factor+offset  
         end do
      end do
   endif
   call HealpixMap_Free(Hmap)

end subroutine HealpixMap_getpatchfromfile


!  subroutine HealpixMap_coaddpatch(MData,Hmap,nmap,maplist,npix,pixlist,factor)
!    real, dimension(:,:) :: MData
!    type(HealpixMap) Hmap
!    integer nmap
!    integer *8 npix
!    integer*8, dimension(npix) :: pixlist
!    integer, dimension(nmap) :: maplist
!    real  :: factor
   
!    integer*8 ss,pp
   
!    do ss=1, nmap
!       do pp=1, npix
!          MData(ss,pp)=MData(ss,pp)+Hmap%TQU(pixlist(pp),maplist(ss))*factor
!       end do
!    end do
   
!  end subroutine HealpixMap_coaddpatch

end module
