program hitmap_to_regionfile

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      use pix_tools
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad

      type(HealpixMap)   :: mask_highres
      type(HealpixMap)   :: mask_lowres
      type(HealpixMap)   :: hitmap_highres
      type(HealpixMap)   :: hitmap_nextres
      type(HealpixMap)   :: hitmap_lowres
      type(HealpixMap)   :: regions_highres
      type(HealpixMap)   :: regions_nextres
      type(HealpixMap)   :: regions_lowres

      real undef      
      integer pp,N_iter,ii,rr,nregion_upgrade,counter
      integer*8, allocatable :: pixlist(:)
      integer*8, allocatable :: region_upgrade(:)


      !Parameters
      logical verbose
      character(LEN=100) :: infile
      character(LEN=100) :: outfile
      real nhits_max

      integer nside_lowres,nside_highres,n_region

      character(len=*), parameter :: CODE = "HITMAP_TO_REGIONFILE"
      
      
      undef= -1.6375e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Divide a region defined by a hitmap into sub-regions'
         print *,CODE,': with a user defined maximum number of hits.'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose = T or F (default = F).'
         print *,'infile = Input Healpix hitmap filename.'
         print *,'outfile = Output Healpix regionfile filename.'
         print *,'nhits_max = Minimum number of hits on each region.'

         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose = Ini_Read_logical('verbose',.false.)
      infile = Ini_Read_String('infile')
      outfile = Ini_Read_String('outfile')
      nhits_max = Ini_Read_real('nhits_max',10000.)

      if(verbose) print *,trim(CODE),': running'


      if(verbose) print *,trim(CODE),': reading ',trim(infile)
      call HealpixMap_read(hitmap_highres,trim(infile))
      where( hitmap_highres%TQU(:,1) .lt. undef)  hitmap_highres%TQU(:,1)=0.


      nside_highres = hitmap_highres%nside
      N_iter        = 5
      nside_lowres  = nside_highres/2**(N_iter)

      !Allocate high res mask to define original high res hitmap region.
      call HealpixMap_init(mask_highres,nside2npix(nside_highres),1)
      where( hitmap_highres%TQU(:,1) .gt. 0.)  mask_highres%TQU(:,1)=1.
      call HealpixMap_udgrade(mask_highres, mask_lowres, nside_lowres)
      where( mask_lowres%TQU(:,1) .lt. 1.)  mask_lowres%TQU(:,1)=0.

      !Allocate pixel list
      allocate(pixlist(nside2npix(nside_highres)))

      !Accumulate hitmap to nside_lowres = nside_highres/(2^N_iter)
      call HealpixMap_udgrade(hitmap_highres, hitmap_lowres, nside_lowres)
      hitmap_lowres%TQU(:,1)=hitmap_lowres%TQU(:,1)*&
           (float(nside_highres)/float(nside_lowres))**2
      where( mask_lowres%TQU(:,1) .eq. 0.)  hitmap_lowres%TQU(:,1)=0.



      !Allocate region file at low res and number pixels 1 to N_region.
      call HealpixMap_init(regions_lowres,nside2npix(nside_lowres),1)
      where( hitmap_lowres%TQU(:,1) .gt. 0.)  regions_lowres%TQU(:,1)=1.
      n_region=0
      do pp=0,regions_lowres%npix-1
         if(hitmap_lowres%TQU(pp,1) .gt. 0) then
            n_region                   = n_region+1
            regions_lowres%TQU(pp,1)   = n_region
            pixlist(n_region)          = pp
         endif
      end do
      

      !Main iteration loop
      do ii = 1, N_iter
         if(verbose) print *,trim(CODE),': Iteration ',ii
         allocate(region_upgrade(n_region))
         nregion_upgrade=0

         do rr = 1, N_region
            !if N_hits(rr) .gt. nhits_max then keep a note of that region number
            if(hitmap_lowres%tqu(pixlist(rr),1) .gt. nhits_max) then
               nregion_upgrade                 = nregion_upgrade+1
               region_upgrade(nregion_upgrade) = rr               
            endif
         end do
         
         !If all regions have less that nhits_max then break out of do loop
         if (nregion_upgrade .eq. 0 ) goto 200
         
         !upgrade region file and hitmap file by one nside.
         call HealpixMap_udgrade(hitmap_lowres, hitmap_nextres, &
              hitmap_lowres%nside*2)
         hitmap_nextres%TQU(:,1) = hitmap_nextres%TQU(:,1)/4.
         call HealpixMap_udgrade(regions_lowres, regions_nextres, &
              regions_lowres%nside*2)
         
         !Work out the 4 pixels where region=rr and increase the region
         !number of three pixels - increasing N_region accordingly.
         do rr = 1, Nregion_upgrade
            counter=0
            do pp = 0, regions_nextres%npix-1
               if (regions_nextres%tqu(pp,1) .eq. region_upgrade(rr)) then
                  counter = counter + 1
                  if (counter .gt. 1) then
                     n_region                  = n_region+1
                     regions_nextres%tqu(pp,1) = n_region
                  end if
               end if
            end do            
         enddo


         !Deallocate low res maps
         call HealpixMap_free(hitmap_lowres)
         call HealpixMap_free(regions_lowres)

         !Reinitialise low res maps at next resolution
         call HealpixMap_init(hitmap_lowres,nside2npix(hitmap_nextres%nside),1)
         call HealpixMap_init(regions_lowres,nside2npix(regions_nextres%nside),1)
         
         hitmap_lowres%tqu(:,:)  = hitmap_nextres%tqu(:,:)
         regions_lowres%tqu(:,:) = regions_nextres%tqu(:,:)

         !Deallocate
         call HealpixMap_free(hitmap_nextres)
         call HealpixMap_free(regions_nextres)

         deallocate(region_upgrade)
      end do   

200   continue
      
      !Upgrade region file to original resolution
      if(regions_lowres%nside .lt. hitmap_highres%nside) then
         call HealpixMap_udgrade(regions_lowres, regions_highres, &
              hitmap_highres%nside)
         
      endif



      !Multiply by mask




      !Ud grade mask
!      if(verbose) print *,trim(CODE),': ud_grading mask from nside ',&
!           hitmap%nside,' to ',nside
!      call HealpixMap_udgrade(hitmap, maskout, nside)
      
      !Clean up mask
!      where( maskout%TQU(:,1) .lt. 1.0)  maskout%TQU(:,1)=0
      
      
      deallocate(pixlist)      

      where( regions_highres%TQU(:,1) .eq. 0)  regions_highres%TQU(:,1)=undef
      if(verbose) print *,trim(CODE),': writing ',trim(outfile)
      call HealpixMap_write(regions_highres,trim(outfile),units=' ')

      

    end program hitmap_to_regionfile
