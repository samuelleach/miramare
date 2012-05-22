module InstrumentData

  type Instrument
     real*8, dimension(:),   pointer :: channel_nu0_ghz    !Central frequency of band
     real*8, dimension(:,:), pointer :: bandpass_ghz       !Frequency abscissa of bandpasses.
     real*8, dimension(:,:), pointer :: bandpass_weight    !Weight of bandpasses.
     real*8, dimension(:,:), pointer :: channel_raw_offset !Average value of data at each channel.
     integer nabscissa
  end type Instrument
  
contains
  
  subroutine Instrument_init(Inst,nchannel,nstokes,nabscissa)
    type(Instrument) Inst
    integer nchannel,nstokes,nabscissa
    
    allocate(Inst%channel_nu0_ghz(nchannel))
    allocate(Inst%bandpass_ghz(nabscissa,nchannel))
    allocate(Inst%bandpass_weight(nabscissa,nchannel))
    allocate(Inst%channel_raw_offset(nchannel,nstokes))
    Inst%nabscissa = nabscissa
  end subroutine Instrument_init
  
end module InstrumentData
