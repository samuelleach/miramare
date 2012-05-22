module mpi_stop

!#ifdef MPIPIX
!include "mpif.h"
!#endif

double precision    :: MPI_StartTime



contains

  subroutine DoStop(S)
    !SL: Warning - need to rework this function.
    character(LEN=*), intent(in), optional :: S
    integer ierror,feedback,mpirank,slow_proposals
    if (present(S)) write (*,*) trim(S)
#ifdef MPI 
    MPI_StartTime = MPI_WTime() - MPI_StartTime 
    if (Feedback > 0 .and. MPIRank==0) then
       
       write (*,*) 'Total time:', nint(MPI_StartTime), &
            '(',MPI_StartTime/(60*60),' hours)'
       
       write (*,*) 'Slow proposals: ', slow_proposals
       
    end if
    call mpi_finalize(ierror)
#endif
    
#ifdef DECONLY

    pause

#endif
    stop
  end subroutine DoStop


  subroutine DoAbort(S)

    character(LEN=*), intent(in), optional :: S
    
    integer ierror
    
    if (present(S)) write (*,*) trim(S)
    
#ifdef MPI 
    
    call MPI_Abort(MPI_COMM_WORLD,ierror,ierror)
    
#endif
    
#ifdef DECONLY
    
    pause

#endif
    
    stop
    
  end subroutine DoAbort
  





end module mpi_stop
