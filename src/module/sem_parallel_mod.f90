module sem_parallel

!======== specification part ==========

  use mpi
  use sem_constants, only: CUSTOM_REAL

  implicit none

!======== implementation part =========
contains

!///////////////////////////////////
subroutine init_mpi()

  implicit none

  integer :: ier

  ! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)
  if (ier /= 0 ) stop 'Error initializing MPI'

end subroutine init_mpi

!///////////////////////////////////
subroutine finalize_mpi()

  implicit none

  integer :: ier

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)
  if (ier /= 0 ) stop 'Error finalizing MPI'

  end subroutine finalize_mpi

!///////////////////////////////////
subroutine abort_mpi()

  use mpi
  
  implicit none
  
  integer :: ier
  
  ! note: MPI_ABORT does not return, and does exit the
  !          program with an error code of 30
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
  
  end subroutine abort_mpi

!///////////////////////////////////
subroutine synchronize_all()

  implicit none

  integer :: ier

  ! synchronizes MPI processes
  call MPI_BARRIER(MPI_COMM_WORLD, ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes'

end subroutine synchronize_all

!///////////////////////////////////
subroutine world_size(sizeval)

  implicit none

  integer,intent(out) :: sizeval

  ! local parameters
  integer :: ier

  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeval,ier)
  if (ier /= 0 ) stop 'Error getting MPI world size'

end subroutine world_size

!///////////////////////////////////
subroutine world_rank(rank)

  implicit none

  integer,intent(out) :: rank

  ! local parameters
  integer :: ier

  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ier)
  if (ier /= 0 ) stop 'Error getting MPI rank'

end subroutine world_rank

!///////////////////////////////////
subroutine sum_all_dp(sendbuf, recvbuf)

  implicit none

  double precision :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ier)

end subroutine sum_all_dp


!///////////////////////////////////
subroutine sum_all_cr(sendbuf, recvbuf)

  implicit none
  
  include 'precision.h'
  
  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
  integer :: ier
  
  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

end subroutine sum_all_cr

!///////////////////////////////////
subroutine max_all_dp(sendbuf, recvbuf)

  use mpi
  
  implicit none

  double precision :: sendbuf, recvbuf
  integer :: ier
  
  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ier)

end subroutine max_all_dp

!///////////////////////////////////
subroutine bcast_all_singledp(buffer)

  use mpi

  implicit none

  double precision :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

end subroutine bcast_all_singledp

!///////////////////////////////////
subroutine bcast_all_singlecr(buffer)

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

end subroutine bcast_all_singlecr


!///////////////////////////////////
subroutine send_i(sendbuf, sendcount, dest, sendtag)

  implicit none

  integer :: dest, sendtag
  integer :: sendcount
  integer, dimension(sendcount) :: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf, sendcount, MPI_INTEGER, dest, sendtag, &
    MPI_COMM_WORLD, ier)

end subroutine send_i


subroutine recv_i(recvbuf, recvcount, dest, recvtag)

  implicit none

  integer :: dest, recvtag
  integer :: recvcount
  integer, dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf, recvcount, MPI_INTEGER, dest, recvtag, &
    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

end subroutine recv_i


!///////////////////////////////////
subroutine send_dp(sendbuf, sendcount, dest, sendtag)

  implicit none
  
  integer :: dest, sendtag
  integer :: sendcount
  double precision, dimension(sendcount) :: sendbuf
  integer :: ier
  
  call MPI_SEND(sendbuf, sendcount, MPI_DOUBLE_PRECISION, dest, sendtag, &
    MPI_COMM_WORLD, ier)

end subroutine send_dp

subroutine send2_dp(sendbuf, sendcount, dest, sendtag)

  implicit none
  
  integer :: dest, sendtag
  integer :: sendcount
  double precision, dimension(:,:) :: sendbuf
  integer :: ier
  
  call MPI_SEND(sendbuf, sendcount, MPI_DOUBLE_PRECISION, dest, sendtag, &
    MPI_COMM_WORLD, ier)

end subroutine send2_dp

subroutine recv_dp(recvbuf, recvcount, dest, recvtag)

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  double precision, dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf, recvcount, MPI_DOUBLE_PRECISION, dest, recvtag, &
    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

end subroutine recv_dp

subroutine recv2_dp(recvbuf, recvcount, dest, recvtag)

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  double precision, dimension(:,:) :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf, recvcount, MPI_DOUBLE_PRECISION, dest, recvtag, &
    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

end subroutine recv2_dp


end module sem_parallel