module parallel

!======== specification part ==========

  use mpi

  implicit none

contains

!======== implementation part =========

subroutine init_mpi()

  implicit none

  integer :: ier

  ! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)
  if (ier /= 0 ) stop 'Error initializing MPI'

end subroutine init_mpi

!
!-----------------------------------
!

subroutine finalize_mpi()

  implicit none

  integer :: ier

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)
  if (ier /= 0 ) stop 'Error finalizing MPI'

  end subroutine finalize_mpi

!
!---------------------------------
!

subroutine synchronize_all()

  implicit none

  integer :: ier

  ! synchronizes MPI processes
  call MPI_BARRIER(MPI_COMM_WORLD, ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes'

end subroutine synchronize_all

!
!---------------------------------
!

subroutine world_size(sizeval)

  implicit none

  integer,intent(out) :: sizeval

  ! local parameters
  integer :: ier

  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeval,ier)
  if (ier /= 0 ) stop 'Error getting MPI world size'

end subroutine world_size

!
!---------------------------------
!

subroutine world_rank(rank)

  implicit none

  integer,intent(out) :: rank

  ! local parameters
  integer :: ier

  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ier)
  if (ier /= 0 ) stop 'Error getting MPI rank'

end subroutine world_rank

!
!---------------------------------
!

subroutine sum_all_dp(sendbuf, recvbuf)

  implicit none

  double precision :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ier)

end subroutine sum_all_dp


end module parallel