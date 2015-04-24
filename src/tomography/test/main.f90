program main

use test 

implicit none

integer :: nglob

nglob = 10

type(mesh(nglob)) :: data
!type (mesh2) :: data2

data%x(:) = 11
!data2%m = 1

print *, data%x

!integer, parameter :: IIN = 40
!integer, parameter :: MAX_STRING_LEN = 512
!character(len=*), parameter :: dmodel_dkernel_list = 'info'
!
!character(len=1024) :: dummy
!integer :: i, ier, nstep_lbfgs
!character(len=MAX_STRING_LEN), dimension(:), allocatable :: dmodel_dir, dkernel_dir
!double precision, dimension(:), allocatable :: step_length, dm_dot_dg, dg_dot_dg
!
!open(unit=IIN, file=trim(dmodel_dkernel_list), status='old', action='read', iostat=ier)
!if (ier /= 0) then
!  write(*,*) 'Error open file: ', trim(dmodel_dkernel_list)
!  stop
!endif
!nstep_lbfgs = 0 ! get number of lines
!do
!  read(IIN,'(a)',iostat=ier) dummy
!  if (ier /= 0) exit
!  print *, 'dummy=', trim(dummy)
!  nstep_lbfgs = nstep_lbfgs + 1
!enddo
!
!allocate(dmodel_dir(nstep_lbfgs), dkernel_dir(nstep_lbfgs))
!allocate(step_length(nstep_lbfgs), dm_dot_dg(nstep_lbfgs), dg_dot_dg(nstep_lbfgs))
!
!rewind(IIN)
!do i = 1, nstep_lbfgs
!  read(IIN,*) dmodel_dir(i), dkernel_dir(i), step_length(i), dm_dot_dg(i), dg_dot_dg(i)
!enddo
!close(IIN)
!
!print *, 'nstep_lbfgs=', nstep_lbfgs
!print *, dmodel_dir
!print *, dkernel_dir
!print *, step_length
!print *, dm_dot_dg 
!print *, dg_dot_dg 

!write(*,'(a)') &
!
!'Name', &
!'', &
!' main', &
!'', &
!'Description', &
!'', &
!' haha'

!real :: x, y
!
!x = -1234.56789
!y = 1234.56789
!
!write(*,'(a,SP,2ES15.7)') 'x,y = ', x, y

!integer, parameter :: n = 10000
!
!real, dimension(n,n) :: array
!real :: sum_array
!
!integer :: m
!real :: smallvar

!m = size(array)
!smallvar = 1.0/m
!array = smallvar
!
!call Kahan_sum(array, sum_array)
!
!print *, 'number of smallvar to sum=', m
!print *, 'smallval=', smallvar
!print *, 'Kahan_summation=', sum_array
!print *, 'plain summation=', sum(array)

!character(len=10) :: arg1 
!real :: n
!
!call get_command_argument(1, arg1)
!
!
!read(arg1,*) n
!
!print *, 'n=', n

!call show_vars()
!call tryit()

contains

subroutine Kahan_sum(array, sum_array)

  real, dimension(:,:), intent(in) :: array
  real, intent(out) :: sum_array

  integer :: i
  real :: t, y, c
  
  c = 0.0
  sum_array = 0.0
  do i = 1, size(array)
    y = array(i,1) - c
    t = sum_array + y
    c = (t - sum_array) - y
    sum_array = t
  end do

end subroutine

end program