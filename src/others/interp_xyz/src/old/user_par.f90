! set these parameters correctly before compiling

! real kind
integer, parameter :: RK=SIZE_DOUBLE

! GLL points (!!! NGLLX=NGLLY=NGLLZ=5)
real(RK), dimension(NGLLX), parameter :: XGLL=(/-1.0,-sqrt(21.0)/7,0.0,sqrt(21.0)/7,1.0/)

! interpolation points
integer, parameter :: NPTS_MAX = 20000

!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!
! directories
character(MAX_STRING_LEN), parameter :: DATDIR='DATABASES_MPI/'
character(MAX_STRING_LEN), parameter :: OUTDIR='profile.out/'
character(MAX_STRING_LEN), parameter :: OUTFNAME='interp'

! model parameter names
integer, parameter :: NPARS=6
character(MAX_STRING_LEN), dimension(npars), parameter :: PARNAME = &
  (/character(MAX_STRING_LEN) :: "vpv","vph","vsv","vsh","eta","rho"/)

! region id: must be crust_mantle=1
integer, parameter :: IREG=1

!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!
! source mesh dimensions(get from values_from_mesher.h: NSPEC_CRUST_MANTLE)
!integer, parameter :: NPROC=144, NSPEC_MAX=7236, NGLOB_MAX=485401
integer, parameter :: NPROC=576, NSPEC_MAX=26240, NGLOB_MAX=1742485

!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!
! for locating grid(1.5 typical element size of old mesh)
!real(RK), parameter :: MARGIN=80.0/288*3.1415926/180*1.5
real(RK), parameter :: MARGIN=80.0/768*3.1415926/180*2.0 ! larger to avoid host element not-found 

integer, dimension(NSPEC_MAX), parameter :: IND_MAX = (/(II,II=1,NSPEC_MAX)/)
