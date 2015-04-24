! set these parameters correctly before compiling

! real kind
integer, parameter :: RK=SIZE_DOUBLE

! mesh layers
integer, parameter :: NLAYER_CRUST=3
integer, dimension(8), parameter :: ID_MESH_LAYER=(/10,11,12,2,3,40,41,5/)

! GLL points (NGLLX=5)
real(RK), dimension(NGLLX), parameter :: XGLL=(/-1.0,-sqrt(21.0)/7,0.0,sqrt(21.0)/7,1.0/)

! directories
character(MAX_STRING_LEN), parameter :: DATDIR1='../DATABASES_MPI_OLD/'
character(MAX_STRING_LEN), parameter :: DATDIR2='../DATABASES_MPI_NEW/'
character(MAX_STRING_LEN), parameter :: OUTDIR2='GLL_interp/'

! model parameter names
integer, parameter :: NPARS=6
character(MAX_STRING_LEN), dimension(npars), parameter :: PARNAME = &
  (/character(MAX_STRING_LEN) :: "vpv","vph","vsv","vsh","eta","rho"/)

! region id: must be crust_mantle=1
integer, parameter :: IREG=1

! old mesh dimensions(get from values_from_mesher.h: NSPEC_CRUST_MANTLE)
integer, parameter :: NPROC1=2, NSPEC1_MAX=7236, NGLOB1_MAX=485401

! new mesh dimensions
integer, parameter :: NPROC2=2, NSPEC2_MAX=26240, NGLOB2_MAX=1742485

! for locating grid(1.5 typical element size of old mesh)
real(RK), parameter :: MARGIN=80.0/288*3.1415926/180*2.0

integer, dimension(NSPEC1_MAX), parameter :: IND1_MAX = (/(II,II=1,NSPEC1_MAX)/)
integer, dimension(NSPEC2_MAX), parameter :: IND2_MAX = (/(II,II=1,NSPEC2_MAX)/)
