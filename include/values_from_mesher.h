
 !
 ! this is the parameter file for static compilation of the solver
 !
 ! mesh statistics:
 ! ---------------
 !
 !
 ! number of chunks =            1
 !
 ! these statistics do not include the central cube
 !
 ! number of processors =            4
 !
 ! maximum number of points per region =      1506509
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      4519527
 !
 ! total elements per slice =        25075
 ! total points per slice =      1663015
 !
 ! the time step of the solver will be DT =   0.152439177    
 !
 ! total for full 1-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !    100300.00000000000     
 ! approximate total number of points in entire mesh = 
 !    6652060.0000000000     
 ! approximate total number of degrees of freedom in entire mesh = 
 !    18764108.000000000     
 !
 ! position of the mesh chunk at the surface:
 ! -----------------------------------------
 !
 ! angular size in first direction in degrees =    20.0000000    
 ! angular size in second direction in degrees =    20.0000000    
 !
 ! longitude of center in degrees =    135.983994    
 ! latitude of center in degrees =    37.6590004    
 !
 ! angle of rotation of the first chunk =   -6.69999981    
 !
 ! corner            1
 ! longitude in degrees =    123.52994859379010     
 ! latitude in degrees =    28.460309587776248     
 !
 ! corner            2
 ! longitude in degrees =    145.61020836088025     
 ! latitude in degrees =    26.414880872256532     
 !
 ! corner            3
 ! longitude in degrees =    123.03384133337009     
 ! latitude in degrees =    48.188913751221484     
 !
 ! corner            4
 ! longitude in degrees =    151.66195832275645     
 ! latitude in degrees =    45.550442233913806     
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =         1440
 ! GLL points along a great circle =         5760
 ! average distance between points in degrees =    1.09083077E-03
 ! average distance between points in km =    6.94968271    
 ! average size of a spectral element in km =    27.7987309    
 !

 ! approximate static memory needed by the solver:
 ! ----------------------------------------------
 !
 ! (lower bound, usually the real amount used is 5% to 10% higher)
 !
 ! (you can get a more precise estimate of the size used per MPI process
 !  by typing "size -d bin/xspecfem3D"
 !  after compiling the code with the DATA/Par_file you plan to use)
 !
 ! size of static arrays per slice =    723.24887200000001       MB
 !                                 =    689.74387359619141       MiB
 !                                 =   0.72324887199999999       GB
 !                                 =   0.67357800155878067       GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    2892.9954880000000       MB
 !                                      =    2758.9754943847656       MiB
 !                                      =    2.8929954879999999       GB
 !                                      =    2.6943120062351227       GiB
 !                                      =    2.8929954880000000E-003  TB
 !                                      =    2.6311640685889870E-003  TiB
 !

 integer, parameter :: NEX_XI_VAL =           80
 integer, parameter :: NEX_ETA_VAL =           80

 integer, parameter :: NSPEC_CRUST_MANTLE =        22800
 integer, parameter :: NSPEC_OUTER_CORE =         2175
 integer, parameter :: NSPEC_INNER_CORE =          100

 integer, parameter :: NGLOB_CRUST_MANTLE =      1506509
 integer, parameter :: NGLOB_OUTER_CORE =       149009
 integer, parameter :: NGLOB_INNER_CORE =         7497

 integer, parameter :: NSPECMAX_ANISO_IC =            1

 integer, parameter :: NSPECMAX_ISO_MANTLE =        22800
 integer, parameter :: NSPECMAX_TISO_MANTLE =        22800
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            1

 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =        22800
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =          100

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =        22800
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =          100

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_AND_ATT =            1
 integer, parameter :: NSPEC_INNER_CORE_STR_AND_ATT =            1

 integer, parameter :: NSPEC_CRUST_MANTLE_STRAIN_ONLY =            1
 integer, parameter :: NSPEC_INNER_CORE_STRAIN_ONLY =            1

 integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT =            1
 integer, parameter :: NSPEC_OUTER_CORE_ADJOINT =            1
 integer, parameter :: NSPEC_INNER_CORE_ADJOINT =            1
 integer, parameter :: NGLOB_CRUST_MANTLE_ADJOINT =            1
 integer, parameter :: NGLOB_OUTER_CORE_ADJOINT =            1
 integer, parameter :: NGLOB_INNER_CORE_ADJOINT =            1
 integer, parameter :: NSPEC_OUTER_CORE_ROT_ADJOINT =            1

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =        22800
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =         2175

 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =      1506509

 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.

 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.

 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.

 logical, parameter :: ATTENUATION_VAL = .true.

 logical, parameter :: ATTENUATION_3D_VAL = .false.

 logical, parameter :: ELLIPTICITY_VAL = .true.

 logical, parameter :: GRAVITY_VAL = .false.

 logical, parameter :: OCEANS_VAL = .true.

 integer, parameter :: NX_BATHY_VAL =         5400
 integer, parameter :: NY_BATHY_VAL =         2700

 logical, parameter :: ROTATION_VAL = .false.
 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =            1

 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .false.

 integer, parameter :: NPROC_XI_VAL =            2
 integer, parameter :: NPROC_ETA_VAL =            2
 integer, parameter :: NCHUNKS_VAL =            1
 integer, parameter :: NPROCTOT_VAL =            4

 integer, parameter :: ATT1_VAL =            5
 integer, parameter :: ATT2_VAL =            5
 integer, parameter :: ATT3_VAL =            5
 integer, parameter :: ATT4_VAL =        22800
 integer, parameter :: ATT5_VAL =          100

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =         1040
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =         1040
 integer, parameter :: NSPEC2D_BOTTOM_CM =          100
 integer, parameter :: NSPEC2D_TOP_CM =         1600
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =           20
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =           20
 integer, parameter :: NSPEC2D_BOTTOM_IC =           25
 integer, parameter :: NSPEC2D_TOP_IC =           25
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =          270
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =          275
 integer, parameter :: NSPEC2D_BOTTOM_OC =           25
 integer, parameter :: NSPEC2D_TOP_OC =          100
 integer, parameter :: NSPEC2D_MOHO =            1
 integer, parameter :: NSPEC2D_400 =            1
 integer, parameter :: NSPEC2D_670 =            1
 integer, parameter :: NSPEC2D_CMB =            1
 integer, parameter :: NSPEC2D_ICB =            1

 logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.
 integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 1
 integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 1

 integer, parameter :: NSPEC_OUTER_CORE_3DMOVIE = 1
 integer, parameter :: NM_KL_REG_PTS_VAL = 1

 integer, parameter :: NGLOB_XY_CM =      1506509
 integer, parameter :: NGLOB_XY_IC =            1

 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.

 logical, parameter :: FORCE_VECTORIZATION_VAL = .false.

 integer, parameter :: NT_DUMP_ATTENUATION =    100000000

 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    20.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    20.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =    37.659000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =   135.984000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =    -6.700000

