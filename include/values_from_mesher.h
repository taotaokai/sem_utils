 
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
 ! number of processors =          144
 !
 ! maximum number of points per region =       485401
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      1456203
 !
 ! total elements per slice =         7992
 ! total points per slice =       539563
 !
 ! the time step of the solver will be DT =   0.1693769    
 !
 ! total for full 1-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !    1150848.00000000     
 ! approximate total number of points in entire mesh = 
 !    77697072.0000000     
 ! approximate total number of degrees of freedom in entire mesh = 
 !    218319984.000000     
 !
 ! position of the mesh chunk at the surface:
 ! -----------------------------------------
 !
 ! angular size in first direction in degrees =    80.00000    
 ! angular size in second direction in degrees =    80.00000    
 !
 ! longitude of center in degrees =    106.0000    
 ! latitude of center in degrees =    32.00000    
 !
 ! angle of rotation of the first chunk =    60.00000    
 !
 ! corner            1
 ! longitude in degrees =    117.915778257699     
 ! latitude in degrees =   -16.6586412980043     
 !
 ! corner            2
 ! longitude in degrees =    165.126222399968     
 ! latitude in degrees =    30.7867888632000     
 !
 ! corner            3
 ! longitude in degrees =    57.4074304719772     
 ! latitude in degrees =    10.0650870551211     
 !
 ! corner            4
 ! longitude in degrees =    54.0787269587550     
 ! latitude in degrees =    75.5308108875256     
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =         1296
 ! GLL points along a great circle =         5184
 ! average distance between points in degrees =   1.2120343E-03
 ! average distance between points in km =    7.721870    
 ! average size of a spectral element in km =    30.88748    
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
 ! size of static arrays per slice =    119.386312000000       MB
 !                                 =    113.855659484863       MiB
 !                                 =   0.119386312000000       GB
 !                                 =   0.111187167465687       GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    17.1916289280000       GB
 !                                      =    16.0109521150589       GiB
 !                                      =   1.719162892800000E-002  TB
 !                                      =   1.563569542486221E-002  TiB
 !
 
 integer, parameter :: NEX_XI_VAL =          288
 integer, parameter :: NEX_ETA_VAL =          288
 
 integer, parameter :: NSPEC_CRUST_MANTLE =         7236
 integer, parameter :: NSPEC_OUTER_CORE =          720
 integer, parameter :: NSPEC_INNER_CORE =           36
 
 integer, parameter :: NGLOB_CRUST_MANTLE =       485401
 integer, parameter :: NGLOB_OUTER_CORE =        51289
 integer, parameter :: NGLOB_INNER_CORE =         2873
 
 integer, parameter :: NSPECMAX_ANISO_IC =            1
 
 integer, parameter :: NSPECMAX_ISO_MANTLE =         7236
 integer, parameter :: NSPECMAX_TISO_MANTLE =         7236
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            1
 
 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =            1
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =            1
 
 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =            1
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =            1
 
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
 
 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =         7236
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =          720
 
 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =            1
 
 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.
 
 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.
 
 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.
 
 logical, parameter :: ATTENUATION_VAL = .false.
 
 logical, parameter :: ATTENUATION_3D_VAL = .false.
 
 logical, parameter :: ELLIPTICITY_VAL = .false.
 
 logical, parameter :: GRAVITY_VAL = .false.
 
 logical, parameter :: OCEANS_VAL = .false.
 
 integer, parameter :: NX_BATHY_VAL = 1
 integer, parameter :: NY_BATHY_VAL = 1
 
 logical, parameter :: ROTATION_VAL = .false.
 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =            1
 
 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .false.
 
 integer, parameter :: NPROC_XI_VAL =           12
 integer, parameter :: NPROC_ETA_VAL =           12
 integer, parameter :: NCHUNKS_VAL =            1
 integer, parameter :: NPROCTOT_VAL =          144
 
 integer, parameter :: ATT1_VAL =            1
 integer, parameter :: ATT2_VAL =            1
 integer, parameter :: ATT3_VAL =            1
 integer, parameter :: ATT4_VAL =            1
 integer, parameter :: ATT5_VAL =            1
 
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          534
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          534
 integer, parameter :: NSPEC2D_BOTTOM_CM =           36
 integer, parameter :: NSPEC2D_TOP_CM =          576
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =           12
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =           12
 integer, parameter :: NSPEC2D_BOTTOM_IC =            9
 integer, parameter :: NSPEC2D_TOP_IC =            9
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =          147
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =          150
 integer, parameter :: NSPEC2D_BOTTOM_OC =            9
 integer, parameter :: NSPEC2D_TOP_OC =           36
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
 
 integer, parameter :: NGLOB_XY_CM =       485401
 integer, parameter :: NGLOB_XY_IC =            1
 
 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.
 
 logical, parameter :: FORCE_VECTORIZATION_VAL = .false.
 
 integer, parameter :: NT_DUMP_ATTENUATION =    100000000
 
 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    80.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    80.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =    32.000000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =   106.000000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =    60.000000
 
