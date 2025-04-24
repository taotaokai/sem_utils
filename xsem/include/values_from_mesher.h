 
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
 ! number of processors =          100
 !
 ! maximum number of points per region =       438161
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      1314483
 !
 ! total elements per slice =         7245
 ! total points per slice =       490499
 !
 ! the time step of the solver will be DT =   0.1150000      (s)
 ! the (approximate) minimum period resolved will be =    14.00232      (s)
 !
 ! total for full 1-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !    724500.000000000     
 ! approximate total number of points in entire mesh = 
 !    49049900.0000000     
 ! approximate total number of degrees of freedom in entire mesh = 
 !    137391900.000000     
 !
 ! position of the mesh chunk at the surface:
 ! -----------------------------------------
 !
 ! angular size in first direction in degrees =    68.00000    
 ! angular size in second direction in degrees =    68.00000    
 !
 ! longitude of center in degrees =    7.000000    
 ! latitude of center in degrees =    52.00000    
 !
 ! angle of rotation of the first chunk =    45.00000    
 !
 ! corner            1
 ! longitude in degrees =    6.99999999999999     
 ! latitude in degrees =    8.40680071311899     
 !
 ! corner            2
 ! longitude in degrees =    64.1611658328499     
 ! latitude in degrees =    34.9434448363920     
 !
 ! corner            3
 ! longitude in degrees =   -50.1611658328499     
 ! latitude in degrees =    34.9434448363920     
 !
 ! corner            4
 ! longitude in degrees =   -173.000000000000     
 ! latitude in degrees =    84.3889812060761     
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =         1270
 ! GLL points along a great circle =         5082
 ! average distance between points in degrees =   7.0833333E-02
 ! average distance between points in km =    7.876307    
 ! average size of a spectral element in km =    31.50523    
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
 ! size of static arrays per slice =    335.947664000000       MB
 !                                 =    320.384658813477       MiB
 !                                 =   0.335947664000000       GB
 !                                 =   0.312875643372536       GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    33.5947664000000       GB
 !                                      =    31.2875643372536       GiB
 !                                      =   3.359476640000000E-002  TB
 !                                      =   3.055426204809919E-002  TiB
 !
 
 integer, parameter :: NEX_XI_VAL =          240
 integer, parameter :: NEX_ETA_VAL =          240
 
 integer, parameter :: NSPEC_CRUST_MANTLE =         6516
 integer, parameter :: NSPEC_OUTER_CORE =          684
 integer, parameter :: NSPEC_INNER_CORE =           45
 
 integer, parameter :: NGLOB_CRUST_MANTLE =       438161
 integer, parameter :: NGLOB_OUTER_CORE =        48789
 integer, parameter :: NGLOB_INNER_CORE =         3549
 
 integer, parameter :: NSPECMAX_ANISO_IC =            0
 
 integer, parameter :: NSPECMAX_ISO_MANTLE =         6516
 integer, parameter :: NSPECMAX_TISO_MANTLE =         6516
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            0
 
 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =         6516
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =           45
 
 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =         6516
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =           45
 
 integer, parameter :: NSPEC_CRUST_MANTLE_STR_AND_ATT =            0
 integer, parameter :: NSPEC_INNER_CORE_STR_AND_ATT =            0
 
 integer, parameter :: NSPEC_CRUST_MANTLE_STRAIN_ONLY =            0
 integer, parameter :: NSPEC_INNER_CORE_STRAIN_ONLY =            0
 
 integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT =            0
 integer, parameter :: NSPEC_OUTER_CORE_ADJOINT =            0
 integer, parameter :: NSPEC_INNER_CORE_ADJOINT =            0
 integer, parameter :: NGLOB_CRUST_MANTLE_ADJOINT =            0
 integer, parameter :: NGLOB_OUTER_CORE_ADJOINT =            0
 integer, parameter :: NGLOB_INNER_CORE_ADJOINT =            0
 integer, parameter :: NSPEC_OUTER_CORE_ROT_ADJOINT =            0
 
 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =         6516
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =          684
 
 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =       438161
 
 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.
 
 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.
 
 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.
 
 logical, parameter :: ATTENUATION_VAL = .true.
 
 logical, parameter :: ATTENUATION_3D_VAL = .false.
 
 logical, parameter :: ELLIPTICITY_VAL = .true.
 
 logical, parameter :: GRAVITY_VAL = .true.
 
 logical, parameter :: OCEANS_VAL = .true.
 
 integer, parameter :: NX_BATHY_VAL =         5400
 integer, parameter :: NY_BATHY_VAL =         2700
 
 logical, parameter :: ROTATION_VAL = .true.
 logical, parameter :: EXACT_MASS_MATRIX_FOR_ROTATION_VAL = .false.
 
 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =          684
 
 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .false.
 
 integer, parameter :: NPROC_XI_VAL =           10
 integer, parameter :: NPROC_ETA_VAL =           10
 integer, parameter :: NCHUNKS_VAL =            1
 integer, parameter :: NPROCTOT_VAL =          100
 
 integer, parameter :: ATT1_VAL =            5
 integer, parameter :: ATT2_VAL =            5
 integer, parameter :: ATT3_VAL =            5
 integer, parameter :: ATT4_VAL =         6516
 integer, parameter :: ATT5_VAL =           45
 
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          498
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          498
 integer, parameter :: NSPEC2D_BOTTOM_CM =           36
 integer, parameter :: NSPEC2D_TOP_CM =          576
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =           15
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =           15
 integer, parameter :: NSPEC2D_BOTTOM_IC =            9
 integer, parameter :: NSPEC2D_TOP_IC =            9
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =          141
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =          144
 integer, parameter :: NSPEC2D_BOTTOM_OC =            9
 integer, parameter :: NSPEC2D_TOP_OC =           36
 integer, parameter :: NSPEC2D_MOHO =            1
 integer, parameter :: NSPEC2D_400 =            1
 integer, parameter :: NSPEC2D_670 =            1
 integer, parameter :: NSPEC2D_CMB =            1
 integer, parameter :: NSPEC2D_ICB =            1
 
 logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.
 integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 0
 integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 0
 
 integer, parameter :: NSPEC_OUTER_CORE_3DMOVIE = 0
 integer, parameter :: NGLOB_XY_CM =       438161
 integer, parameter :: NGLOB_XY_IC =            0
 
 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.
 
 logical, parameter :: FORCE_VECTORIZATION_VAL = .false.
 
 logical, parameter :: UNDO_ATTENUATION_VAL = .true.
 integer, parameter :: NT_DUMP_ATTENUATION_VAL =          550
 
 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    68.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    68.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =    52.000000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =     7.000000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =    45.000000
 
