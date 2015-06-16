
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
 ! maximum number of points per region =       578809
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      1736427
 !
 ! total elements per slice =         9168
 ! total points per slice =       611703
 !
 ! the time step of the solver will be DT =   0.571646929    
 !
 ! total for full 1-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !    36672.000000000000     
 ! approximate total number of points in entire mesh = 
 !    2446812.0000000000     
 ! approximate total number of degrees of freedom in entire mesh = 
 !    7088844.0000000000     
 !
 ! position of the mesh chunk at the surface:
 ! -----------------------------------------
 !
 ! angular size in first direction in degrees =    60.0000000    
 ! angular size in second direction in degrees =    60.0000000    
 !
 ! longitude of center in degrees =    119.500000    
 ! latitude of center in degrees =    40.5000000    
 !
 ! angle of rotation of the first chunk =    50.0000000    
 !
 ! corner            1
 ! longitude in degrees =    122.66077594750357     
 ! latitude in degrees =    1.3825662849418323     
 !
 ! corner            2
 ! longitude in degrees =    168.21552814111038     
 ! latitude in degrees =    33.198170329427846     
 !
 ! corner            3
 ! longitude in degrees =    74.260657455089188     
 ! latitude in degrees =    27.617986690585862     
 !
 ! corner            4
 ! longitude in degrees =    102.45792237405834     
 ! latitude in degrees =    79.229593022618246     
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =          384
 ! GLL points along a great circle =         1536
 ! average distance between points in degrees =    4.09061555E-03
 ! average distance between points in km =    26.0613117    
 ! average size of a spectral element in km =    104.245247    
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
 ! size of static arrays per slice =    247.78929600000001       MB
 !                                 =    236.31028747558594       MiB
 !                                 =   0.24778929599999999       GB
 !                                 =   0.23077176511287689       GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    991.15718400000003       MB
 !                                      =    945.24114990234375       MiB
 !                                      =   0.99115718399999997       GB
 !                                      =   0.92308706045150757       GiB
 !                                      =    9.9115718400000002E-004  TB
 !                                      =    9.0145220747217536E-004  TiB
 !

 integer, parameter :: NEX_XI_VAL =           64
 integer, parameter :: NEX_ETA_VAL =           64

 integer, parameter :: NSPEC_CRUST_MANTLE =         8704
 integer, parameter :: NSPEC_OUTER_CORE =          448
 integer, parameter :: NSPEC_INNER_CORE =           16

 integer, parameter :: NGLOB_CRUST_MANTLE =       578809
 integer, parameter :: NGLOB_OUTER_CORE =        31449
 integer, parameter :: NGLOB_INNER_CORE =         1445

 integer, parameter :: NSPECMAX_ANISO_IC =            1

 integer, parameter :: NSPECMAX_ISO_MANTLE =         8704
 integer, parameter :: NSPECMAX_TISO_MANTLE =         8704
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            1

 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =         8704
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =           16

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =         8704
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =           16

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

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =         8704
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =          448

 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =            1

 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.

 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.

 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.

 logical, parameter :: ATTENUATION_VAL = .true.

 logical, parameter :: ATTENUATION_3D_VAL = .false.

 logical, parameter :: ELLIPTICITY_VAL = .false.

 logical, parameter :: GRAVITY_VAL = .false.

 logical, parameter :: OCEANS_VAL = .false.

 integer, parameter :: NX_BATHY_VAL = 1
 integer, parameter :: NY_BATHY_VAL = 1

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
 integer, parameter :: ATT4_VAL =         8704
 integer, parameter :: ATT5_VAL =           16

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          400
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          400
 integer, parameter :: NSPEC2D_BOTTOM_CM =           64
 integer, parameter :: NSPEC2D_TOP_CM =         1024
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =            4
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =            4
 integer, parameter :: NSPEC2D_BOTTOM_IC =           16
 integer, parameter :: NSPEC2D_TOP_IC =           16
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =           64
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =           64
 integer, parameter :: NSPEC2D_BOTTOM_OC =           16
 integer, parameter :: NSPEC2D_TOP_OC =           64
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

 integer, parameter :: NGLOB_XY_CM =       578809
 integer, parameter :: NGLOB_XY_IC =            1

 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.

 logical, parameter :: FORCE_VECTORIZATION_VAL = .false.

 integer, parameter :: NT_DUMP_ATTENUATION =    100000000

 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    60.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    60.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =    40.500000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =   119.500000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =    50.000000

