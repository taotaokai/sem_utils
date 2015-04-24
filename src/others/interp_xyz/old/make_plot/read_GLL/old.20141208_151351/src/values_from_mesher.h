 
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
 ! number of processors =          200
 !
 ! maximum number of points per region =      2746485
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      8239455
 !
 ! total elements per slice =        47100
 ! total points per slice =      3122911
 !
 ! the time step of the solver will be DT =   4.0650446E-02
 !
 ! total for full 1-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !    9420000.00000000     
 ! approximate total number of points in entire mesh = 
 !    624582200.000000     
 ! approximate total number of degrees of freedom in entire mesh = 
 !    1731819800.00000     
 !
 ! position of the mesh chunk at the surface:
 ! -----------------------------------------
 !
 ! angular size in first direction in degrees =    80.00000    
 ! angular size in second direction in degrees =    40.00000    
 !
 ! longitude of center in degrees =    110.0000    
 ! latitude of center in degrees =    35.00000    
 !
 ! angle of rotation of the first chunk =   0.0000000E+00
 !
 ! corner            1
 ! longitude in degrees =    70.7748650839013     
 ! latitude in degrees =    11.8028684327037     
 !
 ! corner            2
 ! longitude in degrees =    149.225134916099     
 ! latitude in degrees =    11.8028684327037     
 !
 ! corner            3
 ! longitude in degrees =    56.0333135386272     
 ! latitude in degrees =    40.2228269156194     
 !
 ! corner            4
 ! longitude in degrees =    163.966686461373     
 ! latitude in degrees =    40.2228269156194     
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =         3200
 ! GLL points along a great circle =        12800
 ! average distance between points in degrees =   2.8124999E-02
 ! average distance between points in km =    3.127357    
 ! average size of a spectral element in km =    12.50943    
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
 ! size of static arrays per slice =    695.368196000000       MB
 !                                 =    663.154788970947       MiB
 !                                 =   0.695368196000000       GB
 !                                 =   0.647612098604441       GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    139.073639200000       GB
 !                                      =    129.522419720888       GiB
 !                                      =   0.139073639200000       TB
 !                                      =   0.126486738008680       TiB
 !
 
 integer, parameter :: NEX_XI_VAL =          800
 integer, parameter :: NEX_ETA_VAL =          400
 
 integer, parameter :: NSPEC_CRUST_MANTLE =        41600
 integer, parameter :: NSPEC_OUTER_CORE =         5200
 integer, parameter :: NSPEC_INNER_CORE =          300
 
 integer, parameter :: NGLOB_CRUST_MANTLE =      2746485
 integer, parameter :: NGLOB_OUTER_CORE =       354817
 integer, parameter :: NGLOB_INNER_CORE =        21609
 
 integer, parameter :: NSPECMAX_ANISO_IC =            1
 
 integer, parameter :: NSPECMAX_ISO_MANTLE =        41600
 integer, parameter :: NSPECMAX_TISO_MANTLE =            1
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
 
 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =        41600
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =         5200
 
 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =      2746485
 
 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .false.
 
 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.
 
 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.
 
 logical, parameter :: ATTENUATION_VAL = .false.
 
 logical, parameter :: ATTENUATION_3D_VAL = .false.
 
 logical, parameter :: ELLIPTICITY_VAL = .true.
 
 logical, parameter :: GRAVITY_VAL = .false.
 
 logical, parameter :: OCEANS_VAL = .true.
 
 integer, parameter :: NX_BATHY_VAL = NX_BATHY
 integer, parameter :: NY_BATHY_VAL = NY_BATHY
 
 logical, parameter :: ROTATION_VAL = .false.
 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =            1
 
 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .true.
 
 integer, parameter :: NPROC_XI_VAL =           20
 integer, parameter :: NPROC_ETA_VAL =           10
 integer, parameter :: NCHUNKS_VAL =            1
 integer, parameter :: NPROCTOT_VAL =          200
 
 integer, parameter :: ATT1_VAL =            1
 integer, parameter :: ATT2_VAL =            1
 integer, parameter :: ATT3_VAL =            1
 integer, parameter :: ATT4_VAL =            1
 integer, parameter :: ATT5_VAL =            1
 
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =         2180
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =         2180
 integer, parameter :: NSPEC2D_BOTTOM_CM =          100
 integer, parameter :: NSPEC2D_TOP_CM =         1600
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =           60
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =           60
 integer, parameter :: NSPEC2D_BOTTOM_IC =           25
 integer, parameter :: NSPEC2D_TOP_IC =           25
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =          645
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =          650
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
 
 integer, parameter :: NGLOB_XY_CM =      2746485
 integer, parameter :: NGLOB_XY_IC =            1
 
 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.
 
 logical, parameter :: FORCE_VECTORIZATION_VAL = .true.
 
 integer, parameter :: NT_DUMP_ATTENUATION =    100000000
