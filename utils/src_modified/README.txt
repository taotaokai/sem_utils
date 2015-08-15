1. src/specfem3D/locate_source.f90:

  - typical_size = (ANGULAR_WIDTH_XI_IN_DEGREES_VAL*DEGREES_TO_RADIANS/NEX_XI_VAL) * R_UNIT_SPHERE

    To calculate typical correctly according to the header files.
    
  - because I made the source mask file by my own program, the modification on typical_size is irrelevant now. 

2. src/specfem3D/write_output_SAC.f90:

  - correct the sub milli-second error

3. src/shared/get_model_parameters.F90:

!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!  PROBLEM:
!!! 
!!! When I tried to read in GLL model for mesh created using
!!! 1D_(transversely_)isotropic_prem error occurs because the mesh created by setting MODEL = GLL in 
!!! Par_file is incompatible with model created by MODEL = 1D_(transversely_)isotropic_prem.
!!! 
!!! SOLUTION:
!!! 
!!!   1) introduce a new model identifier MODEL = 1D_(transversely_)isotropic_GLL
!!! 
!!! #---- src/shared/get_model_parameters.F90:
!!!  else if (MODEL_ROOT == '1D_isotropic_GLL') then
!!!      HONOR_1D_SPHERICAL_MOHO = .true.
!!!      THREE_D_MODEL = THREE_D_MODEL_GLL ! use this because the package doesn't
!!!                                        ! implement ONE_D_MODEL_GLL 
!!!  else if (MODEL_ROOT == '1D_transversely_isotropic_GLL') then
!!!      HONOR_1D_SPHERICAL_MOHO = .true.
!!!      TRANSVERSE_ISOTROPY = .true.
!!!      THREE_D_MODEL = THREE_D_MODEL_GLL ! use this because the package doesn't
!!!                                        ! implement ONE_D_MODEL_GLL 
!!! #----
!!! 
!!!   set MODEL = 1D_(transversely_)isotropic_GLL in DATA/Par_file as well.
!!! 
!!!   2) because GLL model is intended for 3D model, so we need to avoid any mesh
!!! streching due to internal boundary topography.
!!! 
!!! #---- setup/constants.h.in:
!!! ! to suppress element stretching for 3D moho surface
!!!   logical,parameter :: SUPPRESS_MOHO_STRETCHING = .true.
!!! 
!!! ! to suppress element stretching at 410/660 internal topography
!!! ! (i.e. creates mesh without 410/660 topography for Harvard model
!!! (s362ani,..))
!!!   logical,parameter :: SUPPRESS_INTERNAL_TOPOGRAPHY = .true.
!!! #----
!!! 
!!!   3) don't forget to re-configure before compiling
!!! 
!!!   ./configure FC=ifort MPIFC=mpif90 # this will update the setup/constants.h
!!!   ./make xmeshfem3D xspecfem3D$ 
!!!    
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

4. configure and configure.ac (2015-06-16T10:24:17)

  - change configure and configure.ac: AC_CONFIG_FILES:
    Par_file,CMTSOLUTION,STATIONS -> AC_CONFIG_LINKS. 
    This avoids configure replacing links of these files with the copy, because I usually
    link these control files from another directory.

5. compute_kernels.F90 (2015-06-21T20:15:12)

  - compute_kernels.F90: compute_kernels_hessian()
    use maximum product between forward and backward wavefield amplitudes as
    the approximated illumination amplitude.
  - further modification needed. 

6. about the compiler flags (2015-08-14T19:07:43)
  
  - Simple optimization flag -O 3 is enough. So other flags may confuse the
    compiler and result in less effecient code (e.g. upto 5 X slower).
