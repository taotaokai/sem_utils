# DESCRIPTION

  this package is used for post-processing with specfem3d\_globe

# LIBRARIES 

  fortran90 modules: 

  sem\_constants\_mod, sem\_parallel\_mod

  sem\_io\_mod, sem\_mesh\_mod, sem\_tomography\_mod


# PROGRAMS

  xsem\_measure\_adj: measure adjoint source

  xsem\_interpolate: interpolate SEM mesh data(e.g. make new mesh or create profiles)

# before compilation

  * The header files in include/ are "setup/(constants.h  precision.h) OUTPUT\_FILES/values\_from\_mesher.h",
    which are generated during compilation of specfed3d_globe, 
    you should use the header files for your own application.

  * set correct path to netcdf include and lib directories in Makefile
