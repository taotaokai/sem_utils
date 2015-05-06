DESCRIPTION
----

  this package is used for post-processing with specfem3d\_globe


LIBRARIES 
----

  fortran90 modules: 

  sem\_constants, sem\_parallel: 

  sem\_io, sem\_mesh, sem\_tomography:


PROGRAMS
----

  xsem\_measure\_adj: measure adjoint source

  xsem\_interpolate: interpolate SEM mesh data(e.g. make new mesh or create profiles)

NOTES
----
  * The header files in include/ are setup/\*.h OUTPUT\_FILES/values\_from\_mesher.h
    these files are generated during compilation of specfed3d_globe, you should use the header files for your own application
