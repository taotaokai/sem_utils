# DESCRIPTION

* this package is used for post-processing and managing inversion work flow
    with specfem3d_globe.

#====== Contents

1. include/, src/: Fortran codes for SEM data processing

2. utils/: scripts used to manage the inversion work flow 
    * sem_create.sh, sem_build.sh
    * setup_mesh.sh, setup_event.sh, setup_adjoint.sh, measure_misfit.py
    * update_model.sh, update_kernel.sh
    * make_vtk.sh, make_slice_gcircle.sh, make_slice_sphere.sh
    * plot_slice_*.sh, plot_misfit.sh

#====== Folder Structure for an inversion project

project/
*   sem_utils/ # utility scripts and excutables to run the inversion
|   
*   sem_config/ # configuration files for SEM (model setup,)
|   *   setup/ # these header files are used to build specfem3d_globe
|   |   *   constants.h.in
|   |   *   precision.h (optional)
|   |   *   values_from_mesher.h (optional)
|   *   DATA/
|   |   *   Par_file #used to build specfem3d_globe
|   *   starting_model/
|   ...
|
*   specfem3d_globe/ # build directory of SEM excutables (xmesher3D, xspecfem3D)
|   *   bin/
|   *   ...
|
*   events/ # waveform data
|   *   <event_id>/
|   |   *   data/
|   |   |   *   station.txt # fdsnws-station: text output at channel level  
|   |   |   *   CMTSOLUTION # SEM inputs
|   |   |   *   STATIONS # SEM inputs
|   |   |
|   |   *   disp/
|   |   |   *   <net>.<sta>.<loc>.<cha> # sac file of displacement (after instrument correction )
|   |   |
|   *   <event_id>/
|   |   *   ...
|   ...   
|   
*   iterations/ # iteration database
|   *   iteration.00/
|   |   *   control_file.00 # a shell script contains control parameters
|   |   |
|   |   *   model/
|   |   |   *   DATABASES_MPI/
|   |   |   |   *   prco***_reg1_<model_name>_dmodel.bin # model updates
|   |   |   |   *   prco***_reg1_<v??,rho,eta>.bin # updated model files
|   |   |   |
|   |   *   mesh/
|   |   |   *   DATABASES_MPI/proc*_reg1_solver_data.bin
|   |   |   *   DATA/ # necessary data files, Par_file 
|   |   |   
|   |   *   <event_id>/
|   |   |   *   DATABASES_MPI/ 
|   |   |   *   DATA/ 
|   |   |   *   output_syn/ 
|   |   |   |   *   sac/
|   |   |   *   output_kernel/ 
|   |   |   |   *   kernel/ proc*_reg1_cijkl,rho_kernel.bin
|   |   |   *   output_hess/ 
|   |   |   |   *   kernel/ proc*_reg1_cijkl,rho_kernel.bin
|   |   |   *   output_perturb/ # simulation for perturbed model
|   |   |   |   *   sac/
|   |   |   *   misfit/ 
|   |   |   |   *   CMTSOLUTION.reloc # relocated source parameter
|   |   |   |   *   misfit.pkl # misfit measurments (e.g. CC0, CCmax, ...)
|   |   |   *   adj_kernel/ # adjoint source for kernel calculation 
|   |   |   *   adj_hess/ # random adjoint source for hessian calculation 
|   |   |
|   |   *   <event_id>/
|   |   |   *   ...
|   |   ...
|   |   *   kernel_sum/ # summed kernel (with precondition, smoothing, thresholding etc.)
|   |   |   *   # *_dkernel.bin, *_kernel.bin
|   |   |
|   *   iteration.01/
|   |   ...
|   ...


#====== Project setup

1. setup sem_config/
    
    * DATA/Par_file: define model geometry, mesh slices, simulation properties, etc.

    * setup/constants.h.in: fine tuning of mesh parameters
        - regional_moho, moho_stretching, ...

    * starting_model/DATABASES_MPI: give intial model gll files 
        - NOTE: must has the same mesh geometry as would be created from the above configurations. 
        - proc*_reg1_[vpv,vph,vsv,vsh,eta,rho].bin

2. build specfem3d_globe

    * utils/sem_create.sh

    * utils/sem_build.sh

3. compile utility codes in this directory:

    * put the following SEM header files generated after building the specfem3d_globe:

        - specfem3d_globe/setup/{constants.h, precision.h} 
        - specfem3d_globe/OUTPUT_FILES/values_from_mesher.h,

    into include/. You should use these header files for your own application.
 
    * set correct path to netcdf include and lib directories in Makefile
    
    * make -f Makefile clean all

4. prepare the waveform data directory events/ as described in section:{Folder structure}


#====== Work flow for each iteration

1. setup control parameters for the current iteration.

    * copy and modify utils/control_file into project/iterations/contro_file.<iter>

2. run scripts utils/qsub_iteration to submit jobs for the whole iteration

    > mesh.job
    * utils/update_model.sh -> model/ # check model update direction, step_length
    * utils/setup_mesh.sh -> mesh/
    
    > <event_id>.job
    * utils/setup_event.sh, setup_adjoint -> <event_id>/  (for all events)
    
    > kernel.job
    * utils/update_kernel.sh -> kernel/

3. post-process:

    * plot x-sections: model, kernel

    * plot data misfit: station misfit distribution, waveform profile for each earthquake
