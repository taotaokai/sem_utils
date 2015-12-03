#====== Data structure

project/
*   sem_utils/ # utility scripts and excutables to run the inversion
|   
|   
*   sem_config/ # configuration files for SEM (model setup,)
|   *   setup/
|   |   *   constants.h.in
|   *   DATA/
|   |   *   Par_file
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
|   |   |   *   <net>.<sta>.<loc>.<cha> # sac file of displacement
|   |   |
|   *   <event_id>/
|   |   *   ...
|   ...   
|   
*   iterations/ # iteration database
|   *   iteration.00/
|   |   *   control_file.00 # control parameters
|   |   |
|   |   *   model/
|   |   |   *   DATABASES_MPI/ # *_dmodel.bin, *_v??.bin
|   |   |
|   |   *   mesh/
|   |   |   *   DATABASES_MPI/proc*_reg1_solver_data.bin
|   |   |   *   DATA/ # necessary data files, Par_file 
|   |   |
|   |   *   <event_id>/
|   |   |   *   DATABASES_MPI/ 
|   |   |   *   DATA/ 
|   |   |   *   OUTPUT_forward/ 
|   |   |   *   OUTPUT_adjoint/ 
|   |   |   *   misfit/ 
|   |   |   |   *   CMTSOLUTION.reloc # source parameter
|   |   |   |   *   misfit.json # misfit measured results (e.g. CC0, CCmax, ...)
|   |   |   *   adj/
|   |   |
|   |   *   <event_id>/
|   |   |   *   ...
|   |   ...
|   |   *   kernel/
|   |   |   *   DATABASES_MPI/ # *_dkernel.bin, *_kernel.bin
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

    * sem_utils/utils/sem_create.sh

    * sem_utils/utils/sem_build.sh


#====== Work flow for each iteration

1. define control parameters for the current iteration.

2. run scripts sem_utils/utils/qsub_iteration to submit jobs for the whole iteration

    > one job
    * update_model.sh -> model/ # check model update direction, step_length
    * setup_mesh.sh -> mesh/
    
    > jobs of number of events 
    * setup_event.sh, setup_adjoint -> <event_id>/  (for all events)
    
    > one job
    * update_kernel.sh -> kernel/

3. post-process:

    * plot x-sections: model, kernel

    * plot data misfit: station misfit distribution, waveform profile for each earthquake
