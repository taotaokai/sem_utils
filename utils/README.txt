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
|   |   |   *   DATABASES_MPI/ # link from mesh/
|   |   |   *   DATA/ 
|   |   |   *   green_function/ # forward simulation 
|   |   |   *   source_frechet/ # adjoint simulation for source parameters ()
|   |   |   *   model_frechet/  # adjoint simulation for model parameters
|   |   |   *   misfit/ # data misfit
|   |   |   |   *   CMTSOLUTION.reloc # source parameter
|   |   |   |   *   misfit.json # misfit measured results (e.g. CC0, CCmax, ...)
|   |   |   *   adj/ # adjoint sources, SEM input
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
        - the following parameters must be set properly before compiling. Otherwise 
            re-compilation is needed.

            # mesh geometry
            ANGULAR_WIDTH_XI_IN_DEGREES   = 50.d0      # angular size of a chunk
            ANGULAR_WIDTH_ETA_IN_DEGREES  = 50.d0
            CENTER_LATITUDE_IN_DEGREES    = 38.5d0
            CENTER_LONGITUDE_IN_DEGREES   = 118.0d0
            GAMMA_ROTATION_AZIMUTH        = -10.0d0
            # number of elements at the surface along the two sides of the first chunk
            # (must be multiple of 16 and 8 * multiple of NPROC below)
            NEX_XI                          = 256
            NEX_ETA                         = 256
            # number of MPI processors along the two sides of the first chunk
            NPROC_XI                        = 16
            NPROC_ETA                       = 16
            # parameters describing the Earth model
            OCEANS                          = .true.
            ELLIPTICITY                     = .true.
            TOPOGRAPHY                      = .true.
            GRAVITY                         = .true.
            ROTATION                        = .true.
            ATTENUATION                     = .true.
            # absorbing boundary conditions for a regional simulation
            ABSORBING_CONDITIONS            = .true.
            # to undo attenuation for sensitivity kernel calculations or forward runs with SAVE_FORWARD
            # use one (and only one) of the two flags below. UNDO_ATTENUATION is much better (it is exact)
            # but requires a significant amount of disk space for temporary storage.
            PARTIAL_PHYS_DISPERSION_ONLY    = .false.
            UNDO_ATTENUATION                = .true.
            # mem 
            MEMORY_INSTALLED_PER_CORE_IN_GB = 2.5d0
            PERCENT_OF_MEM_TO_USE_PER_CORE  = 90.d0
            # kernel 
            ANISOTROPIC_KL                  = .true.
            SAVE_TRANSVERSE_KL_ONLY         = .false.
            USE_FULL_TISO_MANTLE            = .true.


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
