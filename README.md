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



## Folder Structure for an inversion project



```shell
project/
    sem_utils/ # utility scripts and excutables to run the inversion
    |   
    sem_config/ # configuration files for SEM (model setup,)
    |   setup/ # these header files are used to build specfem3d_globe
    |   |   constants.h.in
    |   |   precision.h (optional)
    |   |   values_from_mesher.h (optional)
    |   DATA/
    |   |   Par_file #used to build specfem3d_globe
    |   starting_model/
    |   ...
    |
    specfem3d_globe/ # build directory of SEM excutables (xmesher3D, xspecfem3D)
    |   bin/
    |   ...
    |
    events/ # waveform data
    |   <event_id>/
    |   |   channel.txt # fdsnws-station: text output at channel level  
    |   |   CMTSOLUTION # SEM inputs
    |   |   STATIONS # SEM inputs
    |   |   data.h5
    |   ...
    |   
    stage??.[source|structure]/ # iteration database
    |   iter??/
    |   |   control_file.00 # a shell script contains control parameters
    |   |   |
    |   |   model/
    |   |   |   DATABASES_MPI/
    |   |   |   |   prco***_reg1_<model_name>_dmodel.bin # model updates
    |   |   |   |   prco***_reg1_<v??,rho,eta>.bin # updated model files
    |   |   |   |
    |   |   mesh/
    |   |   |   DATABASES_MPI/proc*_reg1_solver_data.bin
    |   |   |   DATA/ # necessary data files, Par_file 
    |   |   |   
    |   |   events/
    |   |   |   <event_id>/
    |   |   |   |   DATABASES_MPI/ 
    |   |   |   |   DATA/ 
    |   |   |   |   output_syn/ 
    |   |   |   |   |   sac/
    |   |   |   |   output_kernel/ 
    |   |   |   |   |   kernel/ proc*_reg1_cijkl,rho_kernel.bin
    |   |   |   |   output_hess/ 
    |   |   |   |   |   kernel/ proc*_reg1_cijkl,rho_kernel.bin
    |   |   |   |   output_perturb/ # simulation for perturbed model
    |   |   |   |   |   sac/
    |   |   |   |   misfit/ 
    |   |   |   |   |   CMTSOLUTION.reloc # relocated source parameter
    |   |   |   |   |   misfit.h5 # misfit measurments (e.g. CC0, CCmax, ...)
    |   |   |   |   SEM/ # adjoint source for kernel calculation 
    |   |   |   ...
    |   |   |
    |   |   kernel_sum/ # summed kernel (with precondition, smoothing, thresholding etc.)
    |   |   |           # *_dkernel.bin, *_kernel.bin
    |   |   |
    |   iter??/
    |   |   ...
```


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



## Work flow of source inversion (one earthquake)

```mermaid
flowchart LR
    classDef GPU fill:#f96
    
    A[initial CMTSOLUTION] --> e1f

    
      direction LR
      e1f[green]:::GPU --> e1m[misfit] --> e1k[srcfrechet]:::GPU --> e1p[perturb] --> es1f[N x forward]:::GPU --> es1d[diff] --> es1s[search]
    
    
  
    es1s --> o[updated CMTSOLUTION]
```

```pseudocode
# Note: [***_job] is a slurm job

[green_job] forward simulation
[misfit_job] measure misfit
[srcfrechet_job] adjoint simulation, get source gradients (e.g. dchi_dxs, dchi_dmt)

[perturb_job] perturb source parameters along some chosen search directions (e.g. dxs, dmt)
[forward_perturb_job] forward simulation for perturbed models

[search_job]
    get waveform differences for each search direction in model space
    while changes in the optimal step lengths greater than a threshold:
      calculate misfit for a range of step lengths, assuming linear relationship between waveform differences and step length
    get optimal step lengths
```



## Work flow of structural inversion

``` mermaid
flowchart LR
    classDef GPU fill:#f96
    
    A[initial model] --> e1f
    A --> edf
    A --> enf

    subgraph e1 [event 1]
        direction LR
      e1f[forward]:::GPU --> e1m[misfit] --> e1k[kernel]:::GPU
    end
    
    subgraph ed [event ...]
      direction LR
      edf[forward] --> edm[misfit] --> edk[kernel]
    end
    
    subgraph en [event N]
        direction LR
      enf[forward] --> enm[misfit] --> enk[kernel]
    end

    e1k --> ks[sum kernel \n preconditioning \n non-linear CG \n search directions]
    edk --> ks
    enk --> ks
    
    ks --> mp[perturb model]

    subgraph es1 [event 1]
        direction LR
        es1f[N x forward]:::GPU --> es1d[diff] --> es1s[search]
    end
    subgraph esd [event ...]
        direction LR
        esdf[forward] --> esdd[diff] --> esds[search]
    end
     subgraph esn [event N]
        direction LR
        esnf[forward] --> esnd[diff] --> esns[search]
    end
    
    mp --> es1f & esdf & esnf
    
    es1s & esds & esns --> o[optimal step lengths]
```

```pseudocode
# Note: [***_job] is a slurm job

for each event:
     [forward_job] forward simulation
     [misfit_job]  measure misfit
     [kernel_job]  kernel simulation, get model gradients (e.g. cijkl_kernel)

[sum_kernel_job] sum model gradients for all events to get averaged model gradients, also with preconditioning

[perturb_model_job] perturb model parameters along some chosen search directions (e.g. dVp, dVs) based on averaged gradients and also nonliear-cg

for each event:
  [forward_perturb_job] forward simulation for perturbed models

[search_job]
  for each event:
    get waveform differences for each search direction in model space
  while changes in the optimal step lengths greater than a threshold:
    for each event:
      calculate misfit for a range of step lengths, assuming linear relationship between waveform differences and step length
    sum search results for all event, get optimal step lengths
```



