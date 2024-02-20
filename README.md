# README file for transporting the CTE-HR fluxes using pre-computed STILT footprints for the PARIS WP6 CO<sub>2</sub> verification games
### Id: README.md, 19-02-2024 D. Kivits $
---
This directory contains all the scripts neccessary to transport PARIS-specific CTE-HR fluxes through the atmosphere using pre-computed STILT influence fields. The current implementation is only limited to forward atmospheric transport, but should form the basis for future iterations of the verification games in which atmospheric inversions are executed. 

The directory is structured as follows:
- the parent directory contains <fp_flux_mult_singlerun.py> (hereafter refered to as 'main script'), and the submit scripts needed to run the main script on the HPC cluster.
- <functions/> contains all the Python functions that are used in the main script. 

This document consists of the following sections:
- Instructions on how to run the fp-flux multiplication pipeline and the command-line argument that the main script takes under '**FP-FLUX MULT PIPELINE WORKFLOW**'.
- A deeper dive into the contents of the fp-flux multiplication pipeline of the main script (and the seperate scripts that compose it) under '**FP-FLUX MULT PIPELINE CONTENTS**'. 
- A more detailed discussion on two options regarding the loading the PARIS-specific CTE-HR fluxes into memory under '**DYNAMIC OR STATIC FLUX LOADING**'.
- Some (staticly explored) optional functionality under '**OPTIONAL EXTRA FUNCTIONALITY**'.
---

## FP-FLUX MULTIPLICATION PIPELINE WORKFLOW
As of the current implementation, the fp-flux multiplication pipeline is executed in the following way: the main script - with its various command line input arguments - is submitted to the HPC cluster using the <PARIS_submit_perturbations.sh script>, that subsequently calls the <fp_flux_mult_submit.sh> script. The former script loops over the different perturbation experiment that need to be executed, and submits the latter script that sets the right directories for each of the perturbation experiments and executes the main script for each of the stations individually. The functions that are used in the main script are located in the <functions/> directory, the functionality of which will be explained in more detail in the section '**FP-FLUX MULT PIPELINE CONTENTS**'.

The <fp_flux_mult_submit.sh> script takes the following command-line arguments, which are given to the script by running the <PARIS_submit_perturbations.sh> script:
- 'sum': Controls whether the script should sum the fluxvariables or not. TRUE: only produce a "mixed" molefraction field. FALSE: transport all fluxvariables individually and store them in ObsPack (bool)
- 'experiment': Controls for which perturbation experiment the fp-flux multiplication script is ran, and sets the input and output directories accordingly (str)

The main script that is ran in the <fp_flux_mult_submit.sh> script takes the following command-line arguments:
- --station <stationcode>: Code of station to run the script for (str)
- --stationsfile <stationsfile>: Path to the (.csv) file that contains all station metadata (str)
- --fluxdir <fluxdir>: Directory where flux files are stored (str)
- --fpdir <fpdir>: Directory where footprint files are stored (str)
- --bgdir <bgdir>: Directory where (TM3) background files are stored (str)
- --outdir <outdir>: Directory where output files are stored (str)
- --stiltdir <stiltdir>: Directory where STILT is executed (str)
- --obspack_path <obspack_path>: Path to the ObsPack collection zip file (str)
- --non_obspack_path <non_obspack_path>: Path to any other extra observation files that are not included in the European Obspack collection zip file (str)
- --start_date <start_date>: Date from when to subset the ObsPacks (str)
- --end_date <end_date>: Date up to where to subset the ObsPacks (str)
- --lat_ll <lat_ll>: Latitude of ll corner to define the grid extent (float)
- --lat_ur <lat_ur>: Latitude of ur corner to define the grid extent (float)
- --lon_ll <lon_ll>: Longitude of ll corner to define the grid extent (float)
- --lon_ur <lon_ur>: Longitude of ur corner to define the grid extent (float)
- --lat_step <lat_step>: Latitude cell size step (float)
- --lon_step <lon_step>: Longitude cell size step (float)
- --sim_len <sim_len>: Simulation length in hours (int)
- --perturbation <perturbation>: Controls for which perturbation experiment the fp-flux multiplication script is ran (str)
- --verbose: Controls whether the script should do logging (bool)
- --sum-variables: Controls whether the script should sum the fluxvariables or not. TRUE: only produce a "mixed" molefraction field. FALSE: transport all fluxvariables individually and store them in ObsPack (bool)

## FP-FLUX MULTIPLICATION PIPELINE CONTENTS
The functionality of the footprint-flux multiplication pipeline given in the main script can be broken down into the following steps:

**Step 1**: the main script reads in the user-defined command-line arguments; defines fluxtype-specific fluxsector names; defines the input grid of the fluxes; defines the list of stations and the path to their respective input obspacks; and runs the 'main' function described in the <functions/fluxfile_functions.py> script. 

**Step 2**: the main function is ran, mainly consists of the following functions:
- **get_3d_station_location**(), which extracts the 3D location of the station from the station metadata file
- **find_fluxfiles**(), which finds all neccessary fluxfiles. This function is executed differently for each user-defined fluxtype, as each fluxtype has its own unique file structure, and therefore requires a different approach to read in the flux files.
- **filter_files_by_date**(), which filters the list of footprint files by date range
- functions that load the fluxfiles into memory, either as a static load into memory or a full lazy-load. T
  - **hard_load_all_fluxes**(), which opens all fluxfiles as xr_mfdataset and appends these fluxes to a dictionairy. Only the variables that are given in the variablelist_vars and 'time' are loaded into this dictionairy. The latter is used for indexing the mole fractions in the ObsPack later.
  - **lazy_load_all_fluxes**(), which opens all fluxfiles as xr_mfdataset. Only the variables that are given in the variablelist_vars are retained in the dataset.
- **create_obs_sim_dict**(), is the most important function of the pipeline. It loops over all footprints files, checks for the presence of both observations in the ObsPack files and footprints in the footprint directory, and creates an 'intermediate' dictionary that contains the modelled timestep as a key and the observed and simulated mole fractions as values. The <create_obs_sim_dict()> contains the following subfunctions:
  - **calculate_mean_bg**(), which calculates the background atmospheric CO<sub>2</sub> concentration at the recorded particle end locations for each set of particles, and then takes the mean of these concentrations. This mean concentration is then used as the background concentration for the footprint.
  - **check_time_in_obspack**(), which checks whether an observation is present in the ObsPack file and skips the current iteration if it is not.
  - **get_latlon_indices**(), convert the sparse indices of the footprint to dense indices that correspond to a given horizontal grid, as defined in the <coordinate_list()> function. These indices are later used in the get_flux_contribution() function. 
  - **find_time_indices**(), converts the time indices in each footprint to the corresponding time indices in the flux files. These indices are later used in the <get_flux_contribution()> function.
  - **get_flux_contribution**(), multiplies the fluxes with the footprints in a non-ortagonal manner
  - **create_intermediate_dict**(), which creates the intermediate dictionary that contains the modelled timestep as a key and the observed and simulated mole fractions as values. If the option 'sum-variables' fed to <fp_flux_mult_singlerun.py> was set to 'True', only the 'mixed' simulated contribution to the atmospheric CO<sub>2</sub> signal is stored, otherwise all simulated sector-specific contributions are stored individually. The intermediate dictionary is later used to write the results to an ObsPack file. 
- **write_dict_to_ObsPack** and **write_dict_to_ObsPack_altsites()**, which write the results of the footprint-flux multiplication to an ObsPack file. The former function is used for the main stations, the latter for the alternative stations (i.e. stations that are not included in the European ObsPack collection). The latter function acts as a starting point for the use of other non-ObsPack structured observation files that need to be included into the fp-flux-multiplcation pipeline. The ObsPacks are subsetted to the user-defined date range and resampled to an hourly time resolution to match the time resolution of the pseudo-observations.

## DYNAMIC OR STATIC FLUX LOADING
As touched on in the previous section, the main script has two options for loading the fluxes into memory: the static flux loading and dynamic flux loading methods. As currently implemented, the dynamic flux loading method is implemented by default. The dynamic flux loading method does a lazy load on all fluxes of the simulation period in the function <lazy_load_all_fluxes()>, creates a dictionairy that loads the fluxes for the first simulation timestep into memory in the function <load_first_flux_timestep()>, and then updates this dictionairy by moving this dictionairy in time in the function <upate_fluxdict()>. The static flux loading method creates a dictionairy and loads all fluxes into memory on the first timestep in the <hard_load_all_fluxes()> function.

The dynamic flux load was developed because we needed a more memory-efficient of handling the fluxes for PARIS. The static flux load seems to load the fluxes into memory in a more time-efficient manner, but is more memory-intensive. For PARIS, loading the fluxing in a static way would require a lot of memory, as the fluxes are then loaded for the entire simulation period for each station individually. To still use the static flux loading method, please source the <fluxfile_functions_partialfluxes.py> script in the main script instead of <fluxfile_functions_dynamicfluxes.py>.

## LEGACY FUNCTIONALITY
For a sensitivity study on the effects of the number of particles on the CTE-HR-STILT simulated mole fractions, I developed some extra (legacy) code that can be found in <functions/stilt_ensemble_functions.py>. The functions in that script were used to compare multiple STILT simulations with various number of particles, and plot the results side-by-side.