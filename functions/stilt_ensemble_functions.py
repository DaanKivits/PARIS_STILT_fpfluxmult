# Daan Kivits, 2024

# This file contains functions used by the FLEXPART flux file creation scripts.
# EDIT 15-09-2023: Added functions to use with the STILT model.

##############################################
########## LOAD NECCESSARY PACKAGES ##########
##############################################
from itertools import groupby
import xarray as xr
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from functions.background_functions import *
from tqdm import tqdm
import logging
from functions.fluxfile_functions_dynamicfluxes import *

def create_obs_sim_dict_ens(flux_dataset, fp_filelist, lats, lons, RDatapath,
                        bgpath, station_lat, station_lon, station_agl, stationname,
                        save_as_pkl=False):
    """ Function to compare the STILT simulation results from simulations with different setups, e.g. a
    sensitivity analysis. Here, the sensitivity analysis is done on the amount of particles released in 
    each of the STILT simulations, with multiple ensemble members to increase statistical significance.
    
    NOTE: Not up to date!
    """
    logging.info('Current station is ' + stationname + ', with lat = ' + str(station_lat) + ', lon = ' + str(station_lon) + ', and agl = ' + str(station_agl))
    logging.info('Starting calculation of flux contributions and pseudo observations.')
    
    # Sort footprint list (required by groupby function)
    fp_filelist.sort()

    # Group subgroup by released number of particles
    fp_filelist_grouped = [list(g) for k, g in groupby(fp_filelist, lambda s: s[-4:-5])]
    print(fp_filelist_grouped)

    # Create empty df to store results
    pseudo_df, mixed_df = (pd.DataFrame() for i in range(2))
    
    for time_specific_list in fp_filelist_grouped:  
        time_specific_multilist = [list(g) for k, g in groupby(time_specific_list, lambda s: s[-9:-6])]
        print(time_specific_multilist)

        for part_specific_list in time_specific_multilist:
            print(part_specific_list)
            
            # Create empty dictionary
            pseudodict, mixeddict = ({} for i in range(2))
            npars = int(part_specific_list[0].split(sep='x')[-2].strip())

            logging.info('Starting up with ' + str(npars) + ' particles.')

            # Loop over all footprint files
            for i in tqdm(range(0, len(part_specific_list))):
                # Get start and end time of footprint
                file = part_specific_list[i]
                with xr.open_dataset(file) as footprint_df:

                    # Get start time and npars of footprint
                    npars = footprint_extract_npars(file)
                    fp_starttime = footprint_extract_time(file)
                    
                    logging.info(f'Current ens_mem_num is {i+1}')
                    
                    # Extract latitude indices from sparse footprint file and insert into xr.DataArray for later indexing
                    lat_indices, lon_indices = get_latlon_indices(footprint_df, lats, lons)
                    hours_into_file = get_time_indices(footprint_df=footprint_df, fp_starttime=fp_starttime,
                                                            fp_filelist = fp_filelist)

                    logging.info(f'This footprint has {len(lat_indices)} values, has {len(np.unique(hours_into_file))} timesteps and goes {len(np.unique(hours_into_file))} hours backwards in time')

                    # Calculate mean background concentration for current footprint
                    bg = calculate_mean_bg(fp_starttime=fp_starttime, RDatapath=RDatapath, bgpath=bgpath,
                                        npars=npars, lat=station_lat, lon=station_lon, agl=station_agl, ens_mem_num=i+1)

                    logging.info(f'Mean background is {bg}')

                    mixed = get_flux_contribution(flux_dataset=flux_dataset, time_index_list=hours_into_file,
                                                lat_index_list=lat_indices, lon_index_list=lon_indices,
                                                footprint_df=footprint_df, infl_varname='Influence')
                    
                    logging.info(f'Mixed flux contribution is {mixed}')

                    # Calculate pseudo observation from bg and flux contribution
                    pseudo_obs = bg + mixed

                    # Save keys and values in summary_dict
                    pseudodict[i+1] = [pseudo_obs]
                    mixeddict[i+1] = [mixed]

            pseudo_subdf = pd.DataFrame(pseudodict, index=[npars])
            pseudo_subdf.index.name = 'Number of particles'
            pseudo_df = pseudo_df.append(pseudo_subdf)

            mixed_subdf = pd.DataFrame(mixeddict, index=[npars])
            mixed_subdf.index.name = 'Number of particles'
            mixed_df = mixed_df.append(mixed_subdf)

            print('Pseudo df: ' + str(pseudo_df))
            print('Mixed df: ' + str(mixed_df))

    logging.info('Finished calculating flux contributions and pseudo observations.')
        
    return pseudo_df, mixed_df
        
def compare_ens_members(filepath, sim_length, fluxdir, fluxfilenamelist, fluxvarnamelist, stationsfile,
        bgfilepath, lats, lons, stationcode):
    
    # Get list of footprint files for station
    sparse_files = sorted(glob.glob(filepath + 'footprint_' + stationcode + '*.nc'))

    # Extract all unique months between start and end time
    mons = footprint_unique_months(sparse_files, sim_length)

    # Only once per station, create list of all CTE-HR flux files
    fluxstring = find_fluxfiles(fluxdir = fluxdir, variablelist_files = fluxfilenamelist, months = mons, fluxtype='CTE-HR')

    # Open all files in fluxstring as xr_mfdataset, and add variables in variablelist
    cte_ds = hard_load_all_fluxes(fluxstring, sim_len=sim_length, variables = fluxvarnamelist)

    # Get 3D station location
    lat,lon,agl = get_3d_station_location(stationsfile, stationcode, colname_lat='lat', colname_lon='lon',
                                            colname_agl='corrected_alt')      
 
    # Loop over all footprints files
    pseudo_df, mixed_df = create_obs_sim_dict_ens(fp_filelist = sparse_files, flux_dataset = cte_ds,
                        lats = lats, lons = lons, RDatapath = filepath, 
                        bgpath = bgfilepath, station_lat = lat, station_lon = lon, 
                        station_agl = agl, stationname = stationcode)
    
    return pseudo_df, mixed_df
    
def plot_ens_members_fromdf(df):
    """ Function to plot the results of the ensemble members from a dataframe in memory. """
    # Create subplots
    fig, axes = plt.subplots(len(df.index), 1, figsize=(15, 5), sharey=True)
    axes[0].set_title('Simulation results for varying amounts of particles')
    
    for i in range(0, len(df.index)):
        npars = df.index[i]
        ax = axes[i]
        df.loc[i].plot(ax=ax, marker='o', label=str(npars) + 'particle simulation')

        ax.set_xlabel('Number of particles')
        ax.set_ylabel('CO2 amount fraction (-)')
        ax.legend()
        ax.grid()
    
    plt.show()

def plot_ens_members_fromcsv(pseudo_csv, mixed_csv, save = True):
    """ Function to visualize the results of STILT simulations with a varying number
    of particles that have been released. The flux contribution is shown next to 
    the total simulated amount fraction. Both of these are read in separately through
    .CSV files. 
    
    """
    # Read in CSV as pd.DataFrame
    pseudo_df = pd.read_csv(pseudo_csv, index_col="Number of particles")
    mixed_df = pd.read_csv(mixed_csv, index_col="Number of particles")

    # create tickslist
    #tickslist = list(range(100, 400, 100))
    tickslist = [100, 200, 300, 400, 500, 600]

    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    ax1 = axes[0]
    ax2 = axes[1]
    ax1.yaxis.get_offset_text().set_fontsize(12)  # Adjust font size if needed

    ax1.set_title('Total simulated mixing ratio \n (pseudo-observations)', pad = 10)
    ax1.set_xlabel('Number of particles')
    ax1.set_ylabel('CO2 amount fraction [ppm]')
    ax1.set_ylim(410,415)
    ax1.set_xticks(tickslist)
    ax1.set_xticklabels(pseudo_df.index)
    ax1.grid()

    ax2.set_title('CTE-HR contribution to total \nCO2 amount fraction', pad = 10)
    ax2.set_ylabel('CTE-HR contribution to total \nCO2 amount fraction [ppm]')
    ax2.set_xlabel('Number of particles')
    ax2.set_ylim(-1,-6)
    ax2.set_xticks(tickslist)
    ax2.set_xticklabels(mixed_df.index)
    ax2.grid()
    ax2.invert_yaxis()

    plt.subplots_adjust(wspace = 0.3)
    
    for i in range(0, len(pseudo_df.index)):
        npars = pseudo_df.index[i]
        ax1.scatter(x=[tickslist[i]] * len(pseudo_df.columns), y=pseudo_df.iloc[i]*1e6, marker='o', c='orange', label=str(npars) + 'particle simulation', s = 10, alpha = 1, zorder=3) 
        ax2.scatter(x=[tickslist[i]] * len(mixed_df.columns), y=mixed_df.iloc[i]*1e6, marker='s', c='red', label=str(npars) + 'particle simulation', s = 10, alpha = 1, zorder = 3)
    
    if (save == True):
        plt.savefig('/projects/0/ctdas/PARIS/Experiments/ens_members_comparison/ens_members_CBW207.png', dpi=300, bbox_inches='tight', overwrite=True)

    plt.show()
    plt.close()