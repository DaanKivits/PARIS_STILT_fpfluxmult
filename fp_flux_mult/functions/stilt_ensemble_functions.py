# Daan Kivits, 2024

# This file contains functions used by the FLEXPART flux file creation scripts.
# EDIT 15-09-2023: Added functions to use with the STILT model.

##############################################
########## LOAD NECCESSARY PACKAGES ##########
##############################################
from dateutil.relativedelta import relativedelta
from itertools import groupby
import csv
import shutil
from datetime import datetime, timedelta, time, date
import os
import xarray as xr
import glob
import time
import matplotlib.pyplot as plt
import netCDF4 as nc
from netCDF4 import stringtochar
import numpy as np
from scipy import stats
import pandas as pd
from functions.background_functions import *
import pickle
from tqdm import tqdm
import logging
import matplotlib.ticker as mtick
import dask
import dask.array as da

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
                    hours_into_file = find_footprint_flux_timeindex(footprint_df=footprint_df, fp_starttime=fp_starttime,
                                                            fp_filelist = fp_filelist)

                    logging.info(f'This footprint has {len(lat_indices)} values, has {len(np.unique(hours_into_file))} timesteps and goes {len(np.unique(hours_into_file))} hours backwards in time')

                    # Calculate mean background concentration for current footprint
                    bg = calculate_mean_bg(fp_starttime=fp_starttime, RDatapath=RDatapath, bgpath=bgpath,
                                        npars=npars, lat=station_lat, lon=station_lon, agl=station_agl, ens_mem_num=i+1)

                    logging.info(f'Mean background is {bg}')

                    mixed = get_flux_contribution_np(flux_dataset=flux_dataset, time_index_list=hours_into_file,
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
    cte_ds = open_fluxfiles(fluxstring, sim_len=sim_length, variables = fluxvarnamelist)

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

def date_range(start, end, intv):
    """ Function to split a certain date range according to a certain interval,
    and yield the start and end date of each interval. To put results in a list, 
    do: list(date_range(start, end, intv)).
    """
    start = datetime.strptime(start,"%Y-%m-%d")
    end = datetime.strptime(end,"%Y-%m-%d")
    diff = (end  - start ) / intv
    for i in range(intv):
        yield (start + diff * i).strftime("%Y-%m-%d")
    yield end.strftime("%Y-%m-%d")


def filter_files_by_date(file_list, start_date, end_date):
    """ Function to filter a list of files by a certain date range. 
    """
    filtered_files = []

    date_format = "%Y-%m-%d"  # adjust the format based on your actual filenames

    start_datetime = datetime.strptime(start_date, date_format)
    end_datetime = datetime.strptime(end_date, date_format)

    for file_name in file_list:
        #date_str = file_name.split("_")[3][:13]  # adjust the index based on your actual filenames
        date_str = os.path.basename(file_name).split("_")[2][:13]  # adjust the index based on your actual filenames
        file_datetime = datetime.strptime(date_str, "%Yx%mx%dx%H")

        if start_datetime <= file_datetime <= end_datetime:
            filtered_files.append(file_name)

    return filtered_files


def parse_datetime(datetime_dataarray):
    """ Function to parse a datetime string from a netcdf file to a datetime object. 
    NOTE: only the first 19 characters are parsed, because the netcdf files often contain
    a string with the format 'YYYY-MM-DDTHH:MM:SS' and some extra unwanted characters.
    
    """
    #return pd.to_datetime(datetime_str[0:19].decode('utf-8'), '%Y-%m-%dT%H:%M:%S')
    datetime_str = datetime_dataarray.astype(str).str[0:19]
    return pd.to_datetime(datetime_str, format='%Y-%m-%dT%H:%M:%S')


def coordinate_list(lat_ll, lat_ur, lon_ll, lon_ur, lat_step, lon_step):
    """ Function to create a list of lons and lats given a 
    range and stepsize.
    """
    #lats = np.array(np.round(np.arange(lat_ll + 0.5 * lat_step, lat_ur + 0.5 * lat_step, lat_step,dtype=np.float32), 2))
    #lons = np.array(np.round(np.arange(lon_ll + 0.5 * lon_step, lon_ur + 0.5 * lon_step, lon_step,dtype=np.float32), 2))
    lats = list(np.round(np.arange(lat_ll + 0.5 * lat_step, lat_ur + 0.5 * lat_step, lat_step,dtype=np.float32), 2))
    lons = list(np.round(np.arange(lon_ll + 0.5 * lon_step, lon_ur + 0.5 * lon_step, lon_step,dtype=np.float32), 2))
    return lats, lons


def footprint_unique_months(fp_filelist, simulation_len):
    """ Function to extract the unique months of the footprint files, to later use to subset CTE-HR flux files (time format YYYYmm). 
    The function accounts for extra possible term in filename, that describes the ensemble run number.
    
    """
    # Define time range
    if len(fp_filelist[0].split(sep='x')) == 8:
        timestr_start = fp_filelist[0][-42:-29]
        timestr_end = fp_filelist[-1][-42:-29]
    elif len(fp_filelist[0].split(sep='x')) == 9:
        timestr_start = fp_filelist[0][-45:-32]
        timestr_end = fp_filelist[-1][-45:-32]
    
    fp_range_start = datetime.strptime(
        timestr_start, '%Yx%mx%dx%H') - timedelta(hours=simulation_len)
    fp_range_end = datetime.strptime(timestr_end, '%Yx%mx%dx%H')

    # Extract all unique months between start and end time
    mons = pd.date_range(fp_range_start, fp_range_end,
                         freq='D').strftime("%Y%m").unique().tolist()
    return mons


def footprint_hours(fp_filelist, simulation_len, shift_forward=False):
    """ Function to extract the hours of a list of footprint files.
        The function accounts for extra possible term in filename, that describes the ensemble run number.
    """
    # Define time string
    if len(fp_filelist[0].split(sep='x')) == 8:
        timestr_start = fp_filelist[0][-42:-29]
        timestr_end = fp_filelist[-1][-42:-29]
    elif len(fp_filelist[0].split(sep='x')) == 9:
        timestr_start = fp_filelist[0][-45:-32]
        timestr_end = fp_filelist[-1][-45:-32]
    
    # Define time range
    if shift_forward:
        fp_range_start = datetime.strptime(
            timestr_start, '%Yx%mx%dx%H') - timedelta(hours=simulation_len)
    else:
        fp_range_start = datetime.strptime(timestr_start, '%Yx%mx%dx%H')
    
    fp_range_end = datetime.strptime(timestr_end, '%Yx%mx%dx%H')

    # Define list of times
    times = pd.date_range(start=fp_range_start, end=fp_range_end, freq='H')

    # Drop times that don't have the same Hour of Day (HOD) as the footprint files
    for time in times:
        if time.hour not in range(fp_range_start.hour, (fp_range_end + timedelta(hours=1)).hour):
            times = times.drop(time)

    return times.tolist()


def footprint_list_extract_starttime(fp_filelist):
    """ Function to extract the start time of A LIST OF footprint files.
    """
    # Define time string
    if len(fp_filelist[0].split(sep='x')) == 8:
        timestr_start = fp_filelist[0][-42:-29]
    elif len(fp_filelist[0].split(sep='x')) == 9:
        timestr_start = fp_filelist[0][-45:-32]

    fp_datetime_start = datetime.strptime(timestr_start, '%Yx%mx%dx%H')

    return fp_datetime_start


def footprint_list_extract_endtime(fp_filelist):
    """ Function to extract the end time of A LIST OF footprint files.
    """
    # Define time string
    if len(fp_filelist[0].split(sep='x')) == 8:
        timestr_end = fp_filelist[-1][-42:-29]
    elif len(fp_filelist[0].split(sep='x')) == 9:
        timestr_end = fp_filelist[-1][-45:-32]
    
    fp_datetime_end = datetime.strptime(timestr_end, '%Yx%mx%dx%H')

    return fp_datetime_end


def footprint_extract_starttime(fp_file):
    """ Function to extract the start time of A SPECIFIC footprint file.
    """
    # Define time string
    if len(fp_file.split(sep='x')) == 8:
        timestr = fp_file[-42:-29]
    elif len(fp_file.split(sep='x')) == 9:
        timestr = fp_file[-45:-32]
    #fp_datetime = pd.to_datetime(timestr, format = '%Yx%mx%dx%H')

    fp_starttime = datetime.strptime(timestr, '%Yx%mx%dx%H')
    
    return fp_starttime


def footprint_extract_endtime(fp_file, simulation_len=240):
    """ Function to extract the time range of A SPECIFIC footprint file.
    """
    # Define time string
    if len(fp_file.split(sep='x')) == 8:
        timestr = fp_file[-42:-29]
    elif len(fp_file.split(sep='x')) == 9:
        timestr = fp_file[-45:-32]
    #fp_datetime = pd.to_datetime(timestr, format = '%Yx%mx%dx%H')
    
    fp_starttime = datetime.strptime(timestr, '%Yx%mx%dx%H')
    fp_endtime = fp_starttime - timedelta(hours=simulation_len) 

    return fp_endtime


def footprint_extract_npars(fp_file):
    """ Function to extract the number of particles in the footprint files.
    """
    # Define time string
    if len(fp_file.split(sep='x')) == 8:
        npars = int(fp_file.split(sep='x')[-1].split(sep='.')[0].strip())
    elif len(fp_file.split(sep='x')) == 9:
        npars = int(fp_file.split(sep='x')[-2].strip())

    return npars


def footprint_extract_time(fp_file):
    """ Function to extract the start time of A SPECIFIC footprint file.
    """
    # Define time string
    if len(fp_file.split(sep='x')) == 8:
        timestr = fp_file[-42:-29]
    elif len(fp_file.split(sep='x')) == 9:
        timestr = fp_file[-45:-32]
    #fp_datetime = pd.to_datetime(timestr, format = '%Yx%mx%dx%H')
    fp_datetime = datetime.strptime(timestr, '%Yx%mx%dx%H')

    return fp_datetime


def find_fluxfiles(fluxdir, fluxtype, perturbationcode=None, variablelist_files=None, months=None):
    """ Function to find all fluxfiles in a fluxdir, given a list of variables and months.
    Choose between multiple fluxtypes, currently still limited to "CTEHR" and "PARIS" (03-11-2023)."""
    fluxstring = []
    
    if fluxtype == "PARIS":
        if perturbationcode == 'BASE_SS':
            fluxstring += sorted(glob.glob(fluxdir + '/BASE/paris_ctehr_*.nc'))
        else:
            fluxstring += sorted(glob.glob(fluxdir + '/' + perturbationcode + '/paris_ctehr_*.nc'))
    
    elif fluxtype == "CTEHR":
        for var in variablelist_files:
            for mon in months:
                fluxstring += sorted(glob.glob(fluxdir + var + '.' + mon + '.nc'))

    elif fluxtype == "Ruben":
        for var in variablelist_files:
            fluxstring += sorted(glob.glob(fluxdir + var + '.2022.nc'))
    else:
        print('Input fluxtype not recognized, please select one that is defined in the find_fluxfiles() function or define one yourself.')

    return fluxstring


def open_fluxfiles(fluxfile_list, sim_len=240, mp_start_date=None, mp_end_date=None, fluxtype=None, variables=None):
    """ Function to open all files in fluxstring as xr_mfdataset, and sum 
    all the variables in variablelist. Choose between multiple fluxtypes, currently still limited to 
    "CTEHR" and "PARIS" (03-11-2023)."""
    fluxdict = {}

    # Add time variable to variables to put the values in the dictionairy later
    variables = variables + ['time']

    if fluxtype == "PARIS":        
        with xr.open_mfdataset(fluxfile_list) as ds:
            
            logging.info('Reading fluxes into memory!')

            if mp_start_date != None and mp_end_date != None:
                
                # Make sure the fluxfiles contain information for the backward simulated hours too!
                start_date_corrected = datetime.strptime(mp_start_date, "%Y-%m-%d") - timedelta(hours=sim_len)
                ds = ds.sel(time=slice(start_date_corrected, mp_end_date))
            
            # Fill dictionairy with all variables in the list of variables
            for var in variables:
                fluxdict[var] = ds[var][:].compute().values

            #return subds[variables].to_dict(data='array')
            #return ds[variables]
            return fluxdict

    elif fluxtype == "CTEHR" or fluxtype == "Ruben":
        with xr.open_mfdataset(fluxfile_list, combine='by_coords') as ds:
                if mp_start_date != None and mp_end_date != None:
                
                    # Make sure the fluxfiles contain information for the backward simulated hours too!
                    start_date_corrected = datetime.strptime(mp_start_date, "%Y-%m-%d") - timedelta(hours=sim_len)
                    ds = ds.sel(time=slice(start_date_corrected, mp_end_date))

                # Fill dictionairy with all variables in variables
                for var in variables:
                    fluxdict[var] = ds[var][:].compute().values
            
                return fluxdict
        
    else:
        print('Input fluxtype not recognized, please select one that is defined in the find_fluxfiles() function or define one yourself.')



def get_3d_station_location(stationsfile, station=str(), colname_lon=str(), colname_lat=str(), colname_agl=str()):
    """ Function to get the lat, lon and agl of a station given the 
    station name and the stationsfile. """
    lat = stationsfile[stationsfile['code'] == station][colname_lat].values[0]
    lon = stationsfile[stationsfile['code'] == station][colname_lon].values[0]
    agl = stationsfile[stationsfile['code'] == station][colname_agl].values[0]

    if (np.isnan(agl)):
        agl = stationsfile[stationsfile['code'] == station]['alt'].values[0]

    return lat, lon, agl


def get_latlon_indices(footprint_df, lats, lons):
    """ Function to get the lat and lon indices of the footprint file,
    given the lat and lon lists of the domain. """
    #lat_indices = [lats.index(i) for i in list(
    #    np.round(footprint_df.variables['Latitude'][:], 2))]
    #lon_indices = [lons.index(i) for i in list(
    #    np.round(footprint_df.variables['Longitude'][:], 2))]

    #lat_indices = xr.DataArray(lat_indices, dims=['pos'])
    #lon_indices = xr.DataArray(lon_indices, dims=['pos'])

    # Find indices of footprint lats and lons in domain-specific lats and lons list + correct R indexing to Python indexing (hence the - 1)
    lat_indices = np.searchsorted(lats, footprint_df.variables['Latitude'][:])
    lon_indices = np.searchsorted(lons, footprint_df.variables['Longitude'][:])

    #lat_indices = np.where(footprint_df.variables['Latitude'].values[:, None] == lats[None, :])[1]
    #lon_indices = np.where(footprint_df.variables['Longitude'].values[:, None] == lons[None, :])[1]

    return lat_indices, lon_indices


def find_footprint_flux_timeindex(flux_ds, footprint_df, fp_starttime):
    """ Function to find the index of the current timestamp in a given flux file. """
    # First calculate a list of hours since the start of the input footprint file
    hours_relative_to_start_fp = np.array(-1*((footprint_df.variables['Time'][:]).astype(int)))
    
    # Calculate time difference between flux dataset and input footprint file
    fp_flux_timediff = flux_ds['time'] - fp_starttime

    # Extract index of the smallest time difference (which is the timestamp 
    # in the flux dataset closest to the start time of the footprint file)
    fp_flux_startindex = np.abs(fp_flux_timediff).argmin().item()

    # Calculate a list of indices of the footprint file relative to the start 
    # of the flux dataset
    fp_flux_diff_indexlist = fp_flux_startindex + hours_relative_to_start_fp
    
    return fp_flux_diff_indexlist

def get_flux_contribution_np(flux_ds, time_index_list, lat_index_list, lon_index_list, footprint_df, infl_varname, variables, fluxtype = None,
                         basepath = None, sum_vars=None, changed_var = None, perturbationcode=None, station=None, fp_starttime = None):
    """ Function to get the flux contribution of a footprint file, with list indexing. 
    Choose between either a mixed or sector-specific flux contribution.
    """
    # Check if observation is in ObsPack or not
    obspack_filepath = glob.glob(basepath + '*' + station.lower() + '*.nc')[0]
    obspack_checktime = check_time_in_obspack(obspack_filepath=obspack_filepath, fp_starttime=fp_starttime)
    
    logging.info(max(lat_index_list))
    logging.info(max(lon_index_list))
    logging.info(min(lat_index_list))
    logging.info(min(lon_index_list))

    if obspack_checktime:
        if sum_vars==True:
            if fluxtype == "PARIS":
                if ((perturbationcode != 'BASE') and (os.path.exists(basepath))):         
                        base_cont_ss = open_obspack(obspack_filepath=obspack_filepath, variables=variables, selection_time=fp_starttime)
                        base_cont = sum_variables(base_cont_ss, listofvars=variables)
                        cont = (flux_ds[changed_var][time_index_list, lat_index_list,
                                                    lon_index_list] * footprint_df.variables[infl_varname][:]).sum()
                        mixed = base_cont + cont
                        return mixed.item()

                elif perturbationcode == 'BASE':
                    base_cont_ss = open_obspack(obspack_filepath=obspack_filepath, variables=variables, selection_time=fp_starttime)
                    base_cont = sum_variables(base_cont_ss, listofvars=variables)
                    
                    return base_cont
            
                elif not os.path.exists(basepath):
                    import inspect
                    logging.error('"BASE" directory not under the parent directory ' + basepath + '! Please put it there or fix the path in the ' + inspect.currentframe().f_code.co_name + ' function.')
                    raise FileNotFoundError
            
            else:
                mixed = sum([(flux_ds[v][time_index_list, lat_index_list, 
                                                        lon_index_list] * footprint_df[infl_varname][:]).sum() 
                                                        for v in variables])
                return mixed.item()
            
        else:
            mixed_list = [(flux_ds[v][time_index_list, lat_index_list, 
                                                        lon_index_list] * footprint_df[infl_varname][:].values).sum() 
                                                        for v in variables]
            return mixed_list

    else:
        return None
    

def create_intermediate_dict(dict_obj, cont, background, key, sum_vars=None, perturbationcode=None):   
    if ((sum_vars==True) & (type(cont) != list)):
        # If sum_vars is True, cont itself is the total contribution
        logging.info(f'Mixed flux contribution is {cont}')           
        
        pseudo_obs = background + cont
        values = [background, cont, pseudo_obs]
    else:
        # If sum_vars is False, cont is a list of values, so we have to sum it to get the total contribution
        contribution_mixed = sum(cont)
        logging.info(f'Mixed flux contribution is {contribution_mixed}')

        pseudo_obs = background + contribution_mixed
        values = cont + [background, contribution_mixed, pseudo_obs]

    # Save keys and values in summary_dict
    dict_obj[key] = values

    return dict_obj
