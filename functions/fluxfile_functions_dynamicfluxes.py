# Daan Kivits, 2023

# This file contains functions used by the FLEXPART flux file creation scripts.
# EDIT 15-09-2023: Added functions to use with the STILT model.

##############################################
########## LOAD NECCESSARY PACKAGES ##########
##############################################
from datetime import datetime, timedelta
import os
import xarray as xr
import glob
import numpy as np
import pandas as pd
from functions.background_functions import *
import pickle
import logging

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
    """ Function to find all fluxfiles in a fluxdir, the method of which is dependent on the fluxtype.
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


def hard_load_all_fluxes(fluxfile_list, sim_len=240, mp_start_date=None, mp_end_date=None, fluxtype=None, variables=None):
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



def lazy_load_all_fluxes(fluxfile_list, sim_len=240, start_date=None, end_date=None, fluxtype=None, variablelist_vars=None):
    """ Function to lazy load all fluxes into an xarray Dataset. These fluxes will be used to update the fluxes for each timestep later.
    Choose between multiple fluxtypes, currently still limited to "CTEHR" and "PARIS" (03-11-2023)."""  
    
    logging.info('Lazy loading fluxes into memory ...')

    # Add time variable to variables to put the values in the dictionairy later
    variablelist_vars = variablelist_vars + ['time']

    if fluxtype == "PARIS" or fluxtype == "Ruben":        
        with xr.open_mfdataset(fluxfile_list) as ds:
            if start_date != None and end_date != None:
               
                # Make sure the fluxfiles contain information for the backward simulated hours too!
                start_date_corrected = datetime.strptime(start_date, "%Y-%m-%d") - timedelta(hours=sim_len)
                ds = ds.sel(time=slice(start_date_corrected, end_date))

            return ds[variablelist_vars]

    elif fluxtype == "CTEHR":
        with xr.open_mfdataset(fluxfile_list, combine='by_coords') as ds:
            
            return ds[variablelist_vars]

    else:
        print('Input fluxtype not recognized, please select one that is defined in the find_fluxfiles() function or define one yourself.')



def load_first_flux_timestep(flux_ds, fp_starttime, fp_endtime, variablelist_vars=None):
    """ Function to load the fluxes corresponding to the first footprint file in the list into memory,
    which can be updated whilst looping over the other footprints using the update_fluxfiles().
    """
    
    # Add time variable to variables to put the values in the dictionairy later
    variablelist_vars = variablelist_vars + ['time']

    # Cut the flux_ds to the time range of the first footprint file
    flux_ds = flux_ds.sel(time=slice(fp_endtime, fp_starttime))
    
    # Load the fluxes into memory in the form of a dictionary
    fluxdict = {}
    for var in variablelist_vars:
        fluxdict[var] = flux_ds[var][:].values
    
    # Sort dictionairy on dates for later modifications to dict based on order
    fluxdict = dict(sorted(fluxdict.items()))    

    return fluxdict


def update_fluxdict(flux_ds, fluxdict, cur_fp_starttime=None, last_fp_starttime=None):
    """ Function to update a flux array by shifting it in time a given amount of hours (timeshift variable).
    """
    timedif = cur_fp_starttime - last_fp_starttime

    # First check if there's a time gap or the files are consecutive
    if cur_fp_starttime == last_fp_starttime:
        return(fluxdict)
    
    if timedif == np.timedelta64(1, 'h'):
        new_fluxds = flux_ds.sel(time=cur_fp_starttime)
        for var in fluxdict.keys():
            fluxdict[var] = fluxdict[var][1:]
            
            if var != 'time':
                fluxdict[var] = np.vstack([fluxdict[var], new_fluxds[var].values[np.newaxis,...]])
            else:
                fluxdict[var] = np.concatenate([fluxdict[var], new_fluxds[var].values[np.newaxis,...]])

    if np.timedelta64(1, 'h') < timedif < np.timedelta64(240, 'h'): 
        skipped_times = pd.date_range(last_fp_starttime + np.timedelta64(1, 'h'), cur_fp_starttime, freq='1H')
        new_fluxds = flux_ds.sel(time=(skipped_times))
        for var in fluxdict.keys():
            fluxdict[var] = fluxdict[var][len(skipped_times):]

            if var != 'time':
                fluxdict[var] = np.vstack([fluxdict[var], new_fluxds[var].values])
            else:
                fluxdict[var] = np.concatenate([fluxdict[var], new_fluxds[var].values])

    if timedif > np.timedelta64(240, 'h'):
        timerange = pd.date_range(cur_fp_starttime - np.timedelta64(239, 'h'), cur_fp_starttime, freq='1H')
        new_fluxds = flux_ds.sel(time=(timerange))
        for var in fluxdict.keys():
            fluxdict[var] = fluxdict[var][len(skipped_times):]

            if var != 'time':
                fluxdict[var] = np.vstack([fluxdict[var], new_fluxds[var].values])
            else:
                fluxdict[var] = np.concatenate([fluxdict[var], new_fluxds[var].values])

    return(fluxdict)
    
    

def check_time_in_obspack(obspack_ds, fp_starttime):
    """ Function to check if a certain time is present in an ObsPack file. """
    if fp_starttime not in obspack_ds['time']:
        return False
    else:
        return True


def open_obspack(obspack_filepath, variables:list, selection_time = None):
    """ Function to open an ObsPack file and extract the variables of interest at a certain time. 
    The function returns a single value that is stored in the obspack at this time! 
    WARNING: Returns array if sum = True, returns dictionairy if sum = False! 
    """
    obspack_dict = {}
    with xr.open_dataset(obspack_filepath) as ds:

        # Select time rangem either a list or a single Timestamp
        if selection_time != None: 
            if type(selection_time) != list:
                ds = ds.sel(time=selection_time)

                # Fill dictionairy with all variables in variables
                for var in variables:
                    obspack_dict[var] = ds[var].compute().values
            
            else:
                ds = ds.sel(time=slice(selection_time[0], selection_time[1]))
                
                # Fill dictionairy with all variables in variables
                for var in variables:
                    obspack_dict[var] = ds[var][:].compute().values

        else:
            raise ValueError('No selection time given! Please give a selection time to the open_obspack() function.')
    
        return obspack_dict


def sum_variables(dset, listofvars):
    """ Function to sum all variables in a dictionairy from a list of given variable names. This works for
    dictionairies that are created by the hard_load_all_fluxes() function, and hence originate from an xarray.Dataset() object."""
    data = 0
    for var in listofvars:
        if type(dset[var]) != xr.DataArray:
            data = data + dset[var]
        else:
            data = data + dset[var][:]

    return data


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


def get_time_indices(flux_ds, footprint_df, fp_starttime):
    """ Function to converts the time indices in a footprint to the corresponding time indices in a given flux file. """
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

def get_flux_contribution(flux_ds, time_index_list, lat_index_list, lon_index_list, footprint_df, infl_varname, variables, ss_obspack=None, 
                             fluxtype = None, sum_vars=None, changed_var = None, perturbationcode=None, fp_starttime = None):
    """ Function to get the flux contribution of a footprint file, with list indexing. 
    Choose between either a mixed or sector-specific flux contribution.
    """

    if sum_vars==True:
        if fluxtype == "PARIS":
            if ((perturbationcode != 'BASE') and (os.path.exists(ss_obspack))):         
                    base_cont_ss = open_obspack(obspack_filepath=ss_obspack, variables=variables, selection_time=fp_starttime)
                    base_cont = sum_variables(base_cont_ss, listofvars=variables)
                    cont = (flux_ds[changed_var][time_index_list, lat_index_list,
                                                lon_index_list] * footprint_df.variables[infl_varname][:]).sum()
                    mixed = base_cont + cont
                    return mixed.item()

            elif perturbationcode == 'BASE':
                base_cont_ss = open_obspack(obspack_filepath=ss_obspack, variables=variables, selection_time=fp_starttime)
                base_cont = sum_variables(base_cont_ss, listofvars=variables)
                
                return base_cont
        
            elif not os.path.exists(ss_obspack):
                import inspect
                logging.error('"BASE_SS" directory not under the parent directory or ' + ss_obspack + ' not in parent directory! Please put it there or fix the path in the ' + inspect.currentframe().f_code.co_name + ' function.')
                raise FileNotFoundError
        
        else:
            mixed = sum([(flux_ds[v][time_index_list, lat_index_list, 
                                                    lon_index_list] * footprint_df[infl_varname][:].values).sum() 
                                                    for v in variables])
            return mixed.item()
        
    else:
        mixed_list = [(flux_ds[v][time_index_list, lat_index_list, 
                                                    lon_index_list] * footprint_df[infl_varname][:].values).sum() 
                                                    for v in variables]
        return mixed_list


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


def create_obs_sim_dict(flux_ds, fp_filelist, lats, lons, RDatapath, sim_len, bgpath, station_lat, station_lon, station_agl, stationname, fluxtype, variables,
                        stilt_rundir, obspack_list, non_obspack_sites=None, start_date=None, end_date=None, outdir = None, sum_vars = None, perturbationcode = None, 
                        changed_var = None, save_as_pkl=False):
    """ Function to create a dictionary with the background, mixed and pseudo observation. 
    """

    logging.info('Current station is ' + stationname + ', with lat = ' + str(station_lat) + ', lon = ' + str(station_lon) + ', and agl = ' + str(station_agl))
    logging.info('Starting calculation of flux contributions and pseudo observations.')

    # Get the previously-computed obspack that contains all sector-specific contributations to atmopsheric concentration signal, 
    # to re-calculate the sum faster in the get_flux_contribtion_np() function
    if sum_vars:
        if fluxtype == "PARIS":
            basepath = os.path.dirname(os.path.dirname(outdir)) + '/BASE_SS/'
        else:
            basepath = outdir
        ss_obspack = glob.glob(basepath + '*' + stationname.lower() + '*.nc')[0]

    # Get the original obspack, to use later as a check whether observations are present or not
    obspack_orig_filepath = [s for s in obspack_list if stationname.lower() in s.lower()][0]

    if stationname in non_obspack_sites:
        obsvarname='co2'
    else:
        obsvarname = 'value'

    # Resample the original obspack to rounded hourly data and open it for the remainder of the for-loop
    with xr.open_dataset(obspack_orig_filepath) as ds:
        logging.info(obspack_orig_filepath)
        logging.info(ds)
        logging.info(ds[obsvarname].sel(time=slice(start_date, end_date)))

        obspack_orig_subds = ds[obsvarname].sel(time=slice(start_date, end_date)).resample({'time':'H'}).mean(dim='time', skipna=True, keep_attrs=True).dropna(dim='time')

        # Import list of missing footprint files
        missing_fp_filestr = stilt_rundir + stationname + '/missing.footprints'
        if os.path.isfile(missing_fp_filestr):
            missing_fp_file = pd.read_csv(missing_fp_filestr)
            missing_fp_filelist = list(missing_fp_file['ident'])

        summary_dict={}

        # Extract the first file from the fp filelist
        file = fp_filelist[0]

        # Extract start and end time of first footprint file
        fp_starttime = footprint_extract_starttime(file)
        fp_endtime = footprint_extract_endtime(file, simulation_len=sim_len)
        
        # Extract amount of particles released per STILT simulation from the fp_filelist
        npars = footprint_extract_npars(file)
        
        # Create timerange of simulation
        timerange_sim = pd.date_range(start=start_date, end=end_date, freq='H').to_list()
        timerange_fp = footprint_hours(fp_filelist=fp_filelist, simulation_len=sim_len, shift_forward=False)

        # Load the fluxes corresponding to the first footprint file in the list into memory
        fluxdict = load_first_flux_timestep(flux_ds=flux_ds, fp_starttime=fp_starttime, fp_endtime=fp_endtime, variablelist_vars=variables)

        last_fp_starttime = np.datetime64(fp_starttime)

        # Loop over all simulation times
        for simtime in timerange_sim:

            # Loop over all footprint files
            if simtime in timerange_fp:            
                
                # Check if observation is in ObsPack or not
                obspack_checktime = check_time_in_obspack(obspack_ds=obspack_orig_subds, fp_starttime=simtime)

                if obspack_checktime:
                    logging.info(f'Footprint found! Working on {simtime} !')

                    # Calculate index of footprint file in footprint list using footprint_hours
                    index_fp = timerange_fp.index(simtime)

                    # Get start and end time of footprint
                    file = fp_filelist[index_fp]
                    ident = os.path.splitext(os.path.basename(file).split(sep='_')[-1])[0]

                    # Create keys and values for summary_dict later
                    key = np.datetime64(simtime)
                    
                    # Calculate mean background concentration for current footprint
                    bg = calculate_mean_bg(fp_starttime=fp_starttime, RDatapath=RDatapath, bgpath=bgpath,
                                npars=npars, lat=station_lat, lon=station_lon, agl=station_agl)

                    logging.info(f'Mean background is {bg}')

                    if os.path.exists(missing_fp_filestr)==False or ident not in missing_fp_filelist:
                        with xr.open_dataset(file) as footprint_df:
                            
                            # First update the old fluxdict so that it shifts an X amount in time
                            fluxdict = update_fluxdict(flux_ds=flux_ds, fluxdict=fluxdict, cur_fp_starttime = key, last_fp_starttime = last_fp_starttime)
                            
                            # Extract latitude indices from sparse footprint file and insert into xr.DataArray for later indexing
                            lat_indices, lon_indices = get_latlon_indices(footprint_df, lats, lons)
                            hours_into_file = get_time_indices(flux_ds=fluxdict, footprint_df=footprint_df, 
                                                                            fp_starttime=key)

                            logging.info(f'This footprint has {len(lat_indices)} values, has {len(np.unique(hours_into_file))} timesteps and goes {len(np.unique(hours_into_file))} hours backwards in time')
                            
                            # Calculate CTE-HR contribution
                            if sum_vars == True:
                                cont = get_flux_contribution(flux_ds=fluxdict, time_index_list=hours_into_file,
                                                                    lat_index_list=lat_indices, lon_index_list=lon_indices,
                                                                    footprint_df=footprint_df, infl_varname='Influence',
                                                                    variables=variables, sum_vars=sum_vars,
                                                                    perturbationcode=perturbationcode, fluxtype=fluxtype,
                                                                    fp_starttime=key, ss_obspack=ss_obspack, changed_var=changed_var)
                            else:
                                cont = get_flux_contribution(flux_ds=fluxdict, time_index_list=hours_into_file,
                                                                    lat_index_list=lat_indices, lon_index_list=lon_indices,
                                                                    footprint_df=footprint_df, infl_varname='Influence',
                                                                    variables=variables, sum_vars=sum_vars,
                                                                    fluxtype=fluxtype, perturbationcode=perturbationcode, fp_starttime=key)

                    elif ident in missing_fp_filelist:
                        logging.info(f'Footprint is empty and no surface exchange with particles is recorded, taking background concentration at station location as pseudo observation.')
                        if sum_vars==True:
                            logging.info('Summing variables!')
                            cont = 0

                        else:
                            logging.info('Not summing variables!')
                            cont = [0] * len(variables)
                    
                    # Calculate pseudo observation from bg and flux contribution            
                    summary_dict = create_intermediate_dict(dict_obj=summary_dict, cont=cont, background=bg, key=key, sum_vars=sum_vars, perturbationcode=perturbationcode)
                
                else:
                    logging.warning(f'No observation was found in the ObsPack for {simtime}, skipping this footprint!')
                    continue

                last_fp_starttime = key

            else:
                logging.info(f'No footprint found for {simtime}, skipping this simulation time!')
                continue 

    logging.info('Finished calculating flux contributions and pseudo observations.')

    if save_as_pkl == True:
        logging.info('Saving dictionary as pickle file')

        # Save dictionary as pickle file
        with open('obs_sim_dict_' + stationname + '_' + fp_starttime.yr + '.pkl', 'wb') as f:
            pickle.dump(summary_dict, f)
    
    elif save_as_pkl == False:
        logging.info('Saving dictionary in Python memory.')

        return summary_dict

def write_dict_to_ObsPack(obspack_list, obs_sim_dict, variables, stationcode, start_date = '2021-01-01', end_date = '2022-01-01', 
                          sum_vars = None, mp_start_date = None, mp_end_date = None,
                          outpath='/projects/0/ctdas/PARIS/DATA/obspacks/'):
    """ Function to write a dictionary (that contains the simulation results, the background,
    and pseudo_observations) to a copy of the original ObsPacks in the ObsPack collection.

    Input variables:
    obspack_list: list of ObsPacks that meet the filtering criteria
    obs_sim_dict: dictionary that contains the simulation results, the background, and observations
    stationcode: stationcode of the station
    outpath: path to the output directory
    sum_vars: boolean that indicates whether the variables in the ObsPack file should be summed
    mp_start_date: start date of the current multiprocessing run
    mp_end_date: end date of the current multiprocessing run
    start_date: start date of the simulation
    end_date: end date of the simulation
    """

    obspack_orig = [s for s in obspack_list if stationcode.lower() in s.lower()][0]

    if obspack_orig is None:
        print(f"No ObsPack found for station code {stationcode}")
        return

    # Check if outpath exists, if not, create it
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    if mp_start_date != None and mp_end_date != None:
        # Put simulation results in seperate run-dependent ObsPack files to combine later
        mp_start_datetime = datetime.strptime(mp_start_date, "%Y-%m-%d")
        mp_end_datetime = datetime.strptime(mp_end_date, "%Y-%m-%d")
        
        # Extract multiprocessing timestr to create ObsPack object to save to later
        timestr = str(mp_start_datetime.strftime("%Y-%m-%d-%H:%M")) + '-' + str(mp_end_datetime.strftime("%Y-%m-%d-%H:%M"))
        obspack = os.path.join(outpath, 'pseudo_' + os.path.basename(os.path.splitext(obspack_orig)[0]) + '_' + timestr + '.nc')
    
    else:
        # Put simulation results in one output ObsPack file
        mp_start_datetime = datetime.strptime(start_date, "%Y-%m-%d")
        mp_end_datetime = datetime.strptime(end_date, "%Y-%m-%d")

        # Put all results into one obspack file (only one run, so no intermediate files needed)
        obspack = os.path.join(outpath, 'pseudo_' + os.path.basename(obspack_orig))

    # Define time since epoch, which is also the start time of the ObsPack file
    obspack_basetime = np.datetime64('1970-01-01T00:00:00')

    obsvarname='value'
    decode=True

    # Open original obspack file
    with xr.open_dataset(obspack_orig, decode_times=decode) as ds:
        ds['time'] = pd.to_datetime(ds['time'], unit='s')

        # Define new variables to fill later
        if sum_vars == True:
            new_variables = ['background', 'mixed', 'pseudo_observation']
        else:
            new_variables = variables + ['background', 'mixed', 'pseudo_observation']

        for new_variable_name in new_variables:
            dummy_var = xr.DataArray(
                name=new_variable_name,
                dims=ds.time.dims,
                coords=ds.coords
            )
            ds[new_variable_name] = dummy_var
    

        # Select only the time range of the simulation
        subds = ds.sel(time=slice(mp_start_datetime, mp_end_datetime))
        #non_numeric_vars = [var for var in subds.variables if subds[var].dtype not in ['float64', 'int64','float32','int32', 'int8']][:]

        # Resample to hourly data
        logging.info('Resampling to hourly data')
        subds_rs = subds[obsvarname].resample({'time':'H'}).mean(dim='time', skipna=True, keep_attrs=True).dropna(dim='time')

        # Change time and value variables to resampled counterparts
        subds['time'] = subds_rs['time']
        subds[obsvarname] = subds_rs
        
        # Manually resample ObsPack datetime element
        subds['datetime'][:] = np.array((parse_datetime(subds['datetime']) - pd.to_timedelta('30T')).round('H').strftime("%Y-%m-%dT%H:%M:%S."), dtype='S20')

        # Calculate time differences only once
        ds_time_seconds = ((subds['time'] - obspack_basetime).values.astype('int64')) / 1e9

        # Loop over all keys and values in dictionary
        logging.info(f'Looping over {len(obs_sim_dict)} footprints')

        logging.info(obs_sim_dict.items())
        for key, values in obs_sim_dict.items():
            if key not in subds['time']:
                logging.info(f'Current footprint time is {key}, for which no observation was found in the ObsPack. Skipping this footprint.')
                continue
            
            else:
                logging.info(f'Current footprint time is {key}')

                # Calculate time in seconds since start of ObsPack file
                obspack_time = int(((key - obspack_basetime) / np.timedelta64(1, 's')).item())
                
                # Calculate time difference between simulation time and ObsPack time
                diff_dt = obspack_time - ds_time_seconds
                diff_index = np.abs(diff_dt).argmin()
                #diff_min = diff_dt[diff_index]

                # Convert values to numpy array
                values = np.array(values)
                logging.info(values)

                # Write values to new variables in new ObsPack file
                for v in new_variables:
                    subds[v][diff_index] = values[new_variables.index(v)]

    # Save and close the new file
    logging.info(f'Writing results to {obspack}')
    subds.to_netcdf(obspack, mode='w')


def write_dict_to_ObsPack_altsites(obspack_list, obs_sim_dict, variables, stationcode, start_date = '2021-01-01', end_date = '2022-01-01', 
                          sum_vars = None, mp_start_date = None, mp_end_date = None,
                          outpath='/projects/0/ctdas/PARIS/DATA/obspacks/'):
    """ Function to write a dictionary (that contains the simulation results, the background,
    and pseudo_observations) to a copy of the original observation file.

    NOTE: This function is adjusted to work with files that are of a non-Obspack format! Please adjust
    to the needs of your own data format if some NetCDF variables are missing or have different names.

    Input variables:
    obspack_list: list of ObsPacks that meet the filtering criteria
    obs_sim_dict: dictionary that contains the simulation results, the background, and observations
    stationcode: stationcode of the station
    outpath: path to the output directory
    sum_vars: boolean that indicates whether the variables in the ObsPack file should be summed
    mp_start_date: start date of the current multiprocessing run
    mp_end_date: end date of the current multiprocessing run
    start_date: start date of the simulation
    end_date: end date of the simulation
    """

    obspack_orig = [s for s in obspack_list if stationcode.lower() in s.lower()][0]

    if obspack_orig is None:
        print(f"No ObsPack found for station code {stationcode}")
        return

    # Check if outpath exists, if not, create it
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    if mp_start_date != None and mp_end_date != None:
        # Put simulation results in seperate run-dependent ObsPack files to combine later
        mp_start_datetime = datetime.strptime(mp_start_date, "%Y-%m-%d")
        mp_end_datetime = datetime.strptime(mp_end_date, "%Y-%m-%d")
        
        # Extract multiprocessing timestr to create ObsPack object to save to later
        timestr = str(mp_start_datetime.strftime("%Y-%m-%d-%H:%M")) + '-' + str(mp_end_datetime.strftime("%Y-%m-%d-%H:%M"))

        # Put all results into one obspack file (only one run, so no intermediate files needed)
        obspack = os.path.join(outpath, 'pseudo_co2_' + stationcode.lower() + '_tower-insitu_xxx_allvalid-' + obspack_orig.split(sep='-')[-2][0:3] + 'magl_' + timestr + '.nc')

    else:
        # Put simulation results in one output ObsPack file
        mp_start_datetime = datetime.strptime(start_date, "%Y-%m-%d")
        mp_end_datetime = datetime.strptime(end_date, "%Y-%m-%d")

        # Put all results into one obspack file (only one run, so no intermediate files needed)
        obspack = os.path.join(outpath, 'pseudo_co2_' + stationcode.lower() + '_tower-insitu_xxx_allvalid-' + obspack_orig.split(sep='-')[-2][0:3] + 'magl.nc')

    # Define time since epoch, which is also the start time of the ObsPack file
    obspack_basetime = np.datetime64('1970-01-01T00:00:00')
    
    obsvarname='co2'
    decode=True

    # Open newfile in 'append' mode
    logging.info('Opening ObsPack file: ' + str(obspack_orig))
    with xr.open_dataset(obspack_orig, decode_times=decode) as ds:
        ds['time'] = pd.to_datetime(ds['time'], unit='s')

        # Select only the time range of the simulation
        subds = ds.sel(time=slice(mp_start_datetime, mp_end_datetime))
        #non_numeric_vars = [var for var in subds.variables if subds[var].dtype not in ['float64', 'int64','float32','int32', 'int8']][:]

        # Resample to hourly data
        #subds_rs = subds[obsvarname].resample({'time':'H'}).mean(dim='time', skipna=True, keep_attrs=True).dropna(dim='time')
        subds_rs = subds.resample({'time':'H'}).mean(dim='time', skipna=True, keep_attrs=True).dropna(dim='time')

        # Change time and value variables to resampled counterparts
        #subds['time'] = subds_rs['time']

        # Correct for units
        subds_rs[obsvarname] = subds_rs[obsvarname] * 1e-6
        subds_rs[obsvarname].attrs['units'] = 'mol mol-1'

        subds_rs['co2_variability'] = subds_rs['co2_variability'] * 1e-6
        subds_rs['co2_variability'].attrs['units'] = 'mol mol-1'

        # Calculate time differences only once
        ds_time_seconds = ((subds_rs['time'] - obspack_basetime).values.astype('int64')) / 1e9

        # Define new variables to fill later
        if sum_vars == True:
            new_variables = ['background', 'mixed', 'pseudo_observation']
        else:
            new_variables = variables + ['background', 'mixed', 'pseudo_observation']

        for new_variable_name in new_variables:
            dummy_var = xr.DataArray(
                name=new_variable_name,
                dims=subds_rs.time.dims,
                coords=subds_rs.coords
            )
            subds_rs[new_variable_name] = dummy_var

        # Loop over all keys and values in dictionary
        logging.info(f'Looping over {len(obs_sim_dict)} footprints')
        logging.info(obs_sim_dict.items())

        for key, values in obs_sim_dict.items():
            if key not in subds_rs['time']:
                logging.info(f'Current footprint time is {key}, for which no observation was found in the ObsPack. Skipping this footprint.')
                continue
            
            else:
                logging.info(f'Current footprint time is {key}')

                # Calculate time in seconds since start of ObsPack file
                obspack_time = int(((key - obspack_basetime) / np.timedelta64(1, 's')).item())
                
                # Calculate time difference between simulation time and ObsPack time
                diff_dt = obspack_time - ds_time_seconds
                diff_index = np.abs(diff_dt).argmin()
                #diff_min = diff_dt[diff_index]

                # Convert values to numpy array
                values = np.array(values)
                logging.info(values)

                # Write values to new variables in new ObsPack file
                for v in new_variables:
                    subds_rs[v][diff_index] = values[new_variables.index(v)]

    # Save and close the new file
    logging.info(f'Writing results to {obspack}')
    subds_rs.to_netcdf(obspack, mode='w')


def main(filepath, sim_length, fluxdir, stilt_rundir, fluxvarnamelist,  stationsfile, stationcode, outpath, bgfilepath, obspack_list,
         lats, lons, non_obspack_sites=None, perturbationcode=None, sum_vars=None, verbose=None, fluxtype=None, start_date=None, end_date=None, mp_start_date = None,
         mp_end_date = None, date_pair=None):
    """ Main function to run the footprint - flux multiplication script. """

    if verbose == True:
        logging.basicConfig(level=logging.INFO, format=' [%(levelname)-7s] (%(asctime)s) py-%(module)-20s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    else:
        logging.basicConfig(level=logging.WARNING, format=' [%(levelname)-7s] (%(asctime)s) py-%(module)-20s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    # Get 3D station location
    lat,lon,agl = get_3d_station_location(stationsfile, stationcode, colname_lat='lat', colname_lon='lon',
                                            colname_agl='corrected_alt')  
    
    # Get list of footprint files for station
    sparse_files = sorted(glob.glob(filepath + 'footprint_' + stationcode + '*.nc'))

    if fluxtype == 'PARIS':
        # Only once per station, create list of all PARIS flux files
        fluxstring = find_fluxfiles(fluxdir = fluxdir, perturbationcode=perturbationcode, fluxtype=fluxtype)

        # Filter list of footprint files by date range (using start_date and end_date, if they exist!)
        mp_start_date, mp_end_date = date_pair
        sparse_files = filter_files_by_date(sparse_files, start_date=mp_start_date, end_date=mp_end_date)
        
        if ((perturbationcode != 'BASE') & (perturbationcode != 'BASE_SS')):
            if perturbationcode=='HGER' or perturbationcode=='ATEN' or perturbationcode=='PTEN' or perturbationcode=='HFRA':
                changed_var = 'flux_ff_exchange_prior'
            elif perturbationcode=='DFIN':
                changed_var = 'flux_bio_exchange_prior'
            fluxvarnamelist.remove(changed_var)
            
            # Open all files in fluxstring as xr_mfdataset, and add variables in variablelist
            #flux_ds = hard_load_all_fluxes(fluxstring, sim_len=sim_length, mp_start_date = mp_start_date, 
            #                        mp_end_date = mp_end_date, fluxtype=fluxtype, variables = [changed_var])
            
            flux_ds = lazy_load_all_fluxes(fluxstring, sim_len=sim_length, start_date = mp_start_date, 
                                    end_date = mp_end_date, variablelist_vars = fluxvarnamelist, 
                                    fluxtype=fluxtype)
            
            # Loop over all footprints files
            obs_sim_dict = create_obs_sim_dict(fp_filelist = sparse_files, flux_ds = flux_ds,
                        lats = lats, lons = lons, RDatapath = filepath, sim_len = sim_length,
                        bgpath = bgfilepath, station_lat = lat, station_lon = lon, start_date = mp_start_date, end_date = mp_end_date,
                        station_agl = agl, stationname = stationcode, variables = fluxvarnamelist, non_obspack_sites=non_obspack_sites,
                        stilt_rundir = stilt_rundir, sum_vars=sum_vars, changed_var=changed_var, perturbationcode=perturbationcode, 
                        outdir = outpath, fluxtype = fluxtype, obspack_list = obspack_list, save_as_pkl=False)
            
        else:
            # Open all files in fluxstring as xr_mfdataset, and add variables in variablelist
            #flux_ds = hard_load_all_fluxes(fluxstring, sim_len=sim_length, mp_start_date = mp_start_date, 
            #                        mp_end_date = mp_end_date, fluxtype=fluxtype, variables = fluxvarnamelist)
            
            flux_ds = lazy_load_all_fluxes(fluxstring, sim_len=sim_length, start_date = mp_start_date, 
                                    end_date = mp_end_date, variablelist_vars = fluxvarnamelist, 
                                    fluxtype=fluxtype)
            
            # Loop over all footprints files
            obs_sim_dict = create_obs_sim_dict(fp_filelist = sparse_files, flux_ds = flux_ds,
                        lats = lats, lons = lons, RDatapath = filepath, sim_len = sim_length,
                        bgpath = bgfilepath, station_lat = lat, station_lon = lon, start_date = mp_start_date, end_date = mp_end_date,
                        station_agl = agl, stationname = stationcode, variables = fluxvarnamelist, non_obspack_sites=non_obspack_sites,
                        stilt_rundir = stilt_rundir, sum_vars=sum_vars, perturbationcode=perturbationcode, 
                        outdir = outpath, fluxtype = fluxtype, obspack_list = obspack_list, save_as_pkl=False)
        
        if stationcode in non_obspack_sites:
            write_dict_to_ObsPack_altsites(obspack_list = obspack_list, obs_sim_dict = obs_sim_dict, sum_vars = sum_vars,
                        variables=fluxvarnamelist, stationcode = stationcode, 
                        mp_start_date = mp_start_date, mp_end_date = mp_end_date,
                        start_date= start_date, end_date = end_date, outpath = outpath)
            
        else:
            write_dict_to_ObsPack(obspack_list = obspack_list, obs_sim_dict = obs_sim_dict, sum_vars = sum_vars,
                            variables=fluxvarnamelist, stationcode = stationcode, 
                            mp_start_date = mp_start_date, mp_end_date = mp_end_date,
                            start_date= start_date, end_date = end_date, outpath = outpath)
 
    elif fluxtype == 'CTEHR':
        # Define which CTE-HR variables to find fluxfiles for
        fluxfilenamelist = ['nep', 'fire', 'ocean', 'anthropogenic']

        # Extract all unique months between start and end time
        mons = footprint_unique_months(sparse_files, sim_length)

        # Only once per station, create list of all CTE-HR flux files
        fluxstring = find_fluxfiles(fluxdir = fluxdir, variablelist_files = fluxfilenamelist, months = mons, fluxtype=fluxtype)
        
        # Open all files in fluxstring as xr_mfdataset, and add variables in variablelist
        #flux_ds = hard_load_all_fluxes(fluxstring, sim_len=sim_length, start_date = start_date, 
        #                        end_date = end_date, variables = fluxvarnamelist, 
        #                        fluxtype=fluxtype)
        
        flux_ds = lazy_load_all_fluxes(fluxstring, sim_len=sim_length, start_date = mp_start_date, 
                        end_date = mp_end_date, variablelist_vars = fluxvarnamelist, 
                        fluxtype=fluxtype)
        
        # Loop over all footprints files
        obs_sim_dict = create_obs_sim_dict(fp_filelist = sparse_files, flux_ds = flux_ds,
                        lats = lats, lons = lons, RDatapath = filepath, sim_len = sim_length,
                        bgpath = bgfilepath, station_lat = lat, station_lon = lon, start_date = mp_start_date, end_date = mp_end_date,
                        station_agl = agl, stationname = stationcode, variables = fluxvarnamelist,
                        stilt_rundir = stilt_rundir, sum_vars=sum_vars, non_obspack_sites=non_obspack_sites,
                        outdir = outpath, fluxtype = fluxtype, obspack_list = obspack_list, save_as_pkl=False)

        if stationcode in non_obspack_sites:
            write_dict_to_ObsPack_altsites(obspack_list = obspack_list, obs_sim_dict = obs_sim_dict, sum_vars = sum_vars,
                        variables=fluxvarnamelist, stationcode = stationcode, 
                        start_date= start_date, end_date = end_date, outpath = outpath)
            
        else:
            write_dict_to_ObsPack(obspack_list = obspack_list, obs_sim_dict = obs_sim_dict, sum_vars = sum_vars,
                            variables=fluxvarnamelist, stationcode = stationcode,
                            start_date= start_date, end_date = end_date, outpath = outpath)
            
    elif fluxtype == 'Ruben':
        # Define which CTE-HR variables to find fluxfiles for
        fluxfilenamelist = ['nep']
        
        # Filter files based on simulation timerange
        sparse_files = filter_files_by_date(sparse_files, start_date=start_date, end_date=end_date)

        # Only once per station, create list of all CTE-HR flux files
        fluxstring = find_fluxfiles(fluxdir = fluxdir, variablelist_files = fluxfilenamelist, fluxtype=fluxtype)
        
        # Open all files in fluxstring as xr_mfdataset, and add variables in variablelist
        #flux_ds = hard_load_all_fluxes(fluxstring, sim_len=sim_length, mp_start_date = mp_start_date, 
        #                            mp_end_date = mp_end_date, fluxtype=fluxtype, variables = fluxvarnamelist)
         
        flux_ds = lazy_load_all_fluxes(fluxstring, sim_len=sim_length, start_date = mp_start_date, 
                        end_date = mp_end_date, variablelist_vars = fluxvarnamelist, 
                        fluxtype=fluxtype)
        
        # Loop over all footprints files
        obs_sim_dict = create_obs_sim_dict(fp_filelist = sparse_files, flux_ds = flux_ds,
                        lats = lats, lons = lons, RDatapath = filepath, sim_len = sim_length,
                        bgpath = bgfilepath, station_lat = lat, station_lon = lon, start_date = mp_start_date, end_date = mp_end_date,
                        station_agl = agl, stationname = stationcode, variables = fluxvarnamelist,
                        stilt_rundir = stilt_rundir, sum_vars=sum_vars, non_obspack_sites=non_obspack_sites,
                        outdir = outpath, fluxtype = fluxtype, obspack_list = obspack_list, save_as_pkl=False)

        if stationcode in non_obspack_sites:
            write_dict_to_ObsPack_altsites(obspack_list = obspack_list, obs_sim_dict = obs_sim_dict, sum_vars = sum_vars,
                        variables=fluxvarnamelist, stationcode = stationcode, 
                        start_date= start_date, end_date = end_date, outpath = outpath)
            
        else:
            write_dict_to_ObsPack(obspack_list = obspack_list, obs_sim_dict = obs_sim_dict, sum_vars = sum_vars,
                            variables=fluxvarnamelist, stationcode = stationcode,
                            start_date= start_date, end_date = end_date, outpath = outpath)
            
    else:
        print("Please provide a valid fluxtype (choose between either 'PARIS' or 'CTEHR' for now)!")

    logging.info(f'Finished writing results to ObsPack file for station {stationcode}')