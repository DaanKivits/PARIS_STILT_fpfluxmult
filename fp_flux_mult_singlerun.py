"""
This script performs the footprint multiplication process using multiple workers.

The script takes input files containing footprints and fluxes, and multiplies them together to calculate the CO2
amount fraction at receptor locations. It uses multiple workers to parallelize the computation and improve performance.

For this workflow, STILT (Stochastic Time-Inverted Lagrangian Transport) influence fields stored in a sparse format are expected as input,
as well as CTE-HR flux fields or PARIS-specific flux fields (which are modified CTE-HR fluxes). 

Usage:
    python fp_flux_mult_MP_multiworkers.py --station <stationcode> --stationsfile <stationsfile> --fluxdir <fluxdir> --fpdir <fpdir> --bgdir <bgdir> --outdir <outdir> 
    --stiltdir <stiltdir> --obspack_path <obspack_path> --non_obspack_path <non_obspack_path> --start_date <start_date> --end_date <end_date> --lat_ll <lat_ll> --lat_ur <lat_ur> 
    --lon_ll <lon_ll> --lon_ur <lon_ur> --lat_step <lat_step> --lon_step <lon_step> --sim_len <sim_len> --npars <npars>
    --perturbation <perturbation> --verbose --sum-variables

Arguments:

    --station <stationcode>: Code of station to run the script for (str)
    --stationsfile <stationsfile>: Path to the (.csv) file that contains all station metadata (str)
    --fluxdir <fluxdir>: Directory where flux files are stored (str)
    --fpdir <fpdir>: Directory where footprint files are stored (str)
    --bgdir <bgdir>: Directory where (TM3) background files are stored (str)
    --outdir <outdir>: Directory where output files are stored (str)
    --stiltdir <stiltdir>: Directory where STILT is executed (str)
    --obspack_path <obspack_path>: Path to the ObsPack collection zip file (str)
    --non_obspack_path <non_obspack_path>: Path to any other extra observation files that are not included in the European Obspack collection zip file (str)
    --start_date <start_date>: Date from when to subset the ObsPacks (str)
    --end_date <end_date>: Date up to where to subset the ObsPacks (str)
    --lat_ll <lat_ll>: Latitude of ll corner to define the grid extent (float)
    --lat_ur <lat_ur>: Latitude of ur corner to define the grid extent (float)
    --lon_ll <lon_ll>: Longitude of ll corner to define the grid extent (float)
    --lon_ur <lon_ur>: Longitude of ur corner to define the grid extent (float)
    --lat_step <lat_step>: Latitude cell size step (float)
    --lon_step <lon_step>: Longitude cell size step (float)
    --sim_len <sim_len>: Simulation length in hours (int)
    --perturbation <perturbation>: Controls for which perturbation experiment the fp-flux multiplication script is ran (str)
    --verbose: Controls whether the script should do logging (bool)
    --sum-variables: Controls whether the script should sum the fluxvariables or not. TRUE: only produce a "mixed" molefraction field. FALSE: transport all fluxvariables individually and store them in ObsPack (bool)

Note:
    This script requires the following modules to be loaded onto the HPC or locally:

        - xarray
        - numpy
        - pandas
        - netCDF4 (for background concentration lookup functionality)

    The rest of the modules should be available in the standard Python 3.11.0 environment.

"""

# Import necessary modules
import pandas as pd
from functions.fluxfile_functions_dynamicfluxes import *
from functions.background_functions import *
from functions.obspack_identification_url import *
import argparse
import logging

# Set logging options
logging.basicConfig(level=logging.INFO, format=' [%(levelname)-7s] (%(asctime)s) py-%(module)-20s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# Load command-line arguments
parser = argparse.ArgumentParser(description='Flux - footprint multiplication script')

######################
##### USER INPUT #####
######################
# Define command-line arguments
parser.add_argument('--station', type=str, help='Code of station to run the script for')
parser.add_argument('--stationsfile', type=str, help='Path to the (.csv) file that contains all station metadata')
parser.add_argument('--fluxdir', type=str, help='Directory where flux files are stored')
parser.add_argument('--fpdir', type=str, help='Directory where footprint files are stored')
parser.add_argument('--bgdir', type=str, help='Directory where (TM3) background files are stored')
parser.add_argument('--outdir', type=str, help='Directory where output files are stored')
parser.add_argument('--stiltdir', type=str, help='Directory where STILT is executed')
parser.add_argument('--obspack_path', type=str, help='Path to the ObsPack collection zip file')
parser.add_argument('--non_obspack_path', type=str, help='Path to any other extra observation files that are not included in the European Obspack collection zip file')
parser.add_argument('--start_date', type=str, help='Date from when to subset the ObsPacks')
parser.add_argument('--end_date', type=str, help='Date up to where to subset the ObsPacks')
parser.add_argument('--lat_ll', type=float, help='Latitude of ll corner to define the grid extent (float)')
parser.add_argument('--lat_ur', type=float, help='Latitude of ur corner to define the grid extent (float)')
parser.add_argument('--lon_ll', type=float, help='Longitude of ll corner to define the grid extent (float)')
parser.add_argument('--lon_ur', type=float, help='Longitude of ur corner to define the grid extent (float)')
parser.add_argument('--lat_step', type=float, help='Latitude cell size step (float)')
parser.add_argument('--lon_step', type=float, help='Longitude cell size step (float)')
parser.add_argument('--sim_len', type=int, help='Simulation length in hours (int)')
parser.add_argument('--fluxtype', type=str, help='A string that defines the type of fluxes that need to be atmospherically transported (str)')
parser.add_argument('--perturbation', nargs='?', type=str, help='If fluxtpye == PARIS, this string controls for which perturbation experiment the fp-flux multiplication script is ran (str)')
parser.add_argument('--verbose', action="store_true", help='Controls whether the script should do logging (bool)')
parser.add_argument('--sum-variables', action="store_true", help='Controls whether the script should sum the fluxvariables or not. TRUE: only produce a "mixed" molefraction field. FALSE: transport all fluxvariables individually and store them in ObsPack (bool)')

# Parse the command-line arguments
args = parser.parse_args()

# Access the values of the arguments
stationcode = args.station
stationsfile = args.stationsfile
fluxdir = args.fluxdir
filepath = args.fpdir
bgfilepath = args.bgdir
outpath = args.outdir
stilt_rundir = args.stiltdir
obspack_filepath = args.obspack_path
non_obspack_filepath = args.non_obspack_path
start_date = args.start_date
end_date = args.end_date
fluxtype=args.fluxtype
perturbation=args.perturbation
verbose=args.verbose
sum_vars=args.sum_variables

# Provide extent
ll_lat = args.lat_ll
ur_lat = args.lat_ur
ll_lon = args.lon_ll
ur_lon = args.lon_ur

# Provide lat/lon step
lat_step = args.lat_step
lon_step = args.lon_step

# Define simulation length
sim_length = args.sim_len

# Define stationlist
stationsfile = pd.read_csv(stationsfile, header=0)
stationslist = stationsfile['code']

# Create a variable list to loop over later, depending on the given fluxtype
if fluxtype == 'PARIS':
    fluxvarnamelist = ['flux_ff_exchange_prior', 'flux_ocean_exchange_prior', 'flux_fire_exchange_prior', 'flux_bio_exchange_prior']
if fluxtype == 'PARIS_SS':
    fluxvarnamelist = ['A_Public_power', 'B_Industry', 'C_Other_stationary_combustion_consumer', 'F_On-road', 'H_Aviation', 'I_Off-road',
 'G_Shipping', 'cement', 'combustion', 'flux_ff_exchange_prior', 'flux_ocean_exchange_prior', 'flux_fire_exchange_prior', 'flux_bio_exchange_prior']
elif fluxtype == 'CTEHR':
    fluxvarnamelist = ['nep', 'fire', 'ocean', 'combustion']
elif fluxtype == 'Ruben':
    fluxvarnamelist = ['nep']

########################
##### DON'T CHANGE #####
########################
# Define lat/lon grid
# here, the extent is given in corner coordinates, but the lats and lons are defined in the center of the grid cells!
ctehr_bbox = [[ll_lon, ur_lon], [ll_lat, ur_lat]]
lats, lons = coordinate_list(lat_ll=ll_lat, lat_ur=ur_lat, lon_ll=ll_lon, lon_ur=ur_lon, lat_step=lat_step, lon_step=lon_step)

# Get dict of ObsPack files that are within the CTE-HR domain, have the highest inlet height, and contain data within the given time range.
# This contains both the filenames of the ObsPacks as well as the Dobj objects.
obspack_list = filter_obspacks(zip_filepath = obspack_filepath, extent = ctehr_bbox, start_date = start_date, end_date = end_date)
non_obspack_sites = unique_stationcodes(non_obspack_filepath=non_obspack_filepath, case='upper')
complete_filelist = add_to_obspacklist(obspack_list=obspack_list, non_obspack_filepath=non_obspack_filepath, extent = ctehr_bbox, start_date = start_date, end_date = end_date)

if __name__ == "__main__":
    if (start_date != None and end_date != None):

        # Create pairs of consecutive dates
        date_pair = [start_date, end_date]
        print("Time range over which script is ran: " + str(date_pair))

        main(
            filepath=filepath,
            fluxdir=fluxdir,
            bgfilepath=bgfilepath,
            outpath=outpath,
            stilt_rundir=stilt_rundir,
            perturbationcode = perturbation,
            fluxvarnamelist=fluxvarnamelist,
            lats=lats,
            lons=lons,
            sim_length=sim_length,
            obspack_list=complete_filelist,
            stationsfile=stationsfile,
            stationcode=stationcode,
            non_obspack_sites=non_obspack_sites,
            fluxtype=fluxtype,
            sum_vars=sum_vars,
            start_date=start_date,
            end_date=end_date,
            date_pair=date_pair,
            verbose=verbose
            )
        
    else:
        print("Please provide a start and end date in the format 'YYYY-MM-DD'")

    # Log total time
    print("Flux multiplication finished for station " + stationcode + "!")