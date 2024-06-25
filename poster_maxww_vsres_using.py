# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 10:50:48 2024

@author: klb236
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from datetime import datetime

# Constants for normalization
wtheta = 300 / (1 * 1004)
wb = (10 * wtheta) / 300
wstar = (wb * 1000) ** (1 / 3)

# Base file paths for UM datasets
um_basepath_conventional = 'C:/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/UM_MONK/um_data_hires/'
um_basepath_standard = 'C:/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/UM_MONK/um_data_hires/diff_sgs/'

# Directories for MONC data
monc_dir_conventional = 'C:/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/UM_MONK/monc_data/sul_pat/'
monc_dir_standard = 'C:/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/UM_MONK/monc_data/sul_pat/diff_sgs/'

# Director for DNS
flf_dir_500m = '/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/data/cbl_500m_flf/run000001columndiags.nc'
flf_dir_250m = '/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/data/cbl_250m_flf/run000001columndiags.nc'
flf_dir_1000m = '/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/data/cbl_1000m_flf/run000001columndiags.nc'

ds_flf_500m = xr.open_dataset(flf_dir_500m)
ds_flf_250m = xr.open_dataset(flf_dir_250m)
ds_flf_1000m = xr.open_dataset(flf_dir_1000m)

# Extract wvar (vertical velocity variance) for both datasets
#wvar_flf_500m = ds_flf_500m['wvar'].isel(time=slice(None, None,10))
wvar_flf_500m = ds_flf_500m['wvar']
wvar_flf_250m = ds_flf_250m['wvar']
wvar_flf_1000m = ds_flf_1000m['wvar']

# Extract time coordinates
time_flf_500m = ds_flf_500m['time']
time_flf_250m = ds_flf_250m['time']
time_flf_1000m = ds_flf_1000m['time']

# Calculate the maximum wvar for each time step across all vertical levels
wvar_flf_max_500m = wvar_flf_500m.max(dim='zp')
wvar_flf_max_250m = wvar_flf_250m.max(dim='zp')
wvar_flf_max_1000m = wvar_flf_1000m.max(dim='zp')

wvar_flf_max_500m_normalised = wvar_flf_max_500m / (wstar ** 2)
wvar_flf_max_250m_normalised = wvar_flf_max_250m / (wstar ** 2)
wvar_flf_max_1000m_normalised = wvar_flf_max_1000m / (wstar ** 2)

flf_res = [200, 500, 1000] #need to change to 250m!!! but wont align otherwise
flf_res_norm = [0.2, 0.5, 1]
wvar_last = []
wvar_last.append(wvar_flf_max_250m_normalised[-1])
wvar_last.append(wvar_flf_max_500m_normalised[-1])
wvar_last.append(wvar_flf_max_1000m_normalised[-1])


# Function to load UM data
def load_um_data(filepath):
    um_dataset = xr.open_dataset(filepath)
    wvar_um = um_dataset['STASH_m01s00i150']
    time_um = um_dataset['min15T0_0']
    return wvar_um, time_um

# Function to load MONC data for the last time step (hour 4)
def load_monc_data_last_time(directory, expt):
    nc_file = f"{directory}/{expt}/cbl_14400.nc"
    with xr.open_dataset(nc_file) as nc:
        w_data = nc['w']
        w = w_data.sel(time_series_300_300 = 14400, method='nearest')
        w = w.rename({'time_series_300_300':'time'})
        print(w.time)
        ww_monc_all = (w**2).mean(dim=['x', 'y'])# Calculate wvar from w
    return ww_monc_all



# Function to convert cftime to datetime
def convert_to_datetime(cftime_obj):
    return datetime(cftime_obj.year, cftime_obj.month, cftime_obj.day, cftime_obj.hour, cftime_obj.minute, cftime_obj.second)

# Define the resolutions to loop through
resolutions = ['50m', '100m', '200m', '250m', '500m', '1km']
#nx_values = [1024, 512, 256, 128, 64]  # Example grid sizes corresponding to the resolutions
nx_values = [50, 100, 200, 250, 500, 1000]
nx_values_norm = [0.05, 0.1, 0.2, 0.25, 0.5, 1]
# Initialize lists to store results
max_ww_monc_conventional_list = []
max_ww_monc_standard_list = []
max_ww_um_conventional_list = []
max_ww_um_standard_list = []

# Loop through each resolution
for resolution, nx in zip(resolutions, nx_values):
    # Load MONC data for the last time step
    ww_monc_conventional_last = load_monc_data_last_time(monc_dir_conventional, resolution)
    ww_monc_standard_last = load_monc_data_last_time(monc_dir_standard, resolution)

    # Normalize MONC data with wstar
    ww_normalized_monc_conventional_last = ww_monc_conventional_last / (wstar ** 2)
    ww_normalized_monc_standard_last = ww_monc_standard_last / (wstar ** 2)

    print(ww_normalized_monc_conventional_last.dims)
    # Extract the maximum ww at the last time point (hour 4)
    max_ww_monc_conventional = ww_normalized_monc_conventional_last.max(dim='z')
    max_ww_monc_standard = ww_normalized_monc_standard_last.max(dim='z')

    # Append to the lists
    max_ww_monc_conventional_list.append(max_ww_monc_conventional)
    max_ww_monc_standard_list.append(max_ww_monc_standard)

    # Load UM data
    um_filepath_conventional = f'{um_basepath_conventional}da249_{resolution}a_pr000.nc'
    um_filepath_standard = f'{um_basepath_standard}da249_{resolution}a_pr000.nc'

    wvar_um_conventional, time_um_conventional = load_um_data(um_filepath_conventional)
    wvar_um_standard, time_um_standard = load_um_data(um_filepath_standard)

    # Calculate ww_normalized for UM conventional
    ww_um_conventional = np.zeros_like(wvar_um_conventional.values)
    num_times = wvar_um_conventional.shape[0]  # Get the number of time steps
    num_levels = wvar_um_conventional.shape[1]  # Get the number of vertical levels
    for t in range(num_times):  # Loop over time steps
        for k in range(num_levels):  # Loop over vertical levels
            w_section = wvar_um_conventional[t, k, :, :].values  # Extract w at each level and time step
            ww_um_conventional[t, k] = np.mean(np.square(w_section))  # Calculate mean squared w for each time and level

    # Normalize ww with wstar
    ww_normalized_um_conventional = ww_um_conventional / (wstar ** 2)

    # Calculate ww_normalized for UM standard
    ww_um_standard = np.zeros_like(wvar_um_standard.values)
    num_times = wvar_um_standard.shape[0]  # Get the number of time steps
    num_levels = wvar_um_standard.shape[1]  # Get the number of vertical levels
    for t in range(num_times):  # Loop over time steps
        for k in range(num_levels):  # Loop over vertical levels
            w_section = wvar_um_standard[t, k, :, :].values  # Extract w at each level and time step
            ww_um_standard[t, k] = np.mean(np.square(w_section))  # Calculate mean squared w for each time and level

    # Normalize ww with wstar
    ww_normalized_um_standard = ww_um_standard / (wstar ** 2)

    # Convert times to datetime
    time_um_conventional_dt = [convert_to_datetime(t) for t in time_um_conventional.values]
    time_um_standard_dt = [convert_to_datetime(t) for t in time_um_standard.values]

    # Normalize times by setting the start time to 0 hours
    start_time_um_conventional = time_um_conventional_dt[0]
    start_time_um_standard = time_um_standard_dt[0]

    time_um_conventional_normalized = [(t - start_time_um_conventional).total_seconds() for t in time_um_conventional_dt]
    time_um_standard_normalized = [(t - start_time_um_standard).total_seconds() for t in time_um_standard_dt]

    # Convert times to hours
    time_um_conventional_hours = np.array(time_um_conventional_normalized) / 3600
    time_um_standard_hours = np.array(time_um_standard_normalized) / 3600

    # Extract the maximum ww at the last time point (hour 4) for UM models
    max_ww_um_conventional = np.interp(4, time_um_conventional_hours, ww_normalized_um_conventional.max(axis=(1, 2, 3)))
    max_ww_um_standard = np.interp(4, time_um_standard_hours, ww_normalized_um_standard.max(axis=(1, 2, 3)))

    # Append to the lists
    max_ww_um_conventional_list.append(max_ww_um_conventional)
    max_ww_um_standard_list.append(max_ww_um_standard)
    
print(max_ww_um_conventional_list[3])
print(max_ww_um_standard_list[3])
print(max_ww_monc_conventional_list[3])
print(max_ww_monc_standard_list[3])


# Plot the maximum ww for each resolution
plt.figure(figsize=(10, 6))
plt.plot(nx_values_norm, max_ww_monc_conventional_list, color='red', marker = 'o', linestyle='-', label='MONC Conventional')
plt.plot(nx_values_norm, max_ww_monc_standard_list, color='red', marker = 's',  linestyle='--', label='MONC Standard')
plt.plot(nx_values_norm, max_ww_um_conventional_list, color='blue', marker = '^', linestyle='-', label='UM Conventional')
plt.plot(nx_values_norm, max_ww_um_standard_list, color='blue', marker = 'd', linestyle='--', label='UM Standard')
plt.scatter(flf_res_norm, wvar_last, color='green', marker='x', label='ILES')


plt.xlabel('Horizontal Resolution ($\Delta x / z_i$)', fontsize=18)
plt.ylabel('$(\max \; w^{2})/(w^*)^2$', fontsize=18)
plt.title('Maximum Vertical Velocity Variance (Time = 4h)', fontsize=18)
plt.xticks((0.05, 0.25, 0.5, 1), fontsize=18)
plt.yticks(fontsize=18)
plt.legend(loc='lower left', fontsize=18)
plt.tight_layout()
plt.savefig('C:/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/UM_MONK/hires_plots/max_wvar_resolution_comparison.png', dpi=300)
plt.show()
