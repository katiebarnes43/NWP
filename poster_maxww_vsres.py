# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 04:16:07 2024

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
resolutions = ['50m', '100m', '200m', '500m', '1km']
nx_values = [1024, 512, 256, 128, 64]  # Example grid sizes corresponding to the resolutions

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
plt.plot(resolutions, max_ww_monc_conventional_list, color='red', linestyle='-', label='MONC Conventional')
plt.plot(resolutions, max_ww_monc_standard_list, color='red', linestyle='--', label='MONC Standard')
plt.plot(resolutions, max_ww_um_conventional_list, color='blue', linestyle='-', label='UM Conventional')
plt.plot(resolutions, max_ww_um_standard_list, color='blue', linestyle='--', label='UM Standard')

plt.xlabel(' Horizontal Resolution (m)')
plt.ylabel('Max Vertical Velocity Variance (m$^{2}$s$^{-2}$/(w$^*)^2$)')
plt.title('Maximum Vertical Velocity Variance at Hour 4 for Each Resolution')
plt.xticks(rotation=45)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig('C:/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/UM_MONK/hires_plots/max_wvar_resolution_comparison.pdf')
plt.show()
