# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 10:40:50 2024

@author: klb236
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from datetime import datetime
import cftime

# For Normalising
wtheta = 300 / (1 * 1004)
wb = (10 * wtheta) / 300
wstar = (wb * 1000) ** (1 / 3)

# File paths for UM datasets
um_filepath_conventional = '/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/UM_MONK/um_data_hires/da249_500ma_pr000.nc'
um_filepath_standard = '/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/UM_MONK/um_data_hires/diff_sgs/da249_500ma_pr000.nc'

# Directories for MONC data
monc_dir_conventional = '/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/UM_MONK/monc_data/sul_pat/'
monc_dir_standard = '/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/UM_MONK/monc_data/sul_pat/diff_sgs/'

# Director for DNS
flf_dir_500m = '/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/data/cbl_500m_flf/run000001columndiags.nc'

# Function to load UM data
def load_um_data(filepath):
    um_dataset = xr.open_dataset(filepath)
    wvar_um = um_dataset['STASH_m01s00i150']
    time_um = um_dataset['min15T0_0']
    return wvar_um, time_um

# Function to load MONC data
def load_monc_data(directory, expt, times, nx):
    ww_monc_all = []
    for time_index, time in enumerate(times):
        nc_file = f"{directory}/{expt}/cbl_{time}.nc"
        with xr.open_dataset(nc_file) as nc:
            w_data = nc['w']  # Shape (time, z, y, x)
            #w_mo = np.swapaxes(w_data, 1, 3)  # Swap axis to get shape (time, x, y, z)
            #w_mo = np.swapaxes(w_data, 1, 3)[0]
            # Calculate mean squared w for all time steps and average them
            ww_monc_all.append((w_data**2).mean(dim=['x', 'y'])) # Average 
            #ww_monc_all[time_index] = np.mean(np.square(w_mo), axis=(0))
    
    ww_monc_all = xr.concat(ww_monc_all, dim='time_series_300_300')
    return ww_monc_all


# Load UM data
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

# Convert cftime to datetime
def convert_to_datetime(cftime_obj):
    return datetime(cftime_obj.year, cftime_obj.month, cftime_obj.day, cftime_obj.hour, cftime_obj.minute, cftime_obj.second)

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

# Define MONC times and normalize
monc_times_seconds = np.array([600, 1200, 1800, 2400, 3000, 3600, 4200, 4800, 5400, 6000, 6600, 7200, 7800, 8400, 9000, 9600, 10200, 10800, 11400, 12000, 12600, 13200, 13800, 14400])
monc_times_hours = monc_times_seconds / 3600

# Load MONC data
resolution = '500m'  # Only using 500m resolution
nx = 128  # No. of points in x
ww_monc_conventional = load_monc_data(monc_dir_conventional, resolution, monc_times_seconds, nx)
ww_monc_standard = load_monc_data(monc_dir_standard, resolution, monc_times_seconds, nx)

# Normalize MONC data with wstar
ww_normalized_monc_conventional = ww_monc_conventional / (wstar ** 2)
ww_normalized_monc_standard = ww_monc_standard / (wstar ** 2)

# Resample UM data to match MONC time steps using interpolation
um_conventional_resampled = np.interp(monc_times_hours, time_um_conventional_hours, ww_normalized_um_conventional.max(axis=(1, 2, 3)))
um_standard_resampled = np.interp(monc_times_hours, time_um_standard_hours, ww_normalized_um_standard.max(axis=(1, 2, 3)))

# Extract maximum wvar for normalized MONC data

print(ww_normalized_monc_conventional.time_series_300_300)
max_ww_monc_conventional = np.max(ww_normalized_monc_conventional, axis=(1))
max_ww_monc_standard = np.max(ww_normalized_monc_standard, axis=(1))

# Define time points in hours for plotting
time_hours = monc_times_hours  # Time points in hours

ds_flf_500m = xr.open_dataset(flf_dir_500m)
wvar_flf_500m = ds_flf_500m['wvar']
time_flf_500m = ds_flf_500m['time']
wvar_flf_max_500m = wvar_flf_500m.max(dim='zp')
wvar_flf_max_500m_normalised = wvar_flf_max_500m / (wstar **2)


print(um_conventional_resampled[-1])
print(um_standard_resampled[-1])
print(max_ww_monc_conventional[-1])
print(max_ww_monc_standard[-1])


# Plot time series of max vertical velocity variance for UM and MONC data
plt.figure(figsize=(10, 6))

# Plot UM data
plt.plot(time_hours, um_conventional_resampled, color='blue', marker = '^', linestyle='-', label='UM Conventional 500m')
plt.plot(time_hours, um_standard_resampled, color='blue', marker = 'd', linestyle='--', label='UM Standard 500m')

# Plot MONC data
plt.plot(ww_normalized_monc_conventional.time_series_300_300/3600, max_ww_monc_conventional, color='red', marker = 'o', linestyle='-', label=f'MONC Conventional {resolution}')
plt.plot(ww_normalized_monc_standard.time_series_300_300/3600, max_ww_monc_standard, color='red', marker = 's', linestyle='--', label=f'MONC Standard {resolution}')

plt.plot(time_flf_500m/3600, wvar_flf_max_500m_normalised, color='green', label='ILES')
plt.xlabel("Time (hours)", fontsize=18)
plt.ylabel("$(\max \; w^{2})/(w^*)^2$", fontsize=18)
plt.ylim(top=0.8)
plt.title("Maximum Vertical Velocity Variance", fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
legend=plt.legend(loc='upper right', ncol=1, fontsize=18)
legend.get_frame().set_alpha(0.5)
#legend.get_frame().set_facecolor((0, 0, 1, 0.1))
plt.tight_layout()
plt.savefig('/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/UM_MONK/hires_plots/time_series_max_wvar_comparison.png', dpi=300)
plt.show()