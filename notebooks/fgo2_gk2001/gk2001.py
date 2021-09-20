import os
from collections import OrderedDict

import yaml
import csv

from netCDF4 import default_fillvals

import numpy as np
import xarray as xr

from . import grid_data

path_to_here = os.path.dirname(os.path.realpath(__file__))

# Data files available here: http://bluemoon.ucsd.edu/publications/ralph/airseaflux/Data_files
files = {
    "fgo2_ann": ("o2flux_ann_global.dat", "O2 flux (annual)"),
    "fgo2_sea": ("o2flux_sea_global.dat", "O2 flux (monthly anomaly)"),
    "fgo2_thm_ann": ("o2thermal_ann_global.dat", "O2 flux (annual, thermal component)"),
    "fgo2_thm_sea": ("o2thermal_sea_global.dat", "O2 flux (monthly anomaly, thermal component)"),
    "fgo2_bio_ann": ("o2flux_bioann_global.dat", "O2 flux (annual, biological component)"),
    "fgo2_bio_sea": ("o2flux_biosea_global.dat", "O2 flux (monthly anomaly, biological component)"),
}

droot = f"{path_to_here}/garcia_keeling_airseaflux"

dpm  = np.array([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
eom = np.cumsum(dpm)
bom = np.concatenate((np.array([0.]), eom[0:11]))
np.testing.assert_array_equal(eom - bom, dpm)


def dat2nc(dat_file_name, varname, long_name, shift_time=None, scaleby=None):
    """read Garcia & Keeling data file and convert to `xarray.Dataset`.
       Optionally apply temporal shift and flux scaling.
    """
    nt, ny, nx = 12, 160, 320
    dss = grid_data.generate_latlon_grid(nx=nx, ny=ny, lon0=-180.)
       
    time = np.vstack((bom, eom)).mean(axis=0)
    date = np.round(2000 * 10000 + np.arange(1,13,1) * 100. + dpm/2.)
    data = np.loadtxt(dat_file_name).ravel().reshape((nt, ny, nx))
    data = data / dpm[:, None, None] * 365.
    data[data==0.] = np.nan
    
    if shift_time is not None:
        time = time + shift_time
        date = date + shift_time
    
    if scaleby is not None:
        data = data * scaleby
            
    time = xr.DataArray(
        time, dims=("time"),
        attrs={"units": "day of year"},
    )    
    dss = dss.assign_coords(time=time)    
    dss["date"] = xr.DataArray(date, dims=("time"))        
    
    dss[varname] = xr.DataArray(
        data, 
        dims=("time","lat","lon"),        
        attrs={
            "long_name": long_name, 
            "units": "mol/m^2/yr",
            "note": f"GK2001 adjustments applied: time shifted = +{shift_time} days; scaleby = {scaleby}",
        }
    )
    dss[varname].encoding["_FillValue"] = default_fillvals["f8"]
    
    return dss

def open_flux_dataset(shift_time=10., scaleby=0.82, clobber=False):
    """open flux dataset"""
    ds = xr.Dataset()
    for v, (dat_file_in, long_name) in files.items():
        file_in = os.path.join(droot, dat_file_in)
        dsi = dat2nc(file_in, v, long_name, shift_time=shift_time, scaleby=scaleby)
        ds = xr.merge((ds, dsi))
        
    return ds