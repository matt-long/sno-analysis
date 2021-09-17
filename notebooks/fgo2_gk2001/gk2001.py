import os
from collections import OrderedDict

import yaml
import csv

from netCDF4 import default_fillvals

import numpy as np
import xarray as xr

path_to_here = os.path.dirname(os.path.realpath(__file__))

Re = 6.37122e6 # m, radius of Earth

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
bom = np.concatenate((np.array([1]), eom[0:11] + 1))


def generate_latlon_grid(nx, ny, lon0=0., include_vertices=False):
    """
    Generate a regular lat,lon grid file.

    Parameters
    ----------
    nx: number of x points.
    ny: number of y points.

    Returns
    -------

    ds: dataset with grid variables
    """

    deg2rad = np.pi / 180.

    dx = 360. / nx
    dy = 180. / ny

    lat = np.arange(-90.+dy/2.,90.,dy)
    lon = np.arange(lon0+dx/2.,lon0+360.,dx)

    yc = np.broadcast_to(lat[:,None],(ny,nx))
    xc = np.broadcast_to(lon[None,:],(ny,nx))

    yv = np.stack((yc-dy/2.,yc-dy/2.,yc+dy/2.,yc+dy/2.),axis=2)
    xv = np.stack((xc-dx/2.,xc+dx/2.,xc+dx/2.,xc-dx/2.),axis=2)

    y0 = np.sin(yv[:,:,0]*deg2rad) # south
    y1 = np.sin(yv[:,:,3]*deg2rad) # north
    x0 = xv[:,:,0]*deg2rad         # west
    x1 = xv[:,:,1]*deg2rad         # east
    area = (y1 - y0) * (x1 - x0) * Re ** 2.

    ds = xr.Dataset(
        {"lat": xr.DataArray(lat,dims=("lat"),
                            attrs={"units":"degrees_north",
                                   "long_name":"Latitude"}),
         "lon": xr.DataArray(lon,dims=("lon"),
                            attrs={"units":"degrees_east",
                            "long_name":"Longitude"})})
    
    if include_vertices:
        ds["xc"] = xr.DataArray(xc,dims=("lat","lon"),
                                attrs={"units":"degrees_east",
                                       "long_name":"longitude of cell centers"})

        ds["yc"] = xr.DataArray(yc,dims=("lat","lon"),
                                attrs={"units":"degrees_north",
                                       "long_name":"latitude of cell centers"})

        ds["xv"] = xr.DataArray(xv,dims=("lat","lon","nv"),
                                attrs={"units":"degrees_east",
                                       "long_name":"longitude of cell corners"})

        ds["yv"] = xr.DataArray(yv,dims=("lat","lon","nv"),
                                attrs={"units":"degrees_north",
                                       "long_name":"latitude of cell corners"})

    ds["area"] = xr.DataArray(area,dims=("lat","lon"),
                              attrs={"units":"m^2",
                                     "long_name":"area"})

    return ds


def dat2nc(dat_file_name, varname, long_name, shift_time=None, scaleby=None):
    """read Garcia & Keeling data file and convert to `xarray.Dataset`.
       Optionally apply temporal shift and flux scaling.
    """
    nt, ny, nx = 12, 160, 320
    dss = generate_latlon_grid(nx=nx, ny=ny, lon0=-180.)
       
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
            "_FillValue": default_fillvals["f8"],
            "note": f"shift_time = {shift_time} days; scaleby = {scaleby}",
        }
    )
    
    return dss

def open_flux_dataset(shift_time=10., scaleby=0.82, clobber=False):
    """open flux dataset"""
    
    cache_file = f"{path_to_here}/garcia_keeling_airseaflux/GK2001.fluxes.dt={int(shift_time)}.scale={scaleby}.nc"
    if os.path.exists(cache_file) and not clobber:
        with xr.open_dataset(cache_file) as ds:
            return ds
    
    ds = xr.Dataset()
    for v, (dat_file_in, long_name) in files.items():
        file_in = os.path.join(droot, dat_file_in)
        dsi = dat2nc(file_in, v, long_name, shift_time=shift_time, scaleby=scaleby)
        ds = xr.merge((ds, dsi))
        
    ds.to_netcdf(cache_file)
    return ds