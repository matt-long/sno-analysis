import os
from subprocess import Popen, PIPE
import tarfile

from calendar import monthrange

import cftime

from netCDF4 import default_fillvals

import numpy as np
import xarray as xr

import util

grid_1x1 = util.generate_latlon_grid(nx=360, ny=180, lon0=-180.)

mwC = 12.01

def _open_dataset(path, add_time=False):
    try:
        units_time = 'days since 2000-01-01 00:00:00'
        ds = xr.open_dataset(path, chunks={'time': 24}).rename({'n_lon': 'lon', 'n_lat': 'lat'})
        
        year, month, day = ds.itime.values[0, 0], ds.itime.values[0, 1], ds.itime.values[0, 2]
        
        dso = xr.Dataset({'lon': grid_1x1.lon, 'lat': grid_1x1.lat, 'area': ds.surface_area})
       
        time = xr.DataArray([cftime.datetime(year, month, day, 12)], dims=('time'))
        time.attrs['bounds'] = 'time_bnds'      
        time.encoding['units'] = units_time
        time.encoding['dtype'] = np.float64    
        time.encoding['_FillValue'] = None        
        
        dso['SFCO2_FF'] = ds.emission.sum('n_hour').expand_dims('time') / mwC / 86400. / 24
        dso.SFCO2_FF.attrs = {k: v for k, v in ds.emission.attrs.items() if k != 'unit'}
        dso.SFCO2_FF.attrs['long_name'] = 'Fossil fuel flux'
        dso.SFCO2_FF.attrs['units'] = 'mol/m^2/s'                    
        dso.SFCO2_FF.attrs['cell_methods'] = 'time: sum'                            
        dso.SFCO2_FF.encoding['_FillValue'] = default_fillvals['f8']
        
        t0 = cftime.date2num(cftime.datetime(year, month, day, 0), units_time)
        time_bounds_data = cftime.num2date(np.array([[t0, t0+1]]), units_time)        
        dso['time_bnds'] = xr.DataArray(
            time_bounds_data,
            dims=('time', 'd2'), 
            name='time_bnds',
        )
        dso.time_bnds.attrs['long_name'] = 'Bounds for time integration'
        dso.time_bnds.encoding['dtype'] = np.float64
        dso.time_bnds.encoding['_FillValue'] = None
        
        dso['time'] = time
        dso = dso.set_coords('time')  
        
        dso.area.attrs['units'] = 'm^2'
        dso.area.attrs['long_name'] = 'Grid cell area'
        dso.area.encoding['_FillValue'] = None
        del dso.area.attrs['unit']
        
        dso.lat.attrs = grid_1x1.lat.attrs
        dso.lon.attrs = grid_1x1.lon.attrs       
        dso.lat.encoding['_FillValue'] = None
        dso.lon.encoding['_FillValue'] = None                    

    except Exception:
        dso = None
        raise
        
    return dso


def retrieve_datasets(date, version='v2020.1'):
    assert version in ['v2020.1']

    year, month = date.year, date.month
    
    _, nday = monthrange(year, month)
    
    cache_storage = f"/glade/scratch/{os.environ['USER']}/OCO2-MIP-SFCO2_FF_{version}"
    os.makedirs(cache_storage, exist_ok=True)


    asset_tar = f'{year:04d}{month:02d}.tar'
    assets_nc = [
        f'{cache_storage}/{year:04d}/{month:02d}/fossil_fuel_1x1_{year:04d}{month:02d}{day:02d}.nc'
        for day in range(1, nday+1)
    ]
    
    if not all([os.path.exists(p) for p in assets_nc]):
        cwd = os.getcwd()
        if version == 'v2020.1':
            url = f'https://zenodo.org/record/4776925/files/{asset_tar}'

        os.chdir(cache_storage)

        p = Popen(['wget', url], stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()
        if p.returncode:    
            print(stderr.decode('UTF-8'))
            print(stdout.decode('UTF-8'))
            raise OSError('data transfer failed')    

        # untar archive
        assert os.path.isfile(asset_tar), f'missing {asset_tar}'
        tar = tarfile.open(asset_tar, "r")
        tar.extractall()
        tar.close()
        
        os.remove(asset_tar)
        
        os.chdir(cwd)

    return [_open_dataset(path) for path in assets_nc]