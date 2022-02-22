import os
from subprocess import Popen, PIPE, check_call
import tarfile

from calendar import monthrange

import cftime

from netCDF4 import default_fillvals

import numpy as np
import xarray as xr

import config
import regrid_tools
import util

src_grid = regrid_tools.grid('latlon', nx=3600, ny=1800, lon0=-180.0)
dst_grid = regrid_tools.grid("latlon", **config.config_dict["flux-dst-grid-kwargs"])

regrid_obj = regrid_tools.regridder(src_grid, dst_grid, method="conserve", clobber=False)

mwCO2 = 44.0 # g-CO2 / mol-CO2
mwO2 = 32.0 # g-O2 / mol-O2

kgCO2mon_to_molCO2m2s = 1e3 / mwCO2 * 12 / 365 / 86400
kgO2mon_to_molO2m2s = 1e3 / mwO2 * 12 / 365 / 86400

def _open_dataset(path):
    try:
        with xr.open_dataset(path, decode_times=False) as ds:
            
            units_time = ds.time.attrs['units']
            area = util.compute_grid_area(ds)

            for group in ['CO2', 'O2']:
                assert group in ['CO2', 'O2']

                with xr.open_dataset(path, group=group) as dsg:
                    data_vars = list(dsg.data_vars)
                    if group == 'CO2':
                        assert set(data_vars) == {'OIL', 'GAS', 'COAL', 'CEMENT', 'BUNKER'}
                    elif group == 'O2':
                        assert set(data_vars) == {'OIL', 'GAS', 'COAL', 'BUNKER'}
                        
                    # compute sum across flux components
                    sum_vars = xr.full_like(dsg[data_vars[0]], fill_value=0.0)
                    for v in dsg.data_vars:
                        sum_vars += dsg[v]

                    # apply unit conversion
                    units = sum_vars.attrs['units']
                    if group == 'CO2':
                        assert units == 'kg CO2 month-1'
                        sum_vars *= kgCO2mon_to_molCO2m2s
                        long_name = 'Total CO2 emissions from: ' + ', '.join(data_vars)
                    elif group == 'O2': 
                        assert units == 'kg O2 month-1'
                        sum_vars *= kgO2mon_to_molO2m2s
                        long_name = 'Total O2 surface flux from: ' + ', '.join(data_vars)

                    # normalize by grid cell area
                    sum_vars /= area

                    # fix attrs
                    sum_vars.attrs['long_name'] = long_name
                    sum_vars.attrs['units'] = 'mol/m^2/s'

                # assignment
                ds[f'SF{group}_FF'] = sum_vars                         
        
        dso = util.generate_latlon_grid(**config.config_dict["flux-dst-grid-kwargs"])[["area"]]        
        dso['area'].encoding['_FillValue'] = None
        
        dso_dst_data = regrid_obj(ds)
        for v in dso_dst_data.data_vars:
            dso[v] = dso_dst_data[v]
            dso[v].encoding = ds[v].encoding
            
        assert ds.sizes['time'] == 12
        year = np.array([d.year for d in cftime.num2date(ds.time.values, units=units_time)])
        print(year)
        assert all(year[0] == year)
        
        time, time_bnds = util.gen_midmonth_cftime_coord((year[0], year[0]))
        dso['time'] = time
        dso['time_bnds'] = time_bnds

        dso.attrs = ds.attrs
        dso.attrs['history'] += '; aggegrated and regridded to 1x1 by Matt Long (NCAR)'
        
    except Exception:
        dso = None
        raise
        
    return dso


def retrieve_datasets(date, version='v2021.3'):
    assert version in ['v2021.3']

    year, month = date.year, date.month
    
    _, nday = monthrange(year, month)
    
    cache_storage = f"/glade/scratch/{os.environ['USER']}/GCP-GridFED-SFCO2_FFF_{version}"
    os.makedirs(cache_storage, exist_ok=True)


    asset_zip = f'GridFED{version}_{year:04d}.zip'
    asset_nc = f'{cache_storage}/GCP_Global_{year:04d}.nc'
    
    if not os.path.exists(asset_nc):
        cwd = os.getcwd()
        if version == 'v2021.3':
            url = f'https://zenodo.org/record/5956612/files/{asset_zip}'

        os.chdir(cache_storage)

        p = Popen(['wget', url], stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()
        if p.returncode:    
            print(stderr.decode('UTF-8'))
            print(stdout.decode('UTF-8'))
            raise OSError(f'data transfer failed: {url}')    

        # unzip archive
        assert os.path.isfile(asset_zip), f'missing {asset_zip}'
        check_call(['unzip', asset_zip])
        os.remove(asset_zip)
        
        os.chdir(cwd)

    return _open_dataset(asset_nc)