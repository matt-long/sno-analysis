import os
from glob import glob
from subprocess import Popen, PIPE
import tarfile

from calendar import monthrange

import cftime
import xarray as xr


import util

grid_1x1 = util.generate_latlon_grid(nx=360, ny=180, lon0=-180.)


def _open_dataset(path, add_time=False):
    try:
        ds = xr.open_dataset(path, chunks={'time': 24})
        ds = ds.rename({
            "n_lon": "lon", 
            "n_lat": "lat", 
            "n_hour": "time", 
            "surface_area": "area"}
        )
        ds["lon"] = grid_1x1.lon
        ds["lat"] = grid_1x1.lon
        ds["time"] = [cftime.datetime(*d) for d in ds.itime]
        ds = ds.drop(["hourly_total_emission", "itime"])
    except Exception:
        ds = None
        raise
    return ds


def retrieve_datasets(date, version='2021'):
    assert version in ['2021']

    year, month = date.year, date.month
    
    _, nday = monthrange(year, month)
    
    cache_storage = f"/glade/scratch/{os.environ['USER']}/ODIAC_SFCO2_FFF_{version}"
    os.makedirs(cache_storage, exist_ok=True)


    asset_tar = f'{year:04d}{month:02d}.tar'
    assets_nc = [
        f'{cache_storage}/{year:04d}/{month:02d}/fossil_fuel_1x1_{year:04d}{month:02d}{day:02d}.nc'
        for day in range(1, nday+1)
    ]
    
    if not all([os.path.exists(p) for p in assets_nc]):
        cwd = os.getcwd()
        if version == '2021':
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