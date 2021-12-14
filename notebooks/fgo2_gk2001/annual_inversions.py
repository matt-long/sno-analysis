import os
import yaml

import numpy as np
from scipy import interpolate
import xarray as xr

from . import grid_data 

path_to_here = os.path.dirname(os.path.realpath(__file__))

Tmolyr_to_molyr = 1e12
Re = 6.37122e6 # m, radius of Earth


def _get_inversion_region_flux_table(product):
    """read database and return inverse fluxes by region"""
    with open(f"{path_to_here}/_o2-inversion-flux-tables.yml") as fid:
        return yaml.safe_load(fid)[product.lower()]["fluxes"]

    
def _get_inversion_region_flux_table_source(product):
    """read database and return inverse fluxes by region"""
    with open(f"{path_to_here}/_o2-inversion-flux-tables.yml") as fid:
        return yaml.safe_load(fid)[product.lower()]["reference_short"]
    
    
def open_inversion_annual_flux(product, gk_grid=False):
    """return inversion flux product"""
    source_info = _get_inversion_region_flux_table_source(product)
    region_flux = _get_inversion_region_flux_table(product)

    print(f"Global sum {product}: {sum(region_flux.values()):0.4f} Tmol/yr")
    
    ds = grid_data._get_inversion_regions(product)
    nlat, nlon = len(ds.lat), len(ds.lon)
    nrgn = len(ds.region)
       
    if gk_grid:        
        nlat, nlon = 160, 320
        dsGK01 = grid_data.generate_latlon_grid(nx=nlon, ny=nlat, lon0=-180.)
        dsGK01["region"] = ds.region
        
        x, y = np.meshgrid(ds.lon, ds.lat)
        xi, yi = np.meshgrid(dsGK01.lon, dsGK01.lat)        
        data = interpolate.griddata(
            (x.ravel(), y.ravel()), ds.REGION_MASK.values.ravel(), (xi, yi), 
            method='nearest',
        )
        dsGK01['REGION_MASK'] = xr.DataArray(data, dims=('lat','lon'))
        dsGK01['REGION_MASK_3D'] = xr.DataArray(
            np.zeros((nrgn, nlat, nlon)), dims=('region','lat','lon'),
        )
        for i in range(nrgn):
            dsGK01.REGION_MASK_3D.values[i, :, :] = np.where(
                dsGK01.REGION_MASK[:, :] == i + 1, 1, 0
            )    
        ds = dsGK01
        
    data = np.zeros((nlat, nlon))    
    for i, region_name in enumerate(ds.region.values):
        area = ds.area.where(ds.REGION_MASK_3D[i, :, :]==1).fillna(0.).sum()
        this_region_flux = region_flux[region_name] / area * Tmolyr_to_molyr
        data = np.where(
            ds.REGION_MASK_3D[i, :, :] == 1, 
            this_region_flux, 
            data,
        )

    ds["fgo2"] = xr.DataArray(
        data, 
        dims=('lat', 'lon'), 
        attrs={
            "long_name": f"Annual mean O2 flux ({source_info})",
            "units": "mol/m^2/yr",
        }
    )        
    ds['time'] = xr.DataArray(np.array([0.]), dims=('time'))    
    return ds[["fgo2", "area", "REGION_MASK", "REGION_MASK_3D"]]    
        
