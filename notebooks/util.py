import os

import numpy as np
import xarray as xr

import matplotlib.pyplot as plt

# make variable so as to enable system dependence
catalog_json = '/glade/collections/cmip/catalog/intake-esm-datastore/catalogs/glade-cmip6.csv.gz'

# constants
mols_to_Tmolmon = 1e-12 * 86400. * 365. / 12.


def get_gridvar(df, source_id, variable_id):
    """get a grid variable from a source_id"""
    df_sub = df.loc[
        (df.source_id==source_id) & (df.variable_id==variable_id)
    ]   
    if len(df_sub) == 0:
        print(f'{source_id}: missing "{variable_id}"')
        return
    return xr.open_dataset(df_sub.iloc[0].path)


def open_cmip_dataset(df, source_id, variable_id, experiment_id,
                      time_slice=None, nmax_members=None):
    """return a dataset for a particular source_id, variable_id, experiment_id"""

    print('='*50) 
    print(f'{source_id}, {experiment_id}, {variable_id}')
    
    df_sub = df.loc[
        (df.source_id==source_id) 
        & (df.variable_id==variable_id) 
        & (df.experiment_id==experiment_id) 
    ]   
    if len(df_sub) == 0: 
        print('no data')
        return

    member_ids = sorted(df_sub.member_id.unique().tolist())
    print(f'\tfound {len(member_ids)} ensemble members')
    
    # optionally truncate ensemble list 
    if nmax_members is not None:
        if len(member_ids) > nmax_members:
            member_ids = member_ids[:nmax_members]
    
    print(f'\treading {len(member_ids)} members: {member_ids}')
    # loop over ensembles and read the datafiles
    ds_list = []
    for member_id in member_ids:
        paths = sorted(list(
            df_sub.loc[(df.member_id == member_id)].path
        ))
        dsi = xr.open_mfdataset(paths)
        if time_slice is not None:
            dsi = dsi.sel(time=time_slice)
        ds_list.append(dsi)    

    # concatenate along ensemble dimension
    print()
    return xr.concat(
        ds_list, 
        dim=xr.DataArray(member_ids, dims=('member_id'), name='member_id')
    )
    

def get_rmask_dict(grid, mask_definition, plot=False):
    """return a dictionary of masked area DataArray's"""
    # determine the latitude variable name
    lat_varname = None
    for lat_varname_test in ['latitude', 'lat', 'nav_lat']:
        if lat_varname_test in grid:
            lat_varname = lat_varname_test
            break

    if lat_varname is None:
        print('cannot determine latitude variable on this grid:')
        grid.info()
        raise ValueError('cannot determine lat_varname')
        
    # define region mask
    if mask_definition == 'SH_NH':
        rmasks = dict(
            NH=grid.areacello.where(grid[lat_varname] >= 20.).fillna(0.),
            SH=grid.areacello.where(grid[lat_varname] <= -20.).fillna(0.),
        ) 
    else:
        raise ValueError('unknown mask definition')
        
    if plot:
        nregion = len(rmasks.keys())
        nrow = int(np.sqrt(nregion))
        ncol = int(nregion/nrow) + min(1, nregion%nrow)
        figsize=(6, 4)
        fig, axs = plt.subplots(
            nrow, ncol, 
            figsize=(figsize[0]*ncol, figsize[1]*nrow),
            squeeze=False)
        for n, (key, rmask) in enumerate(rmasks.items()):
            i, j = np.unravel_index(n, axs.shape)
            ax = axs[i, j]
            rmask.plot(ax=ax)
            ax.set_title(key)
            if 'source_id' in grid.attrs:
                plt.suptitle(grid.attrs['source_id'])
    
    return rmasks


def compute_regional_integral(ds, variable_id, rmasks):
    """return a DataArray of the regional integral of ds[variable_id]"""
    if variable_id == 'fgo2':
        assert (ds[variable_id].attrs['units'] == 'mol m-2 s-1')
        convert = mols_to_Tmolmon
        units_integral = 'Tmol O$_2$ month$^{-1}$'
    else:
        raise NotImplementedError(f'add "{variable_id}" to integral definitions')

    da_list = []
    regions = []
    
    dims_lateral = tuple(d for d in ds[variable_id].dims if d not in ['time', 'member_id'])
    
    for key, rmask in rmasks.items():        
        assert rmask.dims == dims_lateral, 'dimension mismatch on region mask'
        da = ((-1.0) * ds[variable_id] * rmask).sum(dims_lateral) * convert
        da.attrs['units'] = units_integral
        da_list.append(da)
        regions.append(key)    
  
    var =  xr.concat(
        da_list, 
        dim=xr.DataArray(regions, dims=('region'), name='region'),
    )
    var.name = variable_id
    return var
