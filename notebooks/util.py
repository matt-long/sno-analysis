import os

import numpy as np
import xarray as xr

import intake

import matplotlib.pyplot as plt

# make variable so as to enable system dependence
catalog_csv = '/glade/collections/cmip/catalog/intake-esm-datastore/catalogs/glade-cmip6.csv.gz'
catalog_json = '/glade/collections/cmip/catalog/intake-esm-datastore/catalogs/glade-cmip6.json'

cmip6_catalog = intake.open_esm_datastore(catalog_json)

# constants
T0_Kelvin = 273.15
mols_to_Tmolmon = 1e-12 * 86400. * 365. / 12.
µmolkg_to_mmolm3 = 1026. / 1000. # for volume conserving models, makes sense to use constant density


def get_gridvar(df, source_id, variable_id):
    """get a grid variable from a source_id"""
    df_sub = df.loc[
        (df.source_id==source_id) 
        & (df.variable_id==variable_id)
        & (df.grid_label == 'gn')
    ]   
    if len(df_sub) == 0:
        print(f'{source_id}: missing "{variable_id}"')
        return
    return xr.open_dataset(df_sub.iloc[0].path)


def open_cmip_dataset(source_id, variable_id, experiment_id, table_id,
                      time_slice=None, nmax_members=None):
    """return a dataset for a particular source_id, variable_id, experiment_id"""

    print('='*50) 
    print(f'{source_id}, {experiment_id}, {variable_id}')
    
    cat_sub = cmip6_catalog.search(
        source_id=source_id, 
        variable_id=variable_id, 
        experiment_id=experiment_id, 
        table_id=table_id,
        grid_label='gn',
    )
    df_sub = cat_sub.df    
    
    dsets_dict = cat_sub.to_dataset_dict()    
    if len(dsets_dict) == 0: 
        print('no data')
        return
    
    assert len(dsets_dict) == 1, 'expecting single dataset'    
    key, ds = dsets_dict.popitem()
        
    member_ids = sorted(df_sub.member_id.unique().tolist())
    print(f'\tfound {len(member_ids)} ensemble members')
    
    # optionally truncate ensemble list 
    if nmax_members is not None:
        if len(member_ids) > nmax_members:
            ds = ds.sel(member_id=member_ids[:nmax_members])
        
    if time_slice is not None:
        ds = ds.sel(time=time_slice)
        
    return ds


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
    if mask_definition == 'SET_NET':
        rmasks = dict(
            NET=grid.areacello.where(grid[lat_varname] >= 20.).fillna(0.),
            SET=grid.areacello.where(grid[lat_varname] <= -20.).fillna(0.),
        )
    elif mask_definition == 'SH_NH':
        rmasks = dict(
            NH=grid.areacello.where(grid[lat_varname] >= 0.).fillna(0.),
            SH=grid.areacello.where(grid[lat_varname] < 0.).fillna(0.),
        ) 
    elif mask_definition == 'global':
        rmasks = dict(
            NH=grid.areacello.fillna(0.),
            SH=grid.areacello.fillna(0.),
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


def compute_regional_integral(ds, variable_id, rmasks, flipsign=False):
    """return a DataArray of the regional integral of ds[variable_id]"""
    if variable_id == 'fgo2':
        assert (ds[variable_id].attrs['units'] == 'mol m-2 s-1')
        convert = (-1.0) * mols_to_Tmolmon
        units_integral = 'Tmol O$_2$ month$^{-1}$'
        long_name = 'O$_2$ flux'
        
    elif variable_id == 'fgco2':
        #print(ds.attrs['units'])
        #assert (ds[variable_id].attrs['units'] == 'mol m-2 s-1') ## crashing so temporarily commenting out - print says no 'units' attribute
        convert = (-1.0) * mols_to_Tmolmon
        units_integral = 'Tmol CO$_2$ month$^{-1}$'
        long_name = 'CO$_2$ flux'
    else:
        raise NotImplementedError(f'add "{variable_id}" to integral definitions')

    
    if flipsign:
        convert = convert * (-1.0)
    
    da_list = []
    regions = []
    
    dims_lateral = tuple(d for d in ds[variable_id].dims if d not in ['time', 'member_id'])
    
    for key, rmask in rmasks.items():        
        assert rmask.dims == dims_lateral, 'dimension mismatch on region mask'
        da = (ds[variable_id] * rmask).sum(dims_lateral) * convert
        da.attrs['units'] = units_integral
        da.attrs['long_name'] = long_name
        da_list.append(da)
        regions.append(key)    
  
    var =  xr.concat(
        da_list, 
        dim=xr.DataArray(regions, dims=('region'), name='region'),
    )
    var.name = variable_id
    return var


def compute_fgn2(ds):
    """compute N2 from heat flux""" 
    raise NotImplementedError('fgn2 not implemented')
    return ds
    
    
def compute_fgo2_thermal(ds):
    """compute thermal O2 flux"""   
    raise NotImplementedError('fgo2_thermal not implemented')    
    return ds


def O2sol(S, T):
    """
    Solubility of O2 in sea water
    INPUT:
    S = salinity    [PSS]
    T = temperature [degree C]

    conc = solubility of O2 [µmol/kg]

    REFERENCE:
    Hernan E. Garcia and Louis I. Gordon, 1992.
    "Oxygen solubility in seawater: Better fitting equations"
    Limnology and Oceanography, 37, pp. 1307-1312.
    """

    # constants from Table 4 of Hamme and Emerson 2004
    return _garcia_gordon_polynomial(S, T,
                                     A0=5.80871,
                                     A1=3.20291,
                                     A2=4.17887,
                                     A3=5.10006,
                                     A4=-9.86643e-2,
                                     A5=3.80369,
                                     B0=-7.01577e-3,
                                     B1=-7.70028e-3,
                                     B2=-1.13864e-2,
                                     B3=-9.51519e-3,
                                     C0=-2.75915e-7)


def _garcia_gordon_polynomial(S, T,
                              A0=0., A1 = 0., A2 = 0., A3=0., A4=0., A5=0.,
                              B0=0., B1=0., B2=0., B3=0.,
                              C0=0.):

    T_scaled = np.log((298.15 - T) /(T0_Kelvin + T))
    return np.exp(A0 + A1*T_scaled + A2*T_scaled**2. + A3*T_scaled**3. + A4*T_scaled**4. + A5*T_scaled**5. + \
                  S*(B0 + B1*T_scaled + B2*T_scaled**2. + B3*T_scaled**3.) + C0 * S**2.)

