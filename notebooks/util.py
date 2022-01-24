import os
import sys
import subprocess

import warnings
import yaml

from toolz import curry

from datetime import datetime, timezone
import cftime

import numpy as np
import xarray as xr
import pandas as pd

from netCDF4 import default_fillvals

import intake

import matplotlib.pyplot as plt

import pop_tools

import dask
from dask_jobqueue import PBSCluster
from dask.distributed import Client

import solubility



# make variable so as to enable system dependence
catalog_csv = '/glade/collections/cmip/catalog/intake-esm-datastore/catalogs/glade-cmip6.csv.gz'
catalog_json = '/glade/collections/cmip/catalog/intake-esm-datastore/catalogs/glade-cmip6.json'

cache_catalog_csv = "catalogs/cmip6-process-cache.csv.gz" 
cache_catalog_json = "catalogs/cmip6-process-cache.json" 

cmip6_catalog = intake.open_esm_datastore(catalog_json)

# hardwired to only use native grid
grid_label = 'gn'

# constants
mols_to_Tmolmon = 1e-12 * 86400. * 365. / 12.
µmolkg_to_mmolm3 = 1026. / 1000. # for volume conserving models, makes sense to use constant density
kgCO2s_to_Tmolmon = 1000. / 12. * mols_to_Tmolmon
W_to_PW = 1. / 1E15
Re = 6.37122e6 # m, radius of Earth

X_O2 = 0.2094 
X_N2 = 0.7808


project = "NEOL0004"


git_repo = (subprocess
            .check_output(['git', 'config', '--get', 'remote.origin.url'])
            .strip()
            .decode("utf-8")
            .replace('git@github.com:', 'https://github.com/')
            .replace('.git', '')            
           )


def get_ClusterClient(memory="25GB"):
    """return client and cluster"""
    USER = os.environ['USER']
    
    cluster = PBSCluster(
        cores=1,
        memory=memory,
        processes=1,
        queue='casper',
        local_directory=f'/glade/scratch/{USER}/dask-workers',
        log_directory=f'/glade/scratch/{USER}/dask-workers',
        resource_spec='select=1:ncpus=1:mem=25GB',
        project=project,
        walltime='06:00:00',
        interface='ib0',)

    jupyterhub_server_name = os.environ.get('JUPYTERHUB_SERVER_NAME', None)    
    dashboard_link = 'https://jupyterhub.hpc.ucar.edu/stable/user/{USER}/proxy/{port}/status'    
    if jupyterhub_server_name:
        dashboard_link = (
            'https://jupyterhub.hpc.ucar.edu/stable/user/'
            + '{USER}'
            + f'/{jupyterhub_server_name}/proxy/'
            + '{port}/status'
        )
    dask.config.set({'distributed.dashboard.link': dashboard_link})        
    client = Client(cluster)
    return cluster, client



class track_attrs(object):
    """object for tracking variable attributes"""    
    def __init__(self):
        self._variable_attrs_file = "data/cache/cmip/cmip6.variable_attrs.yml"
        if os.path.exists(self._variable_attrs_file):
            with open(self._variable_attrs_file, "r") as fid:
                self._variable_attrs = yaml.safe_load(fid)
        else:        
            self._variable_attrs = {}

        self._attrs_keys = [
            "long_name", "units", "description", "cell_methods", "cell_measures"
        ]
       
    def update_attrs(self, variable_id, attrs, clobber=False):
        if not variable_id in self._variable_attrs or clobber:
            self._variable_attrs[variable_id] = {
                k: attrs[k] for k in self._attrs_keys if k in attrs
            }
            self.persist()
        
    def persist(self):
        with open(self._variable_attrs_file, "w") as fid:
            yaml.dump(self._variable_attrs, fid)
            
    def __getitem__(self, key):
        return self._variable_attrs[key]
    
    
class missing_data_tracker(object):        
    def __init__(self): 
        """construct object for tracking missing data"""
        self._columns = ['source_id', 'experiment_id', 'table_id', 'variable_id',  'grid_label']
        
        self.missing_data_file = 'data/cache/cmip/missing-data.csv'            
        if os.path.exists(self.missing_data_file):
            self._df = pd.read_csv(self.missing_data_file)
        else:
            self._init_df()
    
    @curry
    def ismissing(self, source_id, experiment_id, table_id, variable_id):
        """determine whether a particular query is missing"""
        return (source_id, experiment_id, table_id, variable_id, grid_label) in self._df.set_index(self._columns).index
        
    def set_missing(self, source_id, experiment_id, table_id, variable_id):
        """set a particular query to missing"""
        if self.ismissing(source_id, experiment_id, table_id, variable_id):
            print('already missing')
            return

        df_new = pd.DataFrame(
            dict(source_id=[source_id], 
                 experiment_id=[experiment_id], 
                 table_id=[table_id], 
                 variable_id=[variable_id], 
                 grid_label=[grid_label]))
        self._df = pd.concat((self._df, df_new), ignore_index=True) #.reset_index()
    
    def persist(self):
        """persist missing data dataframe"""
        self._df.to_csv(self.missing_data_file, index=False)
    
    def _init_df(self):
        self._df = pd.DataFrame({k: [] for k in self._columns})         
        
    def clobber(self):
        self._init_df()
        if os.path.exists(self.missing_data_file):
            os.remove(self.missing_data_file)
    

    
def id_and_search_vars(variable_name):
    if ':' in variable_name:
        search_vars = variable_name.split(':')[-1].split(',')
        variable_id = variable_name.split(':')[0]        
    else:
        search_vars = [variable_name]
        variable_id = variable_name
    return variable_id, search_vars
        
    
def get_gridvar(df, source_id, variable_id):
    """get a grid variable from a source_id"""
    df_sub = df.loc[
        (df.source_id==source_id) 
        & (df.variable_id==variable_id)
        & (df.grid_label == grid_label)
    ]   
    if len(df_sub) == 0:
        if "CESM2" in source_id:
            return pop_tools.get_grid("POP_gx1v7")[["TAREA"]].rename({"TAREA": "areacello"})
        else:
            print(f'{source_id}: missing "{variable_id}"')
            return
    return xr.open_dataset(df_sub.iloc[0].path)


def open_cmip_cached(operator_applied, region_mask):
    variable_attrs = track_attrs()
    cat = intake.open_esm_datastore(cache_catalog_json)
    
    dsets = cat.search(
        operator_applied=operator_applied,
        region_mask=region_mask,
    ).to_dataset_dict()

    for key, ds in dsets.items():
        for var_id in ds.data_vars:
            if not ds[var_id].attrs:
                ds[var_id].attrs = variable_attrs[var_id]
    return dsets



def get_member_id_list(source_id, experiment_id, table_id, require_vars):
    
    for i, variable_id in enumerate(require_vars):
        cat_sub = cmip6_catalog.search(
            source_id=source_id, 
            variable_id=variable_id, 
            experiment_id=experiment_id, 
            table_id=table_id,
            grid_label='gn',
        )
        member_ids_var_i = set(cat_sub.df.member_id.unique())
        member_ids = member_ids_var_i if i == 0 else member_ids.intersection(member_ids_var_i)
    
    return list(member_ids)


def open_cmip_dataset(source_id, variable_id, experiment_id, table_id, member_id,
                      time_slice=None, nmax_members=None, aggregate=True,
                      preprocess=None,):
    """return a dataset for a particular source_id, variable_id, experiment_id"""

    print('='*50) 
    print(f'{source_id}, {experiment_id}, {variable_id}')
    
    if source_id == "UKESM1-0-LL":
        cat_sub = cmip6_catalog.search(
            institution_id="MOHC",
            source_id=source_id, 
            variable_id=variable_id, 
            experiment_id=experiment_id, 
            table_id=table_id,
            member_id=member_id,
            grid_label=grid_label,
        )
    else:
        cat_sub = cmip6_catalog.search(
            source_id=source_id, 
            variable_id=variable_id, 
            experiment_id=experiment_id, 
            table_id=table_id,
            member_id=member_id,            
            grid_label=grid_label,
        )
        
    df_sub = cat_sub.df    
    
    with dask.config.set(**{"array.slicing.split_large_chunks": True}):    
        dsets_dict = cat_sub.to_dataset_dict(
            cdf_kwargs=dict(use_cftime=True), 
            preprocess=preprocess,
            aggregate=aggregate,
        )
        
    if len(dsets_dict) == 0: 
        print('no data')
        return
    
    if not aggregate:
        return dsets_dict
    
    assert len(dsets_dict) == 1, (
        f'expecting single dataset; got:\n {dsets_dict}\n\n' + 
        f'{df_sub.path.to_list()}'
    )
    key, ds = dsets_dict.popitem()
        
    member_ids = sorted(df_sub.member_id.unique().tolist())
    print(f'\tfound {len(member_ids)} ensemble members')
    
    # optionally truncate ensemble list 
    if nmax_members is not None:                     
        if len(member_ids) > nmax_members:
            
            # sort ensemble members to prioritize physics and forcing 1 and 
            # to sort realizations numerically 
            real = np.array(member_ids)
            init = np.array(member_ids)
            phys = np.array(member_ids)
            forc = np.array(member_ids)
            for i in range(len(member_ids)):
                real[i] = member_ids[i].split('r')[1].split('i')[0]
                init[i] = member_ids[i].split('i')[1].split('p')[0]
                phys[i] = member_ids[i].split('p')[1].split('f')[0]
                forc[i] = member_ids[i].split('f')[1]
            real = real.astype(int)
            init = init.astype(int)
            phys = phys.astype(int)
            forc = forc.astype(int)      
            member_ids_sorted = np.array(member_ids)[np.lexsort((real, init, phys, forc))]
            
            ds = ds.sel(member_id=member_ids_sorted[:nmax_members])
        
    if time_slice is not None:       
        ds = ds.sel(time=_time_slice_cftime(time_slice, _get_calendar(ds)))
        
    return ds


def _get_calendar(ds):
    if "calendar" in ds.time.encoding:
        return ds.time.encoding["calendar"]
    elif "calendar" in ds.time.attrs:
        return ds.time.attrs["calendar"]
    else:
        raise ValueError("cannot determine calendar")

def _time_slice_cftime(time_slice, calendar):
    """temporary workaround for bug in xarray-pandas-cftime
    See here: https://zulip2.cloud.ucar.edu/#narrow/stream/10-python-questions/topic/datetime.20index/near/44187
    """
    assert len(time_slice.start) == 4
    assert len(time_slice.stop) == 4
    y1, y2 = int(time_slice.start), int(time_slice.stop)
    dec31 = 30 if calendar == "360_day" else 31
    return slice(
        cftime.datetime(y1, 1, 1, calendar=calendar), 
        cftime.datetime(y2, 12, dec31, calendar=calendar)
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
    if mask_definition == 'SET_NET':
        rmasks = dict(
            NET=grid.areacello.where(grid[lat_varname] >= 20.).fillna(0.),
            SET=grid.areacello.where(grid[lat_varname] <= -20.).fillna(0.),
        )
    elif mask_definition == 'SHL_NHL':
        rmasks = dict(
            NHL=grid.areacello.where(grid[lat_varname] >= 45.).fillna(0.),
            SHL=grid.areacello.where(grid[lat_varname] < -45.).fillna(0.),
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


def compute_regional_integral(ds, variable_id, rmasks):
    """return a DataArray of the regional integral of ds[variable_id]"""
    if variable_id == 'fgo2':
        assert (ds[variable_id].attrs['units'] == 'mol m-2 s-1')
        convert = (-1.0) * mols_to_Tmolmon
        units_integral = 'Tmol O$_2$ month$^{-1}$'
        long_name = 'O$_2$ flux'
        sumormean = 'sum'
        
    elif variable_id == 'fgco2':
        assert (ds[variable_id].attrs['units'] == 'kg m-2 s-1')
        convert = (-1.0) * kgCO2s_to_Tmolmon
        units_integral = 'Tmol CO$_2$ month$^{-1}$'
        long_name = 'CO$_2$ flux'
        sumormean = 'sum'
        
    elif variable_id == 'hfds':
        assert (ds[variable_id].attrs['units'] == 'W m-2')
        convert = 1.0 * W_to_PW
        units_integral = 'PW' # convert to 'PJ/month'?
        long_name = 'Heat flux'
        sumormean = 'sum'

    elif variable_id == 'tos':
        assert (ds[variable_id].attrs['units'] == 'degC')
        convert = 1.0
        units_integral = '$^\circ$C'
        long_name = 'Surface temperature'
        sumormean = 'mean'

    elif variable_id == 'sos':
        assert (ds[variable_id].attrs['units'] == '0.001')
        convert = 1.0
        units_integral = 'PSU'
        long_name = 'Surface salinity'
        sumormean = 'mean'
        
    elif variable_id == 'intpp':
        #print(ds[variable_id].attrs['units'])
        assert (ds[variable_id].attrs['units'] == 'mol m-2 s-1')
        convert = 1.0 * mols_to_Tmolmon
        units_integral = 'Tmol C month$^{-1}$' # or say "CO2"?
        long_name = 'NPP'
        sumormean = 'sum'
        
    elif variable_id == 'fgn2':
        convert = 1.0 * mols_to_Tmolmon
        units_integral = 'Tmol N$_2$ month$^{-1}$'
        long_name = 'N$_2$ flux'
        sumormean = 'sum'
        
    elif variable_id == 'fgo2_thermal':
        convert = 1.0 * mols_to_Tmolmon
        units_integral = 'Tmol O$_2$ month$^{-1}$'
        long_name = 'Thermal O$_2$ flux'
        sumormean = 'sum'
        
    elif variable_id == 'fbddtdic':
        convert = 1.0 * mols_to_Tmolmon
        units_integral = 'Tmol C month$^{-1}$' # or say "CO2"?
        long_name = 'Bio DIC change' # Rate of Change of Dissolved Inorganic Carbon Due to Biological Activity
        sumormean = 'sum'
        
    elif variable_id == 'epc100':
        convert = 1.0 * mols_to_Tmolmon
        units_integral = 'Tmol C month$^{-1}$' # or say "CO2"?
        long_name = 'POC export' # Downward Flux of Particulate Organic Carbon
        sumormean = 'sum'
        
    elif variable_id == 'fgapo':
        convert = (-1.0) * mols_to_Tmolmon
        units_integral = 'Tmol APO month$^{-1}$'
        long_name = 'APO flux'
        sumormean = 'sum'
        
    else:
        raise NotImplementedError(f'add "{variable_id}" to integral definitions')

       
    da_list = []
    regions = []
    
    dims_lateral = tuple(d for d in ds[variable_id].dims if d not in ['time', 'member_id'])

    for key, rmask in rmasks.items():
        
        assert rmask.dims == dims_lateral, 'dimension mismatch on region mask'
        
        # ensure that missing values are propagated to the mask
        rmask_v = rmask.where(ds[variable_id].notnull()).fillna(0.)
        
        da = (ds[variable_id] * rmask_v).sum(dims_lateral) * convert

        if sumormean == 'mean':
            da = da / rmask_v.sum(dims_lateral)

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


def compute_fgn2(ds, scaleby=1.):
    """
    compute N2 flux from heat flux and temperature derivative of solubility
    
    using Eq. 2 from Keeling and Shertz, 1992 (and Eq. 19 from Keeling et al., GBC, 1993)
    
    F = -dC/dT * Q/Cp
    
    hfds is in units of W/m^2
    Cp is in units of J/kg/K
    dcdt is in units of umol/kg/K
    """ 
    
    sos = ds['sos']
    tos = ds['tos']
    hfds = ds['hfds']

        
    Cp = 3990.
    dcdt = solubility.N2(sos, tos + 0.5) - solubility.N2(sos, tos - 0.5)
    
    ds['fgn2'] = -1.0 * scaleby *  dcdt * hfds / Cp * 1e-6 # umol/kg/K * W/m^2 / (J/kg/K) ==> mol m-2 s-1 (same as fgo2)
    ds.fgn2.attrs["units"] = "mol m-2 s-1"
    ds.fgn2.attrs["long_name"] = "N2 air-sea flux (computed from heat flux)"
    ds.fgn2.attrs['note'] = f'fluxes computed using F = c * (-dC/dT) * Q/Cp; c = {scaleby:0.4f}'
    
    return ds
        
# old Fortran:
#c       compute solubility and temperature derivative of solubility:
#c
#        call solub(sst(i,j,l),salt(i,j,l),c,dcdt)
#c
#c       compute gas flux:
#c
#        gasf(i,j,l)=-(dcdt/3990.)*heat(i,j,l)*2.628e6/22400.*.028
#c      
#c     3990. is the average heat capacity of seawater in joule/kg/K, see Neumann and Pearson page 47.  This number doesn't vary by more
#c     than about 0.5% with temperature and salinity. dcdt/3990 is then the gas flux / heat flux ratio in ml STP/joule 
#c     2.628e6 converts the flux to per month, 22400. converts the flux from ml to moles, and .028 converts from moles to kgN2.  The
#c     resulting flux is therefore in units of kgN2/gridcell/month, positive upwards. 
    
#      subroutine solub(t,sal,cc,dc)
#c     takes input t (temperature in celsius) and sal (salinity in per mil) and generates output cc (solubility in ml/kg) and dc (temperature derivative of solubility in ml/kg/K).
#      real a(4),b(3),tt,sal,cc,dc,lnc,lnc1,lnc2

#c     The following data for nitrogen in units of ml STP/kg
#      data a/-177.0212, 254.6078, 146.3611, -22.0933/
#      data b/-0.054052, 0.027266, -0.0038430/

#      tt = t+273.15
#      lnc = a(1)+a(2)*(100./tt)+a(3)*log(tt/100.)+a(4)*(tt/100.)+sal*(b(1) + b(2)*(tt/100.)+b(3)*(tt*tt/10000.))
#      tt = tt+0.5
#      lnc1 = a(1)+a(2)*(100./tt)+a(3)*log(tt/100.)+a(4)*(tt/100.)+sal*(b(1) + b(2)*(tt/100.)+b(3)*(tt*tt/10000.))
#      tt = tt-1.0
#      lnc2 = a(1)+a(2)*(100./tt)+a(3)*log(tt/100.)+a(4)*(tt/100.)+sal*(b(1) + b(2)*(tt/100.)+b(3)*(tt*tt/10000.))
#      cc = exp(lnc)
#      dc = exp(lnc1)-exp(lnc2)
#      return
#      end
     
    
def compute_fgo2_thermal(ds):
    """
    compute thermal O2 flux from heat flux and temperature derivative of solubility
    
    using Eq. 2 from Keeling and Shertz, 1992 (and Eq. 19 from Keeling et al., GBC, 1993)
    
    F = -dC/dT * Q/Cp
    
    hfds is in units of W/m^2
    Cp is in units of J/kg/K
    dcdt is in units of umol/kg/K
    """ 
    
    Cp = 3990.
    dcdt = solubility.O2(ds['sos'],ds['tos']+0.5) - solubility.O2(ds['sos'],ds['tos']-0.5)
    
    ds['fgo2_thermal'] = -1. * dcdt * ds['hfds'] / Cp * 1e-6 # umol/kg/K * W/m^2 / (J/kg/K) ==> mol m-2 s-1 (same as fgo2)
    ds.fgo2_thermal.attrs["units"] = "mol m-2 s-1"
    ds.fgo2_thermal.attrs["long_name"] = "O2 flux (thermal component)"
    
    return ds


def compute_fgapo(ds):
    """
    compute APO flux from O2, CO2, and N2 flux
    
    using 
    
    Fapo = Fo2 + 1.1 * Fco2 - X_O2/X_N2 * Fn2 
    """ 
    dsi = compute_fgn2(ds)
    
    ds['fgapo'] = ds['fgo2'] + 1.1 * ds['fgco2'] - X_O2 / X_N2 * dsi['fgn2'] # mol m-2 s-1 (same as fgo2)
    ds.fgapo.attrs["units"] = "mol m-2 s-1"
    ds.fgapo.attrs["long_name"] = "APO flux"
    
    return ds


def generate_latlon_grid(nx, ny, lon0=0.):
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

    dx = 360./nx
    dy = 180./ny

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
        {"lat": xr.DataArray(lat,dims=("lat")),         
         "lon": xr.DataArray(lon,dims=("lon")),
        })

  
    
    ds["xc"] = xr.DataArray(xc,dims=("lat", "lon"),
                            attrs={"units":"degrees_east",
                                   "long_name":"longitude of cell centers"})

    ds["yc"] = xr.DataArray(yc,dims=("lat", "lon"),
                            attrs={"units":"degrees_north",
                                   "long_name":"latitude of cell centers"})

    ds["xv"] = xr.DataArray(xv,dims=("lat", "lon", "nv"),
                            attrs={"units":"degrees_east",
                                   "long_name":"longitude of cell corners"})
    
    ds["yv"] = xr.DataArray(yv,dims=("lat", "lon",  "nv"),
                            attrs={"units":"degrees_north",
                                   "long_name":"latitude of cell corners"})

    ds["area"] = xr.DataArray(area,dims=("lat", "lon"),
                              attrs={"units":"m^2",
                                     "long_name":"area"})

    ds.lon.attrs = {"units": "degrees_east", "long_name":"Longitude"}
    ds.lat.attrs={"units": "degrees_north", "long_name": "Latitude"}
    
    for v in ds.variables:
        ds[v].encoding["_FillValue"] = None

    return ds


def infer_lat_name(ds):
    lat_names = ['latitude', 'lat']
    for n in lat_names:
        if n in ds:
            return n
    raise ValueError('could not determine lat name')    


def infer_lon_name(ds):
    lon_names = ['longitude', 'lon']
    for n in lon_names:
        if n in ds:
            return n
    raise ValueError('could not determine lon name')    
    


def lat_weights_regular_grid(lat):
    """
    Generate latitude weights for equally spaced (regular) global grids.
    Weights are computed as sin(lat+dlat/2)-sin(lat-dlat/2) and sum to 2.0.
    """   
    dlat = np.abs(np.diff(lat))
    np.testing.assert_almost_equal(dlat, dlat[0])
    w = np.abs(np.sin(np.radians(lat + dlat[0] / 2.)) - np.sin(np.radians(lat - dlat[0] / 2.)))

    if np.abs(lat[0]) > 89.9999: 
        w[0] = np.abs(1. - np.sin(np.pi / 2. - np.radians(dlat[0] / 2.)))

    if np.abs(lat[-1]) > 89.9999:
        w[-1] = np.abs(1. - np.sin(np.pi / 2. - np.radians(dlat[0] / 2.)))

    return w


def compute_grid_area(ds, check_total=True):
    """Compute the area of grid cells."""
    
    radius_earth = 6.37122e6 # m, radius of Earth
    area_earth = 4.0 * np.pi * radius_earth**2 # area of earth [m^2]e
    
    lon_name = infer_lon_name(ds)       
    lat_name = infer_lat_name(ds)        
    
    weights = lat_weights_regular_grid(ds[lat_name])
    area = weights + 0.0 * ds[lon_name] # add 'lon' dimension
    area = (area_earth / area.sum(dim=(lat_name, lon_name))) * area
    
    if check_total:
        np.testing.assert_approx_equal(np.sum(area), area_earth)
        
    return xr.DataArray(area, dims=(lat_name, lon_name), attrs={'units': 'm^2', 'long_name': 'area'})  


def yyyymmdd(year, month, day):
    return year * 10000 + month * 100 + day


def gen_midmonth_cftime_coord(year_range, shift_time=0., climatology_year_end=[]):
    
    yr0, yrf = year_range
    s = xr.cftime_range(start=f'{yr0}-01-01', end=f'{yrf}-12-31', freq='MS')
    e = xr.cftime_range(start=f'{yr0}-01-01', end=f'{yrf+1}-01-01', freq='M')
    
    units = f'days since {yr0:04d}-01-01 00:00:00'
    
    time_bounds_data = np.vstack((cftime.date2num(s, units), cftime.date2num(e, units) + 1)).T
    time_data = cftime.num2date(time_bounds_data.mean(axis=1), units)

    time = xr.DataArray(time_data, dims=('time'), name='time')
    time.attrs['shift_time'] = shift_time            
    time.encoding['units'] = units
    time.encoding['dtype'] = np.float64    
    time.encoding['_FillValue'] = None
    
    if climatology_year_end:
        
        for i in range(time_bounds_data.shape[0]):
            da = cftime.num2date(time_bounds_data[i, 1], units)
            time_bounds_data[i, 1] = cftime.date2num(
                cftime.datetime(climatology_year_end, da.month, da.day), units,
            )
            
        time.attrs['climatology'] = 'climatology_bounds'        
            
        time_bnds = xr.DataArray(
            cftime.num2date(time_bounds_data, units),
            dims=('time', 'd2'), 
            coords={'time': time},
            name='climatology_bounds',
        )                    
        
    else:
        time.attrs['bounds'] = 'time_bnds'

        time_bnds = xr.DataArray(
            cftime.num2date(time_bounds_data, units),
            dims=('time', 'd2'), 
            coords={'time': time},
            name='time_bnds',
        )                    


    
    time_bnds.encoding['dtype'] = np.float64
    time_bnds.encoding['_FillValue'] = None
    
    return time, time_bnds


def gen_date_variable(time):
    return xr.DataArray(
        [yyyymmdd(d.year, d.month, d.day) for d in time.values],
        dims=('time'),
        coords={'time': time},
        attrs={'long_name': 'Date', 'units': 'YYYYMMDD'},
        name='date',
    )

def gen_time_components_variable(time, year_offset=0.):
    
    if isinstance(time, xr.DataArray):
        time_data = time.values
    else:
        time_data = time
        
    tc = xr.DataArray(
        [(d.year + year_offset, d.month, d.day, d.hour, d.minute, d.second) for d in time_data],
        dims=('time', 'n_time_components'),
        coords={'time': time},
        attrs={'long_name': 'time components (year, month, day, hour, min, sec)', 'units': 'none'},
        name='time_components',
    )
    tc.encoding['_FillValue'] = None
    tc.encoding['dtype'] = np.int32
    return tc

def gen_daily_cftime_coord(year_range):

    yr0, yrf = year_range
    s = xr.cftime_range(start=f'{yr0}-01-01', end=f'{yrf}-12-31', freq='D')
    e = xr.cftime_range(start=f'{yr0}-01-02', end=f'{yrf+1}-01-01', freq='D')

    units = f'days since {yr0:04d}-01-01 00:00:00'

    time_bounds_data = np.vstack((cftime.date2num(s, units), cftime.date2num(e, units))).T
    time_data = cftime.num2date(time_bounds_data.mean(axis=1), units)

    time = xr.DataArray(time_data, dims=('time'), name='time')
    time.encoding['units'] = units
    time.encoding['dtype'] = np.float64    
    time.encoding['_FillValue'] = None

    time.attrs['bounds'] = 'time_bnds'

    time_bnds = xr.DataArray(
        cftime.num2date(time_bounds_data, units),
        dims=('time', 'd2'), 
        coords={'time': time},
        name='time_bnds',
    )                    

    time_bnds.encoding['dtype'] = np.float64
    time_bnds.encoding['_FillValue'] = None

    return time, time_bnds


def to_netcdf_clean(dset, path, format='NETCDF3_64BIT', **kwargs):
    """wrap to_netcdf method to circumvent some xarray shortcomings"""
    
    dset = dset.copy()
    
    git_sha = subprocess.check_output(['git', 'describe', '--always']).strip().decode("utf-8")
    datestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    provenance_str = f'created by {git_repo}/tree/{git_sha} on {datestamp}'

    if 'history' in dset.attrs:
        dset.attrs['history'] += '; ' + provenance_str
    else:
        dset.attrs['history'] = provenance_str    
    
    
    # ensure _FillValues are not added to coordinates
    for v in dset.coords:
        dset[v].encoding['_FillValue'] = None
    
    for v in dset.data_vars:
        if dset[v].dtype in [np.float32, np.float64]:
            dset[v].encoding['_FillValue'] = default_fillvals['f4']
            dset[v].encoding['dtype'] = np.float32
        
        elif dset[v].dtype in [np.int32, np.int64]:
            dset[v].encoding['_FillValue'] = default_fillvals['i4']
            dset[v].encoding['dtype'] = np.int32            
        elif dset[v].dtype == object:
            pass
        else:
            warnings.warn(f'warning: unrecognized dtype for {v}: {dset[v].dtype}')
    
    sys.stderr.flush()
            
    print('-'*30)
    print(f'Writing {path}')
    dset.to_netcdf(path, format=format, **kwargs)    
    dumps = subprocess.check_output(['ncdump', '-h', path]).strip().decode("utf-8")
    print(dumps)
    dumps = subprocess.check_output(['ncdump', '-k', path]).strip().decode("utf-8")
    print(f'format: {dumps}')    
    print('-'*30)
    
    
class curate_flux_products(object):    
    
    def __init__(self):
        
        self.catalog_file = "catalogs/flux_products-catalog-local.yml"
        if os.path.exists(self.catalog_file):
            with open(self.catalog_file, "r") as fid:
                self.catalog = yaml.safe_load(fid)               
        else:
            self.catalog = yaml.safe_load(
                """
                description: Flux products for transport modeling

                plugins:
                  source:
                    - module: intake_xarray

                sources: {}
                """
            )
    
    def add_source(self, key, urlpath, description, driver='netcdf', xarray_kwargs=None):
        """add a new source to the catalog"""
        if xarray_kwargs is None:
            xarray_kwargs = dict(
                decode_times=False,
            )
        
        self.catalog['sources'][key] = dict(
            driver=driver,
            description=description,
            args=dict(
                urlpath=urlpath,
                xarray_kwargs=xarray_kwargs,
            )
        )
        self.persist()
        
    def persist(self):
        """write the catalog to disk"""        
        with open(self.catalog_file, "w") as fid:
            yaml.dump(self.catalog, fid)    
    
    def open_catalog(self):
        """return as intake catalog"""
        return intake.open_catalog(self.catalog_file)
    
    def __repr__(self):
        return self.catalog.__repr__()
    
    
    
