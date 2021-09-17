import os
import numpy as np

from toolz import curry

import xarray as xr
import pandas as pd

import intake

import matplotlib.pyplot as plt

# make variable so as to enable system dependence
catalog_csv = '/glade/collections/cmip/catalog/intake-esm-datastore/catalogs/glade-cmip6.csv.gz'
catalog_json = '/glade/collections/cmip/catalog/intake-esm-datastore/catalogs/glade-cmip6.json'

cmip6_catalog = intake.open_esm_datastore(catalog_json)

# hardwired to only use native grid
grid_label = 'gn'

# constants
T0_Kelvin = 273.15
mols_to_Tmolmon = 1e-12 * 86400. * 365. / 12.
µmolkg_to_mmolm3 = 1026. / 1000. # for volume conserving models, makes sense to use constant density
kgCO2s_to_Tmolmon = 1000. / 12. * mols_to_Tmolmon
W_to_PW = 1. / 1E15
Re = 6.37122e6 # m, radius of Earth

class missing_data_tracker(object):        
    def __init__(self): 
        """construct object for tracking missing data"""
        self._columns = ['source_id', 'experiment_id', 'table_id', 'variable_id',  'grid_label']
        
        self.missing_data_file = 'data/missing-data.csv'            
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
        os.remove(self.missing_data_file)
    
    
def get_gridvar(df, source_id, variable_id):
    """get a grid variable from a source_id"""
    df_sub = df.loc[
        (df.source_id==source_id) 
        & (df.variable_id==variable_id)
        & (df.grid_label == grid_label)
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
        grid_label=grid_label,
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
            
            # sort ensemble members to prioritize physics and forcing 1, and to sort realizations numerically 
            real = np.array(member_ids)
            init = np.array(member_ids)
            phys = np.array(member_ids)
            forc = np.array(member_ids)
            for i in range(len(member_ids)):
                real[i]=member_ids[i].split('r')[1].split('i')[0]
                init[i]=member_ids[i].split('i')[1].split('p')[0]
                phys[i]=member_ids[i].split('p')[1].split('f')[0]
                forc[i]=member_ids[i].split('f')[1]
            real=real.astype(int)
            init=init.astype(int)
            phys=phys.astype(int)
            forc=forc.astype(int)      
            member_ids_sorted=np.array(member_ids)[np.lexsort((real,init,phys,forc))]
            
            ds = ds.sel(member_id=member_ids_sorted[:nmax_members])
        
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
        if sumormean == 'mean':
            da = da / rmask.sum(dims_lateral)
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
    """
    compute N2 flux from heat flux and temperature derivative of solubility
    
    using Eq. 2 from Keeling and Shertz, 1992 (and Eq. 19 from Keeling et al., GBC, 1993)
    
    F = -dC/dT * Q/Cp
    
    hfds is in units of W/m^2
    Cp is in units of J/kg/K
    dcdt is in units of umol/kg/K
    """ 
    
    Cp = 3990.
    dcdt = _N2sol(ds['sos'],ds['tos']+0.5) - _N2sol(ds['sos'],ds['tos']-0.5)
    
    ds['fgn2'] = -1. * dcdt * ds['hfds'] / Cp * 1e-6 # umol/kg/K * W/m^2 / (J/kg/K) ==> mol m-2 s-1 (same as fgo2)
    
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
    dcdt = _O2sol(ds['sos'],ds['tos']+0.5) - _O2sol(ds['sos'],ds['tos']-0.5)
    
    ds['fgo2_thermal'] = -1. * dcdt * ds['hfds'] / Cp * 1e-6 # umol/kg/K * W/m^2 / (J/kg/K) ==> mol m-2 s-1 (same as fgo2)
    
    return ds


def _O2sol(S, T):
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


def N2solWeiss(S,T):
    '''
    Solubility of N2 in sea water
    INPUT:  
    S = salinity    [PSS]
    T = temperature [degree C]
    
    REFERENCE:
    Weiss, 1970.
    "The solubility of nitrogen, oxygen and argon in water and seawater"
    Deep-sea Research, 17, pp. 721-735.
    
    returns umol/kg
    '''
    
    # T is absolute T
    Tabs = 275.15 + T
    
    rho_ref = 1.026 # g/cm3 (approx. at 15 C)

    # these were coeffs or Bunsen solubility coeff:
    #A1 = -59.6274
    #A2 = 85.7661
    #A3 = 24.3696
    #B1 = -0.051580
    #B2 = 0.026329
    #B3 = -0.0037252

    #N2_sol_an = np.log(A1 + A2*(100.0/T) + S*(B1 + B2*(T/100.0) + B3*((T/100.0)**2)))
    #units_ml_kg__umol_kg = 1.0/0.022391
    #N2_sol = N2_sol_an*units_ml_kg__umol_kg
    #return _umolkg_to_mmolm3(N2_sol)

    # this looks like the equation for Bunsen solubility coeff, but should be np.exp not np.log, also missing A3 term and unit conversion:
    #return np.log(A1 + A2*(100.0/T) + S*(B1 + B2*(T/100.0) + B3*((T/100.0)**2)))
    
    # these are coeffs and equation for ml/kg
    A1 = -177.0212
    A2 = 254.6078
    A3 = 146.3611 
    A4 = -22.0933
    B1 = -0.054052 
    B2 = 0.027266
    B3 = -0.0038430
    
    ml_per_kg_to_mmol_per_m3 = 1 / 22.4 * rho_ref * 1e3 
    ml_to_umol = 1 / 22.4 * 1e3

    return np.exp(A1 + A2*(100.0/Tabs) + A3*np.log(Tabs/100) + A4*(Tabs/100) + S*(B1 + B2*(Tabs/100.0) + B3*(Tabs/100.0)**2)) * ml_to_umol ## * ml_per_kg_to_mmol_per_m3


def N2solHamme(S,T):
    
    '''
    # constants from Table 4 of Hamme and Emerson 2004
    Coef. Ne (nmol/kg) N2 (umol/kg) Ar (umol/kg)
    A0 2.18156 6.42931 2.79150
    A1 1.29108 2.92704 3.17609
    A2 2.12504 4.32531 4.13116
    A3 0 4.69149 4.90379
    B0 -5.94737E-3 -7.44129E-3 -6.96233E-3
    B1 -5.13896E-3 -8.02566E-3 -7.66670E-3
    B2 0 -1.46775E-2 -1.16888E-2
    Check: 7.34121 500.885 13.4622
    check values at temperature of 10 C and salinity of 35 (PSS)
    
    returns umol/kg
    '''
    
    #rho_ref = 1.026 # g/cm3 (approx. at 15 C)

    A0 = 6.42931
    A1 = 2.92704
    A2 = 4.32531
    A3 = 4.69149
    B0 = -7.44129E-3
    B1 = -8.02566E-3
    B2 = -1.46775E-2
    
    #ln C = A0 + A1*Ts + A2*Ts^2 + A3*Ts^3 + S (B0 + B1 * Ts + B2 * Ts^2)
    #Ts = ln((298.15 - t)/(273.15 + t)
    
    T_scaled = np.log((298.15 - T) /(273.15 + T))
    return np.exp(A0 + A1*T_scaled + A2*T_scaled**2. + A3*T_scaled**3. + \
                  S*(B0 + B1*T_scaled + B2*T_scaled**2.)) ## * rho_ref # convert to mmol/m^3/atm


def _N2sol(S, T):
    '''
    Solubility (saturation) of nitrogen (N2) in sea water
    at 1-atm pressure of air including saturated water vapor
   
    INPUT:  (if S and T are not singular they must have same dimensions)
    S = salinity    [PSS]
    T = temperature [degree C]
   
    OUTPUT:
    conc = solubility of N2  [µmol/kg]
   
    REFERENCE:
    Roberta Hamme and Steve Emerson, 2004.
    "The solubility of neon, nitrogen and argon in distilled water and seawater."
    Deep-Sea Research I, 51(11), p. 1517-1528.
    '''
   
    return _garcia_gordon_polynomial(S, T,
                                     A0=6.42931,
                                     A1=2.92704,
                                     A2=4.32531,
                                     A3=4.69149,
                                     B0=-7.44129e-3,
                                     B1=-8.02566e-3,
                                     B2=-1.46775e-2)


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
        {"lat": xr.DataArray(lat,dims=("lat"),
                            attrs={"units":"degrees_north",
                                   "long_name":"Latitude"}),
         "lon": xr.DataArray(lon,dims=("lon"),
                            attrs={"units":"degrees_east",
                            "long_name":"Longitude"})})

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