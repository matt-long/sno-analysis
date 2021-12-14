import os
from collections import OrderedDict

import csv

import numpy as np
import xarray as xr

Re = 6.37122e6 # m, radius of Earth

path_to_here = os.path.dirname(os.path.realpath(__file__))

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




def open_mask_dataset(res=1, clobber=False):
    """
    ### Make a basin mask on 1x1 grid

    read ascii mask files from NODC and make some regional aggregations

    Mask files downloaded from:
       https://www.nodc.noaa.gov/OC5/woa13/masks13.html
   
    Documentation found here: https://data.nodc.noaa.gov/woa/WOA13/DOC/woa13documentation.pdf
    """
    assert res in [1, 0.25]
    res_code = '01' if res == 1 else '04'
    mask_cache_file = f"{path_to_here}/basin-mask_{res_code}.nc"
    
    if os.path.exists(mask_cache_file) and not clobber:
        return xr.open_dataset(mask_cache_file)

    if res == 1:
        nx = 360
        ny = 180
        ds = generate_latlon_grid(nx=nx,ny=ny,lon0=-180., include_vertices=False)
        basinmask_csvfile = f'{path_to_here}/basinmask_01.msk'
        landsea_csvfile = f'{path_to_here}/landsea_01.msk'    

    elif res == 0.25:
        nx = 1440
        ny = 720
        ds = generate_latlon_grid(nx=nx,ny=ny,lon0=-180., include_vertices=False)
        basinmask_csvfile = f'{path_to_here}/basinmask_04.msk'
        landsea_csvfile = f'{path_to_here}/landsea_04.msk'    
        
    ds['mask'] = xr.DataArray(np.zeros(ds.area.shape).astype(int),dims=('lat','lon'))
    ds['kmt'] = xr.DataArray(np.zeros(ds.area.shape).astype(int),dims=('lat','lon'))
        

    #--- basin mask
    # The basin_XX.msk contains the basin code number defined for each grid square at each
    # standard depth from the surface to 5500m. Each basin is identified by a code number that
    # ranges from 1 to 58. The basin code number in a given quarter-degree and one-degree square
    # may change with increased depth level. Appendix 1 lists the geographic basin names, the
    # code number associated with each basin, and the standard depth level at which the given
    # basin is first encountered.
    lat_ndx = 0
    lon_ndx = 1
    mask_ndx = 2
    with open(basinmask_csvfile, 'rt') as f:
        csvdata = csv.reader(f, delimiter=',',skipinitialspace=True)

        for ir, row in enumerate(csvdata):
            if ir < 2:
                pass
            else:
                j = np.where(ds.lat == float(row[lat_ndx]))[0]
                i = np.where(ds.lon == float(row[lon_ndx]))[0]
                if len(j) == 0 or len(i) == 0:
                    print('ERROR: no match')
                    break
                ds.mask.values[j,i] = int(row[mask_ndx])


    #--- land sea mask
    # The landsea_XX.msk contains the standard depth level number at which the bottom of the
    # ocean is first encountered at each quarter-degree or one-degree square for the entire world.
    # Land will have a value of 1
    lat_ndx = 0
    lon_ndx = 1
    mask_ndx = 2
    with open(landsea_csvfile, 'rt') as f:
        csvdata = csv.reader(f, delimiter=',',skipinitialspace=True)

        for ir, row in enumerate(csvdata):
            if ir < 2:
                pass
            else:
                j = np.where(ds.lat == float(row[lat_ndx]))[0]
                i = np.where(ds.lon == float(row[lon_ndx]))[0]
                ds.kmt.values[j,i] = int(row[mask_ndx])            

    ds['mask_orig'] = ds.mask.copy()

    #-- change GIN Sea to Atlantic
    ds.mask.values = np.where(
        (55<ds.lat) & (ds.lat<77) & 
        (ds.mask_orig==11.) & (ds.lon>-45.) & (ds.lon<20.), 1., ds.mask_orig.values)

    #-- change Lab Sea to Atlantic
    ds.mask.values = np.where(
        (55<ds.lat) & (ds.lat<80) & 
        (ds.mask_orig==11.) & (-80<ds.lon) & (ds.lon<-45.), 1., ds.mask.values)

    #-- change definition of Southern Ocean to encompass region south of Tierra del Fuego
    # atlantic
    ds.mask.values = np.where(
        (-58.<ds.lat) & (ds.mask_orig==10.) & 
        (-70<ds.lon) & (ds.lon<20.), 1., ds.mask.values)

    # pacific
    ds.mask.values = np.where(
        (-58.<ds.lat) & (ds.mask_orig==10.) & 
        ((147<=ds.lon) | (ds.lon<=-70.)), 2., ds.mask.values)

    # indian
    ds.mask.values = np.where(
        (-58.<ds.lat) & (ds.mask_orig==10.) & 
        (20<ds.lon) & (ds.lon<147), 3., ds.mask.values)    
    
    
    ds.to_netcdf(mask_cache_file)
    
    return ds      


def _get_basin_aggregation(product):
    """aggregate region mask into basins"""
    if product.lower() == "resplandy":
        my_regions = OrderedDict([('Global' , [1., 2., 3., 5., 10., 11., 12., 56.]),
                                  ('Pacific' , [2., 12.]),
                                  ('Atlantic' , [1., 5.]),
                                  ('Indian' , [3., 56.]),
                                  ('Arctic', [11.]),
                                  ('Southern_Ocean' , [10.])])

    elif product.lower() == "gruber":
        my_regions = OrderedDict([('Global' , [1., 2., 3., 4., 5., 9., 10., 11., 12., 56.]),
                                  ('Pacific' , [2., 12.]),
                                  ('Atlantic' , [1., 4., 5., 9.]),
                                  ('Indian' , [3., 56.]),
                                  ('Arctic', [11.]),
                                  ('Southern_Ocean' , [10.])])
    else:
        raise ValueError(f"not a known basin aggregation: {aggregation}")
        
        
    basin = open_mask_dataset()
    nlat = len(basin.lat)
    nlon = len(basin.lon)
    nrgn = len(my_regions)

    basin['REGION_MASK_3D'] = xr.DataArray(
        np.zeros((nrgn, nlat, nlon)), 
        dims=('region','lat','lon'),
    )
    

    for i, (name, ndx) in enumerate(my_regions.items()):
        tmp = np.zeros((nlat, nlon))
        for n in ndx:
            tmp = np.where(basin.mask.values == n, 1, tmp)
        basin.REGION_MASK_3D.values[i,:,:] = tmp
        basin.REGION_MASK_3D.values = basin.REGION_MASK_3D.where(basin.kmt!=1)
        basin[name] = basin.REGION_MASK_3D.sel(region=i)    
    
    basin['region'] = xr.DataArray(
        list(my_regions.keys()),
        dims=('region'), 
    )    
    return basin


def _get_inversion_regions(product):
    """make mask for inversion region"""
    basin = _get_basin_aggregation(product)
    nlat = len(basin.lat)
    nlon = len(basin.lon)
    M = xr.DataArray(np.ones((nlat, nlon)), dims=('lat','lon'))
    
    if product.lower() == 'resplandy':
        my_regions = OrderedDict([
            ('1' , M.where((basin.Arctic==1.))),
            ('2' , M.where((basin.lat>49.) & (basin.Atlantic == 1.))),
            ('3' , M.where((35.<basin.lat) & (basin.lat<=49.) & (basin.Atlantic == 1.))),
            ('4' , M.where((20.<basin.lat) & (basin.lat<=35.) & (basin.Atlantic == 1.))),
            ('5' , M.where((0.<basin.lat) & (basin.lat<=20.) & (basin.Atlantic == 1.))),
            ('6' , M.where((-20.<basin.lat) & (basin.lat<=0.) & (basin.Atlantic == 1.))),
            ('7' , M.where((-35.<basin.lat) & (basin.lat<=-20.) & (basin.Atlantic == 1.))),
            ('8' , M.where((-44.<basin.lat) & (basin.lat<=-35.) & (basin.Atlantic == 1.))),
            ('9-25-30' , M.where((-58.<basin.lat) & (basin.lat<=-44.))),
            ('10' , M.where((-90.<basin.lat) & (basin.lat<=-58.))),
            ('11' , M.where((50.<basin.lat) & (basin.lat<=70.) & 
                            ~((basin.lon<0.) & (basin.lon>-170.)) & (basin.Pacific == 1.))),
            ('12' , M.where((40.<basin.lat) & (basin.lat<=70.) & 
                            ((basin.lon<0.) & (basin.lon>-170.)) & (basin.Pacific == 1.))),
            ('13-14-15' , M.where(((20.<basin.lat) & (basin.lat<=50.) & 
                                   ~((basin.lon<0.) & (basin.lon>-170.)) & (basin.Pacific == 1.)) |
                                ((20.<basin.lat) & (basin.lat<=40.) & 
                                 ((basin.lon<0.) & (basin.lon>-170.)) & (basin.Pacific == 1.)))),
            ('16' , M.where((0.<basin.lat) & (basin.lat<=20.) & 
                            ~((basin.lon<0.) & (basin.lon>-160.)) & (basin.Pacific == 1.))),
            ('17' , M.where((0.<basin.lat) & (basin.lat<=20.) & 
                            ((basin.lon<0.) & (basin.lon>-160.)) & (basin.Pacific == 1.))),
            ('18' , M.where((-20.<basin.lat) & (basin.lat<=0.) & 
                            ~((basin.lon<0.) & (basin.lon>-160.)) & (basin.Pacific == 1.))),
            ('19' , M.where((-20.<basin.lat) & (basin.lat<=0.) & 
                            ((basin.lon<0.) & (basin.lon>-160.)) & (basin.Pacific == 1.))),
            ('20-22' , M.where(((-30.<basin.lat) & (basin.lat<=-20.) & 
                                ~((basin.lon<0.) & (basin.lon>-130.)) & (basin.Pacific == 1.)) |
                              ((-44.<basin.lat) & (basin.lat<=-30.) & 
                               ~((basin.lon<0.) & (basin.lon>-112.)) & (basin.Pacific == 1.)))),
            ('21-23-24' , M.where(((-30.<basin.lat) & 
                                   (basin.lat<=-20.) & 
                                   ((basin.lon<0.) & (basin.lon>-130.)) & (basin.Pacific == 1.)) |
                                 ((-44.<basin.lat) & (basin.lat<=-30.) & 
                                  ((basin.lon<0.) & (basin.lon>-112.)) & (basin.Pacific == 1.)))),
            ('26-27' , M.where((-20.<basin.lat) & (basin.lat<=50.) & (basin.Indian == 1.))),
            ('28-29' , M.where((-44.<basin.lat) & (basin.lat<=-20.) & (basin.Indian == 1.)))])        
    
    elif product.lower() == 'gruber':
        #-- region 1
        my_regions = OrderedDict([
            ('1' , M.where((basin.Arctic==1.) |
                          ((basin.lat>53.) & (basin.Atlantic == 1.)))),
            ('2' , M.where((13.<basin.lat) & (basin.lat<=53.) & (basin.Atlantic == 1.))),
            ('3' , M.where((-13.<basin.lat) & (basin.lat<=13.) & (basin.Atlantic == 1.))),
            ('4' , M.where((-36.<basin.lat) & (basin.lat<=-13.) & (basin.Atlantic == 1.))),
            ('5' , M.where((-58.<basin.lat) & (basin.lat<=-36.) & ((basin.Atlantic == 1.)))), 
            ('6' , M.where((36.<basin.lat) & (basin.lat<=70.) & (basin.Pacific == 1.))),
            ('7' , M.where((13.<basin.lat) & (basin.lat<=36.) & (basin.Pacific == 1.))),
            ('8&9' , M.where((-13.<basin.lat) & (basin.lat<=13.) & (basin.Pacific == 1.))),
            ('10' , M.where((-36.<basin.lat) & (basin.lat<=-13.) & (basin.Pacific == 1.))),     
            ('11&14' , M.where((-58.<basin.lat) & (basin.lat<=-36.) & ((basin.Pacific == 1.) |
                                                                       (basin.Indian == 1.)))),
            ('12' , M.where((-13.<basin.lat) & (basin.lat<=70.) & (basin.Indian == 1.))),    
            ('13' , M.where((-36.<basin.lat) & (basin.lat<=-13.) & (basin.Indian == 1.))),
            ('15' , M.where((-90.<basin.lat) & (basin.lat<=-58.)))])
    else:
        raise ValueError(f"not a known basin aggregation: {aggregation}")
        
        
    ds = open_mask_dataset()
    
    nlat = len(ds.lat)
    nlon = len(ds.lon)    
    nrgn = len(my_regions)
    ds['REGION_MASK_3D'] = xr.DataArray(
        np.zeros((nrgn, nlat, nlon)), dims=('region','lat','lon')
    )
    
    for i, mask_logic in enumerate(my_regions.values()):
        ds.REGION_MASK_3D.values[i,:,:] = mask_logic.fillna(0.)
    ds.REGION_MASK_3D.values = ds.REGION_MASK_3D.where(ds.kmt!=1)

    ds['REGION_MASK'] = ds.REGION_MASK_3D.isel(region=0)
    for i in range(nrgn):
        ds.REGION_MASK.values = np.where(ds.REGION_MASK_3D[i,:,:] == 1, i+1, ds.REGION_MASK)
        
    ds["region"] = xr.DataArray(list(my_regions.keys()), dims=("region"))
    return ds


