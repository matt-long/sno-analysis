import numpy as np
import xarray as xr
import ESMF

def esmf_create_grid(xcoords, ycoords, xcorners=False, ycorners=False,
                corners=False, domask=False, doarea=False,
                ctk=ESMF.TypeKind.R8):
    """
    Create a 2 dimensional Grid using the bounds of the x and y coordiantes.
    :param xcoords: The 1st dimension or 'x' coordinates at cell centers, as a Python list or numpy Array
    :param ycoords: The 2nd dimension or 'y' coordinates at cell centers, as a Python list or numpy Array
    :param xcorners: The 1st dimension or 'x' coordinates at cell corners, as a Python list or numpy Array
    :param ycorners: The 2nd dimension or 'y' coordinates at cell corners, as a Python list or numpy Array
    :param domask: boolean to determine whether to set an arbitrary mask or not
    :param doarea: boolean to determine whether to set an arbitrary area values or not
    :param ctk: the coordinate typekind
    :return: grid
    """
    [x, y] = [0, 1]

    # create a grid given the number of grid cells in each dimension, the center stagger location is allocated, the
    # Cartesian coordinate system and type of the coordinates are specified
    max_index = np.array([len(xcoords), len(ycoords)])
    grid = ESMF.Grid(max_index,
                     staggerloc=[ESMF.StaggerLoc.CENTER],
                     coord_sys=ESMF.CoordSys.SPH_DEG,
                     num_peri_dims=1,
                     periodic_dim=x,
                     coord_typekind=ctk)

    # set the grid coordinates using pointer to numpy arrays, parallel case is handled using grid bounds
    gridXCenter = grid.get_coords(x)
    x_par = xcoords[grid.lower_bounds[ESMF.StaggerLoc.CENTER][x]:grid.upper_bounds[ESMF.StaggerLoc.CENTER][x]]

    gridXCenter[...] = x_par.reshape((x_par.size, 1))

    gridYCenter = grid.get_coords(y)
    y_par = ycoords[grid.lower_bounds[ESMF.StaggerLoc.CENTER][y]:grid.upper_bounds[ESMF.StaggerLoc.CENTER][y]]
    gridYCenter[...] = y_par.reshape((1, y_par.size))

    # create grid corners in a slightly different manner to account for the bounds format common in CF-like files
    if corners:
        raise ValueError('not tested')
        grid.add_coords([ESMF.StaggerLoc.CORNER])
        lbx = grid.lower_bounds[ESMF.StaggerLoc.CORNER][x]
        ubx = grid.upper_bounds[ESMF.StaggerLoc.CORNER][x]
        lby = grid.lower_bounds[ESMF.StaggerLoc.CORNER][y]
        uby = grid.upper_bounds[ESMF.StaggerLoc.CORNER][y]

        gridXCorner = grid.get_coords(x, staggerloc=ESMF.StaggerLoc.CORNER)
        for i0 in range(ubx - lbx - 1):
            gridXCorner[i0, :] = xcorners[i0+lbx, 0]
        gridXCorner[i0 + 1, :] = xcorners[i0+lbx, 1]

        gridYCorner = grid.get_coords(y, staggerloc=ESMF.StaggerLoc.CORNER)
        for i1 in range(uby - lby - 1):
            gridYCorner[:, i1] = ycorners[i1+lby, 0]
        gridYCorner[:, i1 + 1] = ycorners[i1+lby, 1]

    # add an arbitrary mask
    if domask:
        raise ValueError('not tested')
        mask = grid.add_item(ESMF.GridItem.MASK)
        mask[:] = 1
        mask[np.where((1.75 <= gridXCenter.any() < 2.25) &
                      (1.75 <= gridYCenter.any() < 2.25))] = 0

    # add arbitrary areas values
    if doarea:
        raise ValueError('not tested')
        area = grid.add_item(ESMF.GridItem.AREA)
        area[:] = 5.0

    return grid


def esmf_create_locstream_spherical(lon, lat, coord_sys=ESMF.CoordSys.SPH_DEG,
                                    mask=None):
    """
    :param coord_sys: the coordinate system of the LocStream
    :param domask: a boolean to tell whether or not to add a mask
    :return: LocStream
    """
    if ESMF.pet_count() is not 1:
        raise ValueError("processor count must be 1 to use this function")

    locstream = ESMF.LocStream(len(lon), coord_sys=coord_sys)

    locstream["ESMF:Lon"] = lon
    locstream["ESMF:Lat"] = lat
    if mask is not None:
        locstream["ESMF:Mask"] = mask.astype(np.int32)

    return locstream


def esmf_interp_points(ds_in, locs_lon, locs_lat, lon_field_name='lon',
                lat_field_name='lat'):
    """Use ESMF toolbox to interpolate grid at points."""

    # generate grid object
    grid = esmf_create_grid(ds_in[lon_field_name].values.astype(np.float),
                            ds_in[lat_field_name].values.astype(np.float))


    # generate location stream object
    locstream = esmf_create_locstream_spherical(locs_lon.values.astype(np.float),
                                                locs_lat.values.astype(np.float))
    
    # generate regridding object
    srcfield = ESMF.Field(grid, name='srcfield')
    dstfield = ESMF.Field(locstream, name='dstfield')
    regrid = ESMF.Regrid(srcfield, dstfield,
                         regrid_method=ESMF.RegridMethod.BILINEAR,
                         unmapped_action=ESMF.UnmappedAction.ERROR)

    # construct output dataset
    coords = {c: locs_lon[c] for c in locs_lon.coords}
    dims_loc = locs_lon.dims
    nlocs = len(locs_lon)
    ds_out = xr.Dataset(coords=coords, attrs=ds_in.attrs)

    for name, da_in in ds_in.data_vars.items():

        # get the dimensions of the input dataset; check if it's spatial
        dims_in = da_in.dims
        if lon_field_name not in dims_in or lat_field_name not in dims_in:
            continue

        # get the dimension/shape of output
        non_lateral_dims = dims_in[:-2]
        dims_out = non_lateral_dims + dims_loc
        shape_out = da_in.shape[:-2] + (nlocs,)

        # create output dataset
        da_out = xr.DataArray((np.ones(shape_out)*np.nan).astype(da_in.dtype),
                              name=name,
                              dims=dims_out,
                              attrs=da_in.attrs,
                              coords={c: da_in.coords[c] for c in da_in.coords
                                      if c in non_lateral_dims})
        dstfield.data[...] = np.nan

        if len(non_lateral_dims) > 0:
            da_in_stack = da_in.stack(non_lateral_dims=non_lateral_dims)
            da_out_stack = xr.full_like(da_out, fill_value=np.nan).stack(non_lateral_dims=non_lateral_dims)

            for i in range(da_in_stack.shape[-1]):
                srcfield.data[...] = da_in_stack.data[:, :, i].T
                dstfield = regrid(srcfield, dstfield, zero_region=ESMF.Region.SELECT)
                da_out_stack.data[:, i] = dstfield.data

            da_out.data = da_out_stack.unstack('non_lateral_dims').transpose(*dims_out).data

        else:
            srcfield.data[...] = da_in.data[:, :].T
            dstfield = regrid(srcfield, dstfield, zero_region=ESMF.Region.SELECT)
            da_out.data = dstfield.data

        ds_out[name] = da_out
        
    return ds_out

