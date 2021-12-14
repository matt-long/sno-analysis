import numpy as np
import xarray as xr
from numba import jit


@jit(nopython=True, parallel=True)
def creep_fill(nlat, nlon, var, fillmask, max_iter):
    """Iterative fill algorithm."""

    done = False
    iter_cnt = 0

    work = np.empty_like(var)

    fill_with = np.empty((4,))
    
    while not done:
        iter_cnt += 1

        # assume bottom row is land, so skip it
        for j in range(0, nlat):
            jm1 = j - 1
            jp1 = j + 1

            for i in range(0, nlon):
                # assume periodic in x
                im1 = (i - 1) % nlon
                ip1 = (i + 1) % nlon

                work[j, i] = var[j, i]

                if not fillmask[j, i]:
                    continue
                    
                # self
                if np.isfinite(var[j, i]):
                    continue
                    
                fill_with[:] = np.nan
                propagate = 0
                
                # East
                if np.isfinite(var[j, ip1]):
                    fill_with[0] = var[j, ip1]
                    propagate = 1
                    
                # West
                if np.isfinite(var[j, im1]):
                    fill_with[1] = var[j, im1]
                    propagate = 1
                    
                # South
                if np.isfinite(var[jm1, i]):
                    fill_with[2] = var[jm1, i]
                    propagate = 1

                # North
                if j < nlat - 1:
                    if np.isfinite(var[jp1, i]):
                        fill_with[3] = var[jp1, i]
                        propagate = 1
                        

                if propagate > 0:
                    for n in range(4):
                        if np.isfinite(fill_with[n]):
                            work[j, i] = fill_with[n]
                            break

        var[:, :] = work[:, :]
        if iter_cnt == max_iter:
            done = True