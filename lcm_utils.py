# -*- coding: utf-8 -*-

"""
Created on Tue Apr 30 16:12:38 CEST 2019
Last modified, Thu May 2 09:19:27 CEST 2019

 Copyright 2019 Didier M. Roche <didier.roche@lsce.ipsl.fr>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

@author: Didier M. Roche a.k.a. dmr
"""

# Changes from version 0.0 : initial import of scattered code around

__version__ = "0.1"


# from http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_closest(A,target) :
    return (np.abs(A-target)).argmin()
#enddef find_closest

def distance(lon1, lat1, lon2, lat2):

    import math

    # http://code.activestate.com/recipes/576779-calculating-distance-between-two-geographic-points/
    # http://en.wikipedia.org/wiki/Haversine_formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    q = math.sin(dlat/2)**2 + (math.cos(lat1) * math.cos(lat2) * (math.sin(dlon/2)**2))
    return 2 * math.atan2(math.sqrt(q), math.sqrt(1-q))
#enddef distance

# Searching for CLIO coordinates without sorting on the native grid ...
# implicitly assumes lon is within (0,360.0)
def find_closest_2D(tagt_lon,tagt_lat,lons_array,lats_array):

    import math
    import numpy as np
    dist_array = np.array([distance(math.radians(tagt_lon),math.radians(tagt_lat),math.radians(lons_array[x,y]),math.radians(lats_array[x,y]))
                           for x in range(lons_array.shape[0]) for y in range(lons_array.shape[1])]).reshape(lons_array.shape)
    return np.unravel_index(dist_array.argmin(), dist_array.shape)
#enddef

def load_CLIO_grid(dummy_f="CLIO3_coordinate_file.nc"):

    import netCDF4           as nc
    import numpy             as np

    # Get the coordinates from helper file
    # ====================================

    f = nc.Dataset(dummy_f)

    lons = f.variables['ulon'][...]
    lats = f.variables['ulat'][...]
    depth = f.variables['tdepth'][...]

    #~ salt = f.variables['salt'][0,...].squeeze()

    f.close()

    return lats, lons, depth  # , salt
# end def load_CLIO_grid

def rand_string(N=6):
    import random
    import string

    return ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(N))
#end def rand_string

def plot_CLIO_2D(var2plot,lats,lons,show=False,proj_typ="ortho", clo=None, cla=None):

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import numpy as np
    import numpy.ma as ma

    min_bounds = ma.min(var2plot)
    max_bounds = ma.max(var2plot)
    nbs_bounds = 30
    fix_bounds = np.linspace(min_bounds,max_bounds,nbs_bounds)

    the_chosen_map = plt.cm.coolwarm

    fig = plt.figure(figsize=(10,10))

    if proj_typ == "PC":
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0.0, globe=None))
    elif proj_typ == "ortho":
        if clo != None:
          centralLon = clo
        else:
          centralLon = -70.0
        #end if
        if cla != None:
          centralLat = cla
        else:
          centralLat = 40.0
        #end if

        ax = plt.axes(projection=ccrs.Orthographic(central_longitude=centralLon, central_latitude=centralLat, globe=None))
    #~ ax = plt.axes(projection=ccrs.Miller(central_longitude=360))
    #~ ax = plt.axes(projection=ccrs.LambertConformal(central_longitude=-180.0, central_latitude=39.0, false_easting=0.0, false_northing=0.0, secant_latitudes=None, standard_parallels=None, globe=None, cutoff=-30))
    #~ ax = plt.axes(projection=ccrs.LambertCylindrical(central_longitude=0.0))
    #~ ax = plt.axes(projection=ccrs.Mercator(central_longitude=0.0, min_latitude=-80.0, max_latitude=84.0, globe=None, latitude_true_scale=0.0))
    #endif
    ax.set_global()
    fig.set_facecolor("grey")

    mesh = ax.pcolormesh(lons, lats, var2plot, cmap=the_chosen_map,
                     transform=ccrs.PlateCarree(), edgecolor="white",linewidth=0.05,
                     vmin=min_bounds, vmax=max_bounds)
    plt.colorbar(mesh, orientation='horizontal', shrink=0.75)
    ax.gridlines()
    ax.coastlines()

    if show:
       plt.show()
    #end if
#end def plot_CLIO

# The End of All Things (op. cit.)
