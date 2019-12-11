# -*- coding: utf-8 -*-

"""
Created on Tue Apr 30 16:12:38 CEST 2019
Last modified, Thu May  2 14:09:11 CEST 2019

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
@author: Jean-Yves Peterschmitt a.k.a. jyp
"""

# Changes from version 0.0 : initial import of scattered code around
# Changes from version 0.1 : added smoothing function, suggestion of flhardy

__version__ = "0.2"


# from http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_closest(A,target) :
    import numpy             as np
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

def load_CLIO_grid(dummy_f=None):

    import netCDF4           as nc
    import numpy             as np


    if dummy_f == None:
       import os
       path = os.path.abspath(__file__)
       dir_path = os.path.dirname(path)
       dummy_f = os.path.join(dir_path,"/CLIO3_coordinate_file.nc")
    #end if

    # Suggestion from jyp ...
    if not os.path.exists(dummy_f):
       print("Gros message d'horreur")
    #endif

    # Get the coordinates from helper file
    # ====================================

    f = nc.Dataset(dummy_f)

    lons = f.variables['ulon'][...]
    lats = f.variables['ulat'][...]
    depth = f.variables['tdepth'][...]

    f.close()

    return lats, lons, depth
# end def load_CLIO_grid

def load_CLIO_var(f_name, var_name):

    import netCDF4           as nc

    f = nc.Dataset(f_name)

    var = f.variables[var_name][...].squeeze()

    f.close()

    return var
# end def load_CLIO_var

def rand_string(N=6):
    import random
    import string

    return ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(N))
#end def rand_string

def read_txt(fichname,sep="\t",cols=3):
    """A function reading a text data file assuming separator (sep) and columns (cols)"""

    from numpy import zeros
    from numpy import float32

    fich = open(fichname,'r')

    lignes = fich.readlines()
    taille = len(lignes)

# to be re-written with a "map" argument for any number of fields
    tab = zeros((cols,taille),float32)

    for i in range(taille):
        tab[:,i] = [float(dat) for dat in lignes[i].strip().split(sep)]
    #end for

    fich.close()

    return tab
# end of function read_txt

# from http://stackoverflow.com/questions/5515720/python-smooth-time-series-data
def smooth(x,window_len=30,window='hanning'):

        import numpy as np

        if x.ndim != 1:
            raise ValueError("smooth only accepts 1 dimension arrays.")
        if x.size < window_len:
            raise ValueError("Input vector needs"
                             +" to be bigger than window size.")
        if window_len<3:
            return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise ValueError("Window is one of 'flat',"
                             +"'hanning', 'hamming', 'bartlett', 'blackman'")

        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]

        if window == 'flat': #moving average
            w=np.ones(window_len,'d')
        else:
            w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]

#enddef smooth

# from http://stackoverflow.com/questions/5515720/python-smooth-time-series-data
def smoothD(x,window_len=30,window='hanning'):

        from numpy import ma

        shapeIn = x.shape
        Re_X = x.reshape((x.shape[0],x[0,...].size))
        Re_Y = ma.zeros(Re_X.shape)
        Re_Y.mask = True
        for i in range(Re_X.shape[-1]):
            if not ma.is_masked(Re_X[0,i]):
               Re_Y[:,i] = smooth(Re_X[:,i],window_len,window=window)
        #end for
        return Re_Y.reshape(shapeIn)
#enddef smoothD


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

def read_var_NC(nc_file,var_list):
    # Import the dataset to be manipulated

    import netCDF4 as nc

    # Opening file for read
    f_ile = nc.Dataset(nc_file, mode='r')

    var_listed = []

    for var_name in var_list:
       try:
           var = f_ile.variables[var_name]
           var_listed.append(var)
       except:
           print("Could not retrieve variable:"+var_name)
    #endfor

    return var_listed, f_ile

#end def read_variables_in_file


class ProgressBar:
    """ Creates a text-based progress bar. Call the object with the `print'
        command to see the progress bar, which looks something like this:

        [=======>        22%                  ]

        You may specify the progress bar's width, min and max values on init.

        Modified after Kelvie Wong version, see http://code.activestate.com/recipes/168639/

        D.M. Roche for present version, added the non re-printing if bar did not change,
         plus the deleting (clearing of screen after bar finishes)
    """

    def __init__(self, minValue = 0, maxValue = 100, totalWidth=80):
        self.progBar = "[]"   # This holds the progress bar string
        self.min = minValue
        self.max = maxValue
        self.span = maxValue - minValue
        self.width = totalWidth
        self.amount = 0       # When amount == max, we are 100% done
        self.updateAmount(0)  # Build progress bar string
        self._old_pbar = ''

    def updateAmount(self, newAmount = 0):
        """ Update the progress bar with the new amount (with min and max
            values set at initialization; if it is over or under, it takes the
            min or max value as a default. """
        if newAmount < self.min: newAmount = self.min
        if newAmount > self.max: newAmount = self.max
        self.amount = newAmount

        # Figure out the new percent done, round to an integer
        diffFromMin = float(self.amount - self.min)
        percentDone = (diffFromMin / float(self.span)) * 100.0
        percentDone = int(round(percentDone))

        # Figure out how many hash bars the percentage should be
        allFull = self.width - 2
        numHashes = (percentDone / 100.0) * allFull
        numHashes = int(round(numHashes))

        # Build a progress bar with an arrow of equal signs; special cases for
        # empty and full
        if numHashes == 0:
            self.progBar = "[>%s]" % (' '*(allFull-1))
        elif numHashes == allFull:
            self.progBar = "[%s]" % ('='*allFull)
        else:
            self.progBar = "[%s>%s]" % ('='*(numHashes-1),
                                        ' '*(allFull-numHashes))

        # figure out where to put the percentage, roughly centered
        percentPlace = (len(self.progBar) / 2) - len(str(percentDone))
        percentString = str(percentDone) + "%"

        # slice the percentage into the bar
        self.progBar = ''.join([self.progBar[0:percentPlace], percentString,
                                self.progBar[percentPlace+len(percentString):]
                                ])

    def __str__(self):
        return str(self.progBar)

    def __call__(self, value):
        """ Updates the amount, and writes to stdout. Prints a carriage return
            first, so it will overwrite the current line in stdout."""

        self.updateAmount(value)

        currstr = str(self)
        if currstr != self._old_pbar:
          self._old_pbar = currstr
          print('\r',)
          sys.stdout.write(str(self))
          sys.stdout.flush()

    def __del__(self):
          print('\r',)
          sys.stdout.write("%s" % (' '*self.width))
          print('\r',)
          sys.stdout.flush()

#

# The End of All Things (op. cit.)
