# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 11:21:51 2020

@author: antho
"""

#%%
import os
os.environ['GDAL_DATA'] = r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\gdal'

#%%
import flopy #import FloPy to develop, run, and analyze the model
from flopy.utils import Raster #plot rasters with FloPy
import math
import matplotlib as mp
import pandas as pd
import pyproj #change between WGS84 and Illimap coordinates
from pyproj import Transformer
import rasterio  #import rasters
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import cartopy
import cartopy.crs as ccrs #import projections
import cartopy.feature as cf #import features
from pykrige.uk import UniversalKriging
import pylab #used as a plotting library for spatial data, make contours
from metpy.plots import USCOUNTIES

#%%
#--------------------------------------------------
# Define model domain in lat/long coordinates
# Coordinates for each corner of the domain
sw_lat = 41.174519     #41.174519 #41.356272 #41.41 #southwest latitude
sw_long = -88.408339     #-88.408339 #-88.311217 #-88.24 #southwest longitude
ne_lat = 41.934977     #41.934977 #41.701295 #41.73 #northeast latitude
ne_long = -87.64     #-87.318569 #-87.831109 #-88.03 #northeast longitude

# Note: for now, ne_long must be less than or equal to roughly -87.64 deg, in order to stay within the bounds of Illinois

illimap = {'proj': 'lcc', # code for Lambert Conformal Conic
     'ellps': 'clrk66',
     'lon_0': -89.5,
     'lat_0': 33,
     'lat_1': 33,
     'lat_2': 45,
     'x_0': 2999994*0.3048,
     'y_0': 0}

# Define the transformation from lat/long (WGS84, or EPSG 4326) to Lambert x/y (Illimap)
transformer = Transformer.from_crs("epsg:4326", illimap)

# Note that WGS84 and EPSG 4326 are equivalent

# Transform the coordinates of the NE and SW corners from WGS84 (lat and long) to Illimap (x and y)
nex, ney = transformer.transform(ne_lat, ne_long)
swx, swy = transformer.transform(sw_lat, sw_long)

# Convert nex, ney, swx, swy from m to ft and round
nex, ney = round(nex/0.3048,-4), round(ney/0.3048,-4)
swx, swy = round(swx/0.3048,-4), round(swy/0.3048,-4)

#--------------------------------------------------
# Define spatial and temporal discretization

# Assign discretization variables
Lx = nex-swx # Width of the model domain in feet
Ly = ney-swy # Height of the model domain in feet
nlay = 10 # Number of model layers
dx = 2000 
dy = 2000
nrow = int(Ly/dy) # Number of rows
ncol = int(Lx/dx) # Number of columns

years = []
for i in range(1950,2020):
    i = str(i)
    years.append(i)

nper = len(years) #specify number of stress periods
steady = [False] #specify if stress period is transient or steady-state

#--------------------------------------------------
# Define river elevations

# Import river stage, lambert x, lambert y from river Excel file
dfriv = pd.read_csv(r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\large_files\rivers_625.csv')

# Trim dataframe with river information to the model domain
dfriv = dfriv.loc[dfriv['lamx']<nex]
dfriv = dfriv.loc[dfriv['lamy']<ney]
dfriv = dfriv.loc[dfriv['lamx']>swx]
dfriv = dfriv.loc[dfriv['lamy']>swy]

# Assign all rivers to the upper layer
dfriv['lay'] = 0 #this actually assigns it to layer 1
# Convert lamx to column and lamy to row
dfriv['row'] = np.trunc((ney-dfriv['lamy'])/dy)
dfriv['col'] = np.trunc((dfriv['lamx']-swx)/dx)
# Define the river stage
dfriv['stage'] = dfriv['rvr_stg']
# Define the conductance
dfriv['cond'] = 90000. #ft^2/d (this value was adjusted during calibration)
# Define the river bottom
dfriv['bot'] = dfriv['stage']-3
# Drop unneeded columns
dfriv = dfriv.drop(['STR_ORD_MI','STR_ORD_MA','SUM_LENGTH','rvr_stg','lamx','lamy'],axis=1)

# Group duplicate rivers within a cell
dfriv = dfriv.groupby(['lay','row','col'],as_index=False).mean()

#--------------------------------------------------
# Define top and bottom elevations

# Now load the raster using FloPy's built in Raster toolbox
illinoisdem = Raster.load(r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\large_files\landsurface_el.tif')
bedrock = Raster.load(r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\large_files\bedrock_el.tif')

# Crop the DEM to the model domain
illinoisdem.crop([(swx,swy),(swx,ney),(nex,ney),(nex,swy)])
bedrock.crop([(swx,swy),(swx,ney),(nex,ney),(nex,swy)])
# Define centroid of the southwestern most cell
startx = swx+dx/2 
starty = swy+dy/2
# Calculate the x and y coordinates for the centroid of each cell 
xc = np.arange(swx+dx/2,nex+dx/2,dx) 
yc = np.arange(swy+dy/2,ney+dy/2,dy)
# Create a grid of the x coordinate of each centroid and the y coordinate
xarr, yarr = np.meshgrid(xc,yc)
# Resample the topo raster to the grid of centroids of the model
topgrid = illinoisdem.resample_to_grid(xarr,yarr,1,method='nearest') 
bedrock = bedrock.resample_to_grid(xarr,yarr,1,method='nearest')

# We built our top elevation upside down, let's flip it
topgrid = np.flipud(topgrid) 
bedrockgrid = np.flipud(bedrock)   

# Create ibound as array of ints (1), indicating all cells are active
# Inactivate cells west of the Mississippi River that were originally not present
# Note that because inactive cells would overlap with the river boundaries, this code pushes inactive cells to the west a bit. Adjust per your model domain
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)

# Set the top of Layer 1 = River Elevation
for index, row in dfriv.iterrows():  
    topgrid[int(row['row']),int(row['col'])]=row['stage'] 
    #print(topgrid[int(row['row']),int(row['col'])],row['bot'])  

# Make sure that all layers (combined) are at least 9 ft thick
diff = topgrid-bedrockgrid

diff[diff<=9.0] = 9.0

# Each layer EXCEPT BEDROCK is the same thickness, we need 9 glacial layers
laythick = diff/9 #glacial layers
bothick = 50. #thickness of bedrock. May need to change?? look at report?

# Calculate the bottom of each layer
lay1bot = topgrid-laythick
lay2bot = topgrid-2*laythick
lay3bot = topgrid-3*laythick
lay4bot = topgrid-4*laythick
lay5bot = topgrid-5*laythick
lay6bot = topgrid-6*laythick
lay7bot = topgrid-7*laythick
lay8bot = topgrid-8*laythick
lay9bot = topgrid-9*laythick
lay10bot = lay9bot-bothick 

botgrids = [lay1bot,lay2bot,lay3bot,lay4bot,lay5bot,lay6bot,lay7bot,lay8bot,lay9bot,lay10bot]

#--------------------------------------------------
# Assign hydraulic conductivity

# Assign hydraulic conductivity in ft/day
kc = 220 #predominantly coarse
kf = .3 #predominantly fine
kb = 20 #bedrock

# Determine how to assign hydraulic conductivity
threshold = 15 #25 #19 #anything above this will be assigned kc; anything below will be assigned kf

def kloader(rastername, kc, kf, threshold):
  percent = Raster.load(rastername) #load raster
  percent.crop([(swx,swy),(swx,ney),(nex,ney),(nex,swy)]) #crop array
  percentgrid = percent.resample_to_grid(xarr,yarr,1,method='nearest') #resample to model grid
  percentgrid = np.flipud(percentgrid) #flip the grid
  maxrow = percentgrid.shape[0]
  maxcol = percentgrid.shape[1]
  for row in np.arange(maxrow,0,-1):
    for col in np.arange(maxcol,0,-1):
      if percentgrid[row-1,col-1] < -10:
        percentgrid[row-1,col-1] = percentgrid[row-1,col]
  percentgrid[percentgrid>=threshold] = kc #assign coarse k value
  percentgrid[percentgrid<threshold] = kf #assign fine k value
  return percentgrid

kl1 = kloader(r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\large_files\percentl1.tif',kc,kf,threshold)
kl2 = kloader(r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\large_files\percentl2.tif',kc,kf,threshold)
kl3 = kloader(r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\large_files\percentl3.tif',kc,kf,threshold)
kl4 = kloader(r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\large_files\percentl4.tif',kc,kf,threshold)
kl5 = kloader(r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\large_files\percentl5.tif',kc,kf,threshold)
kl6 = kloader(r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\large_files\percentl6.tif',kc,kf,threshold)
kl7 = kloader(r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\large_files\percentl7.tif',kc,kf,threshold)
kl8 = kloader(r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\large_files\percentl8.tif',kc,kf,threshold)
kl9 = kloader(r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\large_files\percentl9.tif',kc,kf,threshold)
kl10 = kl9-kl9+kb

khlayers = [kl1,kl2,kl3,kl4,kl5,kl6,kl7,kl8,kl9,kl10]
kvlayers=np.divide(khlayers,10.)

#%%
# Define wells

# Import well data from .csv file
dfwel = pd.read_csv(r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\pumping\transient\transient_pumping_WillCounty.csv')
dfwel = dfwel.set_index('p_num') #assign index as p_number so that other columns can be deleted

# Trim dataframe with well information to the model domain
dfwel = dfwel.loc[dfwel['lam_x']<nex]
dfwel = dfwel.loc[dfwel['lam_y']<ney]
dfwel = dfwel.loc[dfwel['lam_x']>swx]
dfwel = dfwel.loc[dfwel['lam_y']>swy]

# Define the flux for each year as the pumpage data from the imported file 
# and convert from gal/year to ft^3/day to match the units of the model.
# Negative pumpage denotes removal of water from the system.

for year in years:
    dfwel[year] = np.float64(dfwel[year]) * -1 / 2730.39

# Define the location of each well as layer, row, and column
# Assign all wells to bedrock layer (10 layers but starts at 0 so the last layer is 9)
dfwel['lay'] = 9
# Convert lamx to column and lamy to row
dfwel['row'] = np.trunc((ney-dfwel['lam_y'])/dy)
dfwel['col'] = np.trunc((dfwel['lam_x']-swx)/dx)

#dfwel['flux']=dfwel['2002']*-1/2730.39

# Drop unneeded columns
dfwel = dfwel.drop(['isws_facility_id','owner','fac_well_num','total_name',
                    'depth_total_last_known','lam_x','lam_y'], axis=1)

print(dfwel)

#%%

# lay row col Q

welldata = {}

for i in range(nper):
    arwell = []
    for j in dfwel.index:
        arwell.append([dfwel.loc[j,'lay'], dfwel.loc[j,'row'], dfwel.loc[j,'col'], dfwel.loc[j,years[i]]])
    welldata[i] = arwell