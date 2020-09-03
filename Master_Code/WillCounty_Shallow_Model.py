# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 14:59:49 2020

@author: antho
"""

#%%

'''
To-Do:

- Download TIFs, so we don't need to import them from Google Drive
- Download MF2005, MFNWT, and MT3DMS executable files
- Add in code for RMS error
- Ignore deprecation warnings produced when creating the error map
- Remove sections of code for uploading files to Colab
'''


#%% IMPORT PACKAGES

#--------------------------------------------------

import flopy # import FloPy to develop, run, and analyze the model
from flopy.utils import Raster # plot rasters with FloPy
import matplotlib as mp
import pandas as pd
import pyproj # change between WGS84 and Illimap coordinates
import rasterio  # import rasters
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import cartopy
import cartopy.crs as ccrs # import projections
import cartopy.feature as cf # import features
from pykrige.uk import UniversalKriging
import pylab # used as a plotting library for spatial data, make contours
from metpy.plots import USCOUNTIES

#--------------------------------------------------
# Packages related to Google Drive

# the following code authorizes you to access files on Google Drive
from google.colab import drive
from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive
from google.colab import auth, files
from oauth2client.client import GoogleCredentials

#%% IMPORT FILES FROM GOOGLE DRIVE

#--------------------------------------------------
# DEM

# First import the land surface .tif from Google Drive
downloaded = drive.CreateFile({'id':"1389l8sgQ8-tsmIZuZosaqvbqpHY40n6l"}) # ft above msl
downloaded.GetContentFile('landsurface_el.tif')

# First import the bedrock elevation .tif from Google Drive
downloaded = drive.CreateFile({'id':"1EZgZDjjILzvRzvY9nf0Qp0NHmspRq4kP"})   
downloaded.GetContentFile('bedrock_el.tif')

# Read in percent thickness of coarse grain for each model layer
downloaded = drive.CreateFile({'id':"18Kw3O6qCzIJ2L6KrVnRPIhea_F8VwyWn"})   
downloaded.GetContentFile('percentl1.tif')

downloaded = drive.CreateFile({'id':"1oZinFPKrGY-FXoE7Zu0okFpAAOe_bwau"})   
downloaded.GetContentFile('percentl2.tif')

downloaded = drive.CreateFile({'id':"1FqVEr4m_ElUyEZeyfnCMwVGDfUqavJZH"})   
downloaded.GetContentFile('percentl3.tif')

downloaded = drive.CreateFile({'id':"1KiHS9TLSP1GAVTjaaJZS4BAwF6gnUeDu"})   
downloaded.GetContentFile('percentl4.tif')

downloaded = drive.CreateFile({'id':"1Z-9EyaAK1NKnRHAlnyGYkI3suvBFC2I6"})   
downloaded.GetContentFile('percentl5.tif')

downloaded = drive.CreateFile({'id':"1pcB9aJpJGfkXOKz10rhs6MpWkQL1_dqr"})   
downloaded.GetContentFile('percentl6.tif')

downloaded = drive.CreateFile({'id':"1Fnh0HIKbUj7pEtlsUKR_Sr7WwfYzWul5"})   
downloaded.GetContentFile('percentl7.tif')

downloaded = drive.CreateFile({'id':"106JacgpwSA3wVAGcBIzGdc8rDVUB6dh7"})   
downloaded.GetContentFile('percentl8.tif')

downloaded = drive.CreateFile({'id':"1WJjhVJ_KSBhZDrgzY3YteNjxaz5nxBid"})   
downloaded.GetContentFile('percentl9.tif')

#--------------------------------------------------
# Rivers

# First import the Excel file from Google Drive
downloaded = drive.CreateFile({'id':"1JsAiGG4RvcfYrQtfgXRW9ZVfAkQ1yRVu"})
downloaded.GetContentFile('rivers_625.csv')

#%% MODEL SETUP

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

prj = pyproj.Proj(illimap) #define "prj" as the illimap coordinate system

wgs84 = pyproj.Proj("epsg:4326") #define "wgs84" as the EPSG 4326 / WGS84 coordinate system

# Note that WGS84 and EPSG 4326 are equivalent

# Transform the coordinates of the NE and SW corners from WGS84 (lat and long) to Illimap (x and y)
nex, ney = pyproj.transform(wgs84,illimap,ne_lat,ne_long) 
swx, swy = pyproj.transform(wgs84,illimap,sw_lat,sw_long)

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

nper = 1 #specify number of stress periods
steady = [True] #specify if stress period is transient or steady-state

#--------------------------------------------------
# Define river elevations

# Import river stage, lambert x, lambert y from river Excel file
dfriv = pd.read_csv('rivers_625.csv')

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

# Import river stage, lambert x, lambert y from river Excel file
dfriv = pd.read_csv('rivers_625.csv')

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
illinoisdem = Raster.load("landsurface_el.tif")
bedrock = Raster.load("bedrock_el.tif")

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

kl1 = kloader('percentl1.tif',kc,kf,threshold)
kl2 = kloader('percentl2.tif',kc,kf,threshold)
kl3 = kloader('percentl3.tif',kc,kf,threshold)
kl4 = kloader('percentl4.tif',kc,kf,threshold)
kl5 = kloader('percentl5.tif',kc,kf,threshold)
kl6 = kloader('percentl6.tif',kc,kf,threshold)
kl7 = kloader('percentl7.tif',kc,kf,threshold)
kl8 = kloader('percentl8.tif',kc,kf,threshold)
kl9 = kloader('percentl9.tif',kc,kf,threshold)
kl10 = kl9-kl9+kb

khlayers = [kl1,kl2,kl3,kl4,kl5,kl6,kl7,kl8,kl9,kl10]
kvlayers=np.divide(khlayers,10.)

#--------------------------------------------------
# Define wells

# Import well data from .csv file
dfwel = pd.read_csv('https://raw.githubusercontent.com/dbabrams/ShallowDolomite_Group/master/pumping/2002_pumping_V2.csv?token=AOLJKYQ7WAUJVO2JL55VRYK7BYM4C')
dfwel = dfwel.set_index('p_num') #assign index as p_number so that other columns can be deleted

# Trim dataframe with well information to the model domain
dfwel = dfwel.loc[dfwel['lam_x']<nex]
dfwel = dfwel.loc[dfwel['lam_y']<ney]
dfwel = dfwel.loc[dfwel['lam_x']>swx]
dfwel = dfwel.loc[dfwel['lam_y']>swy]

# Put the data into the format required for the well package, with columns for layer, row, column, and flux
# Assign all wells to bedrock layer (10 layers but starts at 0 so the last layer is 9)
dfwel['lay'] = 9
# Convert lamx to column and lamy to row
dfwel['row'] = np.trunc((ney-dfwel['lam_y'])/dy)
dfwel['col'] = np.trunc((dfwel['lam_x']-swx)/dx)
# Define the flux as the pumpage data from the imported file and convert from gal/year to ft^3/day to match the units of the model. Negative pumpage denotes removal of water from the system.
dfwel['flux']=dfwel['2002']*-1/2730.39

# Drop unneeded columns
dfwel = dfwel.drop(['isws_facility_id','owner','fac_well_num','depth_total_last_known','lam_x','lam_y','2002'], axis=1)

#print(dfwel)

#--------------------------------------------------
# Define Drains

# Add cells to the drain package wherever low-k material is at land surface.  These are cells in layer 1 with hydraulic conductivity k = kf.

# Create two arrays, one for rows and one for columns, corresponding to the row and column #s of elements in kl1 where k = kf.
rows,cols = np.nonzero(kl1==kf)

# Create a dataframe of the same length as "rows" called "dfdrn".  "dfdrn" includes a single column of zeros called "lay" which assigns all drain cells to the first layer
dfdrn = pd.DataFrame(np.zeros((len(rows), 1)),columns=['lay'])

# Add in the drain cell row and column #s
dfdrn['row'] = rows 
dfdrn['col'] = cols

# Define the elevation of the drain cells (set to land surface elevation)
for value in dfdrn.index:
  dfdrn.loc[value,'elevation'] = topgrid[dfdrn.loc[value,'row'],dfdrn.loc[value,'col']]

# Define the conductance of the drain cells
dfdrn['cond'] = kf*dx*dy/3 #this is the conductance between the cell and the drain, in ft^2/d
# For more information regarding the equation for conductance, see the "Notes on the conductance of drain cells" comment below.

# Notes on the conductance of drain cells
# C=K*dx*dy/delta ; C=K*dx*dy/b
# K is the conductivity of the soil layer.  K is defined as two orders of magnitude greater than kf; the soil layer is assumed to be more weathered than the underlying material.
# delta or b is the depth of the soil layer
# Out of a handful of randomly selected well logs, the depth of the uppermost soil layer (where present) varied significantly.
# For simplicity, the soil layer is assumed to have a thickness of 3 feet, uniform throughout the model.

#print(dfdrn)

#--------------------------------------------------
# Create the MODFLOW model object

# Create a MODFLOW model object and run with MODFLOW 2005.
modelname = "my_model" # name the model
m = flopy.modflow.Modflow(modelname, version = 'mf2005', exe_name = 'mf2005') # create model object m

#--------------------------------------------------
# Append the discretization package to the model object

# Length and time are feet (1) and days (4).
# See https://water.usgs.gov/ogw/modflow/MODFLOW-2005-Guide/index.html?dis.htm 
dis = flopy.modflow.ModflowDis(model=m, nlay=nlay, nrow=nrow, ncol=ncol, 
                               delr=dx, delc=dy, top=topgrid, botm=botgrids, 
                               itmuni = 4, lenuni = 1, 
                               nper=nper, steady=steady)

#--------------------------------------------------
# Basic package

# Create ibound as array of ints (1), indicating all cells are active
#ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)

# Create FloPy bas object
bas = flopy.modflow.ModflowBas(m, ibound=ibound, strt=topgrid)

#--------------------------------------------------
# LPF package

# Define layer type as convertible (1), must be an integer
# For more information, see https://water.usgs.gov/ogw/modflow/MODFLOW-2005-Guide/index.html?dis.htm
laytyp = 0*np.ones((nlay,), dtype=np.int32)

# create the LPF object
lpf = flopy.modflow.ModflowLpf(model=m, hk=khlayers, vka=kvlayers, laytyp=laytyp, ipakcb=1)

#%% BOUNDARY CONDITIONS

#--------------------------------------------------
# Recharge package

# Assign recharge to the model
# Units=ft/d

'''
# Original values
rch1=0.0014 #fine
rch2=0.0024 #coarse
'''

rch1=0.0007 #fine
rch2=0.0012 #coarse

# Create a recharge array based on values in kl1 (the top layer)
# Cells where kl1 = kf are assigned recharge = rch1
# All other cells (where kl1 = kc) are assigned recharge = rch2
recharge=np.where(kl1<=kf, rch1, rch2) 

# Create the recharge package (RCH)
rch = flopy.modflow.mfrch.ModflowRch(model=m, ipakcb=1, nrchop=3, rech = recharge)
#ipakcb=1 signals that cell-by-cell budget data should be saved
#nrchop=3 assigns recharge to the uppermost active cell in each column

#--------------------------------------------------
# River package

# Put into a format that MODFLOW wants
arriv = dfriv.values
riverdata = {0: arriv} #dictionary assigning river data; key (0) is the stress period

# Create the river package (RIV)
riv = flopy.modflow.mfriv.ModflowRiv(model=m, ipakcb=1, stress_period_data = riverdata)
#ipakcb=1 signals that cell-by-cell budget data should be saved

#--------------------------------------------------
# Drain package

# Note: This code block takes about 3 minutes to run

# This code block ensures that no river cells are included in the drain package

# Create an array from the values in the drain dataframe
# In other words, create "ardrn" from "dfdrn"
ardrn = dfdrn.values
ardrn = ardrn.tolist() #convert the drain array from an array to a list

ardrn_initial=len(ardrn)

# Print the initial length of the drain array
print('The initial length of the drain array is',ardrn_initial)

# Search through the river dataframe and the drain array for matches
# If any river cells are included in the drain array, the following loops will remove them
for i in range(len(dfriv)-1):
  for j in range(len(ardrn)-1):
    if dfriv.loc[i,'row']==np.float64(ardrn[j][1]) and dfriv.loc[i,'col']==np.float64(ardrn[j][2]):
      ardrn.remove(ardrn[j])

ardrn_final=len(ardrn)

# Print the final length of the drain array
print('The final length of the drain array is',ardrn_final)
print(ardrn_initial-ardrn_final,'items were removed')

#--------------------------------------------------

# Put into a format that MODFLOW wants
draindata = {0: ardrn} #dictionary assigning drain data; key (0) is the stress period
# Here, "draindata" functions similarly to "riverdata" for creating the river package

# Create the drain package (DRN)
drn = flopy.modflow.mfdrn.ModflowDrn(model=m, ipakcb=1, stress_period_data = draindata)
#ipakcb=1 signals that cell-by-cell budget data should be saved

#--------------------------------------------------
# Well package

# Put into a format that MODFLOW wants
arwell = dfwel.values
welldata = {0: arwell} #dictionary assigning well data; key (0) is the stress period
# Here, "welldata" functions similarly to "riverdata" and "draindata" in the previous code blocks

# Create the well package (WEL)
wel = flopy.modflow.mfwel.ModflowWel(model=m, ipakcb=1, stress_period_data = welldata)
#ipakcb=1 signals that cell-by-cell budget data should be saved

#--------------------------------------------------
# Define output control

# Create oc stress period data. 
spd = {(0, 0): ['print head', 'print budget', 'save head', 'save budget']}
# Create output control object
oc = flopy.modflow.ModflowOc(model=m, stress_period_data=spd, compact=True)

#--------------------------------------------------
# Solver

# We will start by using the PCG solver with default settings
#pcg = flopy.modflow.ModflowPcg(model=m)
pcg = flopy.modflow.ModflowPcg(model=m,mxiter=200,iter1=50,hclose=1e-03,rclose=1e-03,relax=0.98,damp=0.3)

#%% Plot model inputs (boundary conditions, elevations)

#--------------------------------------------------
''''Plot grid and boundary conditions'''

plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.PlotMapView(model=m, layer=0)
#grid = modelmap.plot_grid()
ib = modelmap.plot_ibound()
rvr = modelmap.plot_bc(ftype='RIV')
#add labels and legend
plt.xlabel('Lx (ft)',fontsize = 14)
plt.ylabel('Ly (ft)',fontsize = 14)
plt.title('Ibound', fontsize = 15, fontweight = 'bold')
plt.legend(handles=[mp.patches.Patch(color='blue',label='Const. Head',ec='black'),
                   mp.patches.Patch(color='white',label='Active Cell',ec='black'),
                   mp.patches.Patch(color='black',label='Inactive Cell',ec='black'),
                   mp.patches.Patch(color='green',label='River',ec='green')],
                   bbox_to_anchor=(1.5,1.0))
plt.show()

#--------------------------------------------------
'''Plot hydraulic conductivity'''

plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.map.PlotMapView(model=m, layer=0) #use plotmapview to attach plot to model
#contour_levels = np.linspace(400,800,41)
#topelevations = modelmap.contour_array(topgrid, levels = contour_levels) #create head contours
#plt.clabel(topelevations, inline=True,fontsize=12,fmt='%1.0f')

# Create colormap of named colors
colors = ["saddlebrown","lightgoldenrodyellow"]
cmap = mp.colors.LinearSegmentedColormap.from_list("", colors)
norm = mp.colors.LogNorm(vmin=kf,vmax=kc)
modelmap.plot_array(khlayers[0],norm = norm,cmap=cmap)
rvr = modelmap.plot_bc(ftype='RIV')
ib = modelmap.plot_ibound()

# Display parameters
plt.xlabel('Lx (ft)',fontsize = 14)
plt.ylabel('Ly (ft)',fontsize = 14)
plt.title('Top Elevation (ft AMSL); Hydraulic Conductivity', fontsize = 15, fontweight = 'bold')

# Create a legend for the hydraulic conductivity
plt.legend(handles=[
  mp.patches.Patch(color='green', label='rivers'),
  mp.patches.Patch(color='saddlebrown', label='kf=0.3 ft/d'),
  mp.patches.Patch(color='lightgoldenrodyellow', label='kc=280 ft/d') ], bbox_to_anchor=(1.7,0.8))
plt.show()

#--------------------------------------------------
'''Plot recharge'''

fig=plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.map.PlotMapView(model=m, layer=0) #use plotmapview to attach plot to model
#contour_levels = np.linspace(400,800,41)
#topelevations = modelmap.contour_array(topgrid, levels = contour_levels) #create head contours
#plt.clabel(topelevations, inline=True,fontsize=12,fmt='%1.0f')

# Create colormap of named colors
colors = ["saddlebrown","lightgoldenrodyellow"]
cmap = mp.colors.LinearSegmentedColormap.from_list("", colors)
norm = mp.colors.LogNorm(vmin=rch1,vmax=rch2)
pcm=modelmap.plot_array(recharge, norm = norm,cmap=cmap)
rvr = modelmap.plot_bc(ftype='RIV')
ib = modelmap.plot_ibound()
#display parameters
plt.xlabel('Lx (ft)',fontsize = 14)
plt.ylabel('Ly (ft)',fontsize = 14)
plt.title('Top Elevation (ft AMSL); Recharge Zones', fontsize = 15, fontweight = 'bold')

# Create a legend for the recharge zones 
plt.legend(handles=[
    mp.patches.Patch(color='green', label='rivers'),
    mp.patches.Patch(color='saddlebrown', label='rch1=0.0014 ft/d'),
    mp.patches.Patch(color='lightgoldenrodyellow', label='rch2=0.0024 ft/d') ], bbox_to_anchor=(1.7,0.8))
plt.show()

#--------------------------------------------------
'''Plot transects'''

plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelxsect = flopy.plot.PlotCrossSection(model = m, line={"row":22}) #use plotmapview to attach plot to model. row indicates west to east 
#modelxsect = flopy.plot.PlotCrossSection(model = m, line={"column":15}) #plots north-south transect

# Create colormap of named colors
colors = ["saddlebrown","gray","lightgoldenrodyellow"]
cmap = mp.colors.LinearSegmentedColormap.from_list("", colors)
norm = mp.colors.LogNorm(vmin=kf,vmax=kc)
#modelxsect.plot_grid()
khlaynp = np.array(khlayers)
lines = modelxsect.plot_array(khlaynp,norm=norm, cmap=cmap)
rvr = modelxsect.plot_bc(ftype='RIV')
modelxsect.plot_ibound()
plt.show()

#--------------------------------------------------
'''Plot Transmissivity'''

# Horizontal T = sum(k*b), for each layer.  k is the hydraulic conductivity, and b is the layer thickness
Thz=sum([k for k in khlayers])*laythick-kl10*laythick+kl10*bothick
#----------------------------------------------------------------------------
plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.map.PlotMapView(model=m, layer=0) #use plotmapview to attach plot to model

# Create colormap of named colors
colors = ["darkred","red","pink","white","lightskyblue","blue","indigo"]
cmap = mp.colors.LinearSegmentedColormap.from_list("", colors)
norm = mp.colors.LogNorm(vmin=np.amin(Thz),vmax=np.amax(Thz))
Tmap=modelmap.plot_array(Thz,norm = norm,cmap=cmap)
rvr = modelmap.plot_bc(ftype='RIV')
#display parameters
plt.xlabel('Lx (ft)',fontsize = 14)
plt.ylabel('Ly (ft)',fontsize = 14)
plt.title('Transmissivity of Aquifer in Model Domain', fontsize = 15, fontweight = 'bold')
cbar=plt.colorbar(Tmap,ticks=[np.amin(Thz),np.amax(Thz)], label='Transmissivity (sq.ft/day)')
cbar.set_ticklabels([int(round(np.amin(Thz),0)),int(round(np.amax(Thz),0))])
plt.legend(handles=[mp.patches.Patch(color='green', label='rivers')], bbox_to_anchor=(1.7,0.8))
plt.show()

#%% WRITE AND RUN THE MODFLOW MODEL

# Write the model input
m.write_input()
# Execute the model run
success, mfoutput = m.run_model(pause=False, report=True)
# Report back if the model did not successfully complete
if not success:
    raise Exception('MODFLOW did not terminate normally.')

#%% PLOT OUTPUT DATA

#--------------------------------------------------
'''Extract binary data from head and flow files'''

# Extract binary data from head file as flopy head object
headobj = flopy.utils.binaryfile.HeadFile(modelname+'.hds')

# Extract head data from head object
head = headobj.get_data(totim=1.0)

#print(head[9]) #returns heads in the bottom layer, where pumping occurs

#--------------------------------------------------
'''Plot results with boundary conditions'''

plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.map.PlotMapView(model=m, layer=0) #use plotmapview to attach plot to model
contour_levels = np.linspace(400,1000,61) #set contour levels for contouring head
head_contours = modelmap.contour_array(head[9], levels=contour_levels) #create head contours
plt.clabel(head_contours, inline=True,fontsize=12,fmt='%1.0f')
rvr_map = modelmap.plot_bc(ftype='RIV',color="blue")
drn_map = modelmap.plot_bc(ftype='DRN', color='peachpuff')
wel_map = modelmap.plot_bc(ftype='WEL',plotAll=True) #plot well locations in the model domain, shown on the top layer for visibility

# Display parameters
plt.xlabel('Lx (ft)',fontsize = 14)
plt.ylabel('Ly (ft)',fontsize = 14)
plt.title('Results with Boundary Conditions', fontsize = 15, fontweight = 'bold')

# Create a legend for the boundary conditions
plt.legend(handles=[
  mp.patches.Patch(color='blue', label='Rivers'),
  mp.patches.Patch(color='peachpuff', label='Drains'),
  mp.patches.Patch(color='red', label='Wells') ], bbox_to_anchor=(1.4,0.9))

plt.show()

#--------------------------------------------------
'''Plot potentiometric surface results'''

plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.map.PlotMapView(model=m, layer=0) #use plotmapview to attach plot to model
contour_levels = np.linspace(400,1000,61) #set contour levels for contouring head with values from 400-1000 and with 61 divisions
head_contours = modelmap.contour_array(head[9], levels=contour_levels) #create head contours
plt.clabel(head_contours, inline=True,fontsize=12,fmt='%1.0f')
rvr_map = modelmap.plot_bc(ftype='RIV',color="blue")
wel_map = modelmap.plot_bc(ftype='WEL',plotAll=True) #plot the well sites onto the model domain and show with the top layer

# Display parameters
plt.xlabel('Lx (ft)',fontsize = 14)
plt.ylabel('Ly (ft)',fontsize = 14)
plt.title('Steady-State Model Head(ft) Results', fontsize = 15, fontweight = 'bold')

plt.show()

#%% CALIBRATION

'''

"Can it wait for a bit?  I'm in the middle of some calibrations..."

'''

#--------------------------------------------------
'''1-to-1 plot'''

# Import the observation well data as a dataframe
# "SB_Potent_Surface_points.csv"
pumping_ob = pd.read_csv('https://raw.githubusercontent.com/dbabrams/ShallowDolomite_Group/master/pumping/SB_Potent_Surface_points.csv?token=AOLJKYVBUPMW5FKZK2TPE427BYM5W')

# Trim the dataframe to the model domain
pumping_ob = pumping_ob.loc[pumping_ob['lambx']<nex]
pumping_ob = pumping_ob.loc[pumping_ob['lamby']<ney]
pumping_ob = pumping_ob[pumping_ob['lambx']>swx]
pumping_ob = pumping_ob[pumping_ob['lamby']>swy]

# Reset the dataframe index to the default
pumping_ob=pumping_ob.reset_index(drop=True)

# Convert lambx to column and lamby to row, and convert to grid format
pumping_ob['row'] = np.trunc((ney-pumping_ob['lamby'])/dy)
pumping_ob['col'] = np.trunc((pumping_ob['lambx']-swx)/dx)

# Convert to integers
pumping_ob['row'] = pumping_ob.row.astype("int64")
pumping_ob['col'] = pumping_ob.col.astype("int64")

# Create a new column, containing head values simulated by the model
pumping_ob['simulated']=head[9,pumping_ob['row'],pumping_ob['col']]

#print(pumping_ob)

# Compare the observed and simulated head data
compare=pumping_ob
compare=compare.set_index('Head_ftAMS')
compare=compare.drop(['lambx','lamby','row','col'],axis=1)
#print(compare)

#--------------------------------------------------

# Calculate the mean error and mean absolute error
pumping_ob['error']=pumping_ob['Head_ftAMS']-pumping_ob['simulated'] #difference between observed and simulated heads
pumping_ob['absolute']=pumping_ob.error.abs() #find absolute error
mean=pumping_ob.error.mean() #find mean error
mean_abs = pumping_ob.absolute.mean()
print("Mean Error:", round(mean,1), 'ft')
print("Mean Absolute Error:", round(mean_abs,1), 'ft')

# Calculate the standard deviation of absolute error
sum1 = 0
for i in pumping_ob.absolute.values:
  sum1 += (i - mean_abs)**2
std = np.sqrt(sum1/len(pumping_ob.absolute.values))
print('Standard Deviation:', std)

#--------------------------------------------------

# Based on standard deviation, create a color coding for high and low errors

# Add a legend for colors
thres = std
col = np.where(pumping_ob['simulated'] + thres < pumping_ob['Head_ftAMS'], 'dodgerblue', np.where(pumping_ob['simulated'] - thres > pumping_ob['Head_ftAMS'], 'darkred', 'k'))
s_head = mp.lines.Line2D([],[],marker='o',color='white',label='small simulated head with high errors',markerfacecolor='dodgerblue')
l_head = mp.lines.Line2D([],[],marker='o',color='white',label='large simulated head with high errors',markerfacecolor='darkred')
m_head = mp.lines.Line2D([],[],marker='o',color='white',label='simulated head with low errors',markerfacecolor='k')

# Add a legend for lines
pm_sign = r'$\pm$'
dash_l = mp.lines.Line2D([],[],linewidth=1,linestyle='--',color='blue',label=pm_sign+' 1*Std')
degree_sign = u"\N{DEGREE SIGN}"
grey_l = mp.lines.Line2D([],[],linewidth=4,color='grey',label='45'+degree_sign+'line')

# Plot the calibration data
fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111)
plt.rcParams["figure.figsize"] = (5,5)
ax.scatter(compare.index,compare.simulated, c= col, marker = 'o', s = 55)
x = np.linspace(440,800)
y = x
ax.plot(x,y,linewidth=4,color = 'grey')
ax.plot(x, x+thres, '--b', linewidth=1)
ax.plot(x, x-thres, '--b', linewidth=1)
ax.set_xlim(450,800)
ax.set_ylim(450,800)
ax.set_xlabel('Observed Head (ft)', fontsize = 14)
ax.set_ylabel('Simulated Head (ft)', fontsize = 14)
ax.set_title('Will County Model Calibration', fontsize = 16)
# add first legend for the data
leg1 = ax.legend(handles=[s_head,m_head,l_head],loc='lower right')
# add second legend for dashed lines and hence the first will be removed
leg2 = ax.legend(handles=[dash_l,grey_l])
# Manually add the first legend back
ax.add_artist(leg1)
plt.show()
#print (compare)

#--------------------------------------------------

matrix = np.zeros((nrow,ncol))

error = pumping_ob['error']

x = pumping_ob['row']
y = pumping_ob['col']

#matrix[x.values[0]-1, y.values[0]-1] = error.values[0]
#matrix[x.values[1]-1, y.values[1]-1] = error.values[1]

for i in range(0,102):
    matrix[x.values[i]-1, y.values[i]-1] = error.values[i]



fig = plt.figure(figsize=(20, 21))
cset = plt.pcolor(matrix, cmap='seismic')
plt.clim(-15, 15) 


#This line generates the colorbar 
fig.colorbar(cset)

plt.show()

#--------------------------------------------------
'''Error map'''

# Define the area over which to plot data

# Create and populate new columns in the "pumping_ob" dataframe for the latitudes and longitudes of each observation well
for value in pumping_ob.index:
  pumping_ob.loc[value,'lat'], pumping_ob.loc[value,'long'] = pyproj.transform(illimap,wgs84,pumping_ob.loc[value,'lambx']*0.3048,pumping_ob.loc[value,'lamby']*0.3048) #must convert Lambert x/y from ft to m

#Conduct the Universal Kriging
UK = UniversalKriging(pumping_ob['long'], pumping_ob['lat'],pumping_ob['error'], variogram_model='spherical',nlags=6)

#Create xpoints and ypoints in space, with 0.01 spacing
xpoints = np.arange(sw_long,ne_long,0.01)
ypoints = np.arange(sw_lat,ne_lat,0.01)

# create a meshgrid with xpoints and ypoints, to be used later in the code
X,Y = np.meshgrid(xpoints,ypoints)

# calculate the interpolated grid and fill values.
z, var = UK.execute('grid', xpoints, ypoints)
z = z.filled(fill_value=None)

'''Small plot'''

# Create a figure
fig=plt.figure(figsize=(10,10)) #create 10 x 10 figure

# Define a projection
ax = plt.axes(projection=ccrs.PlateCarree())

# Define the spatial domain to plot
ax.set_xlim(sw_long,ne_long) #The spatial domain spans from minimum longitude to maximum longitude (x-direction)
ax.set_ylim(sw_lat,ne_lat) #The spatial domain spans from minimum longitude to maximum latitude (y-direction)

# Define a title for the plot
ax.set_title('Calculated error over the model area')

# Using CartoPy, import several features from Natural Earth Data (see https://www.naturalearthdata.com/features/)
# Large rivers:
largerivers = cf.NaturalEarthFeature(
    category='physical',
    name='rivers_lake_centerlines',
    scale='110m', # major rivers
    facecolor='none')
# Small rivers:
smallrivers = cf.NaturalEarthFeature(
    category='physical',
    name='rivers_lake_centerlines_scale_rank',
    scale='10m', # smaller rivers (still considered major by many/most people)
    facecolor='none')
# Smallest rivers:
smallestrivers = cf.NaturalEarthFeature(
    category='physical',
    name='rivers_north_america',
    scale='10m',
    facecolor='none')
# Population centers:
popplaces = cf.NaturalEarthFeature(
    category='cultural',
    name='urban_areas', # plots municipal boundaries
    scale='10m',
    facecolor='plum')
# Major roads:
majorroads = cf.NaturalEarthFeature(
    category='cultural',
    name='roads',
    scale='10m',
    facecolor='none')

# Add the features defined above to the plot, with various colors and markers
#ax.add_feature(popplaces,edgecolor='plum',linewidth=1.0) # Population centers will be light purple areas
ax.add_feature(largerivers,edgecolor='aqua',linewidth=2.0) # Large, small, and smallest rivers will be light blue lines
ax.add_feature(smallrivers,edgecolor='aqua',linewidth=2.0)
ax.add_feature(smallestrivers,edgecolor='aqua',linewidth=2.0)
ax.add_feature(majorroads,edgecolor='gray',linewidth=1.0) # Roads will be gray lines
ax.add_feature(USCOUNTIES.with_scale('5m'),linewidth=1.0) # County lines will be black lines

# Create contours from the interpolation
error_contour_levels = [-50,-40,-30,-20,-10,0,10,20,30,40,50] #set contour levels for contouring error
cset_contour = plt.contour(X,Y,z,error_contour_levels,colors='blue')

# Add in a continuous color flood
cset_fill=plt.imshow(z,vmin=-100,vmax=100,cmap=plt.cm.coolwarm,origin='lower',extent=[X.min(), X.max(), Y.min(), Y.max()])

# The following code could also be used, if a discrete color flood is preferred:
#cset_fill = plt.contourf(X,Y,z,vmin=-100,vmax=100,cmap=plt.cm.coolwarm)

# Label contours (makes use of pylab)
pylab.clabel(cset_contour, inline=1, fontsize=10,fmt='%1.0f',colors='black')

# Plot the points that were measured
points=plt.scatter(pumping_ob['long'], pumping_ob['lat'], marker=".", color="black", label="Data Points")

# Add a legend to the plot
cities_towns=mpatches.Patch(color='plum',label='Population Centers') #create a proxy artist to represent population centers on the legend
rivers=mlines.Line2D([],[],linewidth=2.0,color='aqua',label='Rivers') #proxy artist for large and small rivers
roads=mlines.Line2D([],[],linewidth=1.0,color='gray',label='Roads') #proxy artist for major roads
counties=mlines.Line2D([],[],linewidth=1.0,color='black',label='County Lines') #proxy artist for county lines
ax.legend(loc=[1.05,0],handles=[points,rivers,roads,counties])
#ax.legend(loc=[1.05,0],handles=[cities_towns,points,rivers,roads,counties]) #includes population centers

# Add in a colorbar
cax = fig.add_axes([0.94, 0.39, 0.25, 0.03])
fig.colorbar(cset_fill, cax=cax, orientation='horizontal',label='Error (ft difference)')
plt.show()

'''Large plot'''

# Create a figure
fig=plt.figure(figsize=(20,20)) #create 20 x 20 figure

# Define a projection
ax = plt.axes(projection=ccrs.PlateCarree())

# Define the spatial domain to plot
ax.set_xlim(sw_long,ne_long) #The spatial domain spans from minimum longitude to maximum longitude (x-direction)
ax.set_ylim(sw_lat,ne_lat) #The spatial domain spans from minimum longitude to maximum latitude (y-direction)

# Define a title for the plot
ax.set_title('Calculated error over the model area',fontsize=20)

# Using CartoPy, import several features from Natural Earth Data (see https://www.naturalearthdata.com/features/)
# Large rivers:
largerivers = cf.NaturalEarthFeature(
    category='physical',
    name='rivers_lake_centerlines',
    scale='110m', # major rivers
    facecolor='none')
# Small rivers:
smallrivers = cf.NaturalEarthFeature(
    category='physical',
    name='rivers_lake_centerlines_scale_rank',
    scale='10m', # smaller rivers (still considered major by many/most people)
    facecolor='none')
# Smallest rivers:
smallestrivers = cf.NaturalEarthFeature(
    category='physical',
    name='rivers_north_america',
    scale='10m',
    facecolor='none')
# Population centers:
popplaces = cf.NaturalEarthFeature(
    category='cultural',
    name='urban_areas', # plots municipal boundaries
    scale='10m',
    facecolor='plum')
# Major roads:
majorroads = cf.NaturalEarthFeature(
    category='cultural',
    name='roads',
    scale='10m',
    facecolor='none')

# Add the features defined above to the plot, with various colors and markers
#ax.add_feature(popplaces,edgecolor='plum',linewidth=1.0) # Population centers will be light purple areas
ax.add_feature(largerivers,edgecolor='aqua',linewidth=4.0) # Large, small, and smallest rivers will be light blue lines
ax.add_feature(smallrivers,edgecolor='aqua',linewidth=4.0)
ax.add_feature(smallestrivers,edgecolor='aqua',linewidth=4.0)
ax.add_feature(majorroads,edgecolor='gray',linewidth=3.0) # Roads will be gray lines
ax.add_feature(USCOUNTIES.with_scale('5m'),linestyle = '-', linewidth=5.0) # County lines will be black lines

# Create contours from the interpolation
error_contour_levels = [-50,-40,-30,-20,-10,0,10,20,30,40,50] #set contour levels for contouring error
cset_contour = plt.contour(X,Y,z,error_contour_levels,linewidths=2.5,colors='blue')

# Add in a continuous color flood
cset_fill=plt.imshow(z,vmin=-100,vmax=100,cmap=plt.cm.coolwarm,origin='lower',extent=[X.min(), X.max(), Y.min(), Y.max()])

# The following code could also be used, if a discrete color flood is preferred:
#cset_fill = plt.contourf(X,Y,z,vmin=-100,vmax=100,cmap=plt.cm.coolwarm)

# Label contours (makes use of pylab)
pylab.clabel(cset_contour, inline=1, fontsize=20,fmt='%1.0f',colors='black')

# Plot the points that were measured
points=plt.scatter(pumping_ob['long'], pumping_ob['lat'], marker="o", color="black", label="Data Points",)


# Add a legend to the plot
cities_towns=mpatches.Patch(color='plum',label='Population Centers') #create a proxy artist to represent population centers on the legend
rivers=mlines.Line2D([],[],linewidth=4.0,color='aqua',label='Rivers') #proxy artist for large and small rivers
roads=mlines.Line2D([],[],linewidth=3.0,color='gray',label='Roads') #proxy artist for major roads
counties=mlines.Line2D([],[],linewidth=5.0,color='black', linestyle = '-', label='County Lines') #proxy artist for county lines
ax.legend(loc=[1.01,0],handles=[points,rivers,roads,counties])
#ax.legend(loc=[1.01,0],handles=[cities_towns,points,rivers,roads,counties]) #includes population centers

# Location of communities
lons = [-88.0987, -88.0684, -88.0817, -88.0895, -88.2287]
lats = [41.5548, 41.6986, 41.5250, 41.6475, 41.4295]
labels = ['Crest Hill', 'Bolingbrook', 'Joliet', 'Romeoville', 'Channahon']
plt.scatter(lons, lats, marker="*", color="gold", zorder=10)


# This puts the labels on the map
#for i, txt in enumerate(labels):
  #ax.annotate(txt,(lons[i], lats[i]),  fontsize=18, zorder=5, )


# Add in a colorbar
cax = fig.add_axes([0.91, 0.39, 0.03, 0.3])
fig.colorbar(cset_fill, cax=cax, orientation='vertical',label='Error (ft difference)')
plt.show()


#Display minimum and maximum error, for reference
print(pumping_ob.error.min())
print(pumping_ob.error.max())

#%% ZONE BUDGET

#--------------------------------------------------
'''Extract binary data from budget files'''

# Read the cell budget file, a binary output file of the model
cbc_f = flopy.utils.binaryfile.CellBudgetFile(modelname+'.cbc')
#cbc_f.list_records()

#--------------------------------------------------
''' Create and display a ZoneBudget object based on a zone input file'''

# Import a zone input file
# The zone input file is a text file that assigns a zone to each cell in the model
uploaded = files.upload() #upload the most recent zone input file.  Be sure to match the file name below to the actual name of the file you would like to upload.
zon = flopy.utils.zonbud.read_zbarray('zone_file_65x65.txt') #read_zbarray is used to read in the zone input file as an array

del uploaded

# Create a ZoneBudget object and get the budget record array
# set zone 1 as "upper half"; zone 2 as "lower half"
aliases = {1:'Upper half', 2:'Lower half'}
zb = flopy.utils.ZoneBudget(cbc_f, zon, kstpkper=(0, 0), aliases=aliases) #kstpkper specifies the timestep and stress period of interest
zb.get_budget()
#zb.get_record_names()

# Note:  If running this cell produces an error, try re-running the cell again.

# Note:  The zone input file only needs to be uploaded once.  Check the "Files" tab on the left to see if it has already been uploaded.
# If the zone input file has already been uploaded, you can click "Cancel upload" in the output window below.

#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------