# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 09:31:12 2020

@author: antho
"""

#%%
import os
os.environ['GDAL_DATA'] = r'D:\anaconda3\Library\share\gdal'

#%% IMPORT PACKAGES

#--------------------------------------------------
import sys
import platform
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import flopy #import FloPy to develop, run, and analyze the model
from flopy.utils import Raster #raster toolbox built into FloPy
from pyproj import Transformer #transform different coordinate systems and projections
import pandas as pd
import numpy as np

#%% MODEL SETUP

#--------------------------------------------------
'''Set up a workspace for the model'''

# Model input files and output files will reside in the model workspace
model_name = 'Will_Cty_Shallow'
dir_path = r'D:\Documents\GitHub'
workspace = os.path.join(dir_path, model_name)
if not os.path.exists(workspace):
    os.makedirs(workspace)

data_path = os.path.join(dir_path, 'data', 'test001')
if not os.path.exists(data_path):
    os.makedirs(data_path)
assert os.path.isdir(data_path)

#--------------------------------------------------
'''Create a simulation'''

sim = flopy.mf6.MFSimulation(sim_name=model_name, version='mf6', 
                             exe_name=r'D:\Documents\ISWS\downloads from USGS\mf6.1.1\bin\mf6.exe',
                             sim_ws=workspace)

#--------------------------------------------------
'''Create the MODFLOW model object'''

#The model object is named "gwf"
gwf = flopy.mf6.ModflowGwf(sim, modelname=model_name,
                           model_nam_file='{}.nam'.format(model_name))
gwf.name_file.save_flows = True

# Create iterative model solution and register the gwf model with it
ims = flopy.mf6.ModflowIms(sim, pname='ims', print_option='SUMMARY', 
                           complexity='SIMPLE', outer_hclose=0.0001, 
                           outer_maximum=500, under_relaxation='NONE', 
                           inner_maximum=100, inner_hclose=0.0001, 
                           rcloserecord=0.001, linear_acceleration='CG', 
                           scaling_method='NONE', reordering_method='NONE', 
                           relaxation_factor=0.97)
sim.register_ims_package(ims, [gwf.name])

#--------------------------------------------------
'''Define model domain in lat/long coordinates'''

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
'''Define spatial and temporal discretization'''

# Assign discretization variables
Lx = nex-swx # Width of the model domain in feet
Ly = ney-swy # Height of the model domain in feet
nlay = 10 # Number of model layers
dx = 2000 
dy = 2000
nrow = int(Ly/dy) # Number of rows
ncol = int(Lx/dx) # Number of columns

nper = 1 #specify number of stress periods
tdis_rc = [(1.0, 1, 1.0)]   #input for tdis package: [period length, # of time steps, and multiplier for length of successive time steps].
                            #for steady state, a time step of 1 is recommended; period length does not matter

#--------------------------------------------------
'''Define river elevations'''

# Import river stage, lambert x, lambert y from river Excel file
dfriv = pd.read_csv(os.path.join(dir_path,'ShallowDolomite_Group','large_files','rivers_625.csv'))

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
'''Define top and bottom elevations'''

# Now load the raster using FloPy's built in Raster toolbox
illinoisdem = Raster.load(r'D:\Documents\GitHub\ShallowDolomite_Group\large_files\landsurface_el.tif')
bedrock = Raster.load(r'D:\Documents\GitHub\ShallowDolomite_Group\large_files\bedrock_el.tif')

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
'''Assign hydraulic conductivity'''

# Assign hydraulic conductivity in ft/day
kc = 220 #predominantly coarse
kf = .3 #predominantly fine
kb = 20 #bedrock

# Define the threshold to determine how hydraulic conductivity is assigned
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

kl1 = kloader(r'D:\Documents\GitHub\ShallowDolomite_Group\large_files\percentl1.tif',kc,kf,threshold)
kl2 = kloader(r'D:\Documents\GitHub\ShallowDolomite_Group\large_files\percentl2.tif',kc,kf,threshold)
kl3 = kloader(r'D:\Documents\GitHub\ShallowDolomite_Group\large_files\percentl3.tif',kc,kf,threshold)
kl4 = kloader(r'D:\Documents\GitHub\ShallowDolomite_Group\large_files\percentl4.tif',kc,kf,threshold)
kl5 = kloader(r'D:\Documents\GitHub\ShallowDolomite_Group\large_files\percentl5.tif',kc,kf,threshold)
kl6 = kloader(r'D:\Documents\GitHub\ShallowDolomite_Group\large_files\percentl6.tif',kc,kf,threshold)
kl7 = kloader(r'D:\Documents\GitHub\ShallowDolomite_Group\large_files\percentl7.tif',kc,kf,threshold)
kl8 = kloader(r'D:\Documents\GitHub\ShallowDolomite_Group\large_files\percentl8.tif',kc,kf,threshold)
kl9 = kloader(r'D:\Documents\GitHub\ShallowDolomite_Group\large_files\percentl9.tif',kc,kf,threshold)
kl10 = kl9-kl9+kb

khlayers = [kl1,kl2,kl3,kl4,kl5,kl6,kl7,kl8,kl9,kl10]
kvlayers=np.divide(khlayers,10.)

#--------------------------------------------------
'''Spatial and temporal discretization packages'''

# Create the discretization package
dis = flopy.mf6.ModflowGwfdis(model=gwf, pname='dis', nlay=nlay, nrow=nrow, ncol=ncol, 
                              delr=dx, delc=dy, top=topgrid, 
                              botm=botgrids, 
                              filename='{}.dis'.format(model_name))

# Create tdis package
# Discretization information for time is read from the input file for the Timing Module (TDIS)
tdis = flopy.mf6.ModflowTdis(sim, pname='tdis', time_units='DAYS', 
                             nper=nper, perioddata=tdis_rc)

#--------------------------------------------------
'''Initial conditions'''

ic = flopy.mf6.ModflowGwfic(gwf, pname='ic', strt=50.0,
                            filename='{}.ic'.format(model_name))

#--------------------------------------------------
'''Node property flow'''

npf = flopy.mf6.ModflowGwfnpf(gwf, pname='npf', icelltype=1, k=khlayers,
                              k33=kvlayers, save_flows=True)

# When k_horizontal != k_vertical, k=k_horizontal and k33=k_vertical

#--------------------------------------------------
'''Output control'''

oc = flopy.mf6.ModflowGwfoc(gwf, pname='oc', budget_filerecord='{}.cbb'.format(model_name),
                            head_filerecord='{}.hds'.format(model_name),
                            headprintrecord=[('COLUMNS', 10, 'WIDTH', 15,
                                              'DIGITS', 6, 'GENERAL')],
                            saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')])#,
                            #printrecord=[('HEAD', 'FIRST'), ('HEAD', 'LAST'), 
                                         #('BUDGET', 'LAST')])

#--------------------------------------------------
'''Storage package'''

sy = flopy.mf6.ModflowGwfsto.sy.empty(gwf, layered=True)
for layer in range(0,10):
    sy[layer]['data'] = 0.2
    
ss = flopy.mf6.ModflowGwfsto.ss.empty(gwf, layered=True, default_value=0.000001)

sto = flopy.mf6.ModflowGwfsto(gwf, pname='sto', save_flows=True, iconvert=1, 
                              ss=ss, sy=sy, steady_state={0:True})

#%% BOUNDARY CONDITIONS

#--------------------------------------------------
'''Define wells'''

# Import well data from .csv file
dfwel = pd.read_csv(os.path.join(dir_path,'ShallowDolomite_Group','pumping','2002_pumping_V2.csv'))
dfwel = dfwel.set_index('p_num') #assign index as p_number so that other columns can be deleted

# Trim dataframe with well information to the model domain
dfwel = dfwel.loc[dfwel['lam_x']<nex]
dfwel = dfwel.loc[dfwel['lam_y']<ney]
dfwel = dfwel.loc[dfwel['lam_x']>swx]
dfwel = dfwel.loc[dfwel['lam_y']>swy]

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
'''Put the well data into a form that can be added into the well package'''

maxbound = 1000000 #maximum number of well cells that will be specified in a given stress period

# Create a stress_period_data template
stress_period_data_wel = flopy.mf6.ModflowGwfwel.stress_period_data.empty(gwf, maxbound=maxbound, boundnames=False, 
                                                                      aux_vars=None, stress_periods=[0])

for value in dfwel.index:
    stress_period_data_wel[0].append = ((dfwel.loc[value,'lay'],dfwel.loc[value,'row'],dfwel.loc[value,'col']),dfwel.loc[value,'flux'])

wel = flopy.mf6.ModflowGwfwel(gwf, pname='wel', print_input=True, print_flows=True,
                              auxiliary=None, maxbound=maxbound, stress_period_data=stress_period_data_wel,
                              boundnames=False, save_flows=True)

#--------------------------------------------------
'''Define drains'''

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
'''Delete any drain cells that coincide with river cells'''

# This code block ensures that no river cells are accidentally included in the drain package

# Create an array from the values in the drain dataframe
# In other words, create "ardrn" from "dfdrn"
ardrn = dfdrn.values
ardrn = ardrn.tolist() #convert the drain array from an array to a list

# Print the initial length of the drain array
ardrn_initial=len(ardrn)
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

# Convert the drain array back to a dataframe (ardrn -> dfdrn)
dfdrn = pd.DataFrame(ardrn, columns = ['lay','row','col','elevation','cond'])

#--------------------------------------------------
'''Put the drain data into a form that can be added into the drain package'''

maxbound = 1000000 #maximum number of drain cells that will be specified in a given stress period

# Create a stress_period_data template
stress_period_data_drn = flopy.mf6.ModflowGwfdrn.stress_period_data.empty(gwf, maxbound=maxbound, boundnames=False, 
                                                                      aux_vars=None, stress_periods=[0])

for value in dfdrn.index:
    stress_period_data_drn[0].append = ((0,dfdrn.loc[value,'row'],dfdrn.loc[value,'col']),dfdrn.loc[value,'elevation'],dfdrn.loc[value,'cond'])

drn = flopy.mf6.ModflowGwfdrn(gwf, pname='drn', print_input=True, print_flows=True,
                              auxiliary=None, maxbound=maxbound, stress_period_data=stress_period_data_drn,
                              boundnames=False, save_flows=True)

#--------------------------------------------------
'''Define recharge'''

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

#--------------------------------------------------
'''Put the recharge data into a form that can be added into the recharge package'''

maxbound = 1000000 #maximum number of drain cells that will be specified in a given stress period

# Create an array of stress period data

rch_array={}
rch_array[0] = recharge
#rch = flopy.mf6.ModflowGwfrch(gwf, pname='rch', fixed_cell=True, print_input=True,
                              #maxbound=maxbound, stress_period_data=rch_array)

rch = flopy.mf6.ModflowGwfrcha(gwf, pname='rch', fixed_cell=True, print_input=True,
                               recharge=recharge)

#%% Plot model inputs (boundary conditions, elevations)

#--------------------------------------------------
''''Plot grid and boundary conditions'''

'''
plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.PlotMapView(model=gwf, layer=0)
#grid = modelmap.plot_grid()
ib = modelmap.plot_ibound(ibound=ibound)
rvr = modelmap.plot_bc('RIV')
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
'''

#--------------------------------------------------
'''Plot hydraulic conductivity'''

'''
plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.map.PlotMapView(model=gwf, layer=0) #use plotmapview to attach plot to model
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
plt.legend(handles=[mp.patches.Patch(color='green', label='rivers'),
                    mp.patches.Patch(color='saddlebrown', label='kf=0.3 ft/d'),
                    mp.patches.Patch(color='lightgoldenrodyellow', label='kc=280 ft/d') ],
           bbox_to_anchor=(1.7,0.8))
plt.show()
'''

#--------------------------------------------------
'''Plot recharge'''

'''
fig=plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.map.PlotMapView(model=gwf, layer=0) #use plotmapview to attach plot to model
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
'''

#--------------------------------------------------
'''Plot transects'''

'''
plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelxsect = flopy.plot.PlotCrossSection(model = gwf, line={"row":22}) #use plotmapview to attach plot to model. row indicates west to east 
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
'''

#--------------------------------------------------
'''Plot Transmissivity'''

'''
# Horizontal T = sum(k*b), for each layer.  k is the hydraulic conductivity, and b is the layer thickness
Thz=sum([k for k in khlayers])*laythick-kl10*laythick+kl10*bothick
#----------------------------------------------------------------------------
plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.map.PlotMapView(model=gwf, layer=0) #use plotmapview to attach plot to model

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
'''

#%% RUN THE MODEL

# Change where to save simulation
#sim.simulation_data.mfpath.set_sim_path(run_folder)

# Write simulation to new location
sim.write_simulation()

# Print a list of the files that were created in the model workspace
print(os.listdir(workspace))

# Run the simulation
success, buff = sim.run_simulation()
print('\nSuccess is: ', success)

#%% MODEL RESULTS

#--------------------------------------------------
'''Post-process head results'''

'''
# Read the binary head file and plot the results
# We can use the existing Flopy HeadFile class because
# the format of the headfile for MODFLOW 6 is the same
# as for previous MODFLOW verions
headfile = '{}.hds'.format(model_name)
fname = os.path.join(workspace, headfile)
hds = flopy.utils.binaryfile.HeadFile(fname)
h = hds.get_data()

# We can also use the Flopy PlotMapView capabilities for MODFLOW 6
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, aspect='equal')

# Next we create an instance of the ModelMap class
modelmap = flopy.plot.PlotMapView(model=gwf, ax=ax)

# Then we can use the plot_grid() method to draw the grid
# The return value for this function is a matplotlib LineCollection object,
# which could be manipulated (or used) later if necessary.
#quadmesh = modelmap.plot_ibound(ibound=ibd)
linecollection = modelmap.plot_grid()
contours = modelmap.contour_array(h[0])
'''

#--------------------------------------------------
'''Post-process flow results'''

'''
# Read the binary grid file
fname = os.path.join(workspace, '{}.dis.grb'.format(model_name))
bgf = flopy.utils.mfgrdfile.MfGrdFile(fname)

# Data read from the binary grid file is stored in a dictionary
bgf._datadict
'''