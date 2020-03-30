#----------------------------------------------------------------------------
#
# Created on Sun Mar 29 15:06:36 2020
# @author: dbabrams
# Goal for this week: define model domain and river network
# Unless otherwise stated, all units are in feet and days.
#
#----------------------------------------------------------------------------



'''Import packages.'''
#----------------------------------------------------------------------------
import flopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import pyproj
import pandas as pd
#----------------------------------------------------------------------------



'''Create a MODFLOW model object and run with MODFLOW 2005.'''
#----------------------------------------------------------------------------
modelname = "my_model"
m = flopy.modflow.Modflow(modelname, exe_name = 'mf2005')
#----------------------------------------------------------------------------



'''Create the Discretization package'''
#----------------------------------------------------------------------------
# Define model domain in lat/long coordinates
sw_lat = 39.9727 #southwest latitude
sw_long = -90.537 #southwest longitude
ne_lat = 40.6657 #northeast latitude
ne_long = -89.1371 #northeast longitude

# In Illinois, the Illimap projection is used to minimize distortion
# See https://www.spatialreference.org/ref/sr-org/7772/ for details
# Also see http://library.isgs.illinois.edu/Pubs/pdfs/circulars/c451.pdf
# Values originate from here (https://www.spatialreference.org/ref/sr-org/7772/html/)
D = {'proj': 'lcc', # define projection as Lambert Conformal Conic
        'ellps': 'clrk66', # Use the Clarke 1866 ellipsoid
        'lon_0': -89.5, #Central Meridian
        'lat_0': 33, #Latitude of Origin
        'lat_1': 33, #Standard Parallel 1
        'lat_2': 45, #Standard Parallel 2
        'x_0': 2999994*0.3048006096012192, # starting x coordinate is in feet, Python expects meters
        'y_0': 0} # starting y coordinate}

prj = pyproj.Proj(D) # Create a projection object that will be used to convert lat/long to illimap

#Define the northeastern coordintes, round to nearest 10,000
nex, ney = prj(ne_long, ne_lat) # this will return meters
nex, ney = round(nex/0.3048006096012192,-4), round(ney/0.3048006096012192,-4) # convert to feet

#Define the southwestern coordinates, round to nearest 10,000
swx, swy = prj(sw_long, sw_lat) # this will return meters
swx, swy = round(swx/0.3048006096012192,-4), round(swy/0.3048006096012192,-4) # convert to feet

# Assign Discretization variables
Lx = nex-swx # Width of the model domain
Ly = ney-swy # Height of the model domain
ztop = 0. # Model top elevation
zbot = -50. # Model bottom elevation
nlay = 1 # Number of model layers
dx = 2500 # grid spacing (x-direction)
dy = 2500 # grid spacing (y-direction)
nrow = int(Ly/dy) # Number of rows
ncol = int(Lx/dx) # Number of columns
nper = 1 #specify number of stress periods
steady = [True] #specify if stress period is transient or steady-state

# create flopy discretization object
# length and time are feet (1) and days (4).
# See https://water.usgs.gov/ogw/modflow/MODFLOW-2005-Guide/index.html?dis.htm 
dis = flopy.modflow.ModflowDis(model=m, nlay=nlay, nrow=nrow, ncol=ncol, 
                               delr=dx, delc=dy, top=ztop, botm=zbot, 
                               itmuni = 4, lenuni = 1, 
                               nper=nper, steady=steady)
#----------------------------------------------------------------------------


'''Create the Basic Package, which contains ibound and starting heads'''
#----------------------------------------------------------------------------
# Create ibound as array of ints (1), indicating all cells are active
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
#ibound[:, :, 0] = -1 # Designate left boundary cells as constant head
#ibound[:, :, -1] = -1 # Designate right boundary cells as constant head

# Create starting head array, must be floats.
strt = 500*np.ones((nlay, nrow, ncol), dtype=np.float32) #set every cell to 5.0
#strt[:, :, 0] = 10. #set left side head to 10 ft
#strt[:, :, -1] = 0. #set right side head to 0 ft

#Create flopy bas object
bas = flopy.modflow.ModflowBas(m, ibound=ibound, strt=strt)
#----------------------------------------------------------------------------



'''Create the Layer Property Flow Package, which contains information about
hydruaulic conductivity and other information about how to calculate flow'''
#----------------------------------------------------------------------------
hk = 500.*np.ones((nlay,nrow,ncol), dtype=np.float32) #define horizontal hydraulic conductivity
vk = np.ones((nlay,nrow,ncol), dtype=np.float32) #define vertical hydraulic conductivity

#define layer type as convertible (1), must be an integer
#for more information, see https://water.usgs.gov/ogw/modflow/MODFLOW-2005-Guide/index.html?dis.htm
laytyp = np.ones((nlay,), dtype=np.int32)

# create the LPF object
lpf = flopy.modflow.ModflowLpf(model=m, hk=hk, vka=vk, laytyp=laytyp, ipakcb=1)
#----------------------------------------------------------------------------



'''Create the recharge package'''
#----------------------------------------------------------------------------
# assign a uniform recharge of 0.001 ft/day
rch = flopy.modflow.mfrch.ModflowRch(model=m,nrchop=3,rech=1e-3)
#----------------------------------------------------------------------------



'''Create the river package'''
#----------------------------------------------------------------------------
# import stage, stream order, stream length, and centroid coordinate of rivers
dfriv = pd.read_csv('ModelGrid_2500ft_river_cells_centroids.csv')
# trim dataframe with river information to model domain
dfriv = dfriv.loc[dfriv['lamx']<nex]
dfriv = dfriv.loc[dfriv['lamx']>swx]
dfriv = dfriv.loc[dfriv['lamy']<ney]
dfriv = dfriv.loc[dfriv['lamy']>swy]

#assign all rivers to the upper layer
dfriv['lay'] = 0
#convert lamy to row and lamx to column
#note that rows count from top fo bottom
dfriv['row'] = np.trunc((ney-dfriv['lamy'])/dy)
dfriv['col'] = np.trunc((dfriv['lamx']-swx)/dy)

# define the river stage
dfriv['stage'] = dfriv['rvr_stg']
# define the conductance
dfriv['cond'] = 5000 #ft^2/d
# define the river bottom
dfriv['bot'] = dfriv['stage']-3
# drop unneeded files
dfriv = dfriv.drop(['STR_ORD_MI','STR_ORD_MA','SUM_LENGTH','rvr_stg','lamx','lamy'],axis=1)

# create an array and a dictionary containing that array for stress period zero
arriv = dfriv.values # note that the dataframe was developed in the correct order
riverdata = {0: arriv} # stress period 0 is the key
# write the river object
flopy.modflow.mfriv.ModflowRiv(model=m, ipakcb=None, stress_period_data=riverdata)
#----------------------------------------------------------------------------



'''Create the well package'''
#----------------------------------------------------------------------------
# assign the well data to a dictionary for the first stress period
welldata = {0: [0,10,20,-133000]} # well in layer 1, row 6, col 6, pumping -500 cfd
# Create the Well Object
flopy.modflow.mfwel.ModflowWel(model=m, ipakcb=None, stress_period_data=welldata)
#----------------------------------------------------------------------------



'''Create the Output Control Package'''
#----------------------------------------------------------------------------
#create oc stress period data. 
spd = {(0, 0): ['print head', 'print budget', 'save head', 'save budget']}
#create output control object
oc = flopy.modflow.ModflowOc(model=m, stress_period_data=spd, compact=True)
#----------------------------------------------------------------------------



'''Create the PCG Solver Object'''
#----------------------------------------------------------------------------
# for the time being, we will use default settings with the solver
pcg = flopy.modflow.ModflowPcg(model=m)
#----------------------------------------------------------------------------



'''Write MODFLOW input files.'''
#----------------------------------------------------------------------------
m.write_input()
#----------------------------------------------------------------------------



'''Run the model'''
#----------------------------------------------------------------------------
# Executve the model run
success, mfoutput = m.run_model(pause=False, report=True)
# Report back if the model did not successfully complete
if not success:
    raise Exception('MODFLOW did not terminate normally.')
#----------------------------------------------------------------------------
    
    
    
'''Extract binary data from head and flow files'''
#----------------------------------------------------------------------------
#extract binary data from head file as flopy head object
headobj = flopy.utils.binaryfile.HeadFile(modelname+'.hds')
#extract head data from head object
head = headobj.get_data(totim=1.0)

#extract binary data from budget file as flopy budget object
budgobj = flopy.utils.binaryfile.CellBudgetFile(modelname+'.cbc')
#extract flow data from budget object, define face over which flow is reported
frf = budgobj.get_data(text='flow right face', totim=1.0)
fff = budgobj.get_data(text='flow front face', totim=1.0)
#----------------------------------------------------------------------------



'''Plot grid and boundary conditions'''
#----------------------------------------------------------------------------
plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.PlotMapView(model=m, layer=0)
grid = modelmap.plot_grid()
ib = modelmap.plot_ibound()
wll = modelmap.plot_bc(ftype='WEL')
rvr = modelmap.plot_bc(ftype='RIV')
#add labels and legend
plt.xlabel('Lx (ft)',fontsize = 14)
plt.ylabel('Ly (ft)',fontsize = 14)
plt.title('Ibound', fontsize = 15, fontweight = 'bold')
plt.legend(handles=[mp.patches.Patch(color='blue',label='Const. Head',ec='black'),
                   mp.patches.Patch(color='white',label='Active Cell',ec='black'),
                   mp.patches.Patch(color='black',label='Inactive Cell',ec='black')],
                   bbox_to_anchor=(1.5,1.0))
#----------------------------------------------------------------------------



'''Plot results'''
#----------------------------------------------------------------------------
plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.map.PlotMapView(model=m, layer=0) #use plotmapview to attach plot to model
#grid = modelmap.plot_grid() #plot model grid
#contour_levels = np.linspace(head[0].min(),head[0].max(),11) #set contour levels for contouring head
contour_levels = np.linspace(450,650,21) #set contour levels for contouring head
head_contours = modelmap.contour_array(head, levels=contour_levels) #create head contours
plt.clabel(head_contours, inline=True,fontsize=12,fmt='%1.0f')
#flows = modelmap.plot_discharge(frf[0], fff[0], head=head) #create discharge arrows
rvr = modelmap.plot_bc(ftype='RIV')

#display parameters
plt.xlabel('Lx (ft)',fontsize = 14)
plt.ylabel('Ly (ft)',fontsize = 14)
plt.title('Steady-State Model, Flow(ft^3/d) and Head(ft) Results', fontsize = 15, fontweight = 'bold')

#----------------------------------------------------------------------------