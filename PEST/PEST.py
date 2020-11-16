# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 13:19:40 2020

@author: antho
"""

#%%
import os
os.environ['GDAL_DATA'] = r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\gdal'

#%%
import platform
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import flopy
import pyemu

#%%
os.chdir(r'D:\Documents\GitHub\ShallowDolomite_Group\PEST')

nam_file = 'my_model.nam'
org_model_ws = r'D:\Documents\GitHub'
exe_name = r'\\pri-fs1.ad.uillinois.edu\SWSGWmodeling\FloPy_Models\shallow_model\executables\mf2005.exe' #define the file path for the mf2005 executable
new_model_ws = 'new_model_ws'

m = flopy.modflow.Modflow.load(nam_file, exe_name=exe_name, model_ws=org_model_ws)
m.write_input()

#%%
m.run_model()

print(m.get_package_list())

#%% PARAMETERIZATION

# Pilot point multiplier parameters
#pp_props = [["upw.sy",0], ["rch.rech",None]]
pp_props = [["rch.rech",None]]

# Constant parameters (uniform value multipliers)
const_props = []
#for iper in range(m.nper): # recharge for past and future
    #const_props.append(['rch.rech', iper])
#for k in range(m.nlay):
    #const_props.append(["upw.hk",k])
    #const_props.append(["upw.ss",k])

# Grid-scale parameters
grid_props = []
#grid_props = [["upw.sy",0],["upw.vka",1]]
#for k in range(m.nlay):
    #grid_props.append(["upw.hk",k])

# Zones
# The following example file defines zones 1 and 2 as the northern and 
# southern halves of the model domain, respectively
#zn_array = np.loadtxt(r'D:\Documents\GitHub\ShallowDolomite_Group\PEST\example.zones')
#zone_props = [["upw.ss",0], ["rch.rech",0],["rch.rech",1]]
#k_zone_dict = {k:zn_array for k in range(m.nlay)}]

zone_props = []

# Boundary conditions
bc_props = []
#for iper in range(m.nper):
    #bc_props.append(['wel.flux', iper])

#%% OBSERVATIONS

# Here were are building a list of stress period, layer pairs (zero-based) 
# that we will use to setup obserations from every active model cell for a 
# given pair

hds_kperk = []
for iper in range(m.nper):
    for k in range(m.nlay):
        hds_kperk.append([iper,k])

#%% CONSTRUCT PEST INTERFACE
#mfp_boss = pyemu.helpers.PstFromFlopyModel(nam_file,new_model_ws,org_model_ws=temp_model_ws, 
                                           #pp_props=pp_props,spatial_list_props=bc_props, 
                                           #zone_props=zone_props,grid_props=grid_props, 
                                           #const_props=const_props,k_zone_dict=k_zone_dict, 
                                           #remove_existing=True,pp_space=4,sfr_pars=True, 
                                           #sfr_obs=True,hds_kperk=hds_kperk)

mfp = pyemu.helpers.PstFromFlopyModel(model=nam_file, new_model_ws=new_model_ws, org_model_ws=org_model_ws, 
                                      pp_props=pp_props, const_props=const_props, spatial_list_props=bc_props, 
                                      zone_props=zone_props, grid_props=grid_props, remove_existing=True, 
                                      sfr_pars=False, sfr_obs=False, hds_kperk=hds_kperk, 
                                      use_pp_zones=True)