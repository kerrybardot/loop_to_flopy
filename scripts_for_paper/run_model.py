# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 10:55:53 2024

@author: 19904604
"""
# run_model.py
import pandas as pd
import numpy as np
import subprocess
import sys
import os

subprocess.run(['python', 'loop_showcase_functions.py'])
subprocess.run(['python', '../../MODFLOW_Tools/modelling_routines.py'])
subprocess.run(['python', '../../MODFLOW_Tools/meshing_routines.py'])
subprocess.run(['python', '../../MODFLOW_Tools/geomodel_routines.py'])
subprocess.run(['python', '../../MODFLOW_Tools/loop2flopy.py'])

# ------ CREATE PROJECT CLASS -------#
sys.path.append('../../MODFLOW_Tools')
from loop2flopy import Project

P = Project('showcase', boundingbox = [0, 6000, 0, 6000, -500, 0]) # (x0, x1, y0, y1, z0, z1)

# ------------- DISCRETISATION PARAMETERS -----------------------#
P.triExeName = '../exe/triangle.exe'
P.workspace = '../modelfiles/'
P.crs = "epsg:28350" # Put in coordinate reference
P.xypumpbores = [(2000, 2000), (2500, 2000), (2000, 2500), (2500, 2500),] 
P.idpumpbores = ['P0','P1','P2','P3'] 
P.xyobsbores = [(1500,1000), (2500, 5000), (4000, 3000), (1500, 3500), (4500,4500)] 
P.idobsbores = ['OB1', 'OB2', 'OB3', 'OB4', 'OB5'] 
P.nobs = len(P.xyobsbores)
P.nzobs = 3
P.npump = len(P.xypumpbores)

P.r = 40        # refining factor for model boundary. High r has refined edges
P.w = 100      # Interior boundary
P.boundmaxtri = 30000 
P.modelmaxtri = 30000 # 10000 for ref
P.angle = 34   # minimum triangle angles
P.radius1 = 100 # radius of inner circle around pumping bores
P.radius2 = 200 # radius of outer circle around pumping bores
P.boremaxtri = 100000

# ------------- FLOW MODEL PARAMETERS -----------------------#
P.mfexe_name = '../exe/mf6.exe'
P.nlg = 4    # number of geological layers
P.nls = 2    # Number sublayers for conformable
P.geo_pl = 2 # Which geological layer pumping from (zero-based)
P.res = 2    # vertical resolution upon which voxel grid created to pick lithology bottoms

P.rch = 0.4/365 # 0.0027 m/d
P.strt = 0.

P.chfunc = lambda x,z: -(0.005*x) - (z * 0.02)-20 # horizontal gradient of 0.01 and vertical gradient of 0.02
P.xt3d = True

past_years = 2
P.nts_past = past_years * 6
P.tdis_past = [(past_years * 365, P.nts_past, 1.1)] # period length, number of timesteps, tsmult
P.qwell_past = -2000 #m3/d 
future_years = 5
P.nts_future = future_years * 6
P.tdis_future = [(future_years * 365, P.nts_future, 1.1)] # period length, number of timesteps, tsmult
P.qwell_future = -5000 #m3/d 
  
if not os.path.isdir('../results'): os.makedirs('../results', exist_ok=True)
if not os.path.isdir('../modelfiles'): os.makedirs('../modelfiles', exist_ok=True)

# UNPICKLE ZOBS AND WELL SCREENS
import pickle
pickleoff = open('../results/zobs.pkl','rb') 
P.zobs = pickle.load(pickleoff)
pickleoff.close()
pickleoff = open('../results/wel_screens.pkl','rb') 
P.wel_screens = pickle.load(pickleoff)
pickleoff.close()

# MESHING IN 2D
sys.path.append('../../MODFLOW_Tools')
from meshing_routines import createcell2d

P.cell2dvor, P.xcycvor, P.verticesvor, P.vor, vornodes = createcell2d(P, grid = 'vor', fault = False)     

# READIN GEOMODEL DATA
data = pd.read_excel("../data/loop_showcase_data.xls",sheet_name = "pinchout_example")
strat = pd.read_excel("../data/loop_showcase_data.xls",sheet_name = "strat")

pdf = pd.read_csv('../modelfiles/parameters.txt').set_index('parname')#, sep = 'None')#, sep = '/s+')
#print(pdf)#.parval)#['value'])#.to_list())
print(pdf.loc['hk0'].squeeze())

P.hk = [pdf.loc['hk0'].squeeze(), pdf.loc['hk1'].squeeze(), pdf.loc['hk2'].squeeze(), pdf.loc['hk3'].squeeze()]
P.vk = [pdf.loc['vk0'].squeeze(), pdf.loc['vk1'].squeeze(), pdf.loc['vk2'].squeeze(), pdf.loc['vk3'].squeeze()]
P.ss = [0.00009, pdf.loc['ss1'].squeeze(), pdf.loc['ss2'].squeeze(), pdf.loc['ss3'].squeeze()]
P.sy = [pdf.loc['sy0'].squeeze(), 0.1, 0.1,0.1]
P.control_points = (['CP1', 'control', 3000, 3000, -pdf.loc['cp'].squeeze(), -200, 'c', 'lower', 0., 0., 1.],) 

# TRUTH
#hk0, hk1, hk2, hk3, vk0, vk1, vk2, vk3 = 1.7, 0.07, 8.2, 0.05, 0.12,0.007,0.51,0.005
#ss1, ss2, ss3, sy0, cp = 0.00007, 0.00002, 0.00008, 0.12, 200
# PEST
#hk0, hk1, hk2, hk3 = 4.15535, 0.0634193, 1, 0.01
#vk0, vk1, vk2, vk3 = 0.16729,0.029757,2.61737,0.014614
#ss1, ss2, ss3, sy0, cp = 0.0001,6.01E-05,0.0001,0.2,8.66717
#P.hk = [hk0, hk1, hk2, hk3]
#P.vk = [vk0, vk1, vk2, vk3]
#P.ss = [0.00009, ss1, ss2, ss3]
#P.sy = [sy0, 0.1, 0.1,0.1]
#P.control_points = (['CP1', 'control', 3000, 3000, -cp, -200, 'c', 'lower', 0., 0., 1.],) 

from loop_showcase_functions import prepare_geomodel_loopshowcase
from loop_showcase_functions import create_geomodel_loopshowcase

P.data, P.strat = prepare_geomodel_loopshowcase(P, data, strat, include_fault = False)   # Prepare geomodel inputs
P.geomodel = create_geomodel_loopshowcase(P, include_fault = False) # Make geomodel

sys.path.append('../../MODFLOW_Tools')
from loop2flopy import Model
M = Model('run', P, plan = 'vor', transect = 'con') # Create flow model 
M.create_lith_dis_arrays(P)                                 # Create lith and dis arrays
M.create_prop_arrays(P)                                     # Create K arrays
M.create_flow_package_arrays(P)                             # Create flow packages

results = M.write_run_model(P, period = 'Steady', ic_array = P.strt)   
M.gwf, M.head_ss, M.obs_ss = results[0], results[1], results[2]
results = M.write_run_model(P, period = 'Past', ic_array = M.head_ss)   
M.gwf, M.head_past, M.obs_past = results[0], results[1], results[2]

from loop_showcase_functions import process_obs_past
hobs_past = process_obs_past(P, M)
np.savetxt('../modelfiles/heads.txt', hobs_past.flatten())

import pickle
pickle.dump(hobs_past, open('../results/pest_results.pkl','wb'))
    