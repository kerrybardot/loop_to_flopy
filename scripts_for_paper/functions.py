#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Functions for LOOP showcase

def prepare_geomodel_loopshowcase(P, data, strat, include_fault):
    import pandas as pd
    import numpy as np
    import numbers
    df = data.copy()
    units = list(df.columns.values[4:])   # Make a list of formations  
    df.easting = pd.to_numeric(df.easting)    # Make sure Eastings and Northings are numeric values
    df.northing = pd.to_numeric(df.northing)
    df.ground = pd.to_numeric(df.ground)


    # ---------- Prepare borehole data ----------------
    data_list = df.values.tolist()  # Turn data into a list of lists
    formatted_data = []
    for i in range(len(data_list)): #iterate for each row

        #val = strat.val[0] # Add data for groundlevel
        #gx, gy, gz = 0,0,1 # flat
        #formatted_data.append([boreid, easting, northing, groundlevel, val, 'ground', 'ground', gx, gy, gz]) 
        
        boreid = data_list[i][0]
        easting, northing = data_list[i][2], data_list[i][3]
        groundlevel = data_list[i][4]  
        count = 0  # Add data row for each lithology
        for j in range(5,df.shape[1]-1): #iterate through each formation 
            if isinstance(data_list[i][j], numbers.Number) == True:  # Add lithology  
                bottom    = groundlevel + float(data_list[i][j])  # Ground surface - formation bottom (mbgl)
                val       = strat.val[count]                   # designated isovalue
                unit      = strat.unit[count]                  # unit 
                feature   = strat.sequence[count]              # sequence
                gx, gy, gz = 0,0,1                             # normal vector to surface (flat)
                formatted_data.append([boreid, 'raw', easting, northing, bottom, val, unit, feature, gx, gy, gz])    
                current_bottom = np.copy(bottom)            
            count+=1

    # ---------- Add control points ----------------   
   
    for cp in P.control_points:
            formatted_data.append(cp)
    
    data = pd.DataFrame(formatted_data)
    data.columns =['ID','data_type','X','Y','Z','val','lithcode','feature_name', 'gx', 'gy', 'gz']
    
    # ---------- Prepare fault details ----------------   
    if include_fault:
        # Make cloud of points along fault plane
        nh = 10 # points  in x/y plane
        x_array, y_array = [], [] # arrays to create points along fault
        x_array.append(P.fx1)
        y_array.append(P.fy1)
        for i in range(nh-2):
            x_array.append(P.fx1 + (i+1) * (P.fx2-P.fx1)/(nh-1))
            y_array.append(P.fy1 + (i+1) * (P.fy2-P.fy1)/(nh-1))
        x_array.append(P.fx2)
        y_array.append(P.fy2)

        z_array = np.arange(-200, -100, 20) 
        nv = len(z_array) # points in z plane   

        fault_azimuth = np.degrees(np.arcsin((P.fx2-P.fx1)/(P.fy2-P.fy1)))  
        strike, dip = fault_azimuth-180, 90                                                                           

        from LoopStructural.utils import strikedip2vector # array of arrays
        nx, ny, nz = strikedip2vector([strike], [dip])[0]
        #from LoopStructural.utils.helper import strike_dip_vector
        #nx, ny, nz = strike_dip_vector([strike], [dip])[0]

        fault_plane_3d = []
        for v in range(nv):# vertical points 
            for h in range(nh): # horizontal points
                x, y, z = x_array[h], y_array[h], z_array[v]
                fault_plane_3d.append((x,y,z))
                df_new_row = pd.DataFrame.from_records({'ID':['fault'],'data_type':['fault_surface'],
                                                        'X':[x], 'Y':[y], 'Z':[z], 'val':[0.], 
                                                        'feature_name':['Fault'], 'nx': [nx], 'ny': [ny], 'nz': [nz]})
                data = pd.concat([data, df_new_row], ignore_index = True)
            
    return(data, strat)


def create_geomodel_loopshowcase(P, include_fault):
    import numpy as np
    origin  = (P.x0, P.y0, P.z0)
    maximum = (P.x1, P.y1, P.z1)
    
    from LoopStructural import GeologicalModel
    geomodel = GeologicalModel(origin, maximum)
    geomodel.set_model_data(P.data)  
    
    Upper      = geomodel.create_and_add_foliation("upper")
    UpperUC    = geomodel.add_unconformity(Upper, -50) 


    if include_fault:
        print('Fault included!')
        faultfunction = create_faultfunction()
        Fault = geomodel.create_and_add_fault('Fault', 
                                              displacement      = P.fault_max_disp,
                                              fault_center      = P.fault_center,
                                              minor_axis        = P.minor_axis, # fault_influence
                                              #fault_slip_vector = P.fault_slip_vector,
                                              #major_axis        = P.major_axis, # fault_extent
                                              #intermediate_axis = P.intermediate_axis, # fault_vertical_radius
                                              #fault_dip_anisotropy=0.0,
                                              #fault_trace_anisotropy=0.0,
                                              #faultfunction     = faultfunction, #'BaseFault', 
                                              #nelements=4000, 
                                              #steps=4, 
                                              #interpolatortype="FDI", 
                                              #buffer=0.3, 
                                              #solver='cg',
                                              )    
    Lower      = geomodel.create_and_add_foliation("lower") 
    #print(geomodel.feature_name_index)
    
    # Add Strat Column
    stratigraphic_column = {}
    stratigraphic_column["upper"] = {}
    stratigraphic_column["upper"]["a"] = {"min": -50, "max": np.inf, "id": 0}
    stratigraphic_column["lower"] = {}
    stratigraphic_column["lower"]["b"] = {"min": -100, "max": np.inf, "id": 1}
    stratigraphic_column["lower"]["c"] = {"min": -200, "max": -100, "id": 2}
    stratigraphic_column["lower"]["d"] = {"min": -np.inf, "max": -200, "id": 3}


    geomodel.set_stratigraphic_column(stratigraphic_column)
       
    return(geomodel)

def process_obs_steady(P, M):
    import numpy as np
    import os
    import pandas as pd
    M.hobs_steady = np.zeros((P.nobs, P.nzobs, 1), dtype = float)
    fname = str(M.modelname + "_steady.csv")
    csv_file = os.path.join(P.workspace, fname)
    data_set = pd.read_csv(csv_file)#, header=1)
    df = pd.DataFrame(data_set)
    a = df.to_numpy()
    hobs = a[0][1:(P.nobs*P.nzobs+1)]
    for ob_bore in range(P.nobs):
        for z in range(P.nzobs):
            n = ob_bore*P.nzobs + z
            M.hobs_steady[ob_bore, z] = hobs[n]
    return(M.hobs_steady)

def process_obs_past(P, M):
    import numpy as np
    import os
    import pandas as pd
    M.hobs_past = np.zeros((P.nobs, P.nzobs, P.nts_past), dtype = float) #P.nts_past-1)
    fname = str(M.modelname + "_past.csv")
    csv_file = os.path.join(P.workspace, fname)
    data_set = pd.read_csv(csv_file, header=0)
    data_frames = pd.DataFrame(data_set)
    hobs = np.array(data_frames.values)
    hobs = hobs[:,1:]
    hobs = np.swapaxes(hobs, 0, 1)
    for ob_bore in range(P.nobs):
        for z in range(P.nzobs):
            n = ob_bore*P.nzobs + z
            M.hobs_past[ob_bore, z, :] = hobs[n,:]
    return(M.hobs_past)

def process_obs_future(P, M):
    import numpy as np
    import os
    import pandas as pd
    M.hobs_future = np.zeros((P.nobs, P.nzobs, P.nts_future), dtype = float) # P.nts_future-1)
    fname = str(M.modelname + "_future.csv")
    csv_file = os.path.join(P.workspace, fname)
    data_set = pd.read_csv(csv_file, header=0)
    data_frames = pd.DataFrame(data_set)
    hobs = np.array(data_frames.values)
    hobs = hobs[:,1:]
    
    hobs = np.swapaxes(hobs, 0, 1)
    for ob_bore in range(P.nobs):
        for z in range(P.nzobs):
            n = ob_bore*P.nzobs + z
            M.hobs_future[ob_bore, z, :] = hobs[n,:]
    return(M.hobs_future)

# Plot observations
def plot_observations(heads, modelnames, ylim = None):
    import matplotlib.pyplot as plt
    colors = ['gray', 'red', 'green', 'purple', 'orange', 'blue']
    n = len(heads)
    fig = plt.figure(figsize = (12,3))
    fig.suptitle('Steady state heads')

    for ob in range(P.nobs): # New figure for each obs bore (OB1, OB2, OB3, OB4)
        ax = plt.subplot(1,5,ob+1)#, aspect = 'equal')
        ax.set_title(P.idobsbores[ob], size = 10) 
        for i in range(n): # for each option
            hobs = heads[i]   # extract obs data
            ax.plot(P.zobs[ob], hobs[ob,:,0], '-o', markersize = 4, c = colors[i], label = modelnames[i])

        ax.set_xlabel('Obs Depth (m)')
        if ylim != None:
            ax.set_ylim(ylim)
        if ob > 0: ax.set_yticks([])
        if ob == 0: ax.set_ylabel('Head (m)')
        #ax.set_xticks([0,1,2,3,4])
        #ax.set_xticklabels(P.zobs)
        #ax.set_xlim([10,30])
        #ax.set_yticks([30, 40, 50, 60])

    from matplotlib.lines import Line2D
    legend_markers = []
    for i in range(len(options)):
        legend_markers.append((Line2D([0], [0], marker='o', markersize = 4, color=colors[i])))
    ax = plt.subplot(1,5,5, aspect = 'equal')
    ax.set_axis_off()
    ax.legend(legend_markers, modelnames, loc="center", fontsize = 9)#, ncols = 3, bbox_to_anchor=[0.5, 0.7])
    plt.tight_layout()
    plt.show()


# In[ ]:


def find_watertable_disu(P, M, layer):
    model = M.gwf
    water_table = flopy.utils.postprocessing.get_water_table(M.gwf.output.head().get_data())
    M.heads_disv = -1e30 * np.ones_like(M.idomain, dtype=float) 
    for i, h in enumerate(water_table):
        if math.isnan(h) == False: 
            M.heads_disv[M.cellid_disu==i] = h        
    return(M.heads_disv[layer])

def plot_head_diff(P, M, heads1, heads2, vmin = None, vmax = None): 
    fig = plt.figure(figsize=(8, 8))
    ax = plt.subplot(111, aspect="equal")
    #M = pinchout_models[0]      
    pmv = flopy.plot.PlotMapView(modelgrid=M.vgrid)
    H = pmv.plot_array(heads1 - heads2, cmap = 'Spectral', alpha = 0.6)#vmin = vmin, vmax = vmax, 

    for j in range(len(P.xyobsbores)):
        ax.plot(P.xyobsbores[j][0], P.xyobsbores[j][1],'o', ms = '4', c = 'black')
        ax.annotate(P.idobsbores[j], (P.xyobsbores[j][0], P.xyobsbores[j][1]+100), c='black', size = 12) #, weight = 'bold')

    for j in range(len(P.xypumpbores)):
        ax.plot(P.xypumpbores[j][0], P.xypumpbores[j][1],'o', ms = '4', c = 'red')
        ax.annotate(P.idpumpbores[j], (P.xypumpbores[j][0], P.xypumpbores[j][1]+100), c='red', size = 12) #, weight = 'bold')

    if M.plan == 'car': P.sg.plot(ax=ax, edgecolor='black', lw = 0.2)
    if M.plan == 'tri': P.tri.plot(ax=ax, edgecolor='black', lw = 0.2)
    if M.plan == 'vor': P.vor.plot(ax=ax, edgecolor='black', lw = 0.2)
    ax.set_title('Head diffence between two most extreme structural models', size = 10)
    plt.colorbar(H, shrink = 0.4)


# In[1]:


# Plot observations
def param_vs_struct(param_obs_heads, pinchout_obs_heads, xlim = None, ylim = None):

    fig = plt.figure(figsize = (12,3))
    fig.suptitle('Steady state heads')

    for ob in range(P.nobs): # New figure for each obs bore (OB1, OB2, OB3, OB4)
        a = np.array(param_obs_heads)[:,ob,:,0] #16,4,5,1
        b = np.array(pinchout_obs_heads)[:,ob,:,0] #6,4,5,1
        a_max = np.nanpercentile(a, 0, axis=0)
        a_min = np.nanpercentile(a, 100, axis=0)
        b_max = np.nanpercentile(b, 0, axis=0)
        b_min = np.nanpercentile(b, 100, axis=0)
    
        ax = plt.subplot(1,5,ob+1)#, aspect = 'equal')
        ax.set_title(P.idobsbores[ob], size = 10) 
        for i in range(16): # for each option
            h = param_obs_heads[i]   # extract obs data
            ax.plot(P.zobs, h[ob,:,0], '-', lw = 0.2, alpha = 0.7, c = 'blue')
        for i in range(6): # for each option
            h = pinchout_obs_heads[i]   # extract obs data
            ax.plot(P.zobs, h[ob,:,0], '-', lw = 0.2, alpha = 0.7, c = 'red')
        ax.fill_between(P.zobs, a_min, a_max, color = 'blue', alpha = 0.1)
        ax.fill_between(P.zobs, b_min, b_max, color = 'red', alpha = 0.1)
        ax.plot(P.zobs, a_min, '-', lw = 2, c = 'blue')
        ax.plot(P.zobs, a_max, '-', lw = 2, c = 'blue')
        ax.plot(P.zobs, b_min, '-', lw = 2, c = 'red')
        ax.plot(P.zobs, b_max, '-', lw = 2, c = 'red')
        
        ax.set_xlabel('Obs Depth (m)')
        if xlim != None: ax.set_xlim(xlim)
        if ylim != None: ax.set_ylim(ylim)
        if ob > 0: ax.set_yticks([])
        if ob == 0: ax.set_ylabel('Head (m)')
        #ax.set_xticks([0,1,2,3,4])
        #ax.set_xticklabels(P.zobs)
        #ax.set_xlim([10,30])
        #ax.set_yticks([30, 40, 50, 60])

    from matplotlib.lines import Line2D
    legend_markers = []
    legend_markers.append((Line2D([0], [0], marker='o', markersize = 4, color='blue')))
    legend_markers.append((Line2D([0], [0], marker='o', markersize = 4, color='red')))
    ax = plt.subplot(1,5,5, aspect = 'equal')
    ax.set_axis_off()
    ax.legend(legend_markers, ['parameter variation', 'structural variation'], loc="center", fontsize = 9)#, ncols = 3, bbox_to_anchor=[0.5, 0.7])
    plt.tight_layout()
    plt.show()


# In[ ]:





# In[ ]:




