# -*- coding: utf-8 -*-
"""
Created on Thu May  8 09:03:30 2025

@author: vpk410
"""
import os
import pandas as pd
my_path = r'c:\Users\vpk410\Documents\GW_allocation_model-main'

from IOM_v1 import IOM
#%%
Run_attributes = dict()
Run_attributes['run_name'] = r'Baseline_1'
Run_attributes['time_steps'] = 1721
Run_attributes['time_unit'] = 'week'
Run_attributes['volume_unit'] = '1000_m3'
Run_attributes['currency'] = '1000 DKK'
Run_attributes['water_users'] = ['HH','Ind','PS','Agri']
Run_attributes['gw_ind_2_in_use'] = False

Run_attributes['catchment_table'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Input data model\catchments.xlsx'
Run_attributes['catchment_shp'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Shapefiles\Catchment_Watersheds_Capital.shp' 

Run_attributes['wellfields_table'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Input data model\wellfields_active.xlsx'
Run_attributes['wellfields_shp'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Shapefiles\Wellfields_active.shp' 

Run_attributes['waterworks_table'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Input data model\waterworks_active.xlsx'
Run_attributes['waterworks_shp'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Shapefiles\Waterworks_active.shp' 

Run_attributes['wsa_table'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Input data model\water_supply_areas_2019.xlsx'
Run_attributes['wsa_shp'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Shapefiles\WSA_capital_2019.shp'

Run_attributes['wf_ww_table'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Input data model\WF_WW_v2_active.xlsx'
Run_attributes['ww_wsa_table'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Input data model\WW_WSA_2019_no_redundancy_active.xlsx'
Run_attributes['ww_ww_table'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Input data model\WW_WW_v2_active.xlsx'


Run_attributes['savepath'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Results'

Run_attributes['solvername'] = 'cplex'
Run_attributes['solverpath_exe'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\cplex.exe'
#%%
IOM_run = IOM(Run_attributes)
#%%
IOM_run = IOM_run.read_data()
#%%
IOM_run = IOM_run.process_data()
#%%
IOM_run = IOM_run.create_model()
#%%
IOM_run = IOM_run.define_model()
#%%
IOM_run = IOM_run.solve_model()
#%%
IOM_run = IOM_run.save_results()
#%%
%matplotlib qt
IOM_run.plot_catch_basemap()
IOM_run.plot_wsa_basemap()
#%%
IOM_run = IOM_run.wf_ww_basemap()
IOM_run = IOM_run.ww_wsa_basemap()
#%%
IOM_run.plot_dvar_ts(['Allocation_Ind_3168'])
IOM_run.plot_SP_ts(['SP_wb_WW_Storage_1'])
#%%
IOM_run.plot_dvar_bar_shp_t_av('Allocation_HH_',IOM_run.nwsa , unit = '1000 m3/week', shapefile = IOM_run.wsa_shp)
#%%
IOM_run.plot_SP_bar_shp_t_av('SP_wb_WW_',IOM_run.nww , unit = 'DKK/m3', shapefile = IOM_run.waterworks_shp)
#%%
#%%
IOM_run.plot_SP_bar_shp_t_av('SP_wb_WW_',IOM_run.nww , unit = 'DKK/m3', shapefile = IOM_run.waterworks_shp)
#%%
IOM_run.plot_SP_bar_shp_t_av('SP_pumping_WF_',IOM_run.nwf , unit = 'DKK/m3', shapefile = IOM_run.wellfields_shp)
#%%
IOM_run.plot_SP_bar_shp_t_av('SP_lin_res_',IOM_run.ncatch , unit = 'DKK/m3', shapefile = IOM_run.catchment_shp)
#%%
#%%
IOM_run.plot_SP_bar_shp_t_av('SP_min_bf_',IOM_run.ncatch , unit = 'DKK/m3', shapefile = IOM_run.catchment_shp)
#%%
IOM_run.plot_spatial_processed(IOM_run.Fulfilment_per_WSA['HH'] , column = 'HH',title = 'Demand fulfilment HH (%)', shapefile = IOM_run.wsa_shp, shpid = 'WSAID')
#%%
IOM_run.plot_spatial_processed(IOM_run.Fulfilment_per_WSA['Ind'] , column = 'Ind',title = 'Demand fulfilment Ind (%)', shapefile = IOM_run.wsa_shp, shpid = 'WSAID')
#%%
IOM_run.plot_spatial_processed(IOM_run.Fulfilment_per_WSA['PS'] , column = 'PS',title = 'Demand fulfilment PS (%)', shapefile = IOM_run.wsa_shp, shpid = 'WSAID')
#%%
base_str = 'SP_wb_WW_Exchange_'
merged_list = [base_str + str(item) for item in IOM_run.nww_ww]
sp = IOM_run.SPs[merged_list]
sp_mean = sp.mean(axis=0)
sp_mean.index=IOM_run.nww_ww
sp_mean_df = pd.DataFrame(sp_mean,columns = ['sp'])
sp_mean_df['WW1'] = [x[0] for x in sp_mean_df.index]
sp_mean_df['WW2'] = [x[1] for x in sp_mean_df.index]
sp_mean_df['unique_key'] = sp_mean_df['WW1'].astype('float').astype('string') + '_' +  sp_mean_df['WW2'].astype('float').astype('string')
sp_mean_df['sorted_key'] = sp_mean_df.index.map(lambda x: tuple(sorted(x)))
max_values = sp_mean_df.groupby('sorted_key').agg({
    'sp': 'max',  # This finds the maximum of the 'sp' column
    'unique_key': 'first'
    })
#%%    
IOM_run.plot_spatial_processed(max_values , column = 'sp',title = 'Transfer capacity shadow price, DKK/m3/week', 
                               shapefile = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Shapefiles\ww_ww_lines_v2.shp', shpid = 'unique_key',joincol = 'unique_key')