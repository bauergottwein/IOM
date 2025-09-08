# -*- coding: utf-8 -*-
"""
Created on Thu May  8 09:03:30 2025

@author: vpk410
"""
import os
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
Run_attributes['WTP'] = dict()

Run_attributes['catchment_table'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Input data model\catchments.xlsx'
Run_attributes['catchment_shp'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Shapefiles\Catchment_Watersheds_Capital.shp' 

Run_attributes['wellfields_table'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Input data model\wellfields.xlsx'
Run_attributes['wellfields_shp'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Shapefiles\Wellfields.shp' 

Run_attributes['waterworks_table'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Input data model\waterworks.xlsx'
Run_attributes['waterworks_shp'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Shapefiles\Waterworks.shp' 

Run_attributes['wsa_table'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Input data model\water_supply_areas.xlsx'
Run_attributes['wsa_shp'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Shapefiles\WSA_capital.shp'

Run_attributes['wf_ww_table'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Input data model\WF_WW.xlsx'
Run_attributes['ww_wsa_table'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Input data model\WW_WSA.xlsx'
Run_attributes['ww_ww_table'] = r'c:\Users\vpk410\Documents\GW_allocation_model-main\Input data model\WW_WW.xlsx'

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
IOM_run.plot_catch_basemap()
IOM_run.plot_wsa_basemap()
#%%
IOM_run.wf_ww_basemap()
#%%
IOM_run.plot_dvar_ts(['Allocation_Ind_3168'])
IOM_run.plot_SP_ts(['SP_wb_WW_Storage_5780'])
#%%
IOM_run.plot_dvar_bar_shp_t_av('Allocation_HH_',IOM_run.nwsa , unit = '1000 m3/week', shapefile = IOM_run.wsa_shp)
#%%
IOM_run.plot_SP_bar_shp_t_av('SP_wb_WW_',IOM_run.nww[0:-1] , unit = 'DKK/m3', shapefile = IOM_run.waterworks_shp)
