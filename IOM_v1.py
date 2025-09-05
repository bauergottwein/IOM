# -*- coding: utf-8 -*-
"""
@author: Peter Bauer-Gottwein
pbg@ign.ku.dk

"""
import pyomo.environ as pyo
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import os
import time
import geopandas as gpd

class IOM:
    def __init__(self, config_dict):
        print("IOM class constructor")
        for key, value in config_dict.items():
            setattr(self, key, value)
        #**********************************************************************
        # Column names for input tables. Change this if input format is modified
        #**********************************************************************
        # Column names of catchment table
        self.catchid_item = 'Catch_ID'
        self.catchment_area_item = 'Area (m2)'
        self.K_item = 'K (weeks)'
        self.Storage_LR_ini_item = 'Storage_LinRes_ini (1000 m3)'
        self.loss_fraction_item = 'Loss fraction (-)'
        self.demand_string = '_water_dem (1000m3)'
        self.Q_natural_string = 'Q_natural_'
        # Column names of wellfield table
        self.wfid_item = 'WFID'
        self.wf_maxpump_item = 'AnlgTillad (1000 m3/yr)'
        # Column names of waterworks table
        self.wwid_item = 'WWID'
        self.ww_storage_cap_item = 'Storage capacity (1000m3)'
        self.ww_storage_ini_item = 'Storage initial (1000m3)'
        # Column names of the water supply area table
        self.wsaid_item = 'WSAID'
        # Column names WF_WW table
        # Column names WW_WSA table
        # Column names WW_WW table
        self.WW1_item = 'WW1'
        self.WW2_item = 'WW2'
        self.WW_WW_cap = 'Capacity (1000m3/week)'
        # Column names catchment WB files
        self.recharge_item = 'MIKE SHE GW recharge (mm)' # Expects mm/time unit, csv files
        self.real_time_item = 'Start time of weekly time step'
        # Output files
        self.dec_outfile = self.savepath + os.sep + self.run_name + '_Optimal_Decisions.xlsx'
        self.sp_outfile =     self.savepath + os.sep + self.run_name + '_Shadow_Prices.xlsx'

    def read_data(self,**kwargs):
        
        self.ntimes = np.arange(1, self.time_steps+1, 1)   # weekly timesteps
        
        print(self.run_name, 'for', self.time_steps, self.time_unit)
        print(time.strftime("%H:%M:%S") + ' Importing data...')
               
        self.catchment_data = pd.read_csv(self.catchment_table,sep = ';')  
        self.WF = pd.read_csv(self.wellfields_table,sep = ';')
        self.WW=pd.read_csv(self.waterworks_table,sep = ';')
        self.WSA=pd.read_csv(self.wsa_table, sep = ';')
        self.WF_WW=pd.read_csv(self.wf_ww_table,sep = ';') 
        self.WW_WSA=pd.read_csv(self.ww_wsa_table,sep = ';')
        self.WW_WW=pd.read_csv(self.ww_ww_table,sep = ';')
        
        self.ncatch = np.array(self.catchment_data[self.catchid_item])
        
        self.inflow = np.empty((len(self.ncatch), len(self.ntimes)))

        for i in range(1,len(self.ncatch)+1):
            self.data=pd.read_csv(self.path_string_WB + str(int(i)) + '.csv')
            self.inflow[i-1,:] = self.data[self.recharge_item][:len(self.ntimes)]/1000 * self.catchment_data[self.catchment_area_item][i-1]/1000   # from mm to 10^3 m3
            
            #remove the negative data in the inflow ! 
            for t in range(1,len(self.ntimes)):
                if self.inflow[i-1,t-1]<0:
                    self.inflow[i-1,t] += self.inflow[i-1,t-1]
                    self.inflow[i-1,t-1] = 0
            if self.inflow[i-1,len(self.ntimes)-1]<0:
                self.inflow[i-1,len(self.ntimes)-1] = 0   
        
        self.real_time = self.data[self.real_time_item][0:self.time_steps] # real time is defined by inflow input from Mike SHE or similar
        return self

   
    def process_data(self,**kwargs):
     
        print(time.strftime("%H:%M:%S") + ' Processing the data...')
        
        # =============================================================================
        # Indeces (ntimes, nwsa, ncatch, nww, nwf)
        # =============================================================================

        self.nyear=np.arange(1, int(len(self.ntimes)//52.18)+2, 1) # yearly index 
        self.nwsa = np.array(self.WSA[self.wsaid_item])
        self.nww = np.array(self.WW[self.wwid_item])
        self.nwf = np.array(self.WF[self.wfid_item])
        self.nusers = self.water_users
        
        # =============================================================================
        # Double indeces for connections
        # =============================================================================
        
        self.nwf_ww = list(self.WF_WW.itertuples(index=False, name=None))
        self.nww_wsa = list(self.WW_WSA.itertuples(index=False, name=None))
        self.nww_ww = list(self.WW_WW.loc[:,[self.WW1_item,self.WW2_item]].itertuples(index=False, name=None))
        
        # =============================================================================
        # Set the WTP in weekly time steps
        # =============================================================================
        
        self.WTP_weekly = dict()

        for w in self.nwsa:
            for t in self.ntimes:
                for u in self.nusers:
                    self.WTP_weekly[u,w,t] = self.WTP[u]
            
        # =============================================================================
        # Linear reservoirs initial parameters
        # =============================================================================
            
        self.Kgwc = dict(zip(self.catchment_data[self.catchid_item], self.catchment_data[self.K_item]))  # K parameter weeks 
        self.Qbaseini = dict() 
        self.Sinigwc = dict()
        for c in self.ncatch:
                self.Qbaseini[c] = self.catchment_data[self.Storage_LR_ini_item][self.catchment_data[self.catchid_item] == c].iloc[0]/self.Kgwc[c]
                self.Sinigwc[c] = self.catchment_data[self.Storage_LR_ini_item][self.catchment_data[self.catchid_item] == c].iloc[0]
        
        # =============================================================================
        # Water Demand and WTP per WSA
        # =============================================================================     
        
        df_lookup = self.WSA.set_index(self.wsaid_item)             
        self.Demand = dict()
        
        for u in self.nusers:
            for w in self.nwsa:
                for t in self.ntimes:
                    self.Demand[u,w,t] = df_lookup.loc[w, u + self.demand_string] / 52.18

                      
        # =============================================================================
        # Inflow dictionary
        # =============================================================================
        
        self.I_inflow = dict()
        for c in self.ncatch:
            for t in self.ntimes:
                self.I_inflow[c,t] = self.inflow[c-1,t-1]
                
        # =============================================================================
        # Max pumping capacity
        # =============================================================================
        
        self.maxpump = dict()
        for wf in self.nwf:
            self.maxpump[wf] = self.WF.loc[self.WF[self.wfid_item] == wf, self.wf_maxpump_item].values[0]   # yearly maxpump 1000m3
        
        # =============================================================================
        # Storage capacity and initial storage of WaterWorks
        # =============================================================================
              
        self.maxstorage = dict()
        self.Storage_WW_ini = dict()
        for ww in self.nww:
            self.maxstorage[ww] = self.WW.loc[self.WW[self.wwid_item] == ww, self.ww_storage_cap_item].values[0]
            self.Storage_WW_ini[ww] = self.WW.loc[self.WW[self.wwid_item] == ww, self.ww_storage_ini_item].values[0]
            
        # =============================================================================
        # Water exchange dictionary 
        # =============================================================================
        
        self.ww_ww_cap = dict(zip(self.nww_ww,self.WW_WW[self.WW_WW_cap]))  # transfer capacity per 1000m3/day converted to 1000m3/week
                          
        # =============================================================================
        # Create the year / week dictionnary
        # =============================================================================
        
        self.week_in_year = dict()
        for t in self.ntimes:
            self.week_in_year[t]=int(t//52.18)+1 
        
        # =============================================================================
        # Create a month / week dictionnary
        # =============================================================================
        
        self.week_in_month = dict()
        for t in self.ntimes:
            self.week_in_month[t] = int((t % 52.18)//(52.18/12)+1)
        
        # =============================================================================
        # Loss fraction wastewater (what goes to the sea)           
        # =============================================================================
        
        self.loss_fraction=dict()
        for c in self.ncatch:
            self.loss_fraction[c] = self.catchment_data.loc[self.catchment_data[self.catchid_item] == c, self.loss_fraction_item].iloc[0]

        # =============================================================================
        # Natural flow data (Q_natural) for comparison with the baseflow from the model
        # =============================================================================
        
        def linres(Sini, deltat, K, Inflow):
            nper = len(Inflow)
            Qout = np.zeros(nper+1)
            Qout[0] = Sini/K
            for i in range(nper):
                Qout[i+1] = Qout[i]*np.exp(-deltat/K) + Inflow[i]*(1-np.exp(-deltat/K))
            Qout = Qout[1:]
            return Qout
        
        self.Q_natural = pd.DataFrame(index=self.ntimes)
        for c in self.ncatch:
            K = self.catchment_data.loc[self.catchment_data[self.catchid_item] == c, self.K_item].values[0]
            Sini = self.catchment_data.loc[self.catchment_data[self.catchid_item] == c, self.Storage_LR_ini_item].values[0]
            self.Q_natural[self.Q_natural_string + str(c)] = linres(Sini,1,K,self.inflow[c-1])  # timestep 1 week for K in weeks
        
        # minbf = dict(zip(ncatch, 0.75*Q_natural.median()))  # new minbf based on median natural baseflow
        self.minbf = dict(zip(self.ncatch, self.min_BF_fraction_of_natural * self.Q_natural.mean()))  # new minbf based on average natural baseflow
        
        # =============================================================================
        # Area dictionnary
        # =============================================================================
            
        self.area=dict()
        for c in self.ncatch:
            self.area[c] = self.catchment_data.loc[self.catchment_data[self.catchid_item]==c, self.catchment_area_item].values[0]
            
        
        return self
 
    def create_model(self,**kwargs):
        
        
        # =============================================================================
        # Create the model
        # =============================================================================
        print(time.strftime("%H:%M:%S") + ' Creating the model...')
        self.model = pyo.ConcreteModel() # define the model
        
        # =============================================================================
        # Define the index 
        # =============================================================================
        
        self.model.nusers = pyo.Set(dimen=1,initialize=self.nusers) # define user index, set values to nusers
        self.model.ntimes = pyo.Set(dimen=1,initialize=self.ntimes) # define time index, set values to ntimes
        self.model.nyear = pyo.Set(dimen=1,initialize=self.nyear)   # define year index
        self.model.nwsa = pyo.Set(dimen=1,initialize=self.nwsa)     # define WSA index, set values to nwsa
        self.model.nww = pyo.Set(dimen=1,initialize=self.nww)     # define WaterWorks index, set values to nww
        self.model.nwf = pyo.Set(dimen=1,initialize=self.nwf)     # define Wellfields index, set values to nwf
        self.model.ncatch = pyo.Set(dimen=1,initialize=self.ncatch) # define catchment index, set values to ncatch
        self.model.nwf_ww = pyo.Set(dimen=2,initialize=self.nwf_ww) # define index for all wf-ww connections
        self.model.nww_wsa = pyo.Set(dimen=2,initialize=self.nww_wsa) # define index for all ww-wsa connections
        self.model.nww_ww = pyo.Set(dimen=2,initialize=self.nww_ww) # define index for all ww-ww connections
        
        # =============================================================================
        # Declare decision variables - decision variable values will be provided by the optimizer
        # =============================================================================
        
        self.model.Allocations  = pyo.Var(self.model.nusers, self.model.nwsa, self.model.ntimes, within=pyo.NonNegativeReals) # Allocation to households, 1000 m3 per weekly time step
        
        self.model.Pump_WF = pyo.Var(self.model.nwf, self.model.ntimes, within=pyo.NonNegativeReals) # Sum of groundwater pumping for each wellfields 1000 m3 per weekly time step
        self.model.Pump_catch = pyo.Var(self.model.ncatch, self.model.ntimes, within=pyo.NonNegativeReals) # Sum of groundwater pumping for each catchment 1000 m3 per weekly time step
        self.model.Pump_GW_to_BF = pyo.Var(self.model.ncatch, self.model.ntimes, within=pyo.NonNegativeReals)  # Pumping to the river to maintain a min BF 
        self.model.Supply_WF_WW = pyo.Var(self.model.nwf_ww, self.model.ntimes, within=pyo.NonNegativeReals) # Supply from WF to WW 1000m3/week
        self.model.Storage_WW = pyo.Var(self.model.nww, self.model.ntimes, within=pyo.NonNegativeReals) # Water storage for each waterworks 1000m3 per weekly time step 
        self.model.Exchange_WW_WW  = pyo.Var(self.model.nww_ww, self.model.ntimes, within=pyo.NonNegativeReals) # Water transfer from 1 anlaeg to another, therefore 2 times nanlaeg, 1000m3 per weekly time step
        self.model.Supply_WW_WSA = pyo.Var(self.model.nww_wsa, self.model.ntimes, within=pyo.NonNegativeReals) # Water Supply distributed by each Waterworks to all the WSA it serves
        
        self.model.Q_base  = pyo.Var(self.model.ncatch, self.model.ntimes, within=pyo.NonNegativeReals) # Base flow from GW catchment, 1000 m3 per weekly time step
        self.model.Storage_LinRes   = pyo.Var(self.model.ncatch, self.model.ntimes, within=pyo.NonNegativeReals) # One end storage per month and per reservoir. 1000 m3 per weekly time step
        
        # =============================================================================
        # Declare parameters
        # =============================================================================
        
        #model.endtime = Param(initialize = ntimes[-1]) # find end time step of the model
        self.model.Demands  = pyo.Param(self.model.nusers, self.model.nwsa, self.model.ntimes,within=pyo.NonNegativeReals,initialize = self.Demand) # Set Water demands for all use categories
        
        self.model.WTP  = pyo.Param(self.model.nusers,self.model.nwsa, self.model.ntimes,within=pyo.NonNegativeReals,initialize = self.WTP_weekly) # Set Willingness To Pay for all use categories
        
        self.model.maxpump_WF = pyo.Param(self.model.nwf, within=pyo.NonNegativeReals,initialize = self.maxpump)  # Abstraction license 1000 m3/week
        self.model.Storage_WW_ini = pyo.Param(self.model.nww, within=pyo.NonNegativeReals,initialize = self.Storage_WW_ini) # Initial storage in Waterworks 1000m3
        self.model.maxstorage_WW = pyo.Param(self.model.nww, within=pyo.NonNegativeReals,initialize = self.maxstorage) # Max storage capacity 1000 m3
        self.model.maxexchange = pyo.Param(self.model.nww_ww, within=pyo.NonNegativeReals,initialize = self.ww_ww_cap)    # Water transfer capacity between waterworks in 1000m3/week
        
        self.model.pumping_cost = pyo.Param(within=pyo.NonNegativeReals, initialize = self.pumping_cost)   # pump cost in DKK/m3 or thousand DKK/1000m3
        self.model.exchange_cost = pyo.Param(within=pyo.NonNegativeReals, initialize = self.exchange_cost)  # water exchange cost DKK/m3 per distance ?????8
        self.model.loss_fraction_waste = pyo.Param(self.model.ncatch, within=pyo.NonNegativeReals, initialize = self.loss_fraction) #loss fraction of wastewater return flow to the river (what goes to the sea....)
        
        self.model.Storage_LinRes_ini = pyo.Param(self.model.ncatch, within=pyo.NonNegativeReals,initialize = self.Sinigwc) # Set initial GW storage for all groundwater catchments
        self.model.Kgwc = pyo.Param(self.model.ncatch, within=pyo.NonNegativeReals,initialize = self.Kgwc) # Set time constant for all groundwater catchments
        self.model.Qbase_ini = pyo.Param(self.model.ncatch, within=pyo.NonNegativeReals,initialize = self.Qbaseini) # Set initial BaseFlow for all catchments
        self.model.inflow = pyo.Param(self.model.ncatch, self.model.ntimes,within=pyo.Reals,initialize = self.I_inflow) # Set inflow to GW to for all catchments (from MIKE SHE model)
        self.model.minbflow = pyo.Param(self.model.ncatch,within=pyo.Reals,initialize = self.minbf) # Set environmental constraint on flow for all catchments
        
        return self
        
    def define_model(self,**kwargs):
        
        print(time.strftime("%H:%M:%S") + ' Defining the model...')
        
        # =============================================================================
        # Objective function
        # =============================================================================
        
        # Maximize the benefits and minimize the costs
        print(time.strftime("%H:%M:%S") + ' Defining objective function...')
        def obj_rule(model):
            User_ben = sum(model.WTP[u,w,t]*model.Allocations[u,w,t] for u in model.nusers for w in model.nwsa for t in model.ntimes)
            Pump_cost = sum(model.pumping_cost*model.Pump_WF[wf,t]  for wf in model.nwf for t in model.ntimes) + sum(model.pumping_cost*model.Pump_GW_to_BF[c,t] for c in model.ncatch for t in model.ntimes)
            Exchange_cost = sum(model.exchange_cost*model.Exchange_WW_WW[ww_ww,t] for ww_ww in model.nww_ww for t in model.ntimes)
            return User_ben - Pump_cost - Exchange_cost 
        self.model.obj = pyo.Objective(rule=obj_rule, sense = pyo.maximize)
        
        # =============================================================================
        # Allocation constraints
        # =============================================================================
        
        print(time.strftime("%H:%M:%S") + ' Defining allocation constraints...')
        # Household allocation does not exceed household demand. Active for every time step and catchment, thus two indices
        def wd_Users_c(model, u, w, t):
            return model.Allocations[u,w,t] <= model.Demands[u,w,t]
        self.model.wd_Users = pyo.Constraint(self.model.nusers, self.model.nwsa, self.model.ntimes, rule=wd_Users_c)            
        
        # =============================================================================
        # Pumping constraints 
        # =============================================================================
        print(time.strftime("%H:%M:%S") + ' Defining pumping constraints...')
        # Pump_catch variable = Total pumping of the WellFields in one catchment
        def pumping_catch_c(model,c,t):
            # get the list of WF in each catchment
            WF_in_catch = self.WF.drop_duplicates(subset='WFID')[self.WF.drop_duplicates(subset=self.wfid_item)[self.catchid_item] == c][self.wfid_item].values
            return  model.Pump_catch[c,t] == sum([model.Pump_WF[wf,t] for wf in WF_in_catch]) + model.Pump_GW_to_BF[c,t]
        self.model.pumping_catch = pyo.Constraint(self.model.ncatch, self.model.ntimes, rule=pumping_catch_c)
        
        # Total pumping for each Wellfields always below the abstraction license each week
        # def pumping_WF_c(model, wf, t):
        #     return model.Pump_WF[wf,t] <= model.maxpump_WF[wf]/52.18
        # model.pumping_WF = Constraint(model.nwf, model.ntimes, rule=pumping_WF_c)
        
        # Total pumping for each Wellfields always below the abstrction license for each year
        def pumping_WF_c(model, wf, y):
            list_weeks = np.array([week[0] for week in self.week_in_year.items() if week[1] == y])  #get the list of the weeks in a given year
            return sum(model.Pump_WF[wf,t] for t in list_weeks) <= model.maxpump_WF[wf]
        self.model.pumping_WF = pyo.Constraint(self.model.nwf, self.model.nyear, rule=pumping_WF_c)
        
        
        # =============================================================================
        # WaterBalance at the Wellfields, WaterWorks and WSA level
        # =============================================================================
        
        print(time.strftime("%H:%M:%S") + ' Defining wellfield water balance constraint...')
        
        # At the WF level
        def wb_WF_c(model,wf,t):
            # At the WF level
            indices_to_sum = [
                idx for idx in model.nwf_ww
                if idx[0] == wf
                ]
            return sum(
                model.Supply_WF_WW[idx, t]
                for idx in indices_to_sum
                ) == model.Pump_WF[wf,t]
        self.model.wb_WF = pyo.Constraint(self.model.nwf, self.model.ntimes, rule=wb_WF_c)
        
        # DeltaStorage = sum(Pumping)  +- Exchange  - Supply_WW_to_WSA
        
        print(time.strftime("%H:%M:%S") + ' Defining waterworks water balance constraint...')
        
        def wb_WW_c(model,ww,t):
            ww_indices_to_sum = [
                idx for idx in model.nwf_ww
                if idx[1] == ww
                ]
            wsa_indices_to_sum = [
                idx for idx in model.nww_wsa
                if idx[0] == ww
                ]
            ex1_indices_to_sum = [
                idx for idx in model.nww_ww
                if idx[0] == ww
                ]
            ex2_indices_to_sum = [
                idx for idx in model.nww_ww
                if idx[1] == ww
                ]
            if t == 1:
                return model.Storage_WW[ww,t] - model.Storage_WW_ini[ww] == sum(
                    model.Supply_WF_WW[idx, t]
                    for idx in ww_indices_to_sum
                ) - sum(
                    model.Supply_WW_WSA[idx, t]
                    for idx in wsa_indices_to_sum
                ) - sum(
                    model.Exchange_WW_WW[idx, t]
                    for idx in ex1_indices_to_sum
                ) + sum(
                    model.Exchange_WW_WW[idx, t]
                    for idx in ex2_indices_to_sum
                )
            else:
                return model.Storage_WW[ww,t] - model.Storage_WW[ww,t-1] == sum(
                    model.Supply_WF_WW[idx, t]
                    for idx in ww_indices_to_sum
                ) - sum(
                    model.Supply_WW_WSA[idx, t]
                    for idx in wsa_indices_to_sum
                ) - sum(
                    model.Exchange_WW_WW[idx, t]
                    for idx in ex1_indices_to_sum
                ) + sum(
                    model.Exchange_WW_WW[idx, t]
                    for idx in ex2_indices_to_sum
                )
          
        self.model.wb_WW = pyo.Constraint(self.model.nww, self.model.ntimes, rule=wb_WW_c)
        
        # Water allocation in each WSA, should be lower than the sum of the water allocated in the WSA
        
        print(time.strftime("%H:%M:%S") + ' Defining supply area water balance constraint...')
        
        def wb_WSA_c(model, wsa, t):
            # get the list of WW delivering water to each WSA
            # WW_in_WSA = WW[WW['WSAID'] == wsa]['WWID'].values # try both with summing for ww in WW_in_WSA or ww in model.nww and compare time
            wsa_indices_to_sum = [
                idx for idx in model.nww_wsa
                if idx[1] == wsa
            ]
            return sum(model.Allocations[u,wsa,t] for u in model.nusers) == sum(
                    model.Supply_WW_WSA[idx, t]
                    for idx in wsa_indices_to_sum
            )
        
        self.model.wb_WSA = pyo.Constraint(self.model.nwsa, self.model.ntimes, rule=wb_WSA_c)
        
        #For each WaterWorks, Storage WW <= maxstorage WW
        def wb_WW_Storage_c(model, ww, t):
            return model.Storage_WW[ww,t] <= model.maxstorage_WW[ww]
        self.model.wb_WW_Storage = pyo.Constraint(self.model.nww, self.model.ntimes, rule=wb_WW_Storage_c)
        
        # Max exchange capacity
        def wb_WW_Exchange_c(model, ww1, ww2,t):
            return model.Exchange_WW_WW[ww1, ww2, t] <= model.maxexchange[ww1,ww2]
        self.model.wb_WW_Exchange = pyo.Constraint(self.model.nww_ww, self.model.ntimes, rule=wb_WW_Exchange_c)
        
        # # Set the value of Supply_WW_WSA to zero when there is no connexion between WW and WSA
        # def wb_max_WW_WSA_c(model,ww,wsa,t):
        #     return model.Supply_WW_WSA[ww,wsa,t] <= model.WW_WSA[ww,wsa]*10000   # multiply matrix by a high enough number so it's like an infinite transfer capacity in the pipes
        # # model.wb_max_WW_WSA = Constraint(model.nww,model.nwsa, model.ntimes, rule=wb_max_WW_WSA_c)
        
        
        # =============================================================================
        # Linear reservoir constraints 
        # =============================================================================
        print(time.strftime("%H:%M:%S") + ' Defining linear reservoir constraints...')
        
        # Linear reservoirs base flow (from 1st order equation)
        def lin_res_c(model,c,t):
            WSA_in_catch = self.WSA[self.WSA[self.catchid_item] == c][self.wsaid_item].values  # list all WSA in the catchment (for wastewater)
            if t == 1:
                return model.Q_base[c,t] == model.Qbase_ini[c]*np.exp(-1/model.Kgwc[c]) + (model.inflow[c,t]-model.Pump_catch[c,t])*(1-np.exp(-1/model.Kgwc[c])) + (
                    1-model.loss_fraction_waste[c])*sum(model.Allocations[u,wsa,t] for u in model.nusers for wsa in WSA_in_catch) + model.Pump_GW_to_BF[c,t]
            else:
                return model.Q_base[c,t] == model.Q_base[c,t-1]*np.exp(-1/model.Kgwc[c]) + (model.inflow[c,t]-model.Pump_catch[c,t])*(1-np.exp(-1/model.Kgwc[c])) + (
                    1-model.loss_fraction_waste[c])*sum(model.Allocations[u,wsa,t] for u in model.nusers for wsa in WSA_in_catch) + model.Pump_GW_to_BF[c,t]
            
        self.model.lin_res = pyo.Constraint(self.model.ncatch, self.model.ntimes, rule=lin_res_c)  
        
        # Linear reservoirs storage (from 1st order equation), not needed, only for plotting
        def lin_res_stor_c(model,c,t):
            if t == 1:
                return model.Storage_LinRes[c,t] == model.Storage_LinRes_ini[c]*np.exp(-1/model.Kgwc[c]) + (model.inflow[c,t]-model.Pump_catch[c,t])*(1-np.exp(-1/model.Kgwc[c]))*model.Kgwc[c]
            else:
                return model.Storage_LinRes[c,t] == model.Storage_LinRes[c,t-1]*np.exp(-1/model.Kgwc[c]) + (model.inflow[c,t]-model.Pump_catch[c,t])*(1-np.exp(-1/model.Kgwc[c]))*model.Kgwc[c]
        
        self.model.lin_res_stor = pyo.Constraint(self.model.ncatch, self.model.ntimes, rule=lin_res_stor_c)
        
        
        # =============================================================================
        # Environmental constraints
        # =============================================================================
        print(time.strftime("%H:%M:%S") + ' Defining environmental flow constraints...')
        # min baseflow
        def min_bf_c(model,c,y):
            list_weeks = [week[0] for week in self.week_in_year.items() if week[1] == y]  #get the list of the weeks in a given year
            return sum(model.Q_base[c,t] for t in list_weeks)/len(list_weeks) >= model.minbflow[c]   # yearly average flow above minBF
        self.model.min_bf = pyo.Constraint(self.model.ncatch, self.model.nyear, rule=min_bf_c) 
        
        
        # Max pumping of aquifer: from "Model and Ensemble Indicator-Guided Assessment of Robust,
        #                              Exploitable Groundwater Resources for Denmark" table 1, indicator 2
        def gw_ind_2_c(model,c,y):
            list_weeks = [week[0] for week in self.week_in_year.items() if week[1] == y]  #get the list of the weeks in a given year
            return sum(model.Pump_catch[c,t] for t in list_weeks) <= sum(model.inflow[c,t] for t in list_weeks)/2

        if self.gw_ind_2_in_use:
            self.model.gw_ind_2 = pyo.Constraint(self.model.ncatch, self.model.nyear, rule=gw_ind_2_c)
        
        return self
    
    def solve_model(self,**kwargs):
        
        # =============================================================================
        # Dual problem
        # =============================================================================
        
        # formulate dual problem to provide shadow prices
        self.model.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT) 
                
        print(time.strftime("%H:%M:%S") + ' Solving the model...')
        
        # =============================================================================
        # Create a solver
        # =============================================================================
        
        opt = pyo.SolverFactory(self.solvername,executable=self.solverpath_exe)
        # opt = SolverFactory(solvername) # if it works without stating the path of the solver
        
        # =============================================================================
        # Solve
        # =============================================================================
        
        self.results = opt.solve(self.model)
        
        # =============================================================================
        # Check status
        # =============================================================================
        
        print(time.strftime("%H:%M:%S") + ' Model solved! \n')
        
        if (self.results.solver.status == pyo.SolverStatus.ok) and (self.results.solver.termination_condition == pyo.TerminationCondition.optimal):
            # Do something when the solution in optimal and feasible
            print('Optimal and Feasible \n')
        elif (self.results.solver.termination_condition == pyo.TerminationCondition.infeasible):
            # Do something when model in infeasible
            print('Infeasible')
        else:
            # Something else is wrong
            print ('Solver Status: ',  self.results.solver.status)
            print ('Solver termination condition: ', self.results.solver.termination_condition)
                    
        # =============================================================================
        # # Objective value
        # =============================================================================
        
        print("Total Benefit in optimal solution: ", round(pyo.value(self.model.obj)/len(self.model.ntimes)), " thousand DKK per week \n")
                
        # =============================================================================
        # Some results
        # =============================================================================
        
        self.Allocations_tot = dict()
        self.All_Allocations_tot = 0.
        for u in self.model.nusers:
            self.Allocations_tot[u] = sum(pyo.value(self.model.Allocations[u,wsa,t]) for wsa in self.model.nwsa for t in self.model.ntimes)
            self.All_Allocations_tot = self.All_Allocations_tot + self.Allocations_tot[u]
        
        self.Demands_tot = dict()
        self.All_Demands_tot = 0.
        for u in self.model.nusers:
            self.Demands_tot[u] = sum(pyo.value(self.model.Demands[u,wsa,t]) for wsa in self.model.nwsa for t in self.model.ntimes)
            self.All_Demands_tot  = self.All_Demands_tot + self.Demands_tot[u]
        

        self.Deficit_per_WSA = dict()
        for u in self.model.nusers:
            self.Deficit_per_WSA[u] = [100*sum(pyo.value(self.model.Allocations[u,wsa,t]) for t in self.model.ntimes)/(sum(pyo.value(self.model.Demands[u,wsa,t]) for t in self.model.ntimes) + 1E-6) for wsa in self.model.nwsa]
        
    
        for u in self.model.nusers:
            print('Allocation ' + str(u) + ' (% of demand) : ',round(100*self.Allocations_tot[u]/(self.Demands_tot[u]+1E-6),1), '%')
        return self
        
        
    def save_results(self,**kwargs):
        
        print(time.strftime("%H:%M:%S") + ' Saving the results...')
        
        # =============================================================================
        # Process decision variabes data
        # =============================================================================
        
        optimal_Allocations = np.zeros((len(self.model.nusers),len(self.model.ntimes),len(self.model.nwsa)))
        
        optimal_Pump_WF = np.zeros((len(self.model.ntimes),len(self.model.nwf)))
        optimal_Pump_catch = np.zeros((len(self.model.ntimes),len(self.model.ncatch)))
        optimal_Pump_GW_to_BF = np.zeros((len(self.model.ntimes),len(self.model.ncatch)))
        optimal_Storage_WW = np.zeros((len(self.model.ntimes),len(self.model.nww)))
        optimal_Q_base = np.zeros((len(self.model.ntimes),len(self.model.ncatch)))
        optimal_Storage_LinRes = np.zeros((len(self.model.ntimes),len(self.model.ncatch)))
        
        optimal_Exchange_WW_WW = np.zeros((len(self.model.ntimes),len(self.model.nww_ww)))
        optimal_Supply_WW_WSA = np.zeros((len(self.model.ntimes),len(self.model.nww_wsa)))
        optimal_Supply_WF_WW = np.zeros((len(self.model.ntimes),len(self.model.nwf_ww)))
        
        for t in self.model.ntimes:   
            for w in range(len(self.nwsa)):
                for u in range(len(self.nusers)):
                    optimal_Allocations[u,t-1,w] = pyo.value(self.model.Allocations[self.model.nusers.at(u+1),self.model.nwsa.at(w+1),t])

        
            for i in range(1,len(self.model.ncatch)+1):
                optimal_Pump_catch[t-1,i-1] = pyo.value(self.model.Pump_catch[self.model.ncatch.at(i),t])
                optimal_Pump_GW_to_BF[t-1,i-1] = pyo.value(self.model.Pump_GW_to_BF[self.model.ncatch.at(i),t])
                optimal_Q_base[t-1,i-1] = pyo.value(self.model.Q_base[self.model.ncatch.at(i),t])
                optimal_Storage_LinRes[t-1,i-1] = pyo.value(self.model.Storage_LinRes[self.model.ncatch.at(i),t])
                
            for i in range(1,len(self.model.nwf)+1):
                optimal_Pump_WF[t-1,i-1] = pyo.value(self.model.Pump_WF[self.model.nwf.at(i),t])
                
            for i in range(1, len(self.model.nww)+1):
                optimal_Storage_WW[t-1,i-1] = pyo.value(self.model.Storage_WW[self.model.nww.at(i),t])
            
            for i in range(1, len(self.model.nww_ww)+1):    
                optimal_Exchange_WW_WW[t-1,i-1] = pyo.value(self.model.Exchange_WW_WW[self.model.nww_ww.at(i),t])
                    
            for i in range(1, len(self.model.nww_wsa)+1):
                optimal_Supply_WW_WSA[t-1,i-1] = pyo.value(self.model.Supply_WW_WSA[self.model.nww_wsa.at(i),t])
                
            for i in range(1, len(self.model.nwf_ww)+1):
                optimal_Supply_WF_WW[t-1,i-1] = pyo.value(self.model.Supply_WF_WW[self.model.nwf_ww.at(i),t])
        
        # =============================================================================
        # Convert to DataFrame and save
        # =============================================================================
        optimal_Allocation_dfs = dict()
        for u in range(len(self.nusers)):
            optimal_Allocation_dfs[u] = pd.DataFrame(optimal_Allocations[u,:,:], index=self.model.ntimes, columns=['Allocation_' + self.nusers[u] + '_' + str(w) for w in self.model.nwsa])
        list_of_alloc_dfs = []
        for df in optimal_Allocation_dfs.values():
            list_of_alloc_dfs.append(df)
        
        optimal_Storage_WW = pd.DataFrame(optimal_Storage_WW, index=self.model.ntimes, columns=['Storage_WW_'+str(w) for w in self.model.nww])
        optimal_Pump_WF = pd.DataFrame(optimal_Pump_WF, index=self.model.ntimes, columns=['Pump_WF_'+str(w) for w in self.model.nwf])
        optimal_Pump_catch = pd.DataFrame(optimal_Pump_catch, index=self.model.ntimes, columns=['Pump_catch_'+str(c) for c in self.model.ncatch])
        optimal_Pump_GW_to_BF = pd.DataFrame(optimal_Pump_GW_to_BF, index=self.model.ntimes, columns=['Pump_GW_to_BF_'+str(c) for c in self.model.ncatch])
        optimal_Q_base = pd.DataFrame(optimal_Q_base, index=self.model.ntimes, columns=['Q_base_'+str(c) for c in self.model.ncatch])
        optimal_Storage_LinRes = pd.DataFrame(optimal_Storage_LinRes, index=self.model.ntimes, columns=['Send_'+str(c) for c in self.model.ncatch])
        optimal_Exchange_WW_WW = pd.DataFrame(optimal_Exchange_WW_WW, index=self.model.ntimes, columns=['Exchange_'+str(c) for c in self.model.nww_ww])
        optimal_Supply_WW_WSA = pd.DataFrame(optimal_Supply_WW_WSA, index=self.model.ntimes, columns=['Supply_'+str(c) for c in self.model.nww_wsa])
        optimal_Supply_WF_WW = pd.DataFrame(optimal_Supply_WF_WW, index=self.model.ntimes, columns=['Supply_'+str(c) for c in self.model.nwf_ww])
        
        optimal_time = pd.DataFrame({'time':[t for t in self.model.ntimes]}, index=self.model.ntimes)
        
        self.optimal_Decisions = pd.concat([optimal_time, optimal_Storage_WW, optimal_Pump_WF, optimal_Pump_catch, optimal_Pump_GW_to_BF, optimal_Q_base, optimal_Storage_LinRes, optimal_Exchange_WW_WW, optimal_Supply_WW_WSA,optimal_Supply_WF_WW]+list_of_alloc_dfs, axis=1) 
        self.optimal_Decisions.to_excel(self.dec_outfile,sheet_name = 'Decision variables')
        

        
        # =============================================================================
        # Process Shadow prices data
        # =============================================================================
                
        SP_wd_Allocations = np.zeros((len(self.model.nusers),len(self.model.ntimes),len(self.model.nwsa)))
        
        SP_pumping_WF = np.zeros((len(self.model.nyear),len(self.model.nwf)))
        SP_wb_WW_Storage = np.zeros((len(self.model.ntimes),len(self.model.nww)))
        SP_wb_WW_Exchange = np.zeros((len(self.model.ntimes),len(self.model.nww_ww)))
        
        # and Supply_WW_WSA ???
        
        SP_lin_res = np.zeros((len(self.model.ntimes),len(self.model.ncatch)))
        SP_min_bf =np.zeros((len(self.model.nyear),len(self.model.ncatch)))
        SP_gw_ind_2 = np.zeros((len(self.model.nyear),len(self.model.ncatch)))
        
        for t in self.model.ntimes:
            
            for u in range(1,len(self.model.nusers)+1):
                for w in range(1,len(self.model.nwsa)+1):      
                    SP_wd_Allocations [u-1, t-1, w-1] = self.model.dual[self.model.wd_Users[self.model.nusers.at(u),self.model.nwsa.at(w),t]]

                
            for i in range(1,len(self.model.ncatch)+1):
                SP_lin_res [t-1,i-1] = self.model.dual[self.model.lin_res[self.model.ncatch.at(i),t]]
                
                
            for i in range(1, len(self.model.nww)+1):
                SP_wb_WW_Storage[t-1,i-1] = self.model.dual[self.model.wb_WW_Storage[self.model.nww.at(i),t]]
                
            for i in range(1, len(self.model.nww_ww)+1):
                SP_wb_WW_Exchange[t-1,i-1] = self.model.dual[self.model.wb_WW_Exchange[self.model.nww_ww.at(i),t]]
        
        for y in self.model.nyear:
            for i in range(1,len(self.model.nwf)+1):
                SP_pumping_WF[y-1,i-1] = self.model.dual[self.model.pumping_WF[self.model.nwf.at(i),y]]
        
            for i in range(1,len(self.model.ncatch)+1):
                SP_min_bf [y-1,i-1] = self.model.dual[self.model.min_bf[self.model.ncatch.at(i),y]]
                if self.gw_ind_2_in_use:
                    SP_gw_ind_2 [y-1,i-1] = self.model.dual[self.model.gw_ind_2[self.model.ncatch.at(i),y]]
                
        # =============================================================================
        # Convert SP to dataframe and save
        # =============================================================================
        
        SP_wd_Allocations_dfs = dict()
        for u in range(len(self.nusers)):
            SP_wd_Allocations_dfs[u] = pd.DataFrame(SP_wd_Allocations[u,:,:], index=self.model.ntimes, columns=['SP_' + self.nusers[u] + '_' + str(w) for w in self.model.nwsa])
        list_of_SP_dfs = []
        for df in SP_wd_Allocations_dfs.values():
            list_of_SP_dfs.append(df)

        SP_wb_WW_Storage = pd.DataFrame(SP_wb_WW_Storage, index=self.model.ntimes, columns=['SP_wb_WW_Storage_'+str(w) for w in self.model.nww])
        SP_wb_WW_Exchange = pd.DataFrame(SP_wb_WW_Storage, index=self.model.ntimes, columns=['SP_wb_WW_Exchange_'+str(w) for w in self.model.nww_ww])
        SP_pumping_WF = pd.DataFrame(SP_pumping_WF, index=self.model.nyear, columns=['SP_pumping_WF_'+str(w) for w in self.model.nwf])
        
        SP_lin_res = pd.DataFrame(SP_lin_res, index=self.model.ntimes, columns=['SP_lin_res_'+str(c) for c in self.model.ncatch])
        SP_min_bf = pd.DataFrame(SP_min_bf, index=self.model.nyear, columns=['SP_min_bf_'+str(c) for c in self.model.ncatch])
        
        if self.gw_ind_2_in_use:
            SP_gw_ind_2 = pd.DataFrame(SP_gw_ind_2, index=self.model.nyear, columns=['SP_gw_ind_2_'+str(c) for c in self.model.ncatch])
        
        SP_time = pd.DataFrame({'time':[t for t in self.model.ntimes]}, index=self.model.ntimes)
        if self.gw_ind_2_in_use:
            self.SPs = pd.concat([SP_time, SP_wb_WW_Storage, SP_wb_WW_Exchange, SP_pumping_WF, SP_lin_res, SP_min_bf, SP_gw_ind_2] + list_of_SP_dfs, axis=1) 
        else:
            self.SPs = pd.concat([SP_time, SP_wb_WW_Storage, SP_wb_WW_Exchange, SP_pumping_WF, SP_lin_res, SP_min_bf] + list_of_SP_dfs, axis=1) 
        self.SPs.to_excel(self.sp_outfile,sheet_name = 'Shadow prices')
        return self
    
    
    def plot_catch_basemap(self,**kwargs):
        gdf_catch = gpd.read_file(self.catchment_shp)
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        gdf_catch.plot(
            facecolor='lightblue',
            edgecolor='black',
            linewidth=0.5,  # Control the thickness of the border
            figsize=(10, 6),
            ax=ax)
        ax.tick_params(axis='both', labelsize=8)
        ax.xaxis.get_offset_text().set_fontsize(8)
        
        for idx, row in gdf_catch.iterrows():
            # Get the coordinates for the label from the centroid
            x, y = row.geometry.centroid.x, row.geometry.centroid.y
    
            # Annotate the plot with the country name
            ax.annotate(
                text=row[self.catchid_item],         # The text to display
                xy=(x, y),                        # The position of the text
                xytext=(3, 3),            # Offset the text slightly
                textcoords="offset points",
                ha='center',              # Horizontal alignment
                fontsize=6,               # Font size
                color='darkblue'
            )
        
        plt.show()
        
    def plot_wsa_basemap(self,**kwargs):
        gdf_wsa = gpd.read_file(self.wsa_shp)
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        gdf_wsa.plot(
            facecolor='lightcoral',
            edgecolor='black',
            linewidth=0.5,  # Control the thickness of the border
            figsize=(10, 6),
            ax=ax)
        ax.tick_params(axis='both', labelsize=8)
        ax.xaxis.get_offset_text().set_fontsize(8)
        ax.set_xlim(675000,735000)
        ax.set_ylim(6138000,6230000)
        
        for idx, row in gdf_wsa.iterrows():
            # Get the coordinates for the label from the centroid
            x, y = row.geometry.centroid.x, row.geometry.centroid.y
    
            # Annotate the plot with the country name
            ax.annotate(
                text=row[self.wsaid_item],         # The text to display
                xy=(x, y),                        # The position of the text
                xytext=(3, 3),            # Offset the text slightly
                textcoords="offset points",
                ha='center',              # Horizontal alignment
                fontsize=6,               # Font size
                color='darkred'
            )
        
        plt.show()
           
    def plot_dvar_ts(self,dvar_keys,**kwargs):
        # dvar_keys: List of decision variables that should be plotted together
        fig, ax = plt.subplots()
        for d in dvar_keys:
            ax.plot(self.real_time, self.optimal_Decisions[d], label=d)
        ax.set_ylabel(self.volume_unit + ' / ' + self.time_unit)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
        ax.tick_params(axis='x', labelrotation=45)
        ax.legend()
        plt.show()
        
    def plot_SP_ts(self,SP_keys,**kwargs):
        # SP_keys: List of shadow prices that should be plotted together
        fig, ax = plt.subplots()
        for d in SP_keys:
            ax.plot(self.real_time, self.SPs[d], label=d)
        ax.set_ylabel(self.currency + ' / ' + self.volume_unit)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
        ax.tick_params(axis='x', labelrotation=45)
        ax.legend()
        plt.show()
    
