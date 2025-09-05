# -*- coding: utf-8 -*-
"""
@author: Clément Franey
clefraney@gmail.com

modified by Peter Bauer-Gottwein

Updated: June 11 2025


Instructions:
    
Download the folder 'GW_all' anywhere in your computer.

Set the path of the folder where the python model is located in the variable "path_Model_folder".
Set the path of the folder where the solver is located in the variable "solverpath_exe".
Set the name of the sover in the variable "solvername".

You can modify the global parameters as you wish.

"""



import pyomo.environ as pyo
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import time

class IOM:
    def __init__(self, config_dict):
        print("IOM class constructor")
        for key, value in config_dict.items():
            setattr(self, key, value)


    def read_data(self,**kwargs):
        
        self.ntimes = np.arange(1, self.time_steps+1, 1)   # weekly timesteps
        
        print(self.run_name, 'for', self.time_steps, self.time_unit)
        print(time.strftime("%H:%M:%S") + ' Importing data...')
               
        self.catchment_data = pd.read_csv(self.catchment_table,sep = ';')
        self.WTP = pd.read_csv(self.wtp_table,sep = ';')     
        self.WF = pd.read_csv(self.wellfields_table,sep = ';')
        self.WW=pd.read_csv(self.waterworks_table,sep = ';')
        self.WSA=pd.read_csv(self.wsa_table, sep = ';')
        self.WF_WW=pd.read_csv(self.wf_ww_table,sep = ';') 
        self.WW_WSA=pd.read_csv(self.ww_wsa_table,sep = ';')
        self.WW_WW=pd.read_csv(self.ww_ww_table,sep = ';')
        
        self.ncatch = np.array(self.catchment_data[self.catchid_item])
        self.K_optim=pd.read_csv(self.K_optim_table)
        
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
        self.nww_ww = list(self.WW_WW.drop('Capacity (1000m3/day)',axis=1).itertuples(index=False, name=None))
        
        # =============================================================================
        # Convert the WTP in weekly time steps
        # =============================================================================
        
        self.WTP_weekly=pd.DataFrame(columns=['Time step','WTP_HH','WTP_Ind','WTP_PS','WTP_Agri'])
        self.WTP_weekly['Time step']=self.data['Start time of weekly time step']
        for i in range(0,len(self.WTP_weekly)):
            self.WTP_weekly.loc[i, 'WTP_' + self.water_users[0]] = float(self.WTP['WTP Households (DKK/m3)'].iloc[-1])
            self.WTP_weekly.loc[i, 'WTP_Ind'] = self.WTP_weekly.loc[i, 'WTP_HH'] * 1/3
            self.WTP_weekly.loc[i, 'WTP_PS'] = self.WTP_weekly.loc[i, 'WTP_HH'] * 1/2
            self.WTP_weekly.loc[i, 'WTP_Agri'] = self.WTP_weekly.loc[i, 'WTP_HH'] * 1/4
        
        self.WTP_HH = dict()
        self.WTP_Ind = dict()
        self.WTP_PS = dict()
        self.WTP_Agri = dict()

        for w in self.nwsa:
            for t in self.ntimes:
                self.WTP_HH[w,t] = self.WTP_weekly['WTP_HH'][t]
                self.WTP_Ind[w,t] = self.WTP_weekly['WTP_Ind'][t]
                self.WTP_PS[w,t] = self.WTP_weekly['WTP_PS'][t]
                self.WTP_Agri[w,t] = self.WTP_weekly['WTP_Agri'][t]
            
        # =============================================================================
        # Linear reservoirs initial parameters
        # =============================================================================
            
        self.K_optim['Sinigwc'] = 1000   # 1000m3
        self.K_optim['minbf'] = 10        # 1000 m3/week # doesn't work with min baseflow higher than 0 because some catchments are dry during summer...
        self.Kgwc = dict(zip(self.K_optim['Catchment'], self.K_optim['K']))  # K parameter weeks
        self.Sinigwc = dict(zip(self.K_optim['Catchment'], self.K_optim['Sinigwc'])) # Storage intitial GW 1000 m3
        self.minbf = dict(zip(self.K_optim['Catchment'], self.K_optim['minbf'])) # min Baseflow 1000 m3§/week 
        self.Qbaseini = dict()  #initial BaseFlow 1000 m3/week
        for c in self.ncatch:
                self.Qbaseini[c] = self.Sinigwc[c]/self.Kgwc[c]  
        
        # =============================================================================
        # Water Demand and WTP per WSA
        # =============================================================================     
        
        df_lookup = self.WSA.set_index('WSAID')
        df_lookup['D_HH'] = df_lookup['Wateruse households (1000m3)'] / 52.18
        df_lookup['D_Ind'] = df_lookup['Wateruse industries (1000m3)'] / 52.18
        df_lookup['D_PS'] = df_lookup['Wateruse services (1000m3)'] / 52.18
        df_lookup['D_Agri'] = df_lookup['Wateruse agriculture (1000m3)'] / 52.18
               
        self.D_HH = dict()
        self.D_Ind = dict()
        self.D_PS = dict()
        self.D_Agri = dict()      
       
        for w in self.nwsa:
            for t in self.ntimes:
                self.D_HH[w,t] = df_lookup.loc[w, 'D_HH']
                self.D_Ind[w,t] = df_lookup.loc[w, 'D_Ind']
                self.D_PS[w,t] = df_lookup.loc[w, 'D_PS']
                self.D_Agri[w,t] = df_lookup.loc[w, 'D_Agri']
                      
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
            # maxpump[wf] = WF.loc[WF['WFID'] == wf, 'AnlgTillad'].values[0]/52.18/1000 # weekly maxpump 1000m3
            self.maxpump[wf] = self.WF.loc[self.WF['WFID'] == wf, 'AnlgTillad'].values[0]/1000   # yearly maxpump 1000m3
        
        # =============================================================================
        # Storage capacity and initial storage of WaterWorks
        # =============================================================================
              
        self.maxstorage = dict()
        self.Storage_WW_ini = dict()
        for ww in self.nww:
            self.maxstorage[ww] = self.WW.loc[self.WW['WWID'] == ww, 'Storage capacity (1000m3)'].values[0]
            self.Storage_WW_ini[ww] = self.WW.loc[self.WW['WWID'] == ww, 'Storage initial (1000m3)'].values[0]
            
        # =============================================================================
        # Water exchange dictionary 
        # =============================================================================
        
        self.ww_ww_cap = dict(zip(self.nww_ww,self.WW_WW['Capacity (1000m3/day)']*7))  # transfer capacity per 1000m3/day converted to 1000m3/week
                          
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
            self.loss_fraction[c] = self.catchment_data.loc[self.catchment_data['Catch_ID'] == c, 'Loss fraction'].iloc[0]

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
            K = self.K_optim['K'][c-1]
            self.Q_natural['Q_base_'+str(c)] = linres(1000,1,K,self.inflow[c-1])  # timestep 1 week for K in weeks
        
        # minbf = dict(zip(ncatch, 0.75*Q_natural.median()))  # new minbf based on median natural baseflow
        self.minbf = dict(zip(self.ncatch, 0.75*self.Q_natural.mean()))  # new minbf based on average natural baseflow
        
        # =============================================================================
        # Area dictionnary
        # =============================================================================
            
        self.area=dict()
        for c in self.ncatch:
            self.area[c] = self.catchment_data.loc[self.catchment_data['Catch_ID']==c, 'Area (m2)'].values[0]
            
        
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
        
        self.model.A_HH  = pyo.Var(self.model.nwsa, self.model.ntimes, within=pyo.NonNegativeReals) # Allocation to households, 1000 m3 per weekly time step
        self.model.A_Ind  = pyo.Var(self.model.nwsa, self.model.ntimes, within=pyo.NonNegativeReals) # Allocation to Industry, 1000 m3 per weekly time step
        self.model.A_PS  = pyo.Var(self.model.nwsa, self.model.ntimes, within=pyo.NonNegativeReals) # Allocation to Public services, 1000 m3 per weekly time step
        self.model.A_Agri = pyo.Var(self.model.nwsa, self.model.ntimes, within=pyo.NonNegativeReals) # Allocation to Agriculture, 1000 m3 per weekly time step
        
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
        self.model.D_HH  = pyo.Param(self.model.nwsa, self.model.ntimes,within=pyo.NonNegativeReals,initialize = self.D_HH) # Set Houshold water demand to observed household water use, 1000 m3 per weekly time step
        self.model.D_Ind = pyo.Param(self.model.nwsa, self.model.ntimes,within=pyo.NonNegativeReals,initialize = self.D_Ind) # Set Industry water demand to observed Industry water use, 1000 m3 per weekly time step
        self.model.D_PS = pyo.Param(self.model.nwsa, self.model.ntimes,within=pyo.NonNegativeReals,initialize = self.D_PS) # Set Public services water demand to observed Energy supply water use, 1000 m3 per weekly time step
        self.model.D_Agri = pyo.Param(self.model.nwsa, self.model.ntimes,within=pyo.NonNegativeReals,initialize = self.D_Agri) # Set Agriculture water demand to observed water supply water use, 1000 m3 per weekly time step
        
        self.model.WTP_HH  = pyo.Param(self.model.nwsa, self.model.ntimes,within=pyo.NonNegativeReals,initialize = self.WTP_HH) # Set Willingness To Pay for the same use categories
        self.model.WTP_Ind = pyo.Param(self.model.nwsa, self.model.ntimes,within=pyo.NonNegativeReals,initialize = self.WTP_Ind) # 
        self.model.WTP_PS = pyo.Param(self.model.nwsa, self.model.ntimes,within=pyo.NonNegativeReals,initialize = self.WTP_PS) # 
        self.model.WTP_Agri = pyo.Param(self.model.nwsa, self.model.ntimes,within=pyo.NonNegativeReals,initialize = self.WTP_Agri) # 
        
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
        
        print(time.strftime("%H:%M:%S") + ' Defining the constraints...')
        
        # =============================================================================
        # Objective function
        # =============================================================================
        
        # Maximize the benefits and minimize the costs
        print(time.strftime("%H:%M:%S") + ' Defining objective function...')
        def obj_rule(model):
            HH_ben = sum(model.WTP_HH[w,t]*model.A_HH[w,t] for w in model.nwsa for t in model.ntimes)
            Ind_ben = sum(model.WTP_Ind[w,t]*model.A_Ind[w,t]  for w in model.nwsa for t in model.ntimes)
            PS_ben = sum(model.WTP_PS[w,t]*model.A_PS[w,t]  for w in model.nwsa for t in model.ntimes)
            Agri_ben = sum(model.WTP_Agri[w,t]*model.A_Agri[w,t]  for w in model.nwsa for t in model.ntimes)
            Pump_cost = sum(model.pumping_cost*model.Pump_WF[wf,t]  for wf in model.nwf for t in model.ntimes) + sum(model.pumping_cost*model.Pump_GW_to_BF[c,t] for c in model.ncatch for t in model.ntimes)
            Exchange_cost = sum(model.exchange_cost*model.Exchange_WW_WW[ww_ww,t] for ww_ww in model.nww_ww for t in model.ntimes)
            return HH_ben + Ind_ben + PS_ben + Agri_ben - Pump_cost - Exchange_cost 
        self.model.obj = pyo.Objective(rule=obj_rule, sense = pyo.maximize)
        
        # =============================================================================
        # Allocation constraints
        # =============================================================================
        
        print(time.strftime("%H:%M:%S") + ' Defining allocation constraints...')
        # Household allocation does not exceed household demand. Active for every time step and catchment, thus two indices
        def wd_HH_c(model, w, t):
            return model.A_HH[w,t] <= model.D_HH[w,t]
        self.model.wd_HH = pyo.Constraint(self.model.nwsa, self.model.ntimes, rule=wd_HH_c)
        
        # Industrial demand constraint per catchment. Active for every time step and catchment, thus two indices
        def wd_Ind_c(model, w, t):
            return model.A_Ind[w, t] <= model.D_Ind[w,t]
        self.model.wd_Ind = pyo.Constraint(self.model.nwsa, self.model.ntimes, rule=wd_Ind_c)
        
        # Public services demand constraint per catchment. Active for every time step and catchment, thus two indices
        def wd_PS_c(model, w, t):
            return model.A_PS[w,t] <= model.D_PS[w,t]
        self.model.wd_PS = pyo.Constraint(self.model.nwsa, self.model.ntimes, rule=wd_PS_c)
        
        # Agriculture demand constraint per catchment. Active for every time step and catchment, thus two indices
        def wd_Agri_c(model, w, t):
            return model.A_Agri[w,t] <= model.D_Agri[w,t]
        self.model.wd_Agri = pyo.Constraint(self.model.nwsa, self.model.ntimes, rule=wd_Agri_c)
        
        
        # =============================================================================
        # Pumping constraints 
        # =============================================================================
        print(time.strftime("%H:%M:%S") + ' Defining pumping constraints...')
        # Pump_catch variable = Total pumping of the WellFields in one catchment
        def pumping_catch_c(model,c,t):
            # get the list of WF in each catchment
            WF_in_catch = self.WF.drop_duplicates(subset='WFID')[self.WF.drop_duplicates(subset='WFID')['Catch_ID'] == c]['WFID'].values
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
            return model.A_HH[wsa,t] + model.A_Ind[wsa,t] + model.A_PS[wsa,t] + model.A_Agri[wsa,t] == sum(
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
            WSA_in_catch = self.WSA[self.WSA['Catch_ID'] == c]['WSAID'].values  # list all WSA in the catchment (for wastewater)
            if t == 1:
                return model.Q_base[c,t] == model.Qbase_ini[c]*np.exp(-1/model.Kgwc[c]) + (model.inflow[c,t]-model.Pump_catch[c,t])*(1-np.exp(-1/model.Kgwc[c])) + (1-model.loss_fraction_waste[c])*sum([model.A_HH[wsa,t] + model.A_Ind[wsa,t] + model.A_PS[wsa,t] + model.A_Agri[wsa,t] for wsa in WSA_in_catch]) + model.Pump_GW_to_BF[c,t]
            else:
                return model.Q_base[c,t] == model.Q_base[c,t-1]*np.exp(-1/model.Kgwc[c]) + (model.inflow[c,t]-model.Pump_catch[c,t])*(1-np.exp(-1/model.Kgwc[c])) + (1-model.loss_fraction_waste[c])*sum([model.A_HH[wsa,t] + model.A_Ind[wsa,t] + model.A_PS[wsa,t] + model.A_Agri[wsa,t] for wsa in WSA_in_catch]) + model.Pump_GW_to_BF[c,t]
            
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
        
        self.A_HH_tot = sum(pyo.value(self.model.A_HH[wsa,t]) for wsa in self.model.nwsa for t in self.model.ntimes)
        self.A_Ind_tot = sum(pyo.value(self.model.A_Ind[wsa,t]) for wsa in self.model.nwsa for t in self.model.ntimes)
        self.A_PS_tot = sum(pyo.value(self.model.A_PS[wsa,t]) for wsa in self.model.nwsa for t in self.model.ntimes)
        self.A_Agri_tot = sum(pyo.value(self.model.A_Agri[wsa,t]) for wsa in self.model.nwsa for t in self.model.ntimes)
        self.A_tot = self.A_HH_tot + self.A_Ind_tot + self.A_PS_tot + self.A_Agri_tot
        
        self.D_HH_tot = sum(pyo.value(self.model.D_HH[wsa,t]) for wsa in self.model.nwsa for t in self.model.ntimes)
        self.D_Ind_tot = sum(pyo.value(self.model.D_Ind[wsa,t]) for wsa in self.model.nwsa for t in self.model.ntimes)
        self.D_PS_tot = sum(pyo.value(self.model.D_PS[wsa,t]) for wsa in self.model.nwsa for t in self.model.ntimes)
        self.D_Agri_tot = sum(pyo.value(self.model.D_Agri[wsa,t]) for wsa in self.model.nwsa for t in self.model.ntimes)
        self.D_tot = self.D_HH_tot + self.D_Ind_tot + self.D_PS_tot + self.D_Agri_tot
        
        print('Allocation HH (%demand) : ',round(100*self.A_HH_tot/self.D_HH_tot,1), '%')
        print('Allocation Ind (%demand) : ',round(100*self.A_Ind_tot/self.D_Ind_tot,1), '%')
        print('Allocation PS (%demand) : ',round(100*self.A_PS_tot/self.D_PS_tot,1), '%')
        print('Allocation total (%demand)',round(100*self.A_tot/self.D_tot,1), '%')
        
        self.Deficit_HH_per_WSA = [100*sum(pyo.value(self.model.A_HH[wsa,t]) for t in self.model.ntimes)/(sum(pyo.value(self.model.D_HH[wsa,t]) for t in self.model.ntimes)+0.000000001) for wsa in self.model.nwsa]
        self.Deficit_Ind_per_WSA = [100*sum(pyo.value(self.model.A_Ind[wsa,t]) for t in self.model.ntimes)/sum(pyo.value(self.model.D_Ind[wsa,t]) for t in self.model.ntimes) for wsa in self.model.nwsa]
        self.Deficit_PS_per_WSA = [100*sum(pyo.value(self.model.A_PS[wsa,t]) for t in self.model.ntimes)/sum(pyo.value(self.model.D_PS[wsa,t]) for t in self.model.ntimes) for wsa in self.model.nwsa]
        self.Deficit_per_WSA = pd.DataFrame(data={'HH':self.Deficit_HH_per_WSA, 'Ind':self.Deficit_Ind_per_WSA, 'PS':self.Deficit_PS_per_WSA})
        return self
        
        
    def save_results(self,**kwargs):
        
        print(time.strftime("%H:%M:%S") + ' Saving the results...')
        

        self.dec_outfile = self.savepath + os.sep + self.run_name + '_Optimal_Decision.xlsx'
        
        # =============================================================================
        # Process decision variabes data
        # =============================================================================
        
        optimal_A_HH = np.zeros((len(self.model.ntimes),len(self.model.nwsa)))
        optimal_A_Ind = np.zeros((len(self.model.ntimes),len(self.model.nwsa)))
        optimal_A_PS = np.zeros((len(self.model.ntimes),len(self.model.nwsa)))
        optimal_A_Agri = np.zeros((len(self.model.ntimes),len(self.model.nwsa)))
        
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
            
            for i in range(1,len(self.model.nwsa)+1):
                optimal_A_HH[t-1,i-1] = pyo.value(self.model.A_HH[self.model.nwsa.at(i),t])
                optimal_A_Ind[t-1,i-1] = pyo.value(self.model.A_Ind[self.model.nwsa.at(i),t])
                optimal_A_PS[t-1,i-1] = pyo.value(self.model.A_PS[self.model.nwsa.at(i),t])
                optimal_A_Agri[t-1,i-1] = pyo.value(self.model.A_Agri[self.model.nwsa.at(i),t])
        
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
        
        optimal_A_HH = pd.DataFrame(optimal_A_HH, index=self.model.ntimes, columns=['A_HH_'+str(w) for w in self.model.nwsa])
        optimal_A_Ind = pd.DataFrame(optimal_A_Ind, index=self.model.ntimes, columns=['A_Ind_'+str(w) for w in self.model.nwsa])
        optimal_A_PS = pd.DataFrame(optimal_A_PS, index=self.model.ntimes, columns=['A_PS_'+str(w) for w in self.model.nwsa])
        optimal_A_Agri = pd.DataFrame(optimal_A_Agri, index=self.model.ntimes, columns=['A_Agri_'+str(w) for w in self.model.nwsa])
        
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
        
        self.optimal_Decision = pd.concat([optimal_time, optimal_A_HH, optimal_A_Ind, optimal_A_PS, optimal_A_Agri, optimal_Storage_WW, optimal_Pump_WF, optimal_Pump_catch, optimal_Pump_GW_to_BF, optimal_Q_base, optimal_Storage_LinRes, optimal_Exchange_WW_WW, optimal_Supply_WW_WSA,optimal_Supply_WF_WW], axis=1) 
        self.optimal_Decision.to_excel(self.dec_outfile,sheet_name = 'Decision variables')
        

        
        self.sp_outfile =     self.savepath + os.sep + self.run_name + '_Shadow_Prices.xlsx'
        
        # =============================================================================
        # Process Shadow prices data
        # =============================================================================
        
        SP_wd_HH = np.zeros((len(self.model.ntimes),len(self.model.nwsa)))
        SP_wd_Ind = np.zeros((len(self.model.ntimes),len(self.model.nwsa)))
        SP_wd_PS = np.zeros((len(self.model.ntimes),len(self.model.nwsa)))
        SP_wd_Agri = np.zeros((len(self.model.ntimes),len(self.model.nwsa)))
        
        SP_pumping_WF = np.zeros((len(self.model.nyear),len(self.model.nwf)))
        SP_wb_WW_Storage = np.zeros((len(self.model.ntimes),len(self.model.nww)))
        SP_wb_WW_Exchange = np.zeros((len(self.model.ntimes),len(self.model.nww_ww)))
        
        # and Supply_WW_WSA ???
        
        SP_lin_res = np.zeros((len(self.model.ntimes),len(self.model.ncatch)))
        SP_min_bf =np.zeros((len(self.model.nyear),len(self.model.ncatch)))
        SP_gw_ind_2 = np.zeros((len(self.model.nyear),len(self.model.ncatch)))
        
        
        for t in self.model.ntimes:
            
            for i in range(1,len(self.model.nwsa)+1):      
                SP_wd_HH [t-1,i-1] = self.model.dual[self.model.wd_HH[self.model.nwsa.at(i),t]]
                SP_wd_Ind [t-1,i-1] = self.model.dual[self.model.wd_Ind[self.model.nwsa.at(i),t]]
                SP_wd_PS [t-1,i-1] = self.model.dual[self.model.wd_PS[self.model.nwsa.at(i),t]]
                SP_wd_Agri [t-1,i-1] = self.model.dual[self.model.wd_Agri[self.model.nwsa.at(i),t]]
                
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
        
        SP_wd_HH = pd.DataFrame(SP_wd_HH, index=self.model.ntimes, columns=['SP_wd_HH_'+str(w) for w in self.model.nwsa])
        SP_wd_Ind = pd.DataFrame(SP_wd_Ind, index=self.model.ntimes, columns=['SP_wd_Ind_'+str(w) for w in self.model.nwsa])
        SP_wd_PS = pd.DataFrame(SP_wd_PS, index=self.model.ntimes, columns=['SP_wd_PS_'+str(w) for w in self.model.nwsa])
        SP_wd_Agri = pd.DataFrame(SP_wd_Agri, index=self.model.ntimes, columns=['SP_wd_Agri_'+str(w) for w in self.model.nwsa])
        
        SP_wb_WW_Storage = pd.DataFrame(SP_wb_WW_Storage, index=self.model.ntimes, columns=['SP_wb_WW_Storage_'+str(w) for w in self.model.nww])
        SP_wb_WW_Exchange = pd.DataFrame(SP_wb_WW_Storage, index=self.model.ntimes, columns=['SP_wb_WW_Exchange_'+str(w) for w in self.model.nww_ww])
        SP_pumping_WF = pd.DataFrame(SP_pumping_WF, index=self.model.nyear, columns=['SP_pumping_WF_'+str(w) for w in self.model.nwf])
        
        SP_lin_res = pd.DataFrame(SP_lin_res, index=self.model.ntimes, columns=['SP_lin_res_'+str(c) for c in self.model.ncatch])
        SP_min_bf = pd.DataFrame(SP_min_bf, index=self.model.nyear, columns=['SP_min_bf_'+str(c) for c in self.model.ncatch])
        
        if self.gw_ind_2_in_use:
            SP_gw_ind_2 = pd.DataFrame(SP_gw_ind_2, index=self.model.nyear, columns=['SP_gw_ind_2_'+str(c) for c in self.model.ncatch])
        
        SP_time = pd.DataFrame({'time':[t for t in self.model.ntimes]}, index=self.model.ntimes)
        if self.gw_ind_2_in_use:
            self.SP = pd.concat([SP_time, SP_wd_HH, SP_wd_Ind, SP_wd_PS, SP_wd_Agri, SP_wb_WW_Storage, SP_wb_WW_Exchange, SP_pumping_WF, SP_lin_res, SP_min_bf, SP_gw_ind_2], axis=1) 
        else:
            self.SP = pd.concat([SP_time, SP_wd_HH, SP_wd_Ind, SP_wd_PS, SP_wd_Agri, SP_wb_WW_Storage, SP_wb_WW_Exchange, SP_pumping_WF, SP_lin_res, SP_min_bf], axis=1) 
        self.SP.to_excel(self.sp_outfile,sheet_name = 'Shadow prices')
        return self
    
    
    
    
