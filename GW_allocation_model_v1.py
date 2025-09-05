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

#%% Import libraries

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
        #These calibration factors are from Geolux, see also UAWOS D2.1
        self.maxrange = 50
        self.extension = 0.3
        self.fft_halfsize = 2048
        self.internal_dist = 0.1744
        self.corrfac = 1.1469

#%% Set-up the correct path

# =============================================================================
# Insert the path to the folder GW_allocation_model below
# =============================================================================

path_Model_folder = r'c:\Users\vpk410\Documents\\GW_allocation_model-main' # CHANGE TO YOUR PATH

# =============================================================================
# Insert the path to the Solver below
# =============================================================================

solvername='cplex'  # 'glpk' 'cplex'
# solverpath_exe = r'F:\Data\s232484\winglpk-4.65\glpk-4.65\\w64\\glpsol'
solverpath_exe = r'c:\Program Files\IBM\ILOG\CPLEX_Studio2211\cplex\bin\x64_win64\cplex.exe'

# =============================================================================
# Savepath
# =============================================================================

savepath = path_Model_folder + r'/Results' 

#%% Set global parameter

# =============================================================================
# Global parameters of the simulation
# =============================================================================

Catch='Watersheds'   # name differently if using different catchents
Geo='Capital'  # name differently if using different catchents
error_func='NSE'     # NSE, RMSE or r2 (the error function used for parametrizing the GW linear reservoir model)
weeks = 1721 # Timesteps, maximum 1749 weeks  # 10 years = 521 weeks, 20 years = 1043 weeks, 30 years = 1565 weeks, 33 years = 1721 weeks

#%% Load the Data

ntimes = np.arange(1, weeks+1, 1)   # weekly timesteps
Catch_Geo = Catch + '_' + Geo

print('    ', Catch, Geo, 'for', weeks, 'weeks')
start = time.time()
print(time.strftime("%H:%M:%S") + ' Importing data...')

# =============================================================================
# Open the tables
# =============================================================================

os.chdir(path_Model_folder + '/Input data model')

Catchments=pd.read_csv('Table_Catchments_'+Catch_Geo+'.csv') # refers to file h:\Students\Clément\GW_allocation_model-main\Shapefiles\Catchment_Watersheds_Capital.shp
WTP=pd.read_csv('Table_WTP_'+Catch_Geo+'.csv')

WF=pd.read_csv('wellfields.csv',sep = ';')# refers to file h:\Students\Clément\GW_allocation_model-main\Shapefiles\wellfields.shp
WW=pd.read_csv('waterworks.csv',sep = ';')# refers to file h:\Students\Clément\GW_allocation_model-main\Shapefiles\waterworks.shp
WSA=pd.read_csv('Table_WSA_'+Catch_Geo+'.csv', sep = ';')# refers to file h:\Students\Clément\GW_allocation_model-main\Shapefiles\WSA.shp


WF_WW=pd.read_csv('WF_WW.csv',sep = ';') # This table should contain all existing connections from wellfields to waterworks. Often wellfields and waterworks are the same.
WW_WSA=pd.read_csv('WW_WSA.csv',sep = ';')# This table should contain all existing connections from waterworks to water supply areas. Currently, there are no capacity limits.
WW_WW=pd.read_csv('WW_WW.csv',sep = ';')# This table should contain all existing connections between waterworks and corresponding capacities. The assumption is that water can flow both ways for now.

# =============================================================================
# Open the Water balance data 
# =============================================================================

os.chdir(path_Model_folder + '/Input data model/WB_' + Catch_Geo)

ncatch=np.array(Catchments['Catch_ID'])
K_optim=pd.read_csv('K_optim_'+error_func+'.csv')

inflow = np.empty((len(ncatch), len(ntimes))) # 1749 weeks of data, approx 33 years
WB_data=[]
for i in range(1,len(ncatch)+1):
    data=pd.read_csv('WB_SZ_'+Catch+'_'+str(int(i))+'.csv')
    WB_data.append(data)
    inflow[i-1,:] = data['MIKE SHE GW recharge (mm)'][:len(ntimes)]/1000*Catchments['Area (m2)'][i-1]/1000   # from mm to 10^3 m3
    
    #remove the negative data in the inflow ! 
    for t in range(1,len(ntimes)):
        if inflow[i-1,t-1]<0:
            inflow[i-1,t] += inflow[i-1,t-1]
            inflow[i-1,t-1] = 0
    if inflow[i-1,len(ntimes)-1]<0:
        inflow[i-1,len(ntimes)-1] = 0   

#%% Process the data

print(time.strftime("%H:%M:%S") + ' Processing the data...')

# =============================================================================
# Indeces (ntimes, nwsa, ncatch, nww, nwf)
# =============================================================================

ncatch=np.array(Catchments['Catch_ID'])
#ntimes already defined
nyear=np.arange(1, int(len(ntimes)//52.18)+2, 1) # yearly index 
nwsa = np.array(WSA['WSAID'])
nww = np.array(WW['WWID'])
nwf = np.array(WF['WFID'])

# =============================================================================
# Double indeces for connections
# =============================================================================

nwf_ww = list(WF_WW.itertuples(index=False, name=None))
nww_wsa = list(WW_WSA.itertuples(index=False, name=None))
nww_ww = list(WW_WW.drop('Capacity (1000m3/day)',axis=1).itertuples(index=False, name=None))

# =============================================================================
# Convert the WTP in weekly time steps
# =============================================================================

WTP_weekly=pd.DataFrame(columns=['Time step','WTP_HH','WTP_Ind','WTP_PS','WTP_Agri'])
WTP_weekly['Time step']=data['Start time of weekly time step']
for i in range(0,len(WTP_weekly)):
    year=int(WTP_weekly['Time step'][i][:4])
    cond=WTP['Year']==year
    WTP_weekly.loc[i, 'WTP_HH'] = float(WTP['WTP Households (DKK/m3)'][cond].iloc[0])
    WTP_weekly.loc[i, 'WTP_Ind'] = float(WTP['WTP Industry (DKK/m3)'][cond].iloc[0])
    WTP_weekly.loc[i, 'WTP_PS'] = float(WTP['WTP Services (DKK/m3)'][cond].iloc[0])
    WTP_weekly.loc[i, 'WTP_Agri'] = float(WTP['WTP Agriculture (DKK/m3)'][cond].iloc[0])
    

# =============================================================================
# Linear reservoirs initial parameters
# =============================================================================

K_optim['Sinigwc'] = 1000   # 1000m3
K_optim['minbf'] = 10        # 1000 m3/week # doesn't work with min baseflow higher than 0 because some catchments are dry during summer...

K_optim['K'] = K_optim['K']
Kgwc=dict(zip(K_optim['Catchment'], K_optim['K']))  # K parameter weeks
if error_func != 'NSE' and Catch!='Test':
    Kgwc[2] = 10.85 # put a more realistic value for catch 2 (the value from the NSE optim)

Sinigwc = dict(zip(K_optim['Catchment'], K_optim['Sinigwc'])) # Storage intitial GW 1000 m3
minbf = dict(zip(K_optim['Catchment'], K_optim['minbf'])) # min Baseflow 1000 m3§/week 

Qbaseini = dict()  #initial BaseFlow 1000 m3/week
for c in ncatch:
    Qbaseini[c] = Sinigwc[c]/Kgwc[c]  

# =============================================================================
# Water Demand and WTP per WSA
# =============================================================================

df_lookup = WSA.set_index('WSAID')
df_lookup['D_HH'] = df_lookup['Wateruse households (1000m3)'] / 52.18
df_lookup['D_Ind'] = df_lookup['Wateruse industries (1000m3)'] / 52.18
df_lookup['D_PS'] = df_lookup['Wateruse services (1000m3)'] / 52.18
df_lookup['D_Agri'] = df_lookup['Wateruse agriculture (1000m3)'] / 52.18

wtp_hh_val = WTP_weekly['WTP_HH'].iloc[-1]
wtp_agri_val = WTP_weekly['WTP_Agri'].iloc[-1]

D_HH = dict()
D_Ind = dict()
D_PS = dict()
D_Agri = dict()

WTP_HH = dict()
WTP_Ind = dict()
WTP_PS = dict()
WTP_Agri = dict()

for w in nwsa:
    for t in ntimes:
        D_HH[w,t] = df_lookup.loc[w, 'D_HH']
        D_Ind[w,t] = df_lookup.loc[w, 'D_Ind']
        D_PS[w,t] = df_lookup.loc[w, 'D_PS']
        D_Agri[w,t] = df_lookup.loc[w, 'D_Agri']

        WTP_HH[w,t] = wtp_hh_val
        WTP_Ind[w,t] = wtp_hh_val * 1/3
        WTP_PS[w,t] = wtp_hh_val * 1/2
        WTP_Agri[w,t] = wtp_agri_val * 1/4

    print(w)
        
        
# =============================================================================
# Inflow dictionary
# =============================================================================

I_inflow = dict()
for c in ncatch:
    for t in ntimes:
        I_inflow[c,t] = inflow[c-1,t-1]

# =============================================================================
# Max pumping capacity
# =============================================================================

maxpump = dict()
for wf in nwf:
    # maxpump[wf] = WF.loc[WF['WFID'] == wf, 'AnlgTillad'].values[0]/52.18/1000 # weekly maxpump 1000m3
    maxpump[wf] = WF.loc[WF['WFID'] == wf, 'AnlgTillad'].values[0]/1000   # yearly maxpump 1000m3
    if Catch=='Test':
        maxpump[wf]=10*maxpump[wf] # to actually see something

# =============================================================================
# Storage capacity and initial storage of WaterWorks
# =============================================================================

# Change Tinghøj reservoir capacity : 250,000 m3 that is CPH demand for 2 days
WW.loc[WW['WWID'] == 1, 'Storage capacity (1000m3)'] = 250

maxstorage = dict()
Storage_ini = dict()
for ww in nww:
    maxstorage[ww] = WW.loc[WW['WWID'] == ww, 'Storage capacity (1000m3)'].values[0]
    Storage_ini[ww] = WW.loc[WW['WWID'] == ww, 'Storage initial (1000m3)'].values[0]
    
# =============================================================================
# Water exchange dictionary 
# =============================================================================

ww_ww_cap = dict(zip(nww_ww,WW_WW['Capacity (1000m3/day)']*7))  # transfer capacity per 1000m3/day converted to 1000m3/week
                  
# =============================================================================
# Create the year / week dictionnary
# =============================================================================

week_in_year = dict()
for t in ntimes:
    week_in_year[t]=int(t//52.18)+1 

# =============================================================================
# Create a month / week dictionnary
# =============================================================================

week_in_month = dict()
for t in ntimes:
    week_in_month[t] = int((t % 52.18)//(52.18/12)+1)

# =============================================================================
# Loss fraction wastewater (what goes to the sea)           
# =============================================================================

loss_fraction=dict()
for c in ncatch:
    loss_fraction[c] = 0.7

if Catch_Geo == 'Watersheds_Capital':
    loss_fraction[16] = 1   # Copenhagen
    loss_fraction[27] = 1   # Amager, taarnby
    
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

Q_natural = pd.DataFrame(index=ntimes)
for c in ncatch:
    K = K_optim['K'][c-1]
    Q_natural['Q_base_'+str(c)] = linres(1000,1,K,inflow[c-1])  # timestep 1 week for K in weeks

# minbf = dict(zip(ncatch, 0.75*Q_natural.median()))  # new minbf based on median natural baseflow
minbf = dict(zip(ncatch, 0.75*Q_natural.mean()))  # new minbf based on average natural baseflow

# =============================================================================
# Area dictionnary
# =============================================================================
    
area=dict()
for c in ncatch:
    area[c] = Catchments.loc[Catchments['Catch_ID']==c, 'Area (m2)'].values[0]


#%% Create Pyomo Model

# =============================================================================
# Create the model
# =============================================================================
print(time.strftime("%H:%M:%S") + ' Creating the model...')
model = pyo.ConcreteModel() # define the model

# =============================================================================
# Define the index 
# =============================================================================

model.ntimes = pyo.Set(dimen=1,initialize=ntimes) # define time index, set values to ntimes
model.nyear = pyo.Set(dimen=1,initialize=nyear)   # define year index
model.nwsa = pyo.Set(dimen=1,initialize=nwsa)     # define WSA index, set values to nwsa
model.nww = pyo.Set(dimen=1,initialize=nww)     # define WaterWorks index, set values to nww
model.nwf = pyo.Set(dimen=1,initialize=nwf)     # define Wellfields index, set values to nwf
model.ncatch = pyo.Set(dimen=1,initialize=ncatch) # define catchment index, set values to ncatch
model.nwf_ww = pyo.Set(dimen=2,initialize=nwf_ww) # define index for all wf-ww connections
model.nww_wsa = pyo.Set(dimen=2,initialize=nww_wsa) # define index for all ww-wsa connections
model.nww_ww = pyo.Set(dimen=2,initialize=nww_ww) # define index for all ww-ww connections

# =============================================================================
# Declare decision variables - decision variable values will be provided by the optimizer
# =============================================================================

model.A_HH  = pyo.Var(model.nwsa, model.ntimes, within=pyo.NonNegativeReals) # Allocation to households, 1000 m3 per weekly time step
model.A_Ind  = pyo.Var(model.nwsa, model.ntimes, within=pyo.NonNegativeReals) # Allocation to Industry, 1000 m3 per weekly time step
model.A_PS  = pyo.Var(model.nwsa, model.ntimes, within=pyo.NonNegativeReals) # Allocation to Public services, 1000 m3 per weekly time step
model.A_Agri = pyo.Var(model.nwsa, model.ntimes, within=pyo.NonNegativeReals) # Allocation to Agriculture, 1000 m3 per weekly time step

model.Pump_WF = pyo.Var(model.nwf, model.ntimes, within=pyo.NonNegativeReals) # Sum of groundwater pumping for each wellfields 1000 m3 per weekly time step
model.Pump_catch = pyo.Var(model.ncatch, model.ntimes, within=pyo.NonNegativeReals) # Sum of groundwater pumping for each catchment 1000 m3 per weekly time step
model.Pump_GW_to_BF = pyo.Var(model.ncatch, model.ntimes, within=pyo.NonNegativeReals)  # Pumping to the river to maintain a min BF 
model.Supply_WF_WW = pyo.Var(model.nwf_ww, model.ntimes, within=pyo.NonNegativeReals) # Supply from WF to WW 1000m3/week
model.Storage_WW = pyo.Var(model.nww, model.ntimes, within=pyo.NonNegativeReals) # Water storage for each waterworks 1000m3 per weekly time step 
model.Exchange_WW_WW  = pyo.Var(model.nww_ww, model.ntimes, within=pyo.NonNegativeReals) # Water transfer from 1 anlaeg to another, therefore 2 times nanlaeg, 1000m3 per weekly time step
model.Supply_WW_WSA = pyo.Var(model.nww_wsa, model.ntimes, within=pyo.NonNegativeReals) # Water Supply distributed by each Waterworks to all the WSA it serves

model.Q_base  = pyo.Var(model.ncatch, model.ntimes, within=pyo.NonNegativeReals) # Base flow from GW catchment, 1000 m3 per weekly time step
model.Storage_LinRes   = pyo.Var(model.ncatch, model.ntimes, within=pyo.NonNegativeReals) # One end storage per month and per reservoir. 1000 m3 per weekly time step

# =============================================================================
# Declare parameters
# =============================================================================

#model.endtime = Param(initialize = ntimes[-1]) # find end time step of the model
model.D_HH  = pyo.Param(model.nwsa, model.ntimes,within=pyo.NonNegativeReals,initialize = D_HH) # Set Houshold water demand to observed household water use, 1000 m3 per weekly time step
model.D_Ind = pyo.Param(model.nwsa, model.ntimes,within=pyo.NonNegativeReals,initialize = D_Ind) # Set Industry water demand to observed Industry water use, 1000 m3 per weekly time step
model.D_PS = pyo.Param(model.nwsa, model.ntimes,within=pyo.NonNegativeReals,initialize = D_PS) # Set Public services water demand to observed Energy supply water use, 1000 m3 per weekly time step
model.D_Agri = pyo.Param(model.nwsa, model.ntimes,within=pyo.NonNegativeReals,initialize = D_Agri) # Set Agriculture water demand to observed water supply water use, 1000 m3 per weekly time step

model.WTP_HH  = pyo.Param(model.nwsa, model.ntimes,within=pyo.NonNegativeReals,initialize = WTP_HH) # Set Willingness To Pay for the same use categories
model.WTP_Ind = pyo.Param(model.nwsa, model.ntimes,within=pyo.NonNegativeReals,initialize = WTP_Ind) # 
model.WTP_PS = pyo.Param(model.nwsa, model.ntimes,within=pyo.NonNegativeReals,initialize = WTP_PS) # 
model.WTP_Agri = pyo.Param(model.nwsa, model.ntimes,within=pyo.NonNegativeReals,initialize = WTP_Agri) # 

model.maxpump_WF = pyo.Param(model.nwf, within=pyo.NonNegativeReals,initialize = maxpump)  # Abstraction license 1000 m3/week
model.Storage_WW_ini = pyo.Param(model.nww, within=pyo.NonNegativeReals,initialize = Storage_ini) # Initial storage in Waterworks 1000m3
model.maxstorage_WW = pyo.Param(model.nww, within=pyo.NonNegativeReals,initialize = maxstorage) # Max storage capacity 1000 m3
model.maxexchange = pyo.Param(model.nww_ww, within=pyo.NonNegativeReals,initialize = ww_ww_cap)    # Water transfer capacity between waterworks in 1000m3/week

model.pumping_cost = pyo.Param(within=pyo.NonNegativeReals, initialize = 1)   # pump cost in DKK/m3 or thousand DKK/1000m3
model.exchange_cost = pyo.Param(within=pyo.NonNegativeReals, initialize = 1)  # water exchange cost DKK/m3 per distance ?????8
model.loss_fraction_waste = pyo.Param(model.ncatch, within=pyo.NonNegativeReals, initialize = loss_fraction) #loss fraction of wastewater return flow to the river (what goes to the sea....)

model.Storage_LinRes_ini = pyo.Param(model.ncatch, within=pyo.NonNegativeReals,initialize = Sinigwc) # Set initial GW storage for all groundwater catchments
model.Kgwc = pyo.Param(model.ncatch, within=pyo.NonNegativeReals,initialize = Kgwc) # Set time constant for all groundwater catchments
model.Qbase_ini = pyo.Param(model.ncatch, within=pyo.NonNegativeReals,initialize = Qbaseini) # Set initial BaseFlow for all catchments
model.inflow = pyo.Param(model.ncatch, model.ntimes,within=pyo.Reals,initialize = I_inflow) # Set inflow to GW to for all catchments (from MIKE SHE model)
model.minbflow = pyo.Param(model.ncatch,within=pyo.Reals,initialize = minbf) # Set environmental constraint on flow for all catchments


#%% Set up the model

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
model.obj = pyo.Objective(rule=obj_rule, sense = pyo.maximize)

# =============================================================================
# Allocation constraints
# =============================================================================

print(time.strftime("%H:%M:%S") + ' Defining allocation constraints...')
# Household allocation does not exceed household demand. Active for every time step and catchment, thus two indices
def wd_HH_c(model, w, t):
    return model.A_HH[w,t] <= model.D_HH[w,t]
model.wd_HH = pyo.Constraint(model.nwsa, model.ntimes, rule=wd_HH_c)

# Industrial demand constraint per catchment. Active for every time step and catchment, thus two indices
def wd_Ind_c(model, w, t):
    return model.A_Ind[w, t] <= model.D_Ind[w,t]
model.wd_Ind = pyo.Constraint(model.nwsa, model.ntimes, rule=wd_Ind_c)

# Public services demand constraint per catchment. Active for every time step and catchment, thus two indices
def wd_PS_c(model, w, t):
    return model.A_PS[w,t] <= model.D_PS[w,t]
model.wd_PS = pyo.Constraint(model.nwsa, model.ntimes, rule=wd_PS_c)

# Agriculture demand constraint per catchment. Active for every time step and catchment, thus two indices
def wd_Agri_c(model, w, t):
    return model.A_Agri[w,t] <= model.D_Agri[w,t]
model.wd_Agri = pyo.Constraint(model.nwsa, model.ntimes, rule=wd_Agri_c)


# =============================================================================
# Pumping constraints 
# =============================================================================
print(time.strftime("%H:%M:%S") + ' Defining pumping constraints...')
# Pump_catch variable = Total pumping of the WellFields in one catchment
def pumping_catch_c(model,c,t):
    # get the list of WF in each catchment
    WF_in_catch = WF.drop_duplicates(subset='WFID')[WF.drop_duplicates(subset='WFID')['Catch_ID'] == c]['WFID'].values
    return  model.Pump_catch[c,t] == sum([model.Pump_WF[wf,t] for wf in WF_in_catch]) + model.Pump_GW_to_BF[c,t]
model.pumping_catch = pyo.Constraint(model.ncatch, model.ntimes, rule=pumping_catch_c)

# Total pumping for each Wellfields always below the abstraction license each week
# def pumping_WF_c(model, wf, t):
#     return model.Pump_WF[wf,t] <= model.maxpump_WF[wf]/52.18
# model.pumping_WF = Constraint(model.nwf, model.ntimes, rule=pumping_WF_c)

# Total pumping for each Wellfields always below the abstrction license for each year
def pumping_WF_c(model, wf, y):
    list_weeks = np.array([week[0] for week in week_in_year.items() if week[1] == y])  #get the list of the weeks in a given year
    return sum(model.Pump_WF[wf,t] for t in list_weeks) <= model.maxpump_WF[wf]
model.pumping_WF = pyo.Constraint(model.nwf, model.nyear, rule=pumping_WF_c)


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
model.wb_WF = pyo.Constraint(model.nwf, model.ntimes, rule=wb_WF_c)

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
  
model.wb_WW = pyo.Constraint(model.nww, model.ntimes, rule=wb_WW_c)

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

model.wb_WSA = pyo.Constraint(model.nwsa, model.ntimes, rule=wb_WSA_c)

#For each WaterWorks, Storage WW <= maxstorage WW
def wb_WW_Storage_c(model, ww, t):
    return model.Storage_WW[ww,t] <= model.maxstorage_WW[ww]
model.wb_WW_Storage = pyo.Constraint(model.nww, model.ntimes, rule=wb_WW_Storage_c)

# Max exchange capacity
def wb_WW_Exchange_c(model, ww1, ww2,t):
    return model.Exchange_WW_WW[ww1, ww2, t] <= model.maxexchange[ww1,ww2]
model.wb_WW_Exchange = pyo.Constraint(model.nww_ww, model.ntimes, rule=wb_WW_Exchange_c)

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
    WSA_in_catch = WSA[WSA['Catch_ID'] == c]['WSAID'].values  # list all WSA in the catchment (for wastewater)
    if t == 1:
        return model.Q_base[c,t] == model.Qbase_ini[c]*np.exp(-1/model.Kgwc[c]) + (model.inflow[c,t]-model.Pump_catch[c,t])*(1-np.exp(-1/model.Kgwc[c])) + (1-model.loss_fraction_waste[c])*sum([model.A_HH[wsa,t] + model.A_Ind[wsa,t] + model.A_PS[wsa,t] + model.A_Agri[wsa,t] for wsa in WSA_in_catch]) + model.Pump_GW_to_BF[c,t]
    else:
        return model.Q_base[c,t] == model.Q_base[c,t-1]*np.exp(-1/model.Kgwc[c]) + (model.inflow[c,t]-model.Pump_catch[c,t])*(1-np.exp(-1/model.Kgwc[c])) + (1-model.loss_fraction_waste[c])*sum([model.A_HH[wsa,t] + model.A_Ind[wsa,t] + model.A_PS[wsa,t] + model.A_Agri[wsa,t] for wsa in WSA_in_catch]) + model.Pump_GW_to_BF[c,t]
    
model.lin_res = pyo.Constraint(model.ncatch, model.ntimes, rule=lin_res_c)  

# Linear reservoirs storage (from 1st order equation), not needed, only for plotting
def lin_res_stor_c(model,c,t):
    if t == 1:
        return model.Storage_LinRes[c,t] == model.Storage_LinRes_ini[c]*np.exp(-1/model.Kgwc[c]) + (model.inflow[c,t]-model.Pump_catch[c,t])*(1-np.exp(-1/model.Kgwc[c]))*Kgwc[c]
    else:
        return model.Storage_LinRes[c,t] == model.Storage_LinRes[c,t-1]*np.exp(-1/model.Kgwc[c]) + (model.inflow[c,t]-model.Pump_catch[c,t])*(1-np.exp(-1/model.Kgwc[c]))*Kgwc[c]

model.lin_res_stor = pyo.Constraint(model.ncatch, model.ntimes, rule=lin_res_stor_c)


# =============================================================================
# Environmental constraints
# =============================================================================
print(time.strftime("%H:%M:%S") + ' Defining environmental flow constraints...')
# min baseflow
def min_bf_c(model,c,y):
    list_weeks = [week[0] for week in week_in_year.items() if week[1] == y]  #get the list of the weeks in a given year
    return sum(model.Q_base[c,t] for t in list_weeks)/len(list_weeks) >= model.minbflow[c]   # yearly average flow above minBF
model.min_bf = pyo.Constraint(model.ncatch, model.nyear, rule=min_bf_c) 


# Max pumping of aquifer: from "Model and Ensemble Indicator-Guided Assessment of Robust,
#                              Exploitable Groundwater Resources for Denmark" table 1, indicator 2
def gw_ind_2_c(model,c,y):
    list_weeks = [week[0] for week in week_in_year.items() if week[1] == y]  #get the list of the weeks in a given year
    return sum(model.Pump_catch[c,t] for t in list_weeks) <= sum(model.inflow[c,t] for t in list_weeks)/2
in_use = False
if in_use:
    model.gw_ind_2 = pyo.Constraint(model.ncatch, model.nyear, rule=gw_ind_2_c)


# =============================================================================
# Dual problem
# =============================================================================

# formulate dual problem to provide shadow prices
model.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT) 

#%% Solve the model

print(time.strftime("%H:%M:%S") + ' Solving the model...')

# =============================================================================
# Create a solver
# =============================================================================

opt = pyo.SolverFactory(solvername,executable=solverpath_exe)
# opt = SolverFactory(solvername) # if it works without stating the path of the solver

# =============================================================================
# Solve
# =============================================================================

results = opt.solve(model)
end = time.time()

# =============================================================================
# Check status
# =============================================================================

print(time.strftime("%H:%M:%S") + ' Model solved! \n', 'Total computation time = ', round((end-start)/60,1), ' minutes')

if (results.solver.status == pyo.SolverStatus.ok) and (results.solver.termination_condition == pyo.TerminationCondition.optimal):
    # Do something when the solution in optimal and feasible
    print('Optimal and Feasible \n')
elif (results.solver.termination_condition == pyo.TerminationCondition.infeasible):
    # Do something when model in infeasible
    print('Infeasible')
else:
    # Something else is wrong
    print ('Solver Status: ',  results.solver.status)
    print ('Solver termination condition: ', results.solver.termination_condition)


#%% Output

# =============================================================================
# # Objective value
# =============================================================================

print("Total Benefit in optimal solution: ", round(pyo.value(model.obj)/len(model.ntimes)), " thousand DKK per week \n")


# =============================================================================
# Some results
# =============================================================================

A_HH_tot = sum(pyo.value(model.A_HH[wsa,t]) for wsa in model.nwsa for t in model.ntimes)
A_Ind_tot = sum(pyo.value(model.A_Ind[wsa,t]) for wsa in model.nwsa for t in model.ntimes)
A_PS_tot = sum(pyo.value(model.A_PS[wsa,t]) for wsa in model.nwsa for t in model.ntimes)
A_Agri_tot = sum(pyo.value(model.A_Agri[wsa,t]) for wsa in model.nwsa for t in model.ntimes)
A_tot = A_HH_tot + A_Ind_tot + A_PS_tot + A_Agri_tot

D_HH_tot = sum(pyo.value(model.D_HH[wsa,t]) for wsa in model.nwsa for t in model.ntimes)
D_Ind_tot = sum(pyo.value(model.D_Ind[wsa,t]) for wsa in model.nwsa for t in model.ntimes)
D_PS_tot = sum(pyo.value(model.D_PS[wsa,t]) for wsa in model.nwsa for t in model.ntimes)
D_Agri_tot = sum(pyo.value(model.D_Agri[wsa,t]) for wsa in model.nwsa for t in model.ntimes)
D_tot = D_HH_tot + D_Ind_tot + D_PS_tot +D_Agri_tot

print('Allocation HH (%demand) : ',round(100*A_HH_tot/D_HH_tot,1), '%')
print('Allocation Ind (%demand) : ',round(100*A_Ind_tot/D_Ind_tot,1), '%')
print('Allocation PS (%demand) : ',round(100*A_PS_tot/D_PS_tot,1), '%')
print('Allocation total (%demand)',round(100*A_tot/D_tot,1), '%')

Deficit_HH_per_WSA = [100*sum(pyo.value(model.A_HH[wsa,t]) for t in model.ntimes)/(sum(pyo.value(model.D_HH[wsa,t]) for t in model.ntimes)+0.000000001) for wsa in model.nwsa]
Deficit_Ind_per_WSA = [100*sum(pyo.value(model.A_Ind[wsa,t]) for t in model.ntimes)/sum(pyo.value(model.D_Ind[wsa,t]) for t in model.ntimes) for wsa in model.nwsa]
Deficit_PS_per_WSA = [100*sum(pyo.value(model.A_PS[wsa,t]) for t in model.ntimes)/sum(pyo.value(model.D_PS[wsa,t]) for t in model.ntimes) for wsa in model.nwsa]
Deficit_per_WSA = pd.DataFrame(data={'HH':Deficit_HH_per_WSA, 'Ind':Deficit_Ind_per_WSA, 'PS':Deficit_PS_per_WSA})


#%% Save optimal decisions

print(time.strftime("%H:%M:%S") + ' Saving the results...')

os.chdir(savepath)
outfile = r'Optimal_Decision_'+Catch_Geo+'.xlsx'

# =============================================================================
# Process decision variabes data
# =============================================================================

optimal_A_HH = np.zeros((len(model.ntimes),len(model.nwsa)))
optimal_A_Ind = np.zeros((len(model.ntimes),len(model.nwsa)))
optimal_A_PS = np.zeros((len(model.ntimes),len(model.nwsa)))
optimal_A_Agri = np.zeros((len(model.ntimes),len(model.nwsa)))

optimal_Pump_WF = np.zeros((len(model.ntimes),len(model.nwf)))
optimal_Pump_catch = np.zeros((len(model.ntimes),len(model.ncatch)))
optimal_Pump_GW_to_BF = np.zeros((len(model.ntimes),len(model.ncatch)))
optimal_Storage_WW = np.zeros((len(model.ntimes),len(model.nww)))
optimal_Q_base = np.zeros((len(model.ntimes),len(model.ncatch)))
optimal_Storage_LinRes = np.zeros((len(model.ntimes),len(model.ncatch)))

optimal_Exchange_WW_WW = np.zeros((len(model.ntimes),len(model.nww_ww)))
optimal_Supply_WW_WSA = np.zeros((len(model.ntimes),len(model.nww_wsa)))
optimal_Supply_WF_WW = np.zeros((len(model.ntimes),len(model.nwf_ww)))

for t in model.ntimes:
    
    for i in range(1,len(model.nwsa)+1):
        optimal_A_HH[t-1,i-1] = pyo.value(model.A_HH[model.nwsa.at(i),t])
        optimal_A_Ind[t-1,i-1] = pyo.value(model.A_Ind[model.nwsa.at(i),t])
        optimal_A_PS[t-1,i-1] = pyo.value(model.A_PS[model.nwsa.at(i),t])
        optimal_A_Agri[t-1,i-1] = pyo.value(model.A_Agri[model.nwsa.at(i),t])

    for i in range(1,len(model.ncatch)+1):
        optimal_Pump_catch[t-1,i-1] = pyo.value(model.Pump_catch[model.ncatch.at(i),t])
        optimal_Pump_GW_to_BF[t-1,i-1] = pyo.value(model.Pump_GW_to_BF[model.ncatch.at(i),t])
        optimal_Q_base[t-1,i-1] = pyo.value(model.Q_base[model.ncatch.at(i),t])
        optimal_Storage_LinRes[t-1,i-1] = pyo.value(model.Storage_LinRes[model.ncatch.at(i),t])
        
    for i in range(1,len(model.nwf)+1):
        optimal_Pump_WF[t-1,i-1] = pyo.value(model.Pump_WF[model.nwf.at(i),t])
        
    for i in range(1, len(model.nww)+1):
        optimal_Storage_WW[t-1,i-1] = pyo.value(model.Storage_WW[model.nww.at(i),t])
    
    for i in range(1, len(model.nww_ww)+1):    
        optimal_Exchange_WW_WW[t-1,i-1] = pyo.value(model.Exchange_WW_WW[model.nww_ww.at(i),t])
            
    for i in range(1, len(model.nww_wsa)+1):
        optimal_Supply_WW_WSA[t-1,i-1] = pyo.value(model.Supply_WW_WSA[model.nww_wsa.at(i),t])
        
    for i in range(1, len(model.nwf_ww)+1):
        optimal_Supply_WF_WW[t-1,i-1] = pyo.value(model.Supply_WF_WW[model.nwf_ww.at(i),t])

# =============================================================================
# Convert to DataFrame and save
# =============================================================================

optimal_A_HH = pd.DataFrame(optimal_A_HH, index=model.ntimes, columns=['A_HH_'+str(w) for w in model.nwsa])
optimal_A_Ind = pd.DataFrame(optimal_A_Ind, index=model.ntimes, columns=['A_Ind_'+str(w) for w in model.nwsa])
optimal_A_PS = pd.DataFrame(optimal_A_PS, index=model.ntimes, columns=['A_PS_'+str(w) for w in model.nwsa])
optimal_A_Agri = pd.DataFrame(optimal_A_Agri, index=model.ntimes, columns=['A_Agri_'+str(w) for w in model.nwsa])

optimal_Storage_WW = pd.DataFrame(optimal_Storage_WW, index=model.ntimes, columns=['Storage_WW_'+str(w) for w in model.nww])
optimal_Pump_WF = pd.DataFrame(optimal_Pump_WF, index=model.ntimes, columns=['Pump_WF_'+str(w) for w in model.nwf])
optimal_Pump_catch = pd.DataFrame(optimal_Pump_catch, index=model.ntimes, columns=['Pump_catch_'+str(c) for c in model.ncatch])
optimal_Pump_GW_to_BF = pd.DataFrame(optimal_Pump_GW_to_BF, index=model.ntimes, columns=['Pump_GW_to_BF_'+str(c) for c in model.ncatch])
optimal_Q_base = pd.DataFrame(optimal_Q_base, index=model.ntimes, columns=['Q_base_'+str(c) for c in model.ncatch])
optimal_Storage_LinRes = pd.DataFrame(optimal_Storage_LinRes, index=model.ntimes, columns=['Send_'+str(c) for c in model.ncatch])
optimal_Exchange_WW_WW = pd.DataFrame(optimal_Exchange_WW_WW, index=model.ntimes, columns=['Exchange_'+str(c) for c in model.nww_ww])
optimal_Supply_WW_WSA = pd.DataFrame(optimal_Supply_WW_WSA, index=model.ntimes, columns=['Supply_'+str(c) for c in model.nww_wsa])
optimal_Supply_WF_WW = pd.DataFrame(optimal_Supply_WF_WW, index=model.ntimes, columns=['Supply_'+str(c) for c in model.nwf_ww])

optimal_time = pd.DataFrame({'time':[t for t in model.ntimes]}, index=model.ntimes)

optimal_Decision = pd.concat([optimal_time, optimal_A_HH, optimal_A_Ind, optimal_A_PS, optimal_A_Agri, optimal_Storage_WW, optimal_Pump_WF, optimal_Pump_catch, optimal_Pump_GW_to_BF, optimal_Q_base, optimal_Storage_LinRes, optimal_Exchange_WW_WW, optimal_Supply_WW_WSA,optimal_Supply_WF_WW], axis=1) 
optimal_Decision.to_excel(outfile,sheet_name = 'Decision variables')

#%% Save shadow prices

os.chdir(savepath)
outfile =     savepath + os.sep + r'Shadow_Prices_'+Catch_Geo+'.xlsx'

# =============================================================================
# Process Shadow prices data
# =============================================================================

SP_wd_HH = np.zeros((len(model.ntimes),len(model.nwsa)))
SP_wd_Ind = np.zeros((len(model.ntimes),len(model.nwsa)))
SP_wd_PS = np.zeros((len(model.ntimes),len(model.nwsa)))
SP_wd_Agri = np.zeros((len(model.ntimes),len(model.nwsa)))

SP_pumping_WF = np.zeros((len(model.nyear),len(model.nwf)))
SP_wb_WW_Storage = np.zeros((len(model.ntimes),len(model.nww)))
SP_wb_WW_Exchange = np.zeros((len(model.ntimes),len(model.nww_ww)))

# and Supply_WW_WSA ???

SP_lin_res = np.zeros((len(model.ntimes),len(model.ncatch)))
SP_min_bf =np.zeros((len(model.nyear),len(model.ncatch)))
SP_gw_ind_2 = np.zeros((len(model.nyear),len(model.ncatch)))


for t in model.ntimes:
    
    for i in range(1,len(model.nwsa)+1):      
        SP_wd_HH [t-1,i-1] = model.dual[model.wd_HH[model.nwsa.at(i),t]]
        SP_wd_Ind [t-1,i-1] = model.dual[model.wd_Ind[model.nwsa.at(i),t]]
        SP_wd_PS [t-1,i-1] = model.dual[model.wd_PS[model.nwsa.at(i),t]]
        SP_wd_Agri [t-1,i-1] = model.dual[model.wd_Agri[model.nwsa.at(i),t]]
        
    for i in range(1,len(model.ncatch)+1):
        SP_lin_res [t-1,i-1] = model.dual[model.lin_res[model.ncatch.at(i),t]]
        
        
    for i in range(1, len(model.nww)+1):
        SP_wb_WW_Storage[t-1,i-1] = model.dual[model.wb_WW_Storage[model.nww.at(i),t]]
        
    for i in range(1, len(model.nww_ww)+1):
        SP_wb_WW_Exchange[t-1,i-1] = model.dual[model.wb_WW_Exchange[model.nww_ww.at(i),t]]

for y in model.nyear:
    for i in range(1,len(model.nwf)+1):
        SP_pumping_WF[y-1,i-1] = model.dual[model.pumping_WF[model.nwf.at(i),y]]

    for i in range(1,len(model.ncatch)+1):
        SP_min_bf [y-1,i-1] = model.dual[model.min_bf[model.ncatch.at(i),y]]
        SP_gw_ind_2 [y-1,i-1] = model.dual[model.gw_ind_2[model.ncatch.at(i),y]]


# =============================================================================
# Convert SP to dataframe and save
# =============================================================================

SP_wd_HH = pd.DataFrame(SP_wd_HH, index=model.ntimes, columns=['SP_wd_HH_'+str(w) for w in model.nwsa])
SP_wd_Ind = pd.DataFrame(SP_wd_Ind, index=model.ntimes, columns=['SP_wd_Ind_'+str(w) for w in model.nwsa])
SP_wd_PS = pd.DataFrame(SP_wd_PS, index=model.ntimes, columns=['SP_wd_PS_'+str(w) for w in model.nwsa])
SP_wd_Agri = pd.DataFrame(SP_wd_Agri, index=model.ntimes, columns=['SP_wd_Agri_'+str(w) for w in model.nwsa])

SP_wb_WW_Storage = pd.DataFrame(SP_wb_WW_Storage, index=model.ntimes, columns=['SP_wb_WW_Storage_'+str(w) for w in model.nww])
SP_wb_WW_Exchange = pd.DataFrame(SP_wb_WW_Storage, index=model.ntimes, columns=['SP_wb_WW_Exchange_'+str(w) for w in model.nww_ww])
SP_pumping_WF = pd.DataFrame(SP_pumping_WF, index=model.nyear, columns=['SP_pumping_WF_'+str(w) for w in model.nwf])

SP_lin_res = pd.DataFrame(SP_lin_res, index=model.ntimes, columns=['SP_lin_res_'+str(c) for c in model.ncatch])
SP_min_bf = pd.DataFrame(SP_min_bf, index=model.nyear, columns=['SP_min_bf_'+str(c) for c in model.ncatch])
SP_gw_ind_2 = pd.DataFrame(SP_gw_ind_2, index=model.nyear, columns=['SP_gw_ind_2_'+str(c) for c in model.ncatch])

SP_time = pd.DataFrame({'time':[t for t in model.ntimes]}, index=model.ntimes)

SP = pd.concat([SP_time, SP_wd_HH, SP_wd_Ind, SP_wd_PS, SP_wd_Agri, SP_wb_WW_Storage, SP_wb_WW_Exchange, SP_pumping_WF, SP_lin_res, SP_min_bf, SP_gw_ind_2], axis=1) 
SP.to_excel(outfile,sheet_name = 'Shadow prices')




