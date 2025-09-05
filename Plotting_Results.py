# -*- coding: utf-8 -*-
"""
@author: Clément Franey
clefraney@gmail.com

Updated: June 11 2025


Instructions:

Set the path of the folder where the python model is located in the variable "path_Model_folder".

Change the value of the variables "weeks" and "scenario" to plot the results for a given number of weeks and a given scenario.    

"""

#%% Import libraries

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import time
import pickle
from math import log10

#%% Set-up the correct path

# =============================================================================
# Insert the path to the folder GW_allocation_model below
# =============================================================================

path_Model_folder = r'c:\Users\vpk410\Documents\\GW_allocation_model-main' # CHANGE TO YOUR PATH


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
Catch_Geo = Catch + '_' + Geo

# =============================================================================
# Scenarii
# =============================================================================

scenario = 0  # 0 = Baseline , 1 = Maximum abstraction capacity 2 , 2 = local water exchanges

ind_2 = False
if scenario == 1:
    ind_2 = True
    
#%% Load the results of the model

with open(savepath+ os.sep + 'Scenario_' + str(scenario) + '_' + str(weeks) + '_weeks'  + '_indices.pkl', 'rb') as file:
    TIME,list_WSA,list_WW,list_WF,list_catch,list_year = pickle.load(file)
    
with open(savepath+ os.sep + 'Scenario_' + str(scenario) + '_' + str(weeks) + '_weeks'  + '_decision_variables.pkl', 'rb') as file:
    optimal_Decision, optimal_time, optimal_A_HH, optimal_A_Ind, optimal_A_PS, optimal_A_Agri, optimal_Storage_WW, optimal_Pump_WF, optimal_Pump_catch, optimal_Pump_GW_to_BF, optimal_Q_base, optimal_Send, optimal_Exchange, optimal_Supply_WW_WSA = pd.read_pickle(file)
    
with open(savepath+ os.sep + 'Scenario_' + str(scenario) + '_' + str(weeks) + '_weeks'  + '_shadow_prices.pkl', 'rb') as file:
    SP, SP_time, SP_wd_HH, SP_wd_Ind, SP_wd_PS, SP_wd_Agri, SP_wb_WW_Storage, SP_pumping_WF, SP_lin_res, SP_min_bf, SP_gw_ind_2, SP_wb_WW_Exchange = pd.read_pickle(file)


os.chdir(path_Model_folder + '/Input data model')
Catchments=pd.read_csv('Table_Catchments_'+Catch_Geo+'.csv')
WF=pd.read_csv('Table_WF_'+Catch_Geo+'.csv')
WSA=pd.read_csv('Table_WSA_'+Catch_Geo+'.csv')

# =============================================================================
# Open the Water balance data 
# =============================================================================

os.chdir(path_Model_folder + '/Input data model/WB_' + Catch_Geo)

ncatch=np.array(Catchments['Catch_ID'])
K_optim=pd.read_csv('K_optim_'+error_func+'.csv')

inflow = np.empty((len(list_catch), weeks)) # 1749 weeks of data, approx 33 years
WB_data=[]
for i in range(1,len(list_catch)+1):
    data=pd.read_csv('WB_SZ_'+Catch+'_'+str(int(i))+'.csv')
    WB_data.append(data)
    inflow[i-1,:] = data['MIKE SHE GW recharge (mm)'][:weeks]/1000*Catchments['Area (m2)'][i-1]/1000   # from mm to 10^3 m3
    
    #remove the negative data in the inflow ! 
    for t in range(1,weeks):
        if inflow[i-1,t-1]<0:
            inflow[i-1,t] += inflow[i-1,t-1]
            inflow[i-1,t-1] = 0
    if inflow[i-1,weeks-1]<0:
        inflow[i-1,weeks-1] = 0   

# =============================================================================
# Create the year / week dictionnary
# =============================================================================

week_in_year = {}
for t in range(1,weeks+1):
    week_in_year[t]=int(t//52.18)+1 

# =============================================================================
# Create a month / week dictionnary
# =============================================================================

week_in_month = {}
for t in range(1,weeks+1):
    week_in_month[t] = int((t % 52.18)//(52.18/12)+1)
    
# =============================================================================
# Area dictionnary
# =============================================================================
    
area={}
for c in list_catch:
    area[c] = Catchments.loc[Catchments['Catch_ID']==c, 'Area (m2)'].values[0]
    
# =============================================================================
# Max pumping capacity
# =============================================================================

maxpump = {}
for wf in list_WF:
    maxpump[wf] = WF.loc[WF['WFID'] == wf, 'AnlgTillad'].values[0]/1000   # yearly maxpump 1000m3

#%% Plotting time series

plt.figure()
for c in list_catch[20:32]:
    plt.plot(TIME, optimal_Decision['Q_base_'+str(c)], label=str(c))
plt.xlabel('Time')
plt.ylabel('BaseFlow (1000 m3 / week)')
# plt.xticks([date.fromisoformat(str(yr)+'-01-01') for yr in range(1989,2023)])
plt.legend()
plt.title('Baseflow')
plt.show()

plt.figure()
for wsa in list_WSA[:10]:
    plt.plot(SP['time'],SP['SP_wd_HH_'+str(wsa)], label=str(wsa))
# plt.gca().set_xticks([SP['time'][0],SP['time'][500],SP['time'][1000],SP['time'][1500]])
plt.xlabel('Time (weeks)')
plt.ylabel('Shadow price ($kr/m^3$)')
plt.legend()
plt.title('Shadow price of household demand constraint')
plt.show()

plt.figure()
for wsa in list_WSA[:10]:
    plt.plot(SP['time'],SP['SP_wd_Ind_'+str(wsa)], label=str(wsa))
# plt.gca().set_xticks([SP['time'][0],SP['time'][500],SP['time'][1000],SP['time'][1500]])
plt.xlabel('Time (weeks)')
plt.ylabel('Shadow price ($kr/m^3$)')
plt.legend()
plt.title('Shadow price of industry demand constraint')
plt.show()

plt.figure()
for wsa in list_WSA[:10]:
    plt.plot(SP['time'],SP['SP_wd_PS_'+str(wsa)], label=str(wsa))
# plt.gca().set_xticks([SP['time'][0],SP['time'][500],SP['time'][1000],SP['time'][1500]])
plt.xlabel('Time (weeks)')
plt.ylabel('Shadow price ($kr/m^3$)')
plt.legend()
plt.title('Shadow price of public services demand constraint')
plt.show()

plt.figure()
for wsa in list_WSA[:10]:
    plt.plot(SP['time'],SP['SP_wd_Agri_'+str(wsa)], label=str(wsa))
# plt.gca().set_xticks([SP['time'][0],SP['time'][500],SP['time'][1000],SP['time'][1500]])
plt.xlabel('Time (weeks)')
plt.ylabel('Shadow price ($kr/m^3$)')
plt.legend()
plt.title('Shadow price of agriculture demand constraint')
plt.show()

plt.figure()
for wf in list_WF[:10]:
    plt.plot(SP['time'],SP['SP_pumping_WF_'+str(wf)], label=str(wf))
# plt.gca().set_xticks([SP['time'][0],SP['time'][500],SP['time'][1000],SP['time'][1500]])
plt.xlabel('Time (weeks)')
plt.ylabel('Shadow price ($kr/(m^3/week)$)')
plt.legend()
plt.title('Shadow price of pumping WF constraint TO MODIFY FOR YEAR INDEX')
plt.show()

plt.figure()
for c in list_catch[11:25]:
    plt.plot(SP['time'],SP['SP_lin_res_'+str(c)], label=str(c))
# plt.gca().set_xticks([SP['time'][0],SP['time'][500],SP['time'][1000],SP['time'][1500]])
plt.xlabel('Time (weeks)')
plt.ylabel('Shadow price ($kr/(m^3/week)$)')
plt.legend()
plt.title('Shadow price of linear reservoir baseflow constraint')
plt.show()

plt.figure()
for c in list_catch[:10]:
    plt.plot(SP['time'],SP['SP_min_bf_'+str(c)], label=str(c))
# plt.gca().set_xticks([SP['time'][0],SP['time'][500],SP['time'][1000],SP['time'][1500]])
plt.xlabel('Time (weeks)')
plt.ylabel('Shadow price ($kr/(m^3/year)$)')
plt.legend()
plt.title('Shadow price of minimum baseflow constraint TO MODIFY WITH YEAR INDEX')
plt.show()


#%% Plotting bar chart

plt.figure()
plt.bar(np.arange(1,len(list_WF)+1),SP_pumping_WF.mean())
plt.xlabel('Wellfield ID')
plt.ylabel('Shadow price ($DKK/(m^3/week)$)')
# plt.legend()
plt.title('Time-averaged wellfield pumping capacity shadow price')
plt.show()

plt.figure()
plt.bar(np.arange(1,len(list_WW)+1),SP_wb_WW_Storage.mean())
plt.xlabel('WW')
plt.ylabel('Shadow price ($kr/(m^3/week)$)')
# plt.legend()
plt.title('Shadow price of Storage capacity of WW')
plt.show()

plt.figure()
plt.bar(np.arange(1,len(list_catch)+1),SP_lin_res.mean())
plt.xlabel('Catch')
plt.ylabel('Shadow price ($kr/(m^3/week)$)')
# plt.legend()
plt.title('Shadow price of linear reservoir constraint')
plt.show()

plt.figure()
plt.bar(np.arange(1,len(list_catch)+1),SP_min_bf.mean())
plt.xlabel('Catch')
plt.ylabel('Shadow price ($kr/(m^3/year)$)')
# plt.legend()
plt.title('Shadow price of min BaseFlow constraint')
plt.show()

plt.figure()
plt.bar(np.arange(1,len(list_catch)+1),SP_gw_ind_2.mean())
plt.xlabel('Watersheds catchments')
plt.ylabel('Shadow price ($kr/(m^3/year)$)')
# ax.set_xticks(np.arange(1, 33, 1),np.arange(1, 33, 1))
# plt.legend()
plt.title('Yearly average shadow price of the groundwater maximum abstraction constraint')
plt.show()

# =============================================================================
# Water Exchanges 
# =============================================================================

# Average SP on water exchange constraint, received by one WW
ww = 104838
ww_index = list_WW.index(ww)
plt.figure()
plt.bar(np.arange(1,len(list_WW)+1),SP_wb_WW_Exchange[ww_index].mean())
plt.xlabel('WW')
plt.ylabel('Shadow price ($kr/(m^3/week)$)')
# plt.legend()
plt.title('Shadow price of water echange, water received by WW ' + str(list_WW[ww_index]))
plt.show()

# Max of average SP on water exchange constraint, received by all WW
plt.figure()
plt.bar(np.arange(1,len(list_WW)+1),np.array([SP_wb_WW_Exchange[ww].mean().max() for ww in range(0, len(list_WW))]))
plt.xlabel('WW')
plt.ylabel('Shadow price ($kr/(m^3/week)$)')
# plt.legend()
plt.title('Shadow price of water echange constraint, for each of the WW')
plt.show()


# optimal decision + SP on water exchange constraint, from one given WW to another given one

ww1 = 106294
ww2 = 2065
ww2_index = list_WW.index(ww2)

fig, ax1 = plt.subplots()

color = 'tab:blue'
ax1.set_xlabel('Time')
ax1.set_ylabel('Water transfer (1000m3/week)', color=color)
ax1.plot(TIME, optimal_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)], color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'tab:red'
ax2.set_ylabel('Shadow prices ($kr/(m^3/week)$)', color=color)  # we already handled the x-label with ax1
ax2.plot(TIME, SP_wb_WW_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title('Optimal water exchange and SP, from WW '+str(ww1)+' to '+str(ww2))
plt.show()

# optimal decision + SP on water exchange constraint, from Sjælsø to Tinghoj

ww1 = 106294
ww2 = 1
ww2_index = list_WW.index(ww2)

fig, ax1 = plt.subplots()

color = 'tab:blue'
ax1.set_xlabel('Time')
ax1.set_ylabel('Water transfer (1000m3/week)', color=color)
ax1.plot(TIME, optimal_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)], color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'tab:red'
ax2.set_ylabel('Transfer capacity shadow price ($DKK/(m^3/week)$)', color=color)  # we already handled the x-label with ax1
ax2.plot(TIME, SP_wb_WW_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title('Transfer from Sjælsø waterworks to Tinghøj storage')
plt.show()

#Tinghoj storage and SP
ww = 1

plt.figure()
plt.bar(TIME,optimal_Storage_WW['Storage_WW_'+str(ww)])
plt.xlabel('Time (weeks)')
plt.ylabel('Storage ($(1000m^3)$)')
# plt.legend()
plt.title('Storage in Tinghøj')
plt.show()

fig, ax1 = plt.subplots()

color = 'tab:blue'
ax1.set_xlabel('Time')
ax1.set_ylabel('Water Storage (1000m3)', color=color)
ax1.plot(TIME,optimal_Storage_WW['Storage_WW_'+str(ww)], color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'tab:red'
ax2.set_ylabel('Shadow prices ($kr/(m^3/week)$)', color=color)  # we already handled the x-label with ax1
ax2.plot(TIME, SP_wb_WW_Storage['SP_wb_WW_Storage_'+str(ww)], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title('Optimal water storage and SP,in Tinghøj')
plt.show()

# optimal decision + SP on water exchange constraint, from Sondersø to Ballerup

ww1 = 2065
ww2 = 106226
ww2_index =  list_WW.index(ww2)

fig, ax1 = plt.subplots()

color = 'tab:blue'
ax1.set_xlabel('Time')
ax1.set_ylabel('Water transfer (1000m3/week)', color=color)
ax1.plot(TIME, optimal_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)], color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'tab:red'
ax2.set_ylabel('Transfer capacity shadow price ($DKK/(m^3/week)$)', color=color)  # we already handled the x-label with ax1
ax2.plot(TIME, SP_wb_WW_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title('Transfer from Søndersø waterworks to Ballerup waterworks')
plt.show()

# optimal decision + SP on water exchange constraint, from  Roskilde to Lejre

ww1 = 104838
ww2 = 28315
ww2_index = list_WW.index(ww2)

fig, ax1 = plt.subplots()

color = 'tab:blue'
ax1.set_xlabel('Time')
ax1.set_ylabel('Water transfer (1000m3/week)', color=color)
ax1.plot(TIME, optimal_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)], color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'tab:red'
ax2.set_ylabel('Shadow prices ($kr/(m^3/week)$)', color=color)  # we already handled the x-label with ax1
ax2.plot(TIME, SP_wb_WW_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title('Optimal water exchange and SP, from Roskilde to Lejre')
plt.show()

# optimal decision + SP on water exchange constraint, from Slangerup to Hillerød

ww1 = 4230
ww2 = 83381
ww2_index = list_WW.index(ww2)

fig, ax1 = plt.subplots()

color = 'tab:blue'
ax1.set_xlabel('Time')
ax1.set_ylabel('Water transfer (1000m3/week)', color=color)
ax1.plot(TIME, optimal_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)], color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'tab:red'
ax2.set_ylabel('Shadow prices ($kr/(m^3/week)$)', color=color)  # we already handled the x-label with ax1
ax2.plot(TIME, SP_wb_WW_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title('Optimal water exchange and SP, fro Slangerup to Hillerød')
plt.show()

# optimal decision + SP on water exchange constraint, from Sjaelso to Bagsvaerd

ww1 = 106294
ww2 = 106311
ww2_index = list_WW.index(ww2)

fig, ax1 = plt.subplots()

color = 'tab:blue'
ax1.set_xlabel('Time')
ax1.set_ylabel('Water transfer (1000m3/week)', color=color)
ax1.plot(TIME, optimal_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)], color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'tab:red'
ax2.set_ylabel('Shadow prices ($kr/(m^3/week)$)', color=color)  # we already handled the x-label with ax1
ax2.plot(TIME, SP_wb_WW_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title('Optimal water exchange and SP, fro Sjælsø to Gladsaxe')
plt.show()

# optimal decision + SP on water exchange constraint, from Sjaelso to Ermelund (Gentofte)

ww1 = 106294
ww2 = 106292
ww2_index = list_WW.index(ww2)

fig, ax1 = plt.subplots()

color = 'tab:blue'
ax1.set_xlabel('Time')
ax1.set_ylabel('Water transfer (1000m3/week)', color=color)
ax1.plot(TIME, optimal_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)], color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'tab:red'
ax2.set_ylabel('Shadow prices ($kr/(m^3/week)$)', color=color)  # we already handled the x-label with ax1
ax2.plot(TIME, SP_wb_WW_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title('Optimal water exchange and SP, from Sjælsø to Ermelund (Gentofte)')
plt.show()

#%% Climatology figure (data per month, averaged on all the year)

months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# =============================================================================
# Tinghøj reservoir
# =============================================================================

ww = 1 # Tinghøj
avg_per_month = []

for m in range(1,13):
    list_weeks = [week[0] for week in week_in_month.items() if week[1] == m]
    avg=0
    for week in list_weeks:
        avg += optimal_Decision[optimal_Decision['time']==week]['Storage_WW_'+str(ww)].values[0]
    avg = avg/len(list_weeks)
    avg_per_month.append(avg)

plt.figure()
plt.bar(months, avg_per_month)
plt.xlabel('months')
plt.ylabel('Storage ($1000m3$)')
# plt.legend()
plt.title('Average monthly storage in Tinghøj reservoir')
plt.show()

ww = 1 # Tinghøj
avg_per_month = []

for m in range(1,13):
    list_weeks = [week[0] for week in week_in_month.items() if week[1] == m]
    avg=0
    for week in list_weeks:
        avg += SP[SP['time']==week]['SP_wb_WW_Storage_'+str(ww)].values[0]
    avg = avg/len(list_weeks)
    avg_per_month.append(avg)


plt.figure()
plt.bar(months, avg_per_month, color ='tab:red')
plt.xlabel('months')
plt.ylabel('Shadow prices ($kr/(m^3/week)$)')
# plt.legend()
plt.title('Average monthly SP for Tinghøj reservoir capacity constraint')
plt.show()

# =============================================================================
# Exchange Sjaelso Tinghoj
# =============================================================================

ww1 = 106294
ww2 = 1
ww2_index = list_WW.index(ww2)

avg_per_month = []
for m in range(1,13):
    list_weeks = [week[0] for week in week_in_month.items() if week[1] == m]
    avg=0
    for week in list_weeks:
        avg += optimal_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)][week]
    avg = avg/len(list_weeks)
    avg_per_month.append(avg)

plt.figure()
plt.bar(months, avg_per_month)
plt.xlabel('months')
plt.ylabel('Water exchange ($1000m3/week$)')
# plt.legend()
plt.title('Average monthly exchange from Sjælsø to Tinghøj')
plt.show()

avg_per_month = []
for m in range(1,13):
    list_weeks = [week[0] for week in week_in_month.items() if week[1] == m]
    avg=0
    for week in list_weeks:
        avg += SP_wb_WW_Exchange[ww2_index]['Exchange_'+str(ww1)+'_to_'+str(ww2)][week]
    avg = avg/len(list_weeks)
    avg_per_month.append(avg)


plt.figure()
plt.bar(months, avg_per_month, color ='tab:red')
plt.xlabel('months')
plt.ylabel('Shadow prices ($kr/(m^3/week)$)')
# plt.legend()
plt.title('Average monthly SP, water exchange capacity constraint, Sjælsø to Tinghøj')
plt.show()

# =============================================================================
# Water demand fulfilment total
# =============================================================================

avg_per_month = []
Total_D = WSA['Wateruse households (1000m3)'].sum()/52.18 + WSA['Wateruse industries (1000m3)'].sum()/52.18 + WSA['Wateruse services (1000m3)'].sum()/52.18+0.000000001
for m in range(1,13):
    list_weeks = [week[0] for week in week_in_month.items() if week[1] == m]
    avg=0
    for week in list_weeks:
        Total_A = sum(optimal_Decision[optimal_Decision['time']==week]['A_HH_'+str(int(wsa))].values[0] + optimal_Decision[optimal_Decision['time']==week]['A_Ind_'+str(int(wsa))].values[0] + optimal_Decision[optimal_Decision['time']==week]['A_PS_'+str(int(wsa))].values[0] for wsa in list_WSA)
        avg += Total_A/Total_D*100  # % of the demand
    avg = avg/len(list_weeks)
    avg_per_month.append(avg)

plt.figure()
plt.bar(months, avg_per_month)
plt.xlabel('months')
plt.ylabel('Water demand fulfilment (% of the demand)')
# plt.legend()
plt.ylim((50,100))
plt.title('Average monthly water demand fulfilment (all categories)')
plt.show()


# =============================================================================
# Water demand fulfilment total (1 value per year)
# =============================================================================

Years = [1988+x for x in list_year]

avg_per_year = []
Total_D = WSA['Wateruse households (1000m3)'].sum() + WSA['Wateruse industries (1000m3)'].sum() + WSA['Wateruse services (1000m3)'].sum()+0.000000001
for y in list_year:
    list_weeks = [week[0] for week in week_in_year.items() if week[1] == y]
    avg=0
    Total_A = sum(optimal_Decision[optimal_Decision['time']==week]['A_HH_'+str(int(wsa))].values[0] + optimal_Decision[optimal_Decision['time']==week]['A_Ind_'+str(int(wsa))].values[0] + optimal_Decision[optimal_Decision['time']==week]['A_PS_'+str(int(wsa))].values[0] for wsa in list_WSA for week in list_weeks)
    avg += Total_A/Total_D*100  # % of the demand
    avg_per_year.append(avg)

plt.figure()
plt.bar(Years, avg_per_year)
plt.xlabel('Year')
plt.ylabel('Water demand fulfilment (% of the demand)')
# plt.legend()
plt.ylim((50,100))
plt.title('Average yearly water demand fulfilment (all categories)')
plt.show()
# =============================================================================
# Compensation pumping + SP of min baseflow constraint  (maybe not useful because it's yearly constraint)
# =============================================================================

#%% Maps 

# to make PoygonPatch work I run
# pip install pyshp
# pip install shapely==1.8.5 # for PolygonPatch
# pip install descartes

import shapefile
from shapely.geometry import Polygon
from descartes.patch import PolygonPatch
import matplotlib.colors as mcolors


# =============================================================================
# Demand fulfillment
# =============================================================================

sf = shapefile.Reader(path_Model_folder + r'/Shapefiles/Catchment_Watersheds_Capital.shp')

# Households

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    WSA_in_catch = WSA[WSA['Catch_ID'] == c]['WSAID'].values
    
    if len(WSA_in_catch) == 0:
        fc = [0, 1, 0] # all green!
    
    else:
        Total_A = sum(optimal_A_HH['A_HH_'+str(int(wsa))].sum() for wsa in WSA_in_catch)
        Total_D = sum(float(WSA['Wateruse households (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_HH)/52.18 for wsa in WSA_in_catch)+0.000000001
        fc = [(1-Total_A/Total_D), Total_A/Total_D, 0]
        if sum(float(WSA['Wateruse households (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_HH)/52.18 for wsa in WSA_in_catch) == 0:
            fc = [0, 1, 0] # all green!

    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_red", [[1,0,0], [0,1,0]])

# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=0, vmax=1)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Water Demand fulfillment Households')

#Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(0, 1, num=5), ax = ax)
cbar.ax.set_yticklabels([f'{x:.0f}' for x in np.linspace(0, 100, num=5)])
cbar.set_label('Demand fulfillment (%)')

plt.show()


# Industry

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    WSA_in_catch = WSA[WSA['Catch_ID'] == c]['WSAID'].values
    
    if len(WSA_in_catch) == 0:
        fc = [0, 1, 0] # all green!
    else:
        Total_A = sum(optimal_A_Ind['A_Ind_'+str(int(wsa))].sum() for wsa in WSA_in_catch)
        Total_D = sum(float(WSA['Wateruse industries (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_Ind)/52.18 for wsa in WSA_in_catch)+0.000000001
        fc = [(1-Total_A/Total_D), Total_A/Total_D, 0]
        if sum(float(WSA['Wateruse industries (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_Ind)/52.18 for wsa in WSA_in_catch) == 0:
            fc = [0, 1, 0] # all green!
    
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_red", [[1,0,0], [0,1,0]])

# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=0, vmax=1)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(0, 1, num=5),ax = ax)
cbar.ax.set_yticklabels([f'{x:.0f}' for x in np.linspace(0, 100, num=5)])
cbar.set_label('Demand fulfillment (%)')


plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Water Demand fulfillment Industries')
plt.show()

# Public services

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    WSA_in_catch = WSA[WSA['Catch_ID'] == c]['WSAID'].values
    
    if len(WSA_in_catch) == 0:
        fc = [0, 1, 0] # all green!
    else:
        Total_A = sum(optimal_A_PS['A_PS_'+str(int(wsa))].sum() for wsa in WSA_in_catch)
        Total_D = sum(float(WSA['Wateruse services (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_PS)/52.18 for wsa in WSA_in_catch)+0.000000001
        fc = [(1-Total_A/Total_D), Total_A/Total_D, 0]
        if sum(float(WSA['Wateruse services (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_PS)/52.18 for wsa in WSA_in_catch) == 0:
            fc = [0, 1, 0] # all green!
            
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_red", [[1,0,0], [0,1,0]])

# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=0, vmax=1)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(0, 1, num=5),ax = ax)
cbar.ax.set_yticklabels([f'{x:.0f}' for x in np.linspace(0, 100, num=5)])
cbar.set_label('Demand fulfillment (%)')


plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Water Demand fulfillment Public Services')
plt.show()


# Total

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    WSA_in_catch = WSA[WSA['Catch_ID'] == c]['WSAID'].values
    
    if len(WSA_in_catch) == 0:
        fc = [0, 1, 0] # all green!
    else:
        Total_A = sum(optimal_A_HH['A_HH_'+str(int(wsa))].sum() + optimal_A_Ind['A_Ind_'+str(int(wsa))].sum() + optimal_A_PS['A_PS_'+str(int(wsa))].sum() for wsa in WSA_in_catch)
        Total_D = sum(float(WSA['Wateruse households (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_HH)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse industries (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_Ind)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse services (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_PS)/52.18 for wsa in WSA_in_catch)+0.000000001
        fc = [1-Total_A/Total_D, Total_A/Total_D, 0]
        if sum(float(WSA['Wateruse households (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_HH)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse industries (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_Ind)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse services (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_PS)/52.18 for wsa in WSA_in_catch) == 0:
            fc = [0, 1, 0] # all green!        
    
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_red", [[1,0,0], [0,1,0]])

# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=0, vmax=1)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(0, 1, num=5),ax = ax)
cbar.ax.set_yticklabels([f'{x:.0f}' for x in np.linspace(0, 100, num=5)])
cbar.set_label('Demand fulfillment (%)')


plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Water Demand fulfillment Total')
plt.show()

# =============================================================================
# Pumping WF
# =============================================================================

pumping_WF_catch = {}
for c in list_catch:
    WF_in_catch = WF.drop_duplicates(subset='WFID')[WF.drop_duplicates(subset='WFID')['Catch_ID'] == c]['WFID'].values
    pumping_WF_catch[c] = sum(optimal_Pump_WF['Pump_WF_'+str(wf)].mean() for wf in WF_in_catch)

max_pumping_WF_catch = max(pumping_WF_catch.values())
int_max_pumping_WF_catch = 100*(int(max_pumping_WF_catch/100)+1)

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    fc = [1, 1-pumping_WF_catch[c]/int_max_pumping_WF_catch, 1] # gw ind 2
    
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_red", [[1,1,1], [1,0,1]])

# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=0, vmax=int_max_pumping_WF_catch)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(0, int_max_pumping_WF_catch, num=5))
cbar.ax.set_yticklabels([f'{x:.0f}' for x in np.linspace(0, int_max_pumping_WF_catch, num=5)])
cbar.set_label('Pumping (1000m3/week)')

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Weekly average pumping of WF')
plt.show()

# =============================================================================
# Compensating pumping
# =============================================================================

compensating_pumping = {}
for c in list_catch:
    compensating_pumping[c] = optimal_Pump_GW_to_BF['Pump_GW_to_BF_'+str(c)].mean()

max_compensating = max(compensating_pumping.values())
int_max_compensating = 10*(int(max_compensating/10)+1)

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    fc = [1, 1-compensating_pumping[c]/int_max_compensating, 1] # gw ind 2
    
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_red", [[1,1,1], [1,0,1]])

# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=0, vmax=int_max_compensating)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(0, int_max_compensating, num=5))
cbar.ax.set_yticklabels([f'{x:.1f}' for x in np.linspace(0, int_max_compensating, num=5)])
cbar.set_label('Pumping (1000m3/week)')

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Weekly average Compensating pumping')
plt.show()



#%% Water demand per category

water_demand_HH = {}
for c in list_catch:
    WSA_in_catch = WSA[WSA['Catch_ID'] == c]['WSAID'].values
    
    if len(WSA_in_catch) == 0:
        water_demand_HH[c] = 0 
    else:
        Total_D = sum(float(WSA['Wateruse households (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_HH)/52.18 for wsa in WSA_in_catch) + 0.000000001
        water_demand_HH[c] = Total_D/(len(TIME)/52.18)

water_demand_Ind = {}
for c in list_catch:
    WSA_in_catch = WSA[WSA['Catch_ID'] == c]['WSAID'].values
    
    if len(WSA_in_catch) == 0:
        water_demand_Ind[c] = 0 
    else:
        Total_D = sum(float(WSA['Wateruse industries (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_HH)/52.18 for wsa in WSA_in_catch) + 0.000000001
        water_demand_Ind[c] = Total_D/(len(TIME)/52.18)

water_demand_PS = {}
for c in list_catch:
    WSA_in_catch = WSA[WSA['Catch_ID'] == c]['WSAID'].values
    
    if len(WSA_in_catch) == 0:
        water_demand_PS[c] = 0 
    else:
        Total_D = sum(float(WSA['Wateruse services (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_HH)/52.18 for wsa in WSA_in_catch) + 0.000000001
        water_demand_PS[c] = Total_D/(len(TIME)/52.18)

# Water demand per cachment HH

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    WSA_in_catch = WSA[WSA['Catch_ID'] == c]['WSAID'].values
    
    if len(WSA_in_catch) == 0:
        fc = [1, 1, 1] # all white, no demand!
    else:
        Yearly_D = water_demand_HH[c]
        fc = [1-Yearly_D/max(water_demand_HH.values()),1-Yearly_D/max(water_demand_HH.values()), 1]
        # fc = [0, 0, 1/(1+exp(-Total_D/max(water_demand_catch.values())))]        
        if sum(float(WSA['Wateruse households (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_HH)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse industries (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_Ind)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse services (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_PS)/52.18 for wsa in WSA_in_catch) == 0:
            fc = [1, 1, 1] # all white!        
    
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Yearly water demand HH (1000m3)')
plt.show()


# Water demand per cachment Ind

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    WSA_in_catch = WSA[WSA['Catch_ID'] == c]['WSAID'].values
    
    if len(WSA_in_catch) == 0:
        fc = [1, 1, 1] # all white, no demand!
    else:
        Yearly_D = water_demand_Ind[c]
        fc = [1-Yearly_D/max(water_demand_Ind.values()),1-Yearly_D/max(water_demand_Ind.values()), 1]
        # fc = [0, 0, 1/(1+exp(-Total_D/max(water_demand_catch.values())))]        
        if sum(float(WSA['Wateruse households (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_HH)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse industries (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_Ind)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse services (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_PS)/52.18 for wsa in WSA_in_catch) == 0:
            fc = [1, 1, 1] # all white!        
    
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Yearly water demand Ind (1000m3)')
plt.show()


# Water demand per cachment HH

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    WSA_in_catch = WSA[WSA['Catch_ID'] == c]['WSAID'].values
    
    if len(WSA_in_catch) == 0:
        fc = [1, 1, 1] # all white, no demand!
    else:
        Yearly_D = water_demand_PS[c]
        fc = [1-Yearly_D/max(water_demand_PS.values()),1-Yearly_D/max(water_demand_PS.values()), 1]
        # fc = [0, 0, 1/(1+exp(-Total_D/max(water_demand_catch.values())))]        
        if sum(float(WSA['Wateruse households (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_HH)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse industries (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_Ind)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse services (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_PS)/52.18 for wsa in WSA_in_catch) == 0:
            fc = [1, 1, 1] # all white!        
    
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Yearly water demand PS (1000m3)')
plt.show()


#%% Data show water demand and abstraction capacity per cactchment

water_demand_catch = {}
for c in list_catch:
    WSA_in_catch = WSA[WSA['Catch_ID'] == c]['WSAID'].values
    
    if len(WSA_in_catch) == 0:
        water_demand_catch[c] = 0 
    else:
        Total_D = sum(float(WSA['Wateruse households (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_HH)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse industries (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_Ind)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse services (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_PS)/52.18 for wsa in WSA_in_catch) + 0.000000001
        water_demand_catch[c] = Total_D/(len(TIME)/52.18)/area[c]*1000000

average_inflow_catch = {}
for c in list_catch:
    average_inflow_catch[c] = sum(inflow[c-1,t-1] for t in range(1,weeks+1))/(len(TIME)/52.18)/area[c]*1000000

abstraction_license_catch = {}
for c in list_catch:
    WF_in_catch = WF.drop_duplicates(subset='WFID')[WF.drop_duplicates(subset='WFID')['Catch_ID'] == c]['WFID'].values
    
    if len(WF_in_catch) == 0:
        abstraction_license_catch[c] = 0 
    else:
        abstraction_license_catch[c] = sum(maxpump[wf] for wf in WF_in_catch)/area[c]*1000000
        
water_scarcity_catch = {}
for c in list_catch:
    water_scarcity_catch[c] = water_demand_catch[c]/average_inflow_catch[c]
    
demand_abstraction_licence = {}    
for c in list_catch:
    demand_abstraction_licence[c] = water_demand_catch[c]/abstraction_license_catch[c]

max_value = max(max(water_demand_catch.values()),max(average_inflow_catch.values()),max(abstraction_license_catch.values()))
max_demand_abstraction_ratio = max(demand_abstraction_licence.values())

# Water demand per cachment

int_max_water_demand = 100*(int(max(water_demand_catch.values())/100)+1)

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    WSA_in_catch = WSA[WSA['Catch_ID'] == c]['WSAID'].values
    
    if len(WSA_in_catch) == 0:
        fc = [1, 1, 1] # all white, no demand!
    else:
        fc = [1-water_demand_catch[c]/int_max_water_demand,1-water_demand_catch[c]/int_max_water_demand, 1]
        # fc = [0, 0, 1/(1+exp(-Total_D/max(water_demand_catch.values())))]        
        if sum(float(WSA['Wateruse households (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_HH)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse industries (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_Ind)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse services (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_PS)/52.18 for wsa in WSA_in_catch) == 0:
            fc = [1, 1, 1] # all white!        
    
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_blue", ["white", "blue"])

# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=0, vmax=int_max_water_demand)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(0, int_max_water_demand, num=5),ax=ax)
cbar.ax.set_yticklabels([f'{x:.0f}' for x in np.linspace(0, int_max_water_demand, num=5)])
cbar.set_label('Water demand (mm)')

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Average yearly water demand')
plt.show()


# Average yearly inflow per cachment
    
plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    fc = [1-average_inflow_catch[c]/max_value, 1-average_inflow_catch[c]/max_value, 1]
        
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_blue", ["white", "blue"])

# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=0, vmax=100)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(0, 100, num=5),ax=ax)
cbar.ax.set_yticklabels([f'{x:.0f}' for x in np.linspace(0, 100, num=5)])
cbar.set_label('Average recharge (% of highest demand)')

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Average yearly groundwater recharge')
plt.show()


# Water demand/inflow

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    # fc = [1, 1-water_scarcity_catch[c]/max(water_scarcity_catch.values()), 1-water_scarcity_catch[c]/max(water_scarcity_catch.values())]
    
    
    if water_demand_catch[c] > average_inflow_catch[c]:
        fc=[1-(water_demand_catch[c]/average_inflow_catch[c]-1)/2*(1-0.4),0,0]
    elif water_demand_catch[c] > average_inflow_catch[c]/100:
        # fc = [1, 1-water_scarcity_catch[c]/max(water_scarcity_catch.values()), 1-water_scarcity_catch[c]/max(water_scarcity_catch.values())]
        fc = [1, 1-water_demand_catch[c]/average_inflow_catch[c], 1-water_demand_catch[c]/average_inflow_catch[c]]
    else:
        fc=[1,1,1]
        
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_red", [[1,1,1], [1,0,0], [0.7,0,0], [0.4,0,0]])
 
# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=0, vmax=3)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(0, 3, num=7),ax=ax)
cbar.ax.set_yticklabels([f'{x:.1f}' for x in np.linspace(0, 3, num=7)])
cbar.set_label('Demand/Recharge')

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Demand to availability ratio')
plt.show()



# Abstraction license per cachment

int_max_abstraction_licence = 100*(int(max(abstraction_license_catch.values())/100)+1)

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    WF_in_catch = WF.drop_duplicates(subset='WFID')[WF.drop_duplicates(subset='WFID')['Catch_ID'] == c]['WFID'].values
    
    if len(WF_in_catch) == 0:
        fc = [1, 1, 1] # all white, no demand!
    else:
        fc = [1-abstraction_license_catch[c]/int_max_abstraction_licence,1-abstraction_license_catch[c]/int_max_abstraction_licence, 1]
        # fc = [1-abstraction_license_catch[c]/max(abstraction_license_catch.values()),1-abstraction_license_catch[c]/max(abstraction_license_catch.values()), 1]
        if sum(float(WSA['Wateruse households (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_HH)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse industries (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_Ind)/52.18 for wsa in WSA_in_catch) + sum(float(WSA['Wateruse services (1000m3)'][WSA['WSAID']==wsa].iloc[0])*len(optimal_A_PS)/52.18 for wsa in WSA_in_catch) == 0:
            fc = [1, 1, 1] # all white!        
    
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_blue", ["white", "blue"])

# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=0, vmax=int_max_abstraction_licence)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(0, int_max_abstraction_licence, num=5),ax=ax)
cbar.ax.set_yticklabels([f'{x:.0f}' for x in np.linspace(0, int_max_abstraction_licence, num=5)])
cbar.set_label('Abstraction licence (mm)')

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Yearly abstraction licence')
plt.show()


# Water demand/abstraction licence

int_max_demand_abstraction_ratio = 10*(int(max_demand_abstraction_ratio/10+1))
c_over_1 = 0

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    # fc = [1, 1-water_scarcity_catch[c]/max(water_scarcity_catch.values()), 1-water_scarcity_catch[c]/max(water_scarcity_catch.values())]
    if water_demand_catch[c]/abstraction_license_catch[c] <=0.1:
        fc =[1,1,1]
    
    if 0.1 < water_demand_catch[c]/abstraction_license_catch[c] <=1:
        fc = [1, 1+log10(water_demand_catch[c]/abstraction_license_catch[c]), 1+log10(water_demand_catch[c]/abstraction_license_catch[c])]
        
    if 1 <  water_demand_catch[c]/abstraction_license_catch[c]:
        fc=[1-log10(water_demand_catch[c]/abstraction_license_catch[c])/2,0,0]
        c_over_1 += 1
    
    # fc = fc = [1, 1-demand_abstraction_licence[c]/int_max_demand_abstraction_ratio, 1-demand_abstraction_licence[c]/int_max_demand_abstraction_ratio]
        
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_red", [[1,1,1], [1,0,0], [0.5,0,0], [0,0,0]])
# cmap = mcolors.LinearSegmentedColormap.from_list("white_red", [[1,1,1], [1,0,0]])
 
# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=-1, vmax=2)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(-1, 2, num=4),ax=ax)
cbar.ax.set_yticklabels([f'{10**(x):.1f}' for x in np.linspace(-1, 2, num=4)])
cbar.set_label('Demand / Abstraction licence')

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Demand to abstraction licence ratio')
plt.show()

# =============================================================================
# Over exploitation
# =============================================================================

sf = shapefile.Reader(path_Model_folder + r'/Shapefiles/Catchment_Watersheds_Capital.shp')

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

sum_inflow = [sum([inflow[c-1,t-1] for t in range(1,weeks+1)])+0.001 for c in list_catch] # +0.001 to avoid errors due to calculations

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
        
    sum_pump = optimal_Pump_catch['Pump_catch_'+str(c)].sum() +0.000001   # avoid negative value but still smaller than 0.001
    fc = [1, (1-sum_pump/(sum_inflow[c-1]/2)), (1-sum_pump/(sum_inflow[c-1]/2))] # gw ind 2
    
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_red", ["white", "red"])

# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=0, vmax=0.5)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(0, 0.5, num=5),ax=ax)
cbar.ax.set_yticklabels([f'{x:.2f}' for x in np.linspace(0, 0.5, num=5)])
cbar.set_label('GW quality indicator 2')

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Over-exploitation of groundwater reservoirs')
plt.show()



#%% Maps with SP

lin_res_catch = {}
for c in list_catch:
    lin_res_catch[c] = SP_lin_res['SP_lin_res_'+str(int(c))].mean() 

min_bf_catch = {}
for c in list_catch:
    min_bf_catch[c] = SP_min_bf['SP_min_bf_'+str(int(c))].mean() 
    
gw_ind_catch = {}
for c in list_catch:
    gw_ind_catch[c] = SP_gw_ind_2['SP_gw_ind_2_'+str(int(c))].mean()     
    
pumping_WF_catch = {}
for c in list_catch:
    WF_in_catch = WF.drop_duplicates(subset='WFID')[WF.drop_duplicates(subset='WFID')['Catch_ID'] == c]['WFID'].values
    
    if len(WF_in_catch) == 0:
        pumping_WF_catch[c] = 0 
    else:
        pumping_WF_catch[c] = np.mean(np.array([SP_pumping_WF['SP_pumping_WF_'+str(int(wf))].mean() for wf in WF_in_catch]))
        

max_lin_res = max(lin_res_catch.values())
max_min_bf = min(min_bf_catch.values())
max_gw_ind = max(gw_ind_catch.values())
max_pumping_WF = max(pumping_WF_catch.values())


# =============================================================================
# lin res SP
# =============================================================================

int_max_lin_res = 10*(int(max_lin_res/10)+1)

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    fc = [1, 1, 1-lin_res_catch[c]/int_max_lin_res]
    
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_blue", [[1,1,1], [1,1,0]])

# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=0, vmax=int_max_lin_res)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(0, int_max_lin_res, num=5),ax=ax)
cbar.ax.set_yticklabels([f'{x:.1f}' for x in np.linspace(0, int_max_lin_res, num=5)])
cbar.set_label('SP (kr/(m3/week))')

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Linear reservoir constraint SP')
plt.show()


# =============================================================================
# min bf SP
# =============================================================================

int_max_min_bf = 10*(int(max_min_bf/10)-1)

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    fc = [1, 1, 1-min_bf_catch[c]/int_max_min_bf]
    
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_blue", [[1,1,0], [1,1,1]])

# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=int_max_min_bf, vmax=0)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(int_max_min_bf, 0, num=5),ax=ax)
cbar.ax.set_yticklabels([f'{x:.0f}' for x in np.linspace(int_max_min_bf, 0, num=5)])
cbar.set_label('SP (kr/(m3/week))')

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Minimum Baseflow constraint SP')
plt.show()


# =============================================================================
# GW indicator 2
# =============================================================================

int_max_gw_ind = 10*(int(max_gw_ind/10)+1)

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    fc = [1, 1, 1-gw_ind_catch[c]/(int_max_gw_ind+0.0000000001)]
    
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_blue", [[1,1,1], [1,1,0]])

# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=0, vmax=int_max_gw_ind)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(0, int_max_gw_ind, num=5),ax=ax)
cbar.ax.set_yticklabels([f'{x:.1f}' for x in np.linspace(0, int_max_gw_ind, num=5)])
cbar.set_label('SP (kr/(m3/week))')

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('GW indicator 2 constraint SP')
plt.show()


# =============================================================================
# pumping WF
# =============================================================================

int_max_pumping_WF = 10*(int(max_pumping_WF/10)+1)

plt.figure()
ax = plt.axes()
ax.set_aspect('equal')

shape_ex = sf.shape(0) # to get intial values for gobal limits
global_limits = shape_ex.bbox

for c in list_catch:
    shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
    
    fc = [1, 1, 1-pumping_WF_catch[c]/int_max_pumping_WF]
    
    polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
    patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
    ax.add_patch(patch)
    limits = shape_ex.bbox
    global_limits[0]= min(limits[0],global_limits[0])
    global_limits[1]= min(limits[1],global_limits[1])
    global_limits[2]= max(limits[2],global_limits[2])
    global_limits[3]= max(limits[3],global_limits[3])

# Create a color map from white to red
cmap = mcolors.LinearSegmentedColormap.from_list("white_blue", [[1,1,1], [1,1,0]])

# Normalize the values between 0 and 1 for the color map
norm = mcolors.Normalize(vmin=0, vmax=int_max_pumping_WF)

# Create a scalar mappable for the color bar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the color bar to the plot with intermediary values
cbar = plt.colorbar(sm, ticks=np.linspace(0, int_max_pumping_WF, num=5),ax=ax)
cbar.ax.set_yticklabels([f'{x:.1f}' for x in np.linspace(0, int_max_pumping_WF, num=5)])
cbar.set_label('SP (DKK/(m3/week))')

plt.xticks([])
plt.yticks([])
plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
plt.title('Average shadow price per catchment')
plt.show()


#%% Pots for scenario analysis

plt.figure()
plt.plot(min_bf_catch.values(),compensating_pumping.values(),'ro')
plt.xlabel('MinBF shadow prices (kr/(1000m3/week))')
plt.ylabel('Comensating pumping (1000m3/week)')
plt.show()

# =============================================================================
# Scenario 1 - GW indicator 2
# =============================================================================


# Water allocation in cathment 10

WSA_in_catch_10 = WSA[WSA['Catch_ID'] == 10]['WSAID'].values

plt.figure()
for wsa in [2851,2860,3341,3359,4676]:
    plt.plot(TIME, optimal_Decision['A_HH_'+str(wsa)], label=str(wsa))
plt.xlabel('Time')
plt.ylabel('Allocation HH (1000 m3 / week)')
plt.legend()
plt.title('Allocation of water to households')
plt.show()

plt.figure()
total_alloc = 0
total_demand = 0
for wsa in WSA_in_catch_10:
    total_alloc += optimal_Decision['A_HH_'+str(wsa)]+optimal_Decision['A_Ind_'+str(wsa)]+optimal_Decision['A_PS_'+str(wsa)]
    total_demand += float(WSA['Wateruse households (1000m3)'][WSA['WSAID']==wsa].iloc[0])/52.18 + float(WSA['Wateruse industries (1000m3)'][WSA['WSAID']==wsa].iloc[0])/52.18 + float(WSA['Wateruse services (1000m3)'][WSA['WSAID']==wsa].iloc[0])/52.18 + float(WSA['Wateruse agriculture (1000m3)'][WSA['WSAID']==wsa].iloc[0])/52.18
plt.plot(TIME, total_alloc, label='total_allocation')
plt.plot(TIME, total_demand*np.ones(len(optimal_Decision)), label='total_demand')
plt.xlabel('Time')
plt.ylabel('Allocation (1000 m3 / week)')
plt.legend()
plt.title('Total allocation of water in catchment 10')
plt.show()



#%% Find location of a Catchment

# =============================================================================
# Catchment
# =============================================================================

# =============================================================================
# for catch_to_show in range(1,33):
# 
#     sf = shapefile.Reader(path_Model_folder + r'/Shapefiles/Catchment_Watersheds_Capital.shp')
#     
#     plt.figure()
#     ax = plt.axes()
#     ax.set_aspect('equal')
#     
#     shape_ex = sf.shape(0) # to get intial values for gobal limits
#     global_limits = shape_ex.bbox
#     
#     for c in list_catch:
#         shape_ex = sf.shape(c-1) # could break if selected shape has multiple polygons. 
#         
#         fc = [0.9, 1,1]
#         if c == catch_to_show:
#             fc = [1,0,0]
#         
#         polygon = Polygon(shape_ex.points)  # build the polygon from exterior points
#         patch = PolygonPatch(polygon, facecolor=fc, edgecolor=[0,0,0], alpha=0.7, zorder=2)  # [0,0,0] is black and [1,1,1] is white (RGB)
#         ax.add_patch(patch)
#         limits = shape_ex.bbox
#         global_limits[0]= min(limits[0],global_limits[0])
#         global_limits[1]= min(limits[1],global_limits[1])
#         global_limits[2]= max(limits[2],global_limits[2])
#         global_limits[3]= max(limits[3],global_limits[3])
#     
#     plt.xticks([])
#     plt.yticks([])
#     plt.xlim(global_limits[0]-1000,global_limits[2]+1000)
#     plt.ylim(global_limits[1]-1000,global_limits[3]+1000)
#     plt.title('Catchment '+str(catch_to_show))
#     plt.show()
# =============================================================================



