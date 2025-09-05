# -*- coding: utf-8 -*-
"""
@author: Clément Franey
clefraney@gmail.com

Updated: June 04 2025


Instructions:

Set the path of the folder where the python model is located in the variable "folder_path".    
Run this script with the variable "process_WB" equal True at least once.

You can run the script again with processs_WB = False if you want to parametrize the groundwater 
linear reservoir with another error function (MSE, RMSE, NSE, r2).
By default the error function that is used is 'NSE'.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import scipy
import os
from math import sqrt
from sklearn.metrics import r2_score
import shutil


#%% Parameters


# =============================================================================
# Insert the path to the folder GW_allocation_model below
# =============================================================================

folder_path = r'F:\Data\s232484\GW_allocation_model' # CHANGE TO YOUR PATH

# =============================================================================
# Set global parameters
# =============================================================================

Catch='Watersheds' # 'ID15' or 'Watersheds' or 'Test'
Geo='Capital'  # 'Capital' or 'Western' or 'Sjaelland'
Catch_Geo = Catch+'_'+Geo
 
error_func='NSE' #MSE,RMSE,NSE,r2(coef of determination)
calibration=1300 #1300 = 25 years of calibration and 8 years of validation 
process_WB=True # process the WB data and export each to a CSV file, also get the min K value (just need to do it once)

# =============================================================================
# Path to the data
# =============================================================================

raw_data_path = folder_path + '\Raw data\WB_Watersheds_Capital'
input_data_path = folder_path + '\Input data model\WB_Watersheds_Capital'


name='WB_SZ_'

catch_ID=pd.read_csv(raw_data_path + os.sep + Catch_Geo+'_ID.csv')
ncatch=np.array(catch_ID['ID']) # get the list of all catchments!

# copy the file containing the catchments ID
shutil.copyfile(raw_data_path + os.sep + Catch_Geo+'_ID.csv', input_data_path + os.sep + Catch_Geo+'_ID.csv')

#%% Define a processing function

def processing(fname, fnameout,error_func,calibration,process_WB):
    
    #open the data
    data = pd.read_csv(fname,skiprows=4,sep='\t')
    
    #create time serie
    time = np.empty((data.shape[0],), dtype=datetime)
    count = 0
    for t in data[data.keys()[0]]:
        time[count] = datetime.strptime(t[0:11],' %Y %m %d').date()
        count=count+1
    time = time[0:-1]
    
    #extract data
    Inflow = np.array(-(data['qrech'] + data['qrechmp'] ))
    Pumping = np.array(data['qszabsex'])
    Baseflow = np.array(data['qszdrtorivin'] + data['qszrivpos'])
    dstor = np.array(data['dszsto'])
    
    # Process WB data and extract them in a CSV file
    if process_WB==True:
        Inflow_p = np.diff(Inflow)
        Pumping_p = np.diff(Pumping)
        Baseflow_p = np.diff(Baseflow)
        dstor_p = np.diff(dstor)
        
        data_out = pd.DataFrame({"Start time of weekly time step": time, "MIKE SHE GW recharge (mm)":Inflow_p, "MIKE SHE Baseflow (mm)":Baseflow_p, "MIKE SHE Storage change (mm)":dstor_p, "MIKE SHE Pumping (mm)":Pumping_p})
        data_out.to_csv(fnameout)
    
    # Run optimization of the linear reservoir (K parameter)
    Inflow_cal = np.diff(Inflow)[0:calibration]
    Pumping_cal = np.diff(Pumping)[0:calibration]
    Baseflow_cal = np.diff(Baseflow)[0:calibration]
    dstor_cal = np.diff(dstor)[0:calibration]
    
    Inflow_val = np.diff(Inflow)[calibration:]
    Pumping_val = np.diff(Pumping)[calibration:]
    Baseflow_val = np.diff(Baseflow)[calibration:]
    dstor_val = np.diff(dstor)[calibration:]
    
    #find the best K value for the catchment
    def MSE(K, scale, obs, Inflow_cal,Pumping_cal):
        sim, dSto =linres(200,7,K*7,Inflow_cal,Pumping_cal)
        l=len(sim)
        MSE=0
        for i in range(0,l):
            err=(obs[i]-sim[i])**2
            MSE=MSE+err
        return(MSE/l)
    def RMSE(K, scale, obs, Inflow_cal,Pumping_cal):
        sim, dSto =linres(200,7,K*7,Inflow_cal,Pumping_cal)
        l=len(sim)
        MSE=0
        for i in range(0,l):
            err=(obs[i]-sim[i])**2
            MSE=MSE+err
        return(sqrt(MSE/l))
    def NSE(K, scale, obs, Inflow_cal,Pumping_cal):
        sim, dSto =linres(200,7,K*7,Inflow_cal,Pumping_cal)
        l=len(sim)
        num=0
        denom=0
        mean_obs=obs.mean()
        for i in range(0,l):
            num=num+(obs[i]-sim[i])**2
            denom=denom+(mean_obs-sim[i])**2
        return(num/denom-1) #opposite of the NSE for the optimizer !!!
    def r2(K, scale, obs, Inflow_cal,Pumping_cal):           #Coef of determination
        sim, dSto =linres(200,7,K*7,Inflow_cal,Pumping_cal)
        l=len(sim)
        r2=r2_score(obs,sim)
        return(-r2) #opposite of the r2 for the optimizer
    
    # Perform optimization
    K0=5 
    scale=1
    if error_func=='MSE':
        res=scipy.optimize.minimize(fun=MSE,x0=K0,args=(scale, Baseflow_cal,Inflow_cal,Pumping_cal), method='Nelder-Mead')
        Best_K=float(res['x'][0])
        err_cal=MSE(Best_K, scale, Baseflow_cal, Inflow_cal,Pumping_cal)
        err_val=MSE(Best_K, scale, Baseflow_val, Inflow_val,Pumping_val)
    elif error_func=='RMSE':
        res=scipy.optimize.minimize(fun=RMSE,x0=K0,args=(scale, Baseflow_cal,Inflow_cal,Pumping_cal), method='Nelder-Mead')
        Best_K=float(res['x'][0])
        err_cal=RMSE(Best_K, scale, Baseflow_cal, Inflow_cal,Pumping_cal)
        err_val=RMSE(Best_K, scale, Baseflow_val, Inflow_val,Pumping_val)
    elif error_func=='NSE':
        res=scipy.optimize.minimize(fun=NSE,x0=K0,args=(scale, Baseflow_cal,Inflow_cal,Pumping_cal), method='Nelder-Mead')
        Best_K=float(res['x'][0])
        err_cal=-NSE(Best_K, scale, Baseflow_cal, Inflow_cal,Pumping_cal)
        err_val=-NSE(Best_K, scale, Baseflow_val, Inflow_val,Pumping_val)
    elif error_func=='r2':
        res=scipy.optimize.minimize(fun=r2,x0=K0,args=(scale, Baseflow_cal,Inflow_cal,Pumping_cal), method='Nelder-Mead')
        Best_K=float(res['x'][0])
        err_cal=-r2(Best_K, scale, Baseflow_cal, Inflow_cal,Pumping_cal)
        err_val=-r2(Best_K, scale, Baseflow_val, Inflow_val,Pumping_val)
    else:
        None 
    
    return (Best_K, err_cal,err_val)


def linres(Sini, deltat, K, Inflow, Pumping):
    nper = Pumping.shape[0]
    Qout = np.zeros(nper+1)
    Sto = np.zeros(nper+1)
    Qout[0] = Sini/K
    Sto[0] = Sini
    for i in range(nper):
        Qout[i+1] = Qout[i]*np.exp(-deltat/K) + (Inflow[i]-Pumping[i])*(1-np.exp(-deltat/K))
        Sto[i+1] = (Sto[i] + Inflow[i]-Pumping[i])*np.exp(-deltat/K)
    Qout = Qout[1:]
    dSto = np.diff(Sto)
    return Qout, dSto


#%% Process and save the data


K_values=[]
cal_values=[]
val_values=[]

for c in ncatch:
    fname = raw_data_path + os.sep + name + Catch + '_' + str(int(c)) + '.txt'
    fnameout = input_data_path + os.sep + name + Catch + '_' + str(int(c)) + '.csv'
    Best_K, err_cal, err_val =processing(fname, fnameout,error_func,calibration,process_WB)
    K_values.append(Best_K)
    cal_values.append(err_cal)
    val_values.append(err_val)


K_values=np.array(K_values)
cal_values=np.array(cal_values)
K_parameter=pd.DataFrame({"Catchment": ncatch,"K": K_values, "error_cal": cal_values, "error_val": val_values})
K_parameter.to_csv(input_data_path + os.sep + 'K_optim_'+error_func+'.csv')



#%% Find min K value so the linear reservoirs are not emptied during the 33 years period

if process_WB == True: # just need to be done once
    
    K_min = []
    
    for c in ncatch:
        #open the data
        data = pd.read_csv(raw_data_path + os.sep + name+Catch+'_'+str(int(c))+'.txt',skiprows=4,sep='\t')
        
        #extract data
        Inflow = np.array(-(data['qrech'] + data['qrechmp'] ))
        Pumping = np.array(data['qszabsex'])
        Baseflow = np.array(data['qszdrtorivin'] + data['qszrivpos'])
        dstor = np.array(data['dszsto'])
        
        Inflow_p = np.diff(Inflow)
        Pumping_p = np.diff(Pumping)
        Baseflow_p = np.diff(Baseflow)
        dstor_p = np.diff(dstor)
        
    
        Sini = 1000
        K=1
        sim, dSto = linres(Sini, 7, 7*K, Inflow_p, Pumping_p)
        while min(sim) < 0:
            K += 0.1
            sim, dSto = linres(Sini, 7, 7*K, Inflow_p, Pumping_p)
        K_min.append(K)
    
    
    K_min=pd.DataFrame({"Catchment": ncatch,"K": K_min})
    K_min.to_csv(input_data_path + os.sep + 'K_min.csv')

