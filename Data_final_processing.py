# -*- coding: utf-8 -*-
"""
@author: Clément Franey
clefraney@gmail.com

Updated: June 04 2025

Instructions:

Set the path of the folder where the python model is located in the variable "folder_path".    
Run this script once.
""" 

import numpy as np
import pandas as pd
import os
from datetime import datetime

# =============================================================================
# Insert the path to the folder GW_allocation_model below
# =============================================================================

folder_path = r'h:\Students\Clément\GW_allocation_model-main'  # CHANGE TO YOUR PATH

# =============================================================================
# Global parameters
# =============================================================================

os.chdir(folder_path + '\Raw data')
fnameout = folder_path + '\Input data model'


Catch='Watersheds' # 'ID15' or 'Watersheds' or 'Test'
Geo='Capital'  # 'Capital' or 'Western' or 'Sjaelland'

Catch_Geo = Catch + '_' + Geo

#%% TABLE CATCHMENTS

Table_Catchments=pd.read_excel(Catch_Geo+'_area.xlsx')
Table_Catchments.to_csv(fnameout+'\Table_Catchments_'+Catch_Geo+'.csv')


#%% TABLES Wellfields and Waterworks 

# =============================================================================
# Wellfields table
# =============================================================================

wells=pd.read_csv('Wells_Lars_Sjaelland.csv')
Wells_Catchments=pd.read_csv('Wells_'+Catch_Geo+'.csv')
Wells_Catchments=Wells_Catchments.dropna(subset=['Catch_ID'])
WSA_PLANT=pd.read_table('WSA_PLANT.txt', sep="	")
WSA_PLANT.rename(columns={'PLANTID':'AnlaegID'},inplace=True)
#Storage capacity data to find!!!

WF=Wells_Catchments.groupby(['AnlaegID'],as_index=False).first()
WF.drop(WF.columns[[1,2,3,4,5,6,7,8,9,10,11,12,13,49,50,52]], axis=1, inplace=True)


# If there is no Abstraction license, but observed pumping, then abstraction licnese is max(pumping of the dataseries)
#also if the actual pumping is higher than the license, then increase the license to a more realistic value

nb_nul_license=0  # DATA QUALITY CHECK
for i in range(0,len(WF)):   
    if WF['AnlgTillad'][i]==0:
        nb_nul_license+=1
        # WF.loc[i,'AnlgTillad'] = WF.iloc[i][1:36].max()
    WF.loc[i,'AnlgTillad'] = max(WF.iloc[i][1:36].max(),WF['AnlgTillad'][i])
print('Number of wellfields with abstraction license = 0 : ',nb_nul_license)

# =============================================================================
# ADD MISSING DATA
# =============================================================================

# Add missing data by hand, municipality of Glostrup
WF = pd.concat([WF,pd.DataFrame({'AnlaegID':[int(106321)],'AnlgTillad':[int(855000)],'Catch_ID':[int(17)]})], ignore_index=True)
WF = pd.concat([WF,pd.DataFrame({'AnlaegID':[int(107060)],'AnlgTillad':[int(695000)],'Catch_ID':[int(19)]})], ignore_index=True)

# Add Tinghøj reservoir as a WF with 0 abstraction license + add WSA Anlaeg connections
WF = pd.concat([WF,pd.DataFrame({'AnlaegID':[int(1)],'AnlgTillad':[int(0)],'Catch_ID':[int(16)]})], ignore_index=True)

WSA_PLANT = pd.concat([WSA_PLANT,pd.DataFrame({'WSAID':[int(3560)],'AnlaegID':[int(1)],'start':[int(2019)],'end':[int(2019)]})], ignore_index=True)
WSA_PLANT = pd.concat([WSA_PLANT,pd.DataFrame({'WSAID':[int(3564)],'AnlaegID':[int(1)],'start':[int(2019)],'end':[int(2019)]})], ignore_index=True)
WSA_PLANT = pd.concat([WSA_PLANT,pd.DataFrame({'WSAID':[int(3565)],'AnlaegID':[int(1)],'start':[int(2019)],'end':[int(2019)]})], ignore_index=True)
WSA_PLANT = pd.concat([WSA_PLANT,pd.DataFrame({'WSAID':[int(3559)],'AnlaegID':[int(1)],'start':[int(2019)],'end':[int(2019)]})], ignore_index=True)


# =============================================================================
# Link between Wellfields and Waterworks (Anlaeg/overanlaeg)
# =============================================================================

# If one anlaeg is its own overanlaeg, set overanlaeg to 0
Anlaeg_Overanlaeg = pd.read_excel('anlaegid_hieraki.xlsx')
for i in range(0,len(Anlaeg_Overanlaeg)):
    if Anlaeg_Overanlaeg['OVERANLAEG'][i]==Anlaeg_Overanlaeg['ANLAEGID'][i]:
        Anlaeg_Overanlaeg.loc[i,'OVERANLAEG']=0


WW_list = [] # list that will contain all WWID that are connected to a WSA 
# WF['WWID'] = WF['AnlaegID'] #initialize connexion between WF and WW (by default, WF = WW)
WF['WWID'] = 0
Table_WF = pd.DataFrame(columns=['WFID','WWID','AnlgTillad','Catch_ID'])

for i in range(1,len(WF)):   
    start_WFID = WF['AnlaegID'][i]
    queue = [start_WFID]
    
    while queue:
        current_WFID = queue.pop(0)
        if sum(WSA_PLANT['AnlaegID']==current_WFID) >= 1:  # do we find the anlaegID in WSA_PLANT data?
            WWID = current_WFID
            WW_list.append(WWID)
            WF.loc[WF['AnlaegID']==start_WFID,'WWID'] = WWID
            
            # then we store the data of the new WW created in the Table WF data base
            new_row = pd.DataFrame({'WFID':[start_WFID],'WWID':[WWID],'AnlgTillad':[WF[WF['AnlaegID'] == start_WFID]['AnlgTillad'].iloc[0]],'Catch_ID':[int(WF[WF['AnlaegID'] == start_WFID]['Catch_ID'].iloc[0])]})
            Table_WF = pd.concat([Table_WF, new_row])
            
            # and we look for a new AnlaegID if there is an overanlaeg
            new_AnlaegID = Anlaeg_Overanlaeg.loc[Anlaeg_Overanlaeg['ANLAEGID']==current_WFID,'OVERANLAEG']
            queue.extend(new_AnlaegID.tolist())
            
        else:
            new_AnlaegID = Anlaeg_Overanlaeg.loc[Anlaeg_Overanlaeg['ANLAEGID']==current_WFID,'OVERANLAEG']
            queue.extend(new_AnlaegID.tolist())
     
     
# Add missing Regnmeark 25132 abstraction capacity (around 6500000 m3/year)
missing_regnemark = 15000000 - WF[WF['WWID']==25132]['AnlgTillad'].sum()
Table_WF = pd.concat([Table_WF,pd.DataFrame({'WFID':[int(2)],'AnlgTillad':[int(missing_regnemark)],'Catch_ID':[int(0)], 'WWID':[int(25132)] })], ignore_index=True)
     
nb_WF = len(WF)  #DATA QUALITY CHECK
nb_WW_to_WF = len(set(WW_list))  #DATA QUALITY CHECK   #Nb of WW = unique value in the list because there are some redundancies
print('Total nb of WF: ', nb_WF)
print('Total number of WF connected to a WSA via a WW: ', len(WW_list), '\n(The WW can be equal to WF)')
print('Total abastraction license of Wf connected to WSA: ', round(WF[WF['WWID'] > 0]['AnlgTillad'].sum()/1000), '(1000m3)')

# =============================================================================
# Waterworks table
# =============================================================================
    
WW = pd.DataFrame({'AnlaegID':WW_list})
WW = pd.merge(WSA_PLANT, WW, on="AnlaegID",how='right') #how = 'inner', 'left, 'right', 'outer'

WW = WW.drop_duplicates()  # some rows are exactly the same for some reasons... remove them to not mess up with the crosstab

WW['Storage capacity (1000m3)'] = 0
WW['Storage initial (1000m3)'] = 0

nb_WW_to_WSA = len(WW.drop_duplicates('AnlaegID'))  #DATA QUALITY CHECK
print('Total number of WW, when merged with WSA: ', nb_WW_to_WSA)

# =============================================================================
# Tinghoj storage capacity
# =============================================================================

WD_CPH = 53000    # 1000m3/year
ww = 1
WW.loc[WW['AnlaegID']==ww, 'Storage capacity (1000m3)'] = WD_CPH/(365.25/2)   #storage capacity of Tinghoj is CPH waterdemand for a week

# =============================================================================
# Rename index
# =============================================================================

WW.rename(columns={'AnlaegID':'WWID'},inplace=True)
WF.rename(columns={'AnlaegID':'WFID'},inplace=True)

# =============================================================================
# Create a matrix WF WW
# =============================================================================

matrix_WF_WW = pd.crosstab(Table_WF['WFID'], Table_WF['WWID'])

# =============================================================================
# Save data
# =============================================================================

WW.to_csv(fnameout+'\Table_WW_'+Catch_Geo+'.csv')
Table_WF.to_csv(fnameout+'\Table_WF_'+Catch_Geo+'.csv')
matrix_WF_WW.to_csv(fnameout+'\Matrix_WF_WW_'+Catch_Geo+'.csv')


#%% TABLE WSA

# Process WSA_Catchments
WSA_Catchments = pd.read_csv('WSA_'+Catch_Geo+'.csv',sep=',')
WSA_Catchments.drop_duplicates('WSAID', inplace=True)   # drop duplicates because we only need one row per WSA
WSA_Catchments.drop('Area (m2)', axis=1, inplace=True)
WSA_Catchments=WSA_Catchments.dropna(subset=['Catch_ID'])
WSA_Catchments.drop(WSA_Catchments.columns[[3]], axis=1, inplace=True)

nb_WSA = len(WSA_Catchments)  #DATA QUALITY CHECK
print('Total nb of WSA in the catchments: ',nb_WSA)

# =============================================================================
# Drops all WSA that are not linked to any WW
# =============================================================================

WSA_list = WSA_Catchments['WSAID'].tolist()
missing_WSA = pd.DataFrame()  # datafrmae that will contain informaiton about the WSA lost during processing

for w in WSA_list:
    if sum(WW['WSAID']==w) ==0:   # if the WSA index is not in WW table, remove the WSA...
        index= WSA_Catchments[WSA_Catchments['WSAID']==w].index
        
        new_row = WSA_Catchments[WSA_Catchments['WSAID']==w]    # build a new dataframe row per row
        missing_WSA = pd.concat([missing_WSA, new_row], ignore_index=True)
                
        WSA_Catchments.drop(index,axis=0,inplace=True)
        
        
nb_WSA_to_WW = len(WSA_Catchments)  #DATA QUALITY CHECK
print('WSA linked to a WW: ', nb_WSA_to_WW)



# Process WSA Municipality
WSA_Municipality = pd.read_csv('WSA_Municipality.csv',sep=',')

# WSA_Municipality=WSA_Municipality.dropna(subset=['nationalco'])
WSA_Municipality.drop(WSA_Municipality.columns[[0,1,3,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20]], axis=1, inplace=True)
WSA_Municipality.rename(columns={'name_gn_sp':'name_municipality'},inplace=True)

#Calculate municipality area s the sum of all WSA in municipality (this is not the real area, but it allows me to do percentage of WSA in each municipality)
Municipality_area=WSA_Municipality.groupby(['name_municipality'], as_index=False).sum()
# Municipality_area.drop(Municipality_area.columns[[1,2,3,4,5,6,7]],axis=1,inplace=True)
Municipality_area = Municipality_area.filter(['name_municipality','WSA_Area (m3)'])
Municipality_area.rename(columns={'WSA_Area (m3)':'Municipality_Area (m3)'},inplace=True)

WSA_Municipality = pd.merge(WSA_Municipality, Municipality_area, on="name_municipality",how='inner')

Total_job_Ind_Sjaelland=WSA_Municipality.drop_duplicates('name_municipality')['Job_distribution_Industries'].sum()
Total_job_PS_Sjaelland=WSA_Municipality.drop_duplicates('name_municipality')['Job_distribution_Services'].sum()
Total_job_Agri_Sjaelland=WSA_Municipality.drop_duplicates('name_municipality')['Job_distribution_Agriculture'].sum()

wateruse=pd.read_excel('DK_Water_consumption_processed.xlsx')
WSA_Municipality['Wateruse agriculture (1000m3)']= float(wateruse[wateruse['Category']=='Agri']['2023'].iloc[0])*(WSA_Municipality['Job_distribution_Agriculture']/Total_job_Agri_Sjaelland)*(WSA_Municipality['WSA_Area (m3)']/WSA_Municipality['Municipality_Area (m3)'])
WSA_Municipality['Wateruse services (1000m3)']= float(wateruse[wateruse['Category']=='PS']['2023'].iloc[0])*(WSA_Municipality['Job_distribution_Services']/Total_job_PS_Sjaelland)*(WSA_Municipality['WSA_Area (m3)']/WSA_Municipality['Municipality_Area (m3)'])
WSA_Municipality['Wateruse industries (1000m3)']= float(wateruse[wateruse['Category']=='Ind']['2023'].iloc[0])*(WSA_Municipality['Job_distribution_Industries']/Total_job_Ind_Sjaelland)*(WSA_Municipality['WSA_Area (m3)']/WSA_Municipality['Municipality_Area (m3)'])

WSA_total_wateruse = pd.DataFrame(WSA_Municipality['WSAID'])
WSA_total_wateruse['Subtotal (1000m3)'] = WSA_Municipality['Wateruse agriculture (1000m3)'] + WSA_Municipality['Wateruse services (1000m3)'] + WSA_Municipality['Wateruse industries (1000m3)']
missing_WSA = pd.merge(missing_WSA,WSA_total_wateruse, on='WSAID', how = 'left')

# Process WSA_Population
WSA_Population = pd.read_csv('WSA_Population.csv',sep=',')

# WSA_Population=WSA_Population.dropna(subset=['indbygge_2', 'WSAID'])
WSA_Population=WSA_Population.groupby(['WSAID'],as_index=False).sum()
# WSA_Population.drop(WSA_Population.columns[[1,2,3,4,6]], axis=1, inplace=True)
WSA_Population = WSA_Population.filter(['WSAID','indbygge_4'])
WSA_Population.rename(columns={'indbygge_4':'Population #'},inplace=True)

Pop_Sjaelland=WSA_Population['Population #'].sum()  #from indbyggertal 2024

WSA_Population['Wateruse households (1000m3)']=WSA_Population['Population #']*float(wateruse[wateruse['Category']=='HH']['2023'].iloc[0])/Pop_Sjaelland
# WSA_Population.drop(WSA_Population.columns[[2,3]], axis=1, inplace=True)

missing_WSA = pd.merge(missing_WSA,WSA_Population, on='WSAID', how = 'left')
missing_WSA['Total demand (1000m3)']=missing_WSA['Subtotal (1000m3)'] + missing_WSA['Wateruse households (1000m3)']

#Gather all data to create TABLE WSA
Table_WSA=pd.merge(WSA_Catchments, WSA_Municipality, on="WSAID",how='inner')
Table_WSA=pd.merge(Table_WSA, WSA_Population, on="WSAID",how='left')   # for not loosing any information on the industry and services... "left"
Table_WSA['Population #'] = Table_WSA['Population #'].fillna(0)   # replace nan values with  (nan we got from the WSA with no population, but still there are industries etc..)
Table_WSA['Wateruse households (1000m3)'] = Table_WSA['Wateruse households (1000m3)'].fillna(0)
# Table_WSA.to_csv(fnameout+'\Table_WSA_'+Catch_Geo+'.csv')

nb_WSA_final = len(Table_WSA)  #DATA QUALITY CHECK
print('Final number of WSA: ',nb_WSA_final)
print('Total water demand of the capital region: ', round(Table_WSA['Wateruse households (1000m3)'].sum()+Table_WSA['Wateruse industries (1000m3)'].sum()+Table_WSA['Wateruse services (1000m3)'].sum() + missing_WSA['Total demand (1000m3)'].sum()), '(1000m3)')
print('Total water demand for only the connected WSA', round(Table_WSA['Wateruse households (1000m3)'].sum()+Table_WSA['Wateruse industries (1000m3)'].sum()+Table_WSA['Wateruse services (1000m3)'].sum()), '(1000m3)')
print('Percentage of lost water demand in the model: ',round(100*missing_WSA['Total demand (1000m3)'].sum()/(Table_WSA['Wateruse households (1000m3)'].sum()+Table_WSA['Wateruse industries (1000m3)'].sum()+Table_WSA['Wateruse services (1000m3)'].sum() + missing_WSA['Total demand (1000m3)'].sum()),1), '%' )

#%% Matrix

# =============================================================================
# WW to WSA matrix
# =============================================================================

matrix = pd.crosstab(WW['WWID'], WW['WSAID'])

WSA_list = Table_WSA['WSAID'].tolist() 
matrix_sized = matrix[WSA_list].copy()
matrix_sized.to_csv(fnameout+'\Matrix_WW_WSA_'+Catch_Geo+'.csv')


#%% TABLE WATER TRANSFER

Table_Water_Transfer=pd.read_excel('Table_Water_Transfer.xlsx')
Table_Water_Transfer.to_csv(fnameout+'\Table_Water_Transfer_'+Catch_Geo+'.csv')



#%% TABLE WTP

Table_WTP=pd.read_csv('AQUASTAT Dissemination System.csv')
Table_WTP.drop(Table_WTP.columns[[0,1,3,7,8]], axis=1, inplace=True)

Industrial=Table_WTP[Table_WTP['Variable']=='SDG 6.4.1. Industrial Water Use Efficiency']
Industrial.rename(columns={'Value':'WTP Industry (US$/m3)'},inplace=True)
Services=Table_WTP[Table_WTP['Variable']=='SDG 6.4.1. Services Water Use Efficiency']
Services.rename(columns={'Value':'WTP Services (US$/m3)'},inplace=True)
Agriculture=Table_WTP[Table_WTP['Variable']=='SDG 6.4.1. Irrigated Agriculture Water Use Efficiency']
Agriculture.rename(columns={'Value':'WTP Agriculture (US$/m3)'},inplace=True)

Table_WTP=pd.merge(Industrial, Services, on='Year',how='outer')
Table_WTP=pd.merge(Table_WTP, Agriculture, on='Year',how='outer')
# Table_WTP.drop(Table_WTP.columns[[0,3,4,6,7,9]],axis=1,inplace=True)
Table_WTP = Table_WTP.filter(['Year', 'WTP Industry (US$/m3)', 'WTP Services (US$/m3)', 'WTP Agriculture (US$/m3)'])

Rate=pd.read_csv('DEXDNUS.csv')  #https://fred.stlouisfed.org/series/DEXDNUS
Rate['year']=0
for i in range(0,len(Rate)):
    Rate.loc[i, 'year']=int(Rate['observation_date'][i][0:4])
Rate.drop(['observation_date'],axis=1,inplace=True)
Rate=Rate.groupby(['year'],as_index=False).mean()
Rate.rename(columns={'DEXDNUS':'Rate'},inplace=True)

Table_WTP['WTP Industry (US$/m3)']=Table_WTP['WTP Industry (US$/m3)']*Rate['Rate'][0:len(Table_WTP)]
Table_WTP['WTP Services (US$/m3)']=Table_WTP['WTP Services (US$/m3)']*Rate['Rate'][0:len(Table_WTP)]
Table_WTP['WTP Agriculture (US$/m3)']=Table_WTP['WTP Agriculture (US$/m3)']*Rate['Rate'][0:len(Table_WTP)]
Table_WTP.rename(columns={'WTP Industry (US$/m3)':'WTP Industry (DKK/m3)','WTP Services (US$/m3)':'WTP Services (DKK/m3)','WTP Agriculture (US$/m3)':'WTP Agriculture (DKK/m3)'},inplace=True)

Table_WTP['WTP Households (DKK/m3)']=77.25+12.48 # from Vand i tal 2024 (DANVA) and Brundtland 

Avg_Agri=Table_WTP['WTP Agriculture (DKK/m3)'].mean()
Table_WTP['WTP Agriculture (DKK/m3)'] = Table_WTP['WTP Agriculture (DKK/m3)'].fillna(Avg_Agri)

new_row = pd.DataFrame({"Year": [2022], "WTP Industry (DKK/m3)": [float(Table_WTP['WTP Industry (DKK/m3)'][Table_WTP['Year']==2021].iloc[0])], "WTP Services (DKK/m3)": [float(Table_WTP['WTP Services (DKK/m3)'][Table_WTP['Year']==2021].iloc[0])], "WTP Agriculture (DKK/m3)": [float(Table_WTP['WTP Agriculture (DKK/m3)'][Table_WTP['Year']==2021].iloc[0])], "WTP Households (DKK/m3)": [float(Table_WTP['WTP Households (DKK/m3)'][Table_WTP['Year']==2021].iloc[0])]})
df = pd.concat([Table_WTP, new_row], ignore_index=True)
# print(df)

df.to_csv(fnameout+'\Table_WTP_'+Catch_Geo+'.csv') 


#%% Adjusting the wateruse data

HH = Table_WSA['Wateruse households (1000m3)'].sum()
Ind = Table_WSA['Wateruse industries (1000m3)'].sum()
PS = Table_WSA['Wateruse services (1000m3)'].sum()
TOTAL = HH + Ind + PS
pop = Table_WSA['Population #'].sum()

avg_HH = 1000*HH/pop
avg_tot = 1000*TOTAL/pop

# data from DANVA Vand i tal 2024
ref_HH = 35.89
ref_tot = 59.4

factor_HH = ref_HH/avg_HH
factor_Ind_PS = (ref_tot-ref_HH)/(avg_tot-avg_HH)

Table_WSA['Wateruse households (1000m3)'] = factor_HH * Table_WSA['Wateruse households (1000m3)']
Table_WSA['Wateruse industries (1000m3)'] = factor_Ind_PS * Table_WSA['Wateruse industries (1000m3)']
Table_WSA['Wateruse services (1000m3)'] = factor_Ind_PS * Table_WSA['Wateruse services (1000m3)']

Table_WSA.to_csv(fnameout+'\Table_WSA_'+Catch_Geo+'.csv')