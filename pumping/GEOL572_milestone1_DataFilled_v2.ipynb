"""
Created on Wed Apr  8 20:30:45 2020

@author: ckpch
"""

#!/usr/bin/env python
# coding: utf-8

# In[188]:


import numpy as np
import pandas as pd
import urllib, json

#This line will have the full data set printed out. 
#Toggle off and clear output if you no longer want to see full dataset
pd.set_option("display.max_rows", None, "display.max_columns", None)


# In[189]:

#raw_data provides all of the data from the CSV file
#df1 drops deep wells and those with unknown locations or depths
#read in the csv
raw_csv = 'Pumpage_Data_WillCounty.csv'

#Create a copy of the data is not edited
#This will be needed for McKaleigh for 3d
raw_data = pd.read_csv(raw_csv)
print(raw_data)

print('size of original dataframe')
print(raw_data.shape)

#This copy of the data we will edit heavily 
df1 = pd.read_csv(raw_csv)

#drop rows if the depth, lamx, or lamy is unknown
df1.dropna(subset=['depth_total_last_known', 'lam_x', 'lam_y'], inplace = True)
    
#drop all of the rows that are deep wells (>400 ft)
df1.drop(df1[df1['depth_total_last_known'] > 400].index, inplace = True) 

print('size of updated dataframe')
print(df1.shape)

#%%
#df2 provides all of the old data but combined by facility; small fry facilities (<0.1mgd) were removed
#reset index to be based on well owner
df2=df1.set_index('owner')
#delete extra rows so that the dataframe will just be owner and the years of pumping
df2=df2.drop(['p_num', 'isws_facility_id', 'fac_well_num','depth_total_last_known', 'lam_x', 'lam_y'], axis=1)
#sum together all pumping for the same well owners so that the dataset is by pumping by facility
df2=df2.groupby(level=0).sum(min_count=1)

#calculate the maximum pumping for each facility
df2['max']=df2[['1981', '1982', '1983',
       '1984', '1985', '1986', '1987', '1988', '1989', '1990', '1991', '1992',
       '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001',
       '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010',
       '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019']].max(axis=1) 
#remove all facilities with nan in max and where maximum well pumping is less than 0.1mgd / "small fry" wells
#df3 drops all rows where the max value is < 0 and wells that don't meet minimum pumpage, then these calculation columns are removed in df4
df3=df2.dropna(subset=['max'],axis=0)
df3['small']=df3['max']-36524250
df4=df3[df3.small > 0]
df4=df4.drop(['small', 'max'], axis=1)

#loop to calculate the average of each row and then turn any columns with values greater than twice the average to nan
for index, row in df4.iterrows():
    ave = df4.mean(axis=1)
    df4['ave']=df4.index.map(ave)
    for columns in df4.columns:
        df4[columns][df4[columns] > 2*df4.ave]=np.nan
        print(df4[columns])

#%%
#
#remove the more recent years
subset = df4[['2013', '2014', '2015', '2016', '2017', '2018', '2019']]

#drop the orignal recent years from the dataframe
modify3 = df4.drop(['2013', '2014', '2015', '2016', '2017', '2018', '2019'], axis='columns')


#to bridge the more recent data, forward fill, back fill, and subbing in 0s
subset.fillna(method='ffill', axis='columns', inplace=True, limit=3, downcast=None)
subset.fillna(method='bfill', axis='columns', inplace=True, limit=3, downcast=None)
subset.fillna(value = 0, axis=1, inplace=True)


#replace updated recent years and forming a new dataframe
df5 = pd.concat([modify3, subset], axis = 1)

#fill in for up to two years of lapsed reporting data, the rest NaNs become 0
df5 = df5.fillna(method='ffill', axis='columns', limit=1)
df5 = df5.fillna(method='bfill', axis='columns', limit=1)
df5 = df5.fillna(0)
