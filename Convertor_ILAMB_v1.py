# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 14:38:06 2016

@author: kang
"""

import numpy as np
import datetime
import os
import shutil
from netCDF4 import Dataset

varname = 'snd';
long_name = 'Snow Depth'
units = 'm';

data0 = np.loadtxt('snd.csv', 
                   skiprows=1, 
                   delimiter=',', 
                   usecols=(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17))

lat0 = data0[:,1]
lon0 = data0[:,2]
year0 = data0[:,3]
name0 = data0[:,0]
#year = np.unique(year)
year_selected = np.arange(2001, 2015)

data000 = data0[:,0:3]
p = np.unique(data0[:,0])

lat_uniq = p + 0;
lon_uniq = p + 0;

site_name_all = '';

for i in np.arange(0,np.size(lat_uniq)):
    id0 = np.where(data0[:,0] == p[i])
    id0 = id0[0]
    id0 = id0[0]
    
    lat_uniq[i] = data0[id0,1]
    lon_uniq[i] = data0[id0,2]
    
    site_name0 = 'STN_'+str(i+1);
    
    if i == 0:
        site_name_all = site_name0;
    else:
        site_name_all = site_name_all+','+site_name0;

time_offset = datetime.date(1850, 1, 1)

for i in year_selected:
    
    i_year = i
    
    data_index_yearly = np.where(year0 == i_year)

    path_out = varname+'/pbs/derived/'+str(i_year);

    if os.path.isdir(path_out):
        shutil.rmtree(path_out)
    
    os.makedirs(path_out);

    for j in np.arange(1,13):    
    
        i_month = j
        date = datetime.date(i_year,i_month,1)
        time = datetime.date.toordinal(date)-datetime.date.toordinal(time_offset)+1
        time = np.double(time)        
        
        name1 = name0[data_index_yearly]
        lat1 = lat0[data_index_yearly]
        lon1 = lon0[data_index_yearly]  
        
        data_1 = data0[data_index_yearly, j+3]
        data_1[np.where(data_1<-100.)] = np.nan
        data_1 = data_1+0;
        data_1[np.where(np.isnan(data_1))] = -999.
#        data_1 = np.float(data_1)
        
        # unifing to same sites:
        z = np.in1d(p,name1).nonzero()
        data_2 = p *0.0 - 999.;    
        data_2[z] = data_1;
        
        outfilename = path_out+'/'+varname+'_'+str(i_year)+str(i_month).zfill(2)+'.nc'
        
        ncid = Dataset(outfilename, 'w', format='NETCDF4')
        
        ncid.site_name = site_name_all
        
        # dimension         
        ncid.createDimension('data', np.size(data_2))

        # variables
        lat = ncid.createVariable('lat','f',('data',), fill_value=-999.0)   
        lon = ncid.createVariable('lon','f',('data',), fill_value=-999.0) 
        tsl = ncid.createVariable(varname,'f',('data',), fill_value=-999.0)
        
        # att
        lat.units = 'degrees_north';
        lat.long_name ='latitude';

        lon.units = 'degrees_east';
        lon.long_name ='longitude'; 
        
        tsl.units = units;
        tsl.long_name = long_name
        tsl.time_unit = 'No. of days since 1850-1-1'
        tsl.time = time;
        
        # data
        lat[:] = lat_uniq;
        lon[:] = lon_uniq;
        tsl[:] = data_2;     
        
        ncid.close()

