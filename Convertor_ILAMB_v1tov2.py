# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 14:38:06 2016

@author: kang
"""

def get_lat_lon_datatype(top_path, time_start, varname):
    
    import numpy as np
    from netCDF4 import Dataset
    import os
    
    file_path_temp = os.listdir(top_path+np.str(np.int(time_start)));
    file_path_temp = top_path+np.str(np.int(time_start))+ '/'+file_path_temp[5]
    
#    file_path_temp = top_path+np.str(np.int(time_start))+ '/' + varname + '_' + np.str(np.int(time_start)) + np.str(np.int(1)).zfill(2) + '.nc';
    
    fh = Dataset(file_path_temp, mode='r')
    lons = fh.variables['lon'][:]; #n_lons = np.size(lons)
    lats = fh.variables['lat'][:]; #n_lats = np.size(lats)
    tmp_info = fh.variables[varname]
    tmp  = tmp_info[:]
    tmp_units = tmp_info.getncattr('units')
    
    source_pathname = top_path.split('/')
    source_pathname = source_pathname[-3]

    print varname + ' + ' + source_pathname
    
    if varname == 'pr' and source_pathname == 'GPCP2':
    
        tmp_miss = -999.0;
    
    else:
            
        tmp_miss = tmp_info.getncattr('_FillValue')
        
        
    tmp_longname = tmp_info.getncattr('long_name')
    
    fh.close()
    
    print "Variable:" + tmp_longname
    print "Original Units:" + tmp_units;

    dim_data = np.shape(tmp)
    dim_data = np.array(dim_data)
    
    print "Dim:" + np.str(dim_data)

    return lons, lats, dim_data, tmp_miss, tmp_longname

def get_all_years(path_to_folder):
    
    import numpy as np
    import os
    
    folder_list  = os.listdir(path_to_folder)
    n_s = np.size(folder_list)
    year0 = np.zeros(n_s)
    for i in np.arange(0, n_s):
    
        names = folder_list[i]
        pd = os.path.isdir(path_to_folder+folder_list[i])
    
        if pd:
            year0[i] = np.double(names)

    year0 = year0[np.where(year0>0)];
              
    onset = np.int(min(year0));
    final = np.int(max(year0));

    return onset, final

def creat_noleap_date(time_start, time_end):
    
    # this function to creat a date list for given onset and final year
    #   1 the date unit is the days since 1850 - 01 - 01 00:00:00
    #   2 only for month
    #   3 noleap calendar

    import numpy as np

    year0 = 1850;
    year1 = 4000;
    
    n_years       = year1 - year0 + 1;
    n_months      = n_years * 12;
    days_in_month = np.array([31,28,31,30,31,30,31,31,30,31,30,31]);
    
    year_all      = np.arange(year0, year1+1)
    
    final = np.zeros(shape=(n_months,6))
    
    for i in np.arange(0,n_years):
        
        idx = np.arange(i*12, (i+1)*12)
        
        final[idx,0] = year_all[i]
    
        final[idx,1] = np.arange(0,12) + 1;
    
        final[idx,2] = days_in_month;
        
    final[:,3] = np.cumsum(final[:,2]) - final[:,2] * 0.5;
    
    final[:,5] = np.cumsum(final[:,2])
    
    final[1:,4] = final[0:-1,5]
    
    idx = np.where((final[:,0]>=time_start) & (final[:,0]<=time_end));
    
    idx = np.array(idx);
    idx = idx[0,:]
    
    final = final[idx,:]
    
    date_1 = final[:,[0,1]]; # year and month
    date_2 = final[:,2]    ; # days in a month
    date_3 = final[:,3]    ; # center date for each month
    date_4 = final[:,[4,5]]; # date boundaries
    
    return date_1, date_2, date_3, date_4

def convert_benchmark_dataset_1_variable_1_source(pathname_in, varname, sourcename, 
                                                  pathname_out, varunit_out, 
                                                  convert_type, convert_factor):
    
    import numpy as np
    import os
    import shutil
    from netCDF4 import Dataset
    
    # 1: set up some parameters:
    
    top_path = pathname_in+'/'+varname+'/'+sourcename+'/'+'derived'+'/';
    
    [time_start, time_end] = get_all_years(top_path);
    [lons, lats, dim_data, missing_value, longname] = get_lat_lon_datatype(top_path, time_start, varname);
    
    print "From "+np.str(time_start)+" to "+np.str(time_end)
    
    n_lats = np.size(lats)
    n_lons = np.size(lons)
    
    # varunit      ='Kg m-2 s';
    
    pathname    = pathname_out;
    varname_out  =varname;
#    varunit_out  ='kg m-2 s-1'; # same
    varlongname_out = longname;
    
    ndim         = np.size(dim_data); # 1: site; 2:spatial grid
    
    # pathname    = '/Users/kangwang/Desktop/ILAMB_CRU_Test/DATA/tas/CRU/'

    if os.path.exists(pathname+'/'+sourcename):
        shutil.rmtree(pathname+'/'+sourcename)
    
    os.mkdir(pathname+'/'+sourcename)
    
    if ndim == 1:
        suffix = ''
    if ndim == 2:
        suffix = '_0.5x0.5'
        
    outfilename = pathname+'/'+sourcename+'/'+varname_out+suffix+'.nc'
        
    #    prepare the noleap date list for given oneset and final year:
    [date0, days_in_month, center_date, date_bnds] = creat_noleap_date(time_start=time_start, time_end=time_end); 
    
    # 2: read data from nc file and compose to a single array
    
    n_years  = np.int(time_end - time_start +1);
    n_months = np.int(n_years * 12);
    
    if ndim ==1:
        out_array_size = np.array([n_months, n_lats])
    
    if ndim ==2:
        out_array_size = np.array([n_months, n_lats, n_lons])
    
    final = np.zeros(out_array_size);
            
    for i in np.arange(0, n_years):
        
    #    print path_to_folder+np.str(np.int(time_start+i))
        year_folder_path = top_path+np.str(np.int(time_start+i));
        
        for j in np.arange(0,12):
            
            my_example_nc_file = year_folder_path + '/' + varname + suffix +'_' + np.str(np.int(time_start+i)) + np.str(np.int(j+1)).zfill(2)+'.nc'
            
            print my_example_nc_file
            
            fh = Dataset(my_example_nc_file, mode='r')
            
    #        lons = fh.variables['lon'][:]; n_lons = np.size(lons)
    #        lats = fh.variables['lat'][:]; n_lats = np.size(lats)
    
            tmp_info = fh.variables[varname]
    #        tmp_units = tmp_info.getncattr('units')
            tmp = tmp_info[:];
            
            fh.close()
            
#            tmp1 = tmp.data;

            u = np.int(i*12+j)
                    
            final[u,:] = tmp;
        
#    idx00 = np.where(final == missing_value)
#        
#    idx00 = np.array(idx00)
#    idx00 = idx00.size
#                    
#    if idx00>0:
#            
#        idx11 = np.where(final == missing_value)
#        
#        final[idx11] = np.nan;
    #
    ## 2.1: convert the data to match standard unit in CFUNITS:
    #
    ##days_in_month = np.repeat(days_in_month, n_lons*n_lats);
    ##secs_in_month = np.reshape(days_in_month,[n_times,n_lats,n_lons]) * 3600.0*24.0; # just for conversion of precp units
    ##
    ##convert_factor = 1./secs_in_month; # if nothing to change, please set it to ONE. ******************* TO be SET
    #
    # convert_factor = convert_factor;

    if convert_type.lower() == '+':
        final = final + convert_factor;
        missing_value = missing_value + convert_factor;
    if convert_type.lower() == '*':
        final = final * convert_factor;
        missing_value = missing_value * convert_factor;
#        
#    if idx00>0: 
#        final[np.where(np.isnan(final))] = missing_value;
        
    #
    #####################################################
    # 2.2 MAKE SURE the time series is noleap calendar:
    #####################################################
    
    times = center_date; # Change original date to noleap calendar
    
    # 3: write out the file:
    
    ncid = Dataset(outfilename, 'w', format='NETCDF4',clobber=False)
    
    if ndim == 1:
        
        # dimension         
        ncid.createDimension('time', n_months)
        ncid.createDimension('data', n_lats)
        
        # variables
        time = ncid.createVariable('time','f',('time')) 
        lat = ncid.createVariable('lat','f',('data'))   
        lon = ncid.createVariable('lon','f',('data')) 
        tas = ncid.createVariable(varname_out,'f',('time','data'), fill_value=missing_value)
        
    if ndim == 2:
        
        # dimension         
        ncid.createDimension('time', n_months)
        ncid.createDimension('lat', n_lats)
        ncid.createDimension('lon', n_lons)
        
        # variables
        time = ncid.createVariable('time','f',('time')) 
        lat = ncid.createVariable('lat','f',('lat'))   
        lon = ncid.createVariable('lon','f',('lon')) 
        tas = ncid.createVariable(varname_out,'f',('time','lat','lon'), fill_value=missing_value)
    
    # att
    lat.units = 'degrees_north';
    lat.long_name ='latitude';
    
    lon.units = 'degrees_east';
    lon.long_name ='longitude'; 
    
    tas.units = varunit_out;
    tas.long_name = varlongname_out
    
    time.units = 'days since 1850-01-01 00:00:00';
    time.long_name = 'time';
    time.standard_name = 'time';
    time.calendar = 'noleap';

    time[:] = times;
    tas[:]  = final; 
    lat[:]  = lats;
    lon[:]  = lons;
    
    ncid.close()
    
# =======================================================
# =======================================================

import os
import numpy as np

# varname         = 'co2';
#varname         = 'burntArea';
#varname         = 'albedo';
#varname         = 'tas';
#varname         = 'rsds';
#varname          = 'tsl';
varname  = 'snd';

#varunit_out     = 'Watt m-2';
#varunit_out     = 'kg m-2 s-1'
#varunit_out     = 'kg m-2'
# varunit_out     = 'ppm'
#varunit_out     = 'K'
#varunit_out     = '%';
#varunit_out     = '1';
#varunit_out     = 'unitless';
varunit_out = 'm';

convert_type    = '*'; # NOT available 
convert_factor  = 1; # NOT available 
pathname_in     = '.'; 

source_list = os.listdir(pathname_in+'/'+varname); nn = np.size(source_list)

for i in np.arange(nn):
    
#    source_list = 'GPCC';
    
    if os.path.isdir(pathname_in+'/'+varname+'/'+source_list[i]):
    
        sourcename      = source_list[i]
        pathname_out    = 'v2/'+varname
        
        if os.path.exists(pathname_out) == False:
            os.mkdir(pathname_out)
        
        print "Starting to convert the data from "+sourcename+"...";
    
        convert_benchmark_dataset_1_variable_1_source(pathname_in, varname, sourcename, 
                                                      pathname_out, varunit_out, 
                                                      convert_type, convert_factor)