import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, DateFormatter
import datetime
from scipy import stats
from scipy.optimize import curve_fit
import math as m
from scipy.stats import linregress
from netCDF4 import Dataset

import scipy.io

def mean(a):
    return np.nansum(a) / np.count_nonzero(~np.isnan(a))

# -------------------------------------------------------------------------------------------
# TEST script to create 5 minute time scan cadance
# --> DON'T average in time together
# (((--> average every 15 minutes together in 1 profile)))
# -------------------------------------------------------------------------------------------
#taken out bc HDF5 '160113',
lidar_dates = ['160105', '160106',  '160114', '160118', '160120', '160121', '160123', '160126', '160127', '160129', '160201', '160215', '160216']

for obs_date in lidar_dates:


    m_file = scipy.io.loadmat('/Users/marinstanev/Dropbox/MISU/Data/esrange_mat_format/16jan/Temp_matrix_'+obs_date+'.mat')

    print obs_date
    #print m_file
    #print m_file['Temp_matrix'][12][0]
    #print len(m_file['Temp_matrix'][2])
    #print len(m_file['Temp_matrix'][12])
    #print len(m_file['date_time'])
    #print len(m_file['Temp_matrix'])
    
    # Get altitudes
    alts = [x[0]/1000. for x in m_file['H']]
    #print len(alts)
    #sys.exit(0)


    # Create NetCDF file
    f = Dataset('/Users/marinstanev/Dropbox/MISU/Data/esrange_mat_format/16jan/Rayleigh_channels_'+obs_date+'.nc', 'w')

    f.createDimension('Altitude', len(alts))
    f.createDimension('Ray_XLow', None)
    f.createDimension('Ray_Low', None)
    f.createDimension('Ray_XHigh', None)
    f.createDimension('Ray_Med', None)
    f.createDimension('Vib_Ram', None)
    f.createDimension('Rot_Ram_p', None)
    f.createDimension('Rot_Ram_X', None)

    nc_alt = f.createVariable('Altitude', 'f', ('Altitude',), zlib=True)
    nc_alt[:] = alts

    # 2==Rayleigh Medium, 4==Rayleigh Low, 6==Rayleigh X-High, 8==Rayleigh X-Low
    # 10==Vibrational Raman, 12==Rotational Raman parallel, 14==Rotational Raman vertical
    Cs = [2,4,6,8,10,12,14]

    for C in Cs:
        
        print 'Lidar channel ', C
        # ----------------------------------------------------------------------------------
        # Insert aprrox. 5 minute cadance
        # --> insert a 'scan' of NaNs if the time difference is larger than 5 minutes
        # ----------------------------------------------------------------------------------

        # Get data of Rayleigh X Low channel
        ray_m = m_file['Temp_matrix'][:,C].tolist()
        print len(ray_m), len(ray_m[0])


        # Convert date and time into python datetime items
        dates = []
        x_old = datetime.datetime(int(m_file['date_time'][0][0]), int(m_file['date_time'][0][1]), int(m_file['date_time'][0][2]), int(m_file['date_time'][0][3]), int(m_file['date_time'][0][4]), int(m_file['date_time'][0][5]))
        dates.append(x_old)
        
        for i in range(len(m_file['date_time'])):
            x = datetime.datetime(int(m_file['date_time'][i][0]), int(m_file['date_time'][i][1]), int(m_file['date_time'][i][2]), int(m_file['date_time'][i][3]), int(m_file['date_time'][i][4]), int(m_file['date_time'][i][5]))
            #print x
            x_orig = x
            if i>1:
                dates.append(x)
                dates_size = len(dates)-1
                # If step between times is larger than 5 minutes, insert a NaN scan
                while (x-x_old)>datetime.timedelta(seconds=5*60):
                    #print dates
                    x_new = x-datetime.timedelta(seconds=5*60)
                    #x_new = x_old+datetime.timedelta(seconds=5*60)
                    print x_new
                    dates.insert(dates_size,x_new)
                    for j in range(len(ray_m)):
                        ray_m[j].insert(i,np.nan)
                    x = x_new
                #dates.append(x_orig)
                x_old=x_orig
                #print dates
                print '---'
                #if i==3:
                #    sys.exit(0)
            elif i==1:
                dates.append(x)
                x_old = x
        # ----------------------------------------------------------------------------------

        # Put dates in string format for saving into NetCDF
        dates2 = []
        for i in range(len(dates)):
            x = dates[i]
            dates2.append(int(datetime.datetime.strftime(x, '%Y%m%d%H%M%S')))
        dates = dates2




        print len(ray_m), len(ray_m[0]), len(dates)
        t,count=0,0
        while t<len(dates)-1:
            #print t, count, round(len(dates)/3.)
            
            # Create once the NetCDF dimensions and variables
            if C==2 and count==0:
                f.createDimension('Date', int(round((len(dates)/1.-1))))#3.))))
                nc_date = f.createVariable('Date', 'i8', ('Date',), zlib=True)
                nc_XLow= f.createVariable('Ray_XLow', 'f', ('Altitude','Date',), zlib=True)
                nc_Low = f.createVariable('Ray_Low', 'f', ('Altitude','Date',), zlib=True)
                nc_XHigh= f.createVariable('Ray_XHigh', 'f', ('Altitude','Date',), zlib=True)
                nc_Med = f.createVariable('Ray_Med', 'f', ('Altitude','Date',), zlib=True)
                nc_vib_R = f.createVariable('Vib_Ram', 'f', ('Altitude','Date',), zlib=True)
                nc_rot_R_p= f.createVariable('Rot_Ram_p', 'f', ('Altitude','Date',), zlib=True)
                nc_rot_R_X = f.createVariable('Rot_Ram_X', 'f', ('Altitude','Date',), zlib=True)
            
            # Get the middle time of the 3 times
            # COULD BE BETTER: AVERAGED OVER THE 3 TIMES ISO TAKING THE MIDDLE ONE !!!
            nc_date[count] = dates[t]#+1]
            
            #o_ray_m = [ray_m[x][t:t+3] for x in range(len(ray_m))]
            #o_ray_m = map(mean,zip(o_ray_m))

            o_ray_m = [ray_m[x][t] for x in range(len(ray_m))]
            o_ray_m = map(mean,zip(o_ray_m))
            
            if C==2:
                nc_Med[:,count] = o_ray_m
            elif C==4:
                nc_Low[:,count] = o_ray_m
            elif C==6:
                nc_XHigh[:,count] = o_ray_m
            elif C==8:
                nc_XLow[:,count] = o_ray_m
            elif C==10:
                nc_vib_R[:,count] = o_ray_m
            elif C==12:
                nc_rot_R_p[:,count] = o_ray_m
            elif C==14:
                nc_rot_R_X[:,count] = o_ray_m
                    
            
            t+=1#3
            count+=1

        print '--'

    print '---'
    f.close()

