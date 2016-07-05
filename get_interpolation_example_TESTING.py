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
# TEST script to perform linear/parabola background fit
# --> check every 10th scan and make sure the altitude ranges are correct for PSCs
# -------------------------------------------------------------------------------------------

dates = ['160105', '160106', '160113', '160114', '160118', '160120', '160121', '160123', '160126', '160127', '160129', '160201', '160215', '160216']

date = dates[0]

# 4==Rayleigh Low, 8 = Rayleigh X-Low
C = 8
m_file = scipy.io.loadmat('./data/Lidar_original/Temp_matrix_'+date+'.mat')

print date

#print len(m_file)
#print len(m_file['Temp_matrix']), len(m_file['Temp_matrix'][0])
#print len(m_file['H'])
#print len(m_file['date_time']), len(m_file['date_time'][0])


# Get altitudes
alts = [x[0]/1000. for x in m_file['H']]

if int(date)==160105:   a_s, a_l, a_m, a_h = 30, 90, 145, 200
# For 2016-01-06
elif int(date)==160106:   a_s, a_l, a_m, a_h = 0, 75, 140, 200
elif int(date)==160113:   a_s, a_l, a_m, a_h = 0, 90, 150, 200
elif int(date)==160114:   a_s, a_l, a_m, a_h = 0, 15, 140, 200
elif int(date)==160118:   a_s, a_l, a_m, a_h = 0, 25, 140, 200
# For 2016-01-20
elif int(date)==160120:   a_s, a_l, a_m, a_h = 0, 45, 140, 200
elif int(date)==160121:   a_s, a_l, a_m, a_h = 32, 50, 140, 200
elif int(date)==160123:   a_s, a_l, a_m, a_h = 0, 80, 140, 200
elif int(date)==160126:   a_s, a_l, a_m, a_h = 0, 95, 160, 240
# For 2016-01-27
elif int(date)==160127:   a_s, a_l, a_m, a_h = 30, 90, 155, 200
elif int(date)==160129:   a_s, a_l, a_m, a_h = 40, 95, 160, 200
elif int(date)==160201:   a_s, a_l, a_m, a_h = 20, 55, 140, 200
elif int(date)==160215:   a_s, a_l, a_m, a_h = 0, 55, 140, 200
elif int(date)==160216:   a_s, a_l, a_m, a_h = 55, 95, 140, 200

else:
    a_s, a_l, a_m, a_h = 55, 95, 140, 170

print 'Altitude regions for parabola fit: ', alts[a_s], alts[a_l], alts[a_m], alts[a_h]
# All the selected altitudes
old_alts = alts[a_s: a_h]
# Select altitude ranges where we think no PSCs are present
x1, x2 = alts[a_s:a_l], alts[a_m:a_h]

new_alts = x2
NEWEST_alts = x1+x2





days = np.arange(0,100,10)
for day in days:
    print day
    
    # Get data of Rayleigh (X) Low channel
    ray_m = m_file['Temp_matrix'][:,C]

    # Select the data and take logarithm
    old_ray_m = ray_m[:,day][a_s: a_h]
    o_ray_m, o_alts, orig_ray_m = [], [], []
    for x,val in enumerate(old_ray_m):
        if val>0:
            orig_ray_m.append(old_ray_m[x])
            o_ray_m.append(m.log(old_ray_m[x]))
            o_alts.append(old_alts[x])

    # Select the data where we think there is no PSC
    y1 = ray_m[a_s:a_l,day].tolist()
    y2 = ray_m[a_m:a_h,day].tolist()
    # Points for the linear fit
    ray_m = y2
    nw_ray_m, nw_alts = [], []
    for x,val in enumerate(ray_m):
        if val>1:
            nw_ray_m.append(m.log(ray_m[x]))
            nw_alts.append(new_alts[x])

    # Points for the parabolic fit
    NEWEST_ray_m = y1+y2
    NW2_ray_m, NW2_alts = [], []
    for x,val in enumerate(NEWEST_ray_m):
        if val>1:
            NW2_ray_m.append(m.log(NEWEST_ray_m[x]))
            NW2_alts.append(NEWEST_alts[x])


    fig, ax = plt.subplots(1)

    ax.plot(o_ray_m,o_alts,c='green')
    
    if len(nw_alts)>0:
        # Fit a straight line
        coeff = np.polyfit(nw_alts, nw_ray_m, 1)
        poly = np.poly1d(coeff)
        new_pol = poly(o_alts)

        ax.scatter(nw_ray_m,nw_alts,c='blue')
        ax.plot(new_pol,o_alts,c='black')
        
        # Fit a parabola
        coeff = np.polyfit(NW2_alts, NW2_ray_m, 2)
        poly = np.poly1d(coeff)
        new_pol = poly(o_alts)
        
        ax.scatter(NW2_ray_m,NW2_alts,c='blue')
        ax.plot(new_pol,o_alts,c='red')


    plt.show()

