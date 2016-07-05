import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, DateFormatter
import datetime
from scipy import stats
from netCDF4 import Dataset
from dateutil.rrule import rrule, HOURLY


# --------------------------------------------------------------------------------
# Script to plot the original data, but where every 5 minutes data is present (data or NaNs)
# --> masks where NaNs are present
# --> save directly in '/Figures/raw_data'
# --------------------------------------------------------------------------------
#taken out bc no nc file bc HDF5 format: '160113',
lidar_dates = ['160105', '160106',  '160114', '160118', '160120', '160121', '160123', '160126', '160127', '160129', '160201', '160215', '160216']

for obs_date in lidar_dates:


    m_file = Dataset('/Users/marinstanev/Dropbox/MISU/Data/esrange_mat_format/16jan/Rayleigh_channels_'+obs_date+'.nc')



    # Get altitudes
    alts = [x for x in m_file['Altitude']]

    # Get data of Rayleigh X Low channel
    ray = m_file['Ray_Low']
    ram = m_file['Rot_Ram_p']

    Xray = m_file['Ray_XLow']
    Xram = m_file['Rot_Ram_X']

    dates = m_file['Date']
    dates = [datetime.datetime.strptime(str(x), '%Y%m%d%H%M%S') for x in dates]
    
    # Get an even hourly spaced list
    h_start = dates[0] - datetime.timedelta(minutes=dates[0].minute, seconds=dates[0].second)
    h_end = dates[-1] - datetime.timedelta(minutes=dates[-1].minute, seconds=dates[-1].second) + datetime.timedelta(minutes=60)
    #print h_start, h_end
    hours = [dt for dt in rrule(HOURLY, dtstart=h_start, until=h_end)]


    print len(ray), len(ray[0]), len(alts), len(dates)
    
    """
    fig,ax=plt.subplots(1)
    ax.plot(ray[:,100],alts)
    ax.set_xscale('log')
    plt.show()
    for i in range(len(dates)):
        print i, np.nanmean(ray[300:500,i])
        #print ray[300:500,i]
    sys.exit(0)
    """
    
    x1, x2 = 30, 200
    print alts[x1], alts[x2]
    alts = alts[x1:x2]
    #sys.exit(0)
    
    ray, ram = ray[x1:x2], ram[x1:x2]
    Xray, Xram = Xray[x1:x2], Xram[x1:x2]
    

    # Correction factor to be applied on the Raman channel
    constant, Xconstant = [], []
    for j in range(len(dates)):   # Loop over scans
        x1,x2 = [], []
        Xx1,Xx2 = [], []
        
        for i in range(len(alts)):
            if alts[i]>26.:
                # You want to compare to measurements at each altitude, with z>=27km (aerosol free)
                if not np.isnan(ray[i,j]) and not np.isnan(ram[i,j]):
                    if ray[i,j]>0 and ram[i,j]>0:
                        x1.append(ray[i,j])
                        x2.append(ram[i,j])
        
                if not np.isnan(Xray[i,j]) and not np.isnan(Xram[i,j]):
                    if Xray[i,j]>0 and Xram[i,j]>0:
                        Xx1.append(Xray[i,j])
                        Xx2.append(Xram[i,j])

        if len(x1)>0 and len(x2)>0: constant.append(np.max(x1)/np.max(x2))
        else: constant.append(np.nan)
        
        if len(Xx1)>0 and len(Xx2)>0: Xconstant.append(np.max(Xx1)/np.max(Xx2))
        else: Xconstant.append(np.nan)
    

    # Calculate the backscatter ratio now, with correction factor applied
    ratio = [[ray[i,j]/(constant[j]*ram[i,j]) for j in range(len(ram[i]))] for i in range(len(ram))]
    ratio = np.asarray(ratio)

    Xratio = [[Xray[i,j]/(Xconstant[j]*Xram[i,j]) for j in range(len(Xram[i]))] for i in range(len(Xram))]
    Xratio = np.asarray(Xratio)


    """
    fig,ax=plt.subplots(1)
    ax.plot(ray[:,100],alts,color='b')
    ax.plot(ram[:,100],alts,color='r')
    ax.scatter(ray[:,100],alts,color='b')
    ax.scatter(ram[:,100],alts,color='r')
    bla = [constant[100]*ram[i,100] for i in range(len(alts))]
    ax.plot(bla,alts,color='g')
    ax.set_xscale('log')
    
    ax.scatter(ratio[:,100],alts)
    ax.plot(ratio[:,200],alts)
    plt.show()
    sys.exit(0)
    """

    # Calculate depolarisation of the aerosol (delta_aero), in percent
    d_aero = np.ndarray(shape=(len(alts),len(dates)))
    for i in range(len(alts)):
        for j in range(len(dates)):
            if ratio[i,j]>1.06 and Xratio[i,j]>1.06 and not np.isinf(Xratio[i,j]) and not np.isinf(ratio[i,j]):
                delta = (Xratio[i,j]-1)/(ratio[i,j]-1)*0.36
            else: delta = np.nan

            d_aero[i,j] = delta


    # Get an hourly profile
    hour_integrated_x = np.ndarray(shape=(len(alts), len(hours)-1))
    hour_integrated_y = np.ndarray(shape=(len(alts), len(hours)-1))

    for i in range(len(alts)):
        for j in range(len(hours)-1):
            x,y = [], []
            for k in range(len(dates)):
                # Check if scan is within this hour range
                if hours[j]<=dates[k] and dates[k]<hours[j+1]:
                    # Check whether ratios are large enough
                    if ratio[i,k]>1.06 and Xratio[i,k]>1.06 and not np.isinf(ratio[i,k]) and not np.isinf(Xratio[i,k]):
                        x.append(ratio[i,k])
                        y.append(Xratio[i,k])
            if len(x)>0 and len(y)>0:
                hour_integrated_x[i,j], hour_integrated_y[i,j] = np.nanmean(x), np.nanmean(y)
            else:
                hour_integrated_x[i,j], hour_integrated_y[i,j] = np.nan, np.nan



    fig, ax = plt.subplots(1)

    X,Y = np.meshgrid(mpl.dates.date2num(dates),alts)

    ax.xaxis.set_major_locator(DayLocator())
    ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_minor_locator(HourLocator())
    ax.autoscale_view()

    # For plotting the backscatter ratio R
    lev_exp = np.arange(1.06,2.1,0.1)
    levels = np.power(10,lev_exp)
    im = ax.contourf(X, Y, Xratio)#,levels=levels)#, norm=mpl.colors.LogNorm())
    cbar = plt.colorbar(im)

    # For plotting the depolarisation of the aerosols
    #levels = np.arange(2,10.01,0.1)
    #im = ax.contourf(X, Y, d_aero, levels=levels,extend='both')
    #cbar = plt.colorbar(im,ticks=np.arange(2,10.1,1))

    ax.set_ylabel('Altitude [km]')
    plt.title(obs_date)
    plt.savefig('/Users/marinstanev/Dropbox/MISU/Plots/PSC/raw_data_'+obs_date+'.png')
    plt.show()
    

    #"""
    fig, ax = plt.subplots(1)
    
    # To plot every single value
    #x,y = [], []
    #for i in range(len(ratio)):
    #    for j in range(len(ratio[i])):
    #        if ratio[i,j]>1.06 and Xratio[i,j]>1.06:
    #            x.append(ratio[i,j])
    #            y.append(Xratio[i,j])
    #ax.scatter(x,y)
    

    R, xR, R_alt, R_time = [], [], [], []

    xx,yy = np.arange(0,1e2,1).tolist(), np.arange(0,1e4,1).tolist()
    print len(xx), len(yy)
    X, Y = np.meshgrid(yy,xx)

    counts = np.ndarray(shape=(len(xx),len(yy)))
    for i in range(len(counts)):
        for j in range(len(counts[i])):
            counts[i,j] = 0

    # Loop over the altitudes and hourly integrated scans
    for i in range(len(hour_integrated_x)):
        for j in range(len(hour_integrated_x[i])):
            
            if not np.isnan(hour_integrated_x[i,j]) and not np.isnan(hour_integrated_y[i,j]):
                
                a = xx.index(int(hour_integrated_x[i,j]))
                b = yy.index(int(hour_integrated_y[i,j]))
                counts[a,b] += 1

                R.append(hour_integrated_x[i,j])
                xR.append(hour_integrated_y[i,j])
                R_alt.append(alts[i])
                R_time.append(j)


    # Plot intensity of the cloud in bins
    #counts = np.ma.masked_where(counts==0, counts)
    #lev_exp = np.arange(0,3.5,0.1)
    #levels = np.power(10,lev_exp)
    #im = ax.contourf(Y,X,counts, norm=mpl.colors.LogNorm(), levels=levels)
    #cbar = plt.colorbar(im)
    
    # Plot backscatter ratios with altitude colour coded values
    cax = ax.scatter(R,xR,c=R_alt)
    cbar = plt.colorbar(cax)

    # Plot hourly integrated backscatter ratios
    #ax.scatter(hour_integrated_x,hour_integrated_y)
    
    ax.set_ylabel('Perpendicular backscatter ratio')
    ax.set_xlabel('Parallel backscatter ratio')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([1e0, 1e2])
    ax.set_ylim([1e0, 1e4])

    # For plotting the depolarisation lines
    pct = [0.3, 0.1, 0.03, 0.007, 0.0035]
    R_range = np.arange(1e0,1e2,0.1)
    for p in pct:
        l_1 = [(x-1.)*p/0.0036+1 for x in R_range]
        ax.plot(R_range,l_1,c='black')
        plt.text(R_range[-350], l_1[-500], str(p*100)+'%')
    plt.grid()
    plt.savefig('/Users/marinstanev/Dropbox/MISU/Plots/PSC/depol_'+obs_date+'.png')
    plt.show()
    
    
    R_array=np.array(R)
    xR_array=np.array(xR)
    R_alt_array=np.array(R_alt)
    R_time_array=np.array(R_time)
    depol_10line = xR_array-((R_array-1.)*0.1/0.0036+1)
    depol_4line= xR_array-((R_array-1.)*0.004/0.0036+1)
    psc_NAT=R_array[np.where((R_array<=2.0) & (depol_10line>0))]
    psc_NAT_alt=R_alt_array[np.where((R_array<=2.0) & (depol_10line>0))]
    psc_NAT_time=R_time_array[np.where((R_array<=2.0) & (depol_10line>0))]
    psc_ICE=R_array[np.where((R_array>2.0) & (depol_10line>0))]
    psc_ICE_alt=R_alt_array[np.where((R_array>2.0) & (depol_10line>0))]
    psc_ICE_time=R_time_array[np.where((R_array>2.0) & (depol_10line>0))]
    psc_STS=R_array[np.where((R_array<5.0) & (depol_4line<0))]
    psc_STS_alt=R_alt_array[np.where((R_array<5.0) & (depol_4line<0))]
    psc_STS_time=R_time_array[np.where((R_array<5.0) & (depol_4line<0))]
    psc_MIX1=R_array[np.where((R_array<5.0) & (depol_4line>0) & (depol_10line<0))]
    psc_MIX1_alt=R_alt_array[np.where((R_array<5.0) & (depol_4line>0) & (depol_10line<0))]
    psc_MIX1_time=R_time_array[np.where((R_array<5.0) & (depol_4line>0) & (depol_10line<0))]
    psc_MIX2=R_array[np.where((R_array>=5.0) & (depol_10line<0))]
    psc_MIX2_alt=R_alt_array[np.where((R_array>=5.0) & (depol_10line<0))]
    psc_MIX2_time=R_time_array[np.where((R_array>=5.0) & (depol_10line<0))]
    
    fig, ax = plt.subplots(1)
    plt.scatter(psc_NAT_time,psc_NAT_alt,c="yellow",marker="_",label="NAT")
    plt.scatter(psc_STS_time,psc_STS_alt,c="red",marker="_",label="STS")
    plt.scatter(psc_ICE_time,psc_ICE_alt,c="blue",marker="_",label="ICE")
    plt.scatter(psc_MIX1_time,psc_MIX1_alt,c="green",marker="_",label="MIX")
    plt.scatter(psc_MIX2_time,psc_MIX2_alt,c="green",marker="_")
    plt.legend(loc=4)
    plt.ylabel('Altitude in km')
    plt.xlabel('Time in h')
    plt.title('Time resolved PSC observation Date:'+obs_date)
    #plt.legend('NAT','STS','ICE','MIX')
    plt.savefig('/Users/marinstanev/Dropbox/MISU/Plots/PSC/type_'+obs_date+'.png')
    #plt.set_ylabel('alt in km')
    #plt.set_xlabel('Parallel backscatter ratio')
    #plt.title(obs_date)
    #plt.show()
    

   # """
