#!/usr/bin/env python

##This script creates a plot of salinity vs salinity difference for saildrone and each SSS platform
##The script also calculates the percentage of SMAP measurements that overestimate salinity 
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import datetime
import sys
import cftime
from decimal import Decimal

#data directory for saildrone and satellite orbital data
data_sat = '~/paper_software/2020_ATOMIC_Salinity/data/sss_collocations_orbital_norepeat/'
data_sat2 = '~/paper_software/2020_ATOMIC_Salinity/data/sss_collocations_8day_norepeat/'

#data directory for HYCOM data
data_dir1 = '~/hycom_files/'
files = glob(data_dir1+'*nc4')
print(files)
hycom=xr.open_mfdataset(data_dir1+'*nc4',concat_dim='time').isel(depth=0)

##this next section removes duplicate timesteps
_, index = np.unique(hycom['time'], return_index=True)
hycom2=hycom.isel(time=index)

#change hycom coordinates to match with saildrone(0-359:-180-179)
hycom_lon=hycom2.assign_coords(longitude=(((hycom2.lon + 180) % 360) - 180))
hycom2=hycom_lon.swap_dims({'lon':'longitude'})

#remove nans from hycom data
filled=hycom2.chunk({'time':-1}).interpolate_na(dim="time",method="nearest",fill_value='extrapolate')
filled2=filled.interpolate_na(dim="lat", method="nearest",fill_value='extrapolate')
filled3=filled2.interpolate_na(dim="lon", method="nearest",fill_value='extrapolate')

#Collocated Saildrone and Satellite data
JPL_0 = glob(data_sat+'*jpl*.nc')
RSS_0 = glob(data_sat+'*rss*.nc')

legend_properties = {'weight': 'semibold', 'size': '12'}

JPL=xr.open_dataset(JPL_0[0],decode_times=False)
RSS=xr.open_dataset(RSS_0[0],decode_times=False)

JPL_1=xr.open_dataset(JPL_0[1],decode_times=False)
RSS_1=xr.open_dataset(RSS_0[1],decode_times=False)

JPL_2=xr.open_dataset(JPL_0[2],decode_times=False)
RSS_2=xr.open_dataset(RSS_0[2],decode_times=False)

# Convert times from nano seconds to seconds
ns = 1e-9
rss_time = RSS.time.values * ns
jpl_time = JPL.time.values * ns

rss_time_1 = RSS_1.time.values * ns
jpl_time_1 = JPL_1.time.values * ns

rss_time_2 = RSS_2.time.values * ns
jpl_time_2 = JPL_2.time.values * ns

test_1=rss_time.astype(np.float)
test2_1=jpl_time.astype(np.float)

test_2=rss_time_1.astype(np.float)
test2_2=jpl_time_1.astype(np.float)

test_3=rss_time_2.astype(np.float)
test2_3=jpl_time_2.astype(np.float)

# Swap Dimensions
RSS = RSS.swap_dims({'ob': 'time'})
JPL = JPL.swap_dims({'ob': 'time'})

RSS_1 = RSS_1.swap_dims({'ob': 'time'})
JPL_1 = JPL_1.swap_dims({'ob': 'time'})

RSS_2 = RSS_2.swap_dims({'ob': 'time'})
JPL_2 = JPL_2.swap_dims({'ob': 'time'})

ss_times = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test_1]
jp_times = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test2_1]

ss_times_1 = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test_2]
jp_times_1 = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test2_2]

ss_times_2 = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test_3]
jp_times_2 = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test2_3]

RSS.assign_coords(time=ss_times)
JPL.assign_coords(time=jp_times)

RSS_1.assign_coords(time=ss_times_1)
JPL_1.assign_coords(time=jp_times_1)

RSS_2.assign_coords(time=ss_times_2)
JPL_2.assign_coords(time=jp_times_2)

num = RSS.trajectory.values
j_num = JPL.trajectory.values

num_1 = RSS_1.trajectory.values
j_num_1 = JPL_1.trajectory.values

num_2 = RSS_2.trajectory.values
j_num_2 = JPL_2.trajectory.values
print(num)

# interp HYCOM data to saildrone dimensions
hysal = filled3.interp(lat=JPL.lat, longitude=JPL.lon, time=jp_times, method='nearest')
hysal2 = hysal.salinity

# interp HYCOM data to saildrone dimensions
hysal_1 = filled3.interp(lat=JPL_1.lat, longitude=JPL_1.lon, time=jp_times_1, method='nearest')
hysal2_1 = hysal_1.salinity

# interp HYCOM data to saildrone dimensions
hysal_2 = filled3.interp(lat=JPL_2.lat, longitude=JPL_2.lon, time=jp_times_2, method='nearest')
hysal2_2 = hysal_2.salinity

# Mean Difference 1026
mean_sbe_JPL = JPL.SAL_CTD_MEAN #.mean       #(dim='trajectory')
mean_JPL = JPL.smap_SSS #.mean       #(dim='trajectory')
diff_JPL= mean_sbe_JPL-mean_JPL

p_JPL=diff_JPL.where(diff_JPL > 0, drop=True)
#print((len(p_JPL)/len(diff_JPL))*100) #51.2

# Mean Difference 1060
mean_sbe_JPL_1 = JPL_1.SAL_CTD_MEAN #.mean       #(dim='trajectory')
mean_JPL_1 = JPL_1.smap_SSS #.mean       #(dim='trajectory')
diff_JPL_1= mean_sbe_JPL_1-mean_JPL_1

p_JPL_1=diff_JPL_1.where(diff_JPL_1 > 0, drop=True)
#print((len(p_JPL_1)/len(diff_JPL_1))*100) #61.9

# Mean Difference 1061
mean_sbe_JPL_2 = JPL_2.SAL_CTD_MEAN #.mean       #(dim='trajectory')
mean_JPL_2 = JPL_2.smap_SSS #.mean       #(dim='trajectory')
diff_JPL_2= mean_sbe_JPL_2-mean_JPL_2

p_JPL_2=diff_JPL_2.where(diff_JPL_2 > 0, drop=True)
#print((len(p_JPL_2)/len(diff_JPL_2))*100) #57.4

##percentage of all JPL values from all three saildrones that are less than 0
JPL_all=len(diff_JPL)+len(diff_JPL_1)+len(diff_JPL_2)
JPL_neg=len(p_JPL)+len(p_JPL_1)+len(p_JPL_2)
print(JPL_neg)
print(JPL_all)

JPL_percent=(JPL_neg/JPL_all)*100
print("JPL")
print(JPL_percent)


# Mean Difference
mean_RSS = RSS.smap_SSS#.mean#(dim='trajectory') #.mean(dim='ob')
mean_sbe_RSS = RSS.SAL_CTD_MEAN
mean_RSS_40km = RSS.smap_SSS_40km#.mean#(dim='trajectory')
diff_RSS= mean_sbe_RSS-mean_RSS
diff_RSS_40=mean_sbe_RSS-mean_RSS_40km

p_RSS=diff_RSS.where(diff_RSS > 0, drop=True)
p_RSS_40=diff_RSS_40.where(diff_RSS_40 > 0, drop=True)

#print((len(p_RSS)/len(diff_RSS))*100) #74.06
#print((len(p_RSS_40)/len(diff_RSS_40))*100) #61.3


# Mean Difference
mean_RSS_1 = RSS_1.smap_SSS#.mean#(dim='trajectory') #.mean(dim='ob')
mean_sbe_RSS_1 = RSS_1.SAL_CTD_MEAN
mean_RSS_40km_1 = RSS_1.smap_SSS_40km#.mean#(dim='trajectory')
diff_RSS_1= mean_sbe_RSS_1-mean_RSS_1
diff_RSS_40_1=mean_sbe_RSS_1-mean_RSS_40km_1

p_RSS_1=diff_RSS_1.where(diff_RSS_1 > 0, drop=True)
p_RSS_40_1=diff_RSS_40_1.where(diff_RSS_40_1 > 0, drop=True)
#print((len(p_RSS_1)/len(diff_RSS_1))*100)#72.7
#print((len(p_RSS_40_1)/len(diff_RSS_40_1))*100)#61.72

# Mean Difference
mean_RSS_2 = RSS_2.smap_SSS#.mean#(dim='trajectory') #.mean(dim='ob')
mean_sbe_RSS_2 = RSS_2.SAL_CTD_MEAN
mean_RSS_40km_2 = RSS_2.smap_SSS_40km#.mean#(dim='trajectory')
diff_RSS_2= mean_sbe_RSS_2-mean_RSS_2
diff_RSS_40_2=mean_sbe_RSS_2-mean_RSS_40km_2

p_RSS_2=diff_RSS_2.where(diff_RSS_2 > 0, drop=True)
p_RSS_40_2=diff_RSS_40_2.where(diff_RSS_40_2 > 0, drop=True)
#print(len(p_RSS_2)/len(diff_RSS_2)*100) #72.38
#print(len(p_RSS_40_2)/len(diff_RSS_40_2)*100) #61.04

##percentage of all JPL values from all three saildrones that are less than 0
RSS_all=len(diff_RSS)+len(diff_RSS_1)+len(diff_RSS_2)
RSS_neg=len(p_RSS)+len(p_RSS_1)+len(p_RSS_2)

RSS_percent=(RSS_neg/RSS_all)*100
print("RSS")
print(RSS_percent)

##percentage of all JPL values from all three saildrones that are less than 0
RSS_40all=len(diff_RSS_40)+len(diff_RSS_40_1)+len(diff_RSS_40_2)
RSS_40neg=len(p_RSS_40)+len(p_RSS_40_1)+len(p_RSS_40_2)

RSS_40_percent=(RSS_40neg/RSS_40all)*100
print("RSS_40")
print(RSS_40_percent)

#Mean Difference HYCOM
mean_HYCOM = hysal2
diff_HYCOM=mean_sbe_JPL-mean_HYCOM


#Mean Difference HYCOM
mean_HYCOM_1 = hysal2_1
diff_HYCOM_1=mean_sbe_JPL_1-mean_HYCOM_1


#Mean Difference HYCOM
mean_HYCOM_2 = hysal2_2
diff_HYCOM_2=mean_sbe_JPL_2-mean_HYCOM_2


p_HYCOM=diff_HYCOM.where(diff_HYCOM > 0, drop=True)
p_HYCOM_1=diff_HYCOM_1.where(diff_HYCOM_1 > 0, drop=True)
p_HYCOM_2=diff_HYCOM_2.where(diff_HYCOM_2 > 0, drop=True)

##percentage of all JPL values from all three saildrones that are less than 0
HYCOM_all=len(diff_HYCOM)+len(diff_HYCOM_1)+len(diff_HYCOM_2)
HYCOM_neg=len(p_HYCOM)+len(p_HYCOM_1)+len(p_HYCOM_2)

HYCOM_percent=(HYCOM_neg/HYCOM_all)*100
print("HYCOM")
print(HYCOM_percent)

fig = plt.subplots(figsize=(18, 8))
ax = plt.subplot(2, 2, 1)
plt.scatter(mean_sbe_RSS.values, diff_RSS.values, s=10, color='g', label=str(num))
plt.scatter(mean_sbe_RSS_1.values, diff_RSS_1.values,s=10, color='b', label=str(num_1))
plt.scatter(mean_sbe_RSS_2.values, diff_RSS_2.values, s=10,color='r', label=str(num_2))


'''
#line of best fit
m,b = np.polyfit(mean_sbe_RSS.values, diff_RSS.values,1)
m1,b1 = np.polyfit(mean_sbe_RSS_1.values, diff_RSS_1.values,1)
m2,b2 = np.polyfit(mean_sbe_RSS_2.values, diff_RSS_2.values,1)
plt.plot(mean_sbe_RSS,m*mean_sbe_RSS +b, 'g')
plt.plot(mean_sbe_RSS_1,m*mean_sbe_RSS_1 +b1, 'b')
plt.plot(mean_sbe_RSS_2,m*mean_sbe_RSS_2 +b2, 'r') '''

#plt.title("In situ - RSS70", fontweight='semibold')
ax.set_ylabel("$\Delta$SSS (psu)", fontsize=15, fontweight='semibold')
#ax.set_xlabel("Salinity (psu)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
ax.set_ylim(-2, 2)
leg = ax.legend(loc='upper left', prop=legend_properties)


for h, t in zip(leg.legendHandles, leg.get_texts()):
    t.set_color(h.get_facecolor()[0])

plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
ax.text(.5,.9,"In situ - RSS70",
    horizontalalignment='center',
    transform=ax.transAxes,fontweight='semibold')
#ax.set_title("(a)", loc="left",y=0.9, x=-0.1,fontweight='semibold')
ax.text(33.7, -1, '(a)', color='k', style='normal',fontsize='15',fontweight='semibold')
#plt.title(figure_title, y=1.08)
a = ax.get_ygridlines()
b = a[2]
b.set_color('black')
b.set_linewidth(1.5)


ax1 = plt.subplot(2, 2, 3)
plt.scatter(mean_sbe_JPL.values, diff_JPL.values, s=10,color='g', label=str(j_num))
plt.scatter(mean_sbe_JPL_1.values, diff_JPL_1.values,s=10, color='b', label=str(j_num_1))
plt.scatter(mean_sbe_JPL_2.values, diff_JPL_2.values, s=10,color='r', label=str(j_num_2))

'''
#line of best fit
m,b = np.polyfit(mean_sbe_JPL.values, diff_JPL.values,1)
m1,b1 = np.polyfit(mean_sbe_JPL_1.values, diff_JPL_1.values,1)
m2,b2 = np.polyfit(mean_sbe_JPL_2.values, diff_JPL_2.values,1)
plt.plot(mean_sbe_JPL,m*mean_sbe_JPL +b, 'g')
plt.plot(mean_sbe_JPL_1,m*mean_sbe_JPL_1 +b1, 'b')
plt.plot(mean_sbe_JPL_2,m*mean_sbe_JPL_2 +b2, 'r')'''

#plt.title("In situ - JPL", fontweight='semibold')
ax1.set_ylabel("$\Delta$SSS (psu)", fontsize=15, fontweight='semibold')
ax1.set_xlabel("Salinity (psu)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
ax1.set_ylim(-2, 2)
#plt.legend(loc='upper left', prop=legend_properties)
ax1.tick_params(axis='y')
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
ax1.text(.5,.9,"In situ - JPL",
    horizontalalignment='center',
    transform=ax1.transAxes,fontweight='semibold')
#ax1.set_title("(b)", loc="left",y=0.9, x=-0.1,fontweight='semibold')
ax1.text(33.7, -1, '(b)', color='k', style='normal',fontsize='15',fontweight='semibold')
a = ax1.get_ygridlines()
b = a[2]
b.set_color('black')
b.set_linewidth(1.5)

ax2 =  plt.subplot(2, 2, 2)
ax2.set_ylabel("$\Delta$SSS (psu)", fontsize=15, fontweight='semibold')
#ax2.set_xlabel("Salinity (psu)", fontsize=15, fontweight='semibold')
plt.scatter(mean_sbe_RSS.values, diff_RSS_40.values, s=10,color='g', label=str(num))
plt.scatter(mean_sbe_RSS_1.values, diff_RSS_40_1.values,s=10, color='b', label=str(num_1))
plt.scatter(mean_sbe_RSS_2.values, diff_RSS_40_2.values, s=10,color='r', label=str(num_2))

'''
#line of best fit
m,b = np.polyfit(mean_sbe_RSS.values, diff_RSS_40.values,1)
m1,b1 = np.polyfit(mean_sbe_RSS_1.values, diff_RSS_40_1.values,1)
m2,b2 = np.polyfit(mean_sbe_RSS_2.values, diff_RSS_40_2.values,1)

plt.plot(mean_sbe_RSS,m*mean_sbe_RSS +b, 'g')
plt.plot(mean_sbe_RSS_1,m*mean_sbe_RSS_1 +b1, 'b')
plt.plot(mean_sbe_RSS_2,m*mean_sbe_RSS_2 +b2, 'r')'''


ax2.tick_params(axis='y')
plt.tick_params(axis='both', which='major', labelsize=15)
#plt.legend(loc='upper left', prop=legend_properties)
#plt.title("In situ - RSS40", fontweight='semibold')
ax2.set_ylim(-2, 2)
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
#ax2.set_title("(c)", loc="left",y=0.9, x=-0.1,fontweight='semibold')
ax2.text(33.7, -1, '(c)', color='k', style='normal',fontsize='15',fontweight='semibold')
ax2.text(.5,.9,"In situ - RSS40",  horizontalalignment='center', transform=ax2.transAxes,fontweight='semibold')

a = ax2.get_ygridlines()
b = a[2]
b.set_color('black')
b.set_linewidth(1.5)

ax3 =  plt.subplot(2, 2, 4)
ax3.set_ylabel("$\Delta$SSS (psu)", fontsize=15, fontweight='semibold')
ax3.set_xlabel("Salinity (psu)", fontsize=15, fontweight='semibold')
plt.scatter(mean_sbe_JPL.values, diff_HYCOM.values, s=10,color='g', label=str(num))
plt.scatter(mean_sbe_JPL_1.values, diff_HYCOM_1.values, s=10,color='b', label=str(num_1))
plt.scatter(mean_sbe_JPL_2.values, diff_HYCOM_2.values, s=10,color='r', label=str(num_2))


#line of best fit
#m,b = np.polyfit(mean_sbe_JPL.values, diff_HYCOM.values,1)
#m1,b1 = np.polyfit(mean_sbe_JPL_1.values, diff_HYCOM_1.values,1)
#m2,b2 = np.polyfit(mean_sbe_JPL_2.values, diff_HYCOM_2.values,1)
#plt.plot(mean_sbe_JPL,m*mean_sbe_JPL +b, 'g')
#plt.plot(mean_sbe_JPL_1,m*mean_sbe_JPL_1 +b1, 'b')
#plt.plot(mean_sbe_JPL_2,m*mean_sbe_JPL_2 +b2, 'r')


ax2.tick_params(axis='y')
ax3.set_ylim(-2, 2)
#plt.title("In situ - HYCOM", fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
#plt.legend(loc='upper left', prop=legend_properties)
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
ax3.text(.5,.9,"In situ - HYCOM",
    horizontalalignment='center',
    transform=ax3.transAxes,fontweight='semibold')
#ax3.set_title("d)", loc="left",y=0.9, x=-0.1,fontweight='semibold')
ax3.text(33.7, -1, '(d)', color='k', style='normal',fontsize='15',fontweight='semibold')
a = ax3.get_ygridlines()
b = a[2]
b.set_color('black')
b.set_linewidth(1.5)
plt.show()
