#!/usr/bin/env python

##This script creates plots that salinity vs salinity difference for saildrone and each SSS platform, Figure 3
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import datetime
import sys
import cftime
from decimal import Decimal

#data directory for saildrone and satellite orbital data
data_sat = 'C:/Users/intern-1/Documents/GitHub/paper_software/2020_ATOMIC_Salinity/data/sss_collocations_orbital_norepeat/'
data_sat2 = 'C:/Users/intern-1/Documents/GitHub/paper_software/2020_ATOMIC_Salinity/data/sss_collocations_8day_norepeat/'

#data directory for HYCOM data
data_dir1 = 'C:/Users/intern-1/Documents/hycom_files/'
files = glob(data_dir1+'*nc4')
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
JPL = glob(data_sat+'*jpl*.nc')
RSS = glob(data_sat+'*rss*.nc')

legend_properties = {'weight': 'semibold', 'size': '12'}

#Saildrone 1026
JPL=xr.open_dataset(JPL[0],decode_times=False)
RSS=xr.open_dataset(RSS[0],decode_times=False)

# Convert times from nano seconds to seconds
ns = 1e-9
rss_time = RSS.time.values * ns
jpl_time = JPL.time.values * ns

test=rss_time.astype(np.float)
test2=jpl_time.astype(np.float)

# Swap Dimensions
RSS = RSS.swap_dims({'ob': 'time'})
JPL = JPL.swap_dims({'ob': 'time'})

ss_times = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test]
jp_times = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test2]
RSS.assign_coords(time=ss_times)
JPL.assign_coords(time=jp_times)

# interp HYCOM data to saildrone dimensions
hysal = filled3.interp(lat=JPL.lat, longitude=JPL.lon, time=jp_times, method='nearest')
hysal2 = hysal.salinity

# Mean Difference
mean_sbe_JPL = JPL.SAL_CTD_MEAN #.mean       #(dim='trajectory')
mean_JPL = JPL.smap_SSS #.mean       #(dim='trajectory')
diff_JPL= mean_sbe_JPL-mean_JPL

# Mean Difference
mean_RSS = RSS.smap_SSS#.mean#(dim='trajectory') #.mean(dim='ob')
mean_sbe_RSS = RSS.SAL_CTD_MEAN
mean_RSS_40km = RSS.smap_SSS_40km#.mean#(dim='trajectory')
diff_RSS= mean_sbe_RSS-mean_RSS
diff_RSS_40=mean_sbe_RSS-mean_RSS_40km

#Mean Difference HYCOM
mean_HYCOM = hysal2
diff_HYCOM=mean_sbe_JPL-mean_HYCOM

#Create plot
fig = plt.subplots(figsize=(18, 10))
ax = plt.subplot(2, 2, 1)
plt.scatter(mean_sbe_RSS.values, diff_RSS.values, color='b', label='SBE37 - RSS70')
ax.set_ylabel("Difference", fontsize=15, fontweight='semibold')
#ax.set_xlabel("Salinity (psu)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
ax.set_ylim(-2, 2)
plt.legend(loc='upper left', prop=legend_properties)
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
#ax.set_title("a)", loc="left",y=0.9, x=-0.1,fontweight='semibold')
ax.text(33.7, 1, 'a)', color='k', style='normal',fontsize='15',fontweight='semibold')
#plt.title(figure_title, y=1.08)
a = ax.get_ygridlines()
b = a[2]
b.set_color('black')
b.set_linewidth(1.5)


ax1 = plt.subplot(2, 2, 3)
plt.scatter(mean_sbe_JPL.values, diff_JPL.values, color='b', label='SBE37 - JPL')
ax1.set_ylabel("Difference", fontsize=15, fontweight='semibold')
ax1.set_xlabel("Salinity (psu)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
ax1.set_ylim(-2, 2)
plt.legend(loc='upper left', prop=legend_properties)
ax1.tick_params(axis='y')
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
#ax1.set_title("b)", loc="left",y=0.9, x=-0.1,fontweight='semibold')
ax1.text(33.7, 1, 'b)', color='k', style='normal',fontsize='15',fontweight='semibold')
a = ax1.get_ygridlines()
b = a[2]
b.set_color('black')
b.set_linewidth(1.5)

ax2 =  plt.subplot(2, 2, 2)
ax2.set_ylabel("Difference", fontsize=15, fontweight='semibold')
#ax2.set_xlabel("Salinity (psu)", fontsize=15, fontweight='semibold')
plt.scatter(mean_sbe_RSS.values, diff_RSS_40.values, color='b', label='SBE37 - RSS40')
ax2.tick_params(axis='y')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.legend(loc='upper left', prop=legend_properties)
ax2.set_ylim(-2, 2)
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
#ax2.set_title("c)", loc="left",y=0.9, x=-0.1,fontweight='semibold')
ax2.text(33.7, 1, 'c)', color='k', style='normal',fontsize='15',fontweight='semibold')
a = ax2.get_ygridlines()
b = a[2]
b.set_color('black')
b.set_linewidth(1.5)

ax3 =  plt.subplot(2, 2, 4)
ax3.set_ylabel("Difference", fontsize=15, fontweight='semibold')
ax3.set_xlabel("Salinity (psu)", fontsize=15, fontweight='semibold')
plt.scatter(mean_sbe_JPL.values, diff_HYCOM.values, color='b', label='SBE37 - HYCOM' )
ax2.tick_params(axis='y')
ax3.set_ylim(-2, 2)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.legend(loc='upper left', prop=legend_properties)
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
#ax3.set_title("d)", loc="left",y=0.9, x=-0.1,fontweight='semibold')
ax3.text(33.7, 1, 'd)', color='k', style='normal',fontsize='15',fontweight='semibold')
a = ax3.get_ygridlines()
b = a[2]
b.set_color('black')
b.set_linewidth(1.5)

plt.show()
