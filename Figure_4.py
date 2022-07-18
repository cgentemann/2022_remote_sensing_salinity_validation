# This script creates the timeseries for salinity and SST observed by the NASA Saildrone(SD1026) and satellite datasets for the entire ATOMIC cruise
import xarray as xr
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import HourLocator, DateFormatter, DayLocator
from matplotlib.ticker import FormatStrFormatter
import os
import datetime
import numpy as np

data_sat= '~/paper_software/2020_ATOMIC_Salinity/data/sss_collocations_orbital_norepeat/'
data_hycom='~/hycom_files/'
PNG = '~/hycom_files/'

# Opening data

#IMPORT HYCOM
files = glob(data_hycom+'*nc4')
hycom=xr.open_mfdataset(data_hycom+'*nc4',concat_dim='time').isel(depth=0)

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
hysal_1=filled3.salinity

##Collocated Data
J_file=glob(data_sat+'*jpl*.nc')
R_file=glob(data_sat+'*rss*.nc')

JPL = xr.open_dataset(J_file[0], decode_times=False)
JPL_1 = xr.open_dataset(J_file[1], decode_times=False)
JPL_2= xr.open_dataset(J_file[2], decode_times=False)

RSS= xr.open_dataset(R_file[0], decode_times=False)
RSS_1 = xr.open_dataset(R_file[1], decode_times=False)
RSS_2 = xr.open_dataset(R_file[2], decode_times=False)

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

ss_times_60 = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test_1]
jp_times_61 = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test2_1]

ss_times_61 = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test_2]
jp_times_26 = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test2_2]

ss_times_26 = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test_3]
jp_times_60 = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test2_3]

RSS.assign_coords(time=ss_times_60)
JPL.assign_coords(time=jp_times_61)

RSS_1.assign_coords(time=ss_times_61)
JPL_1.assign_coords(time=jp_times_26)

RSS_2.assign_coords(time=ss_times_26)
JPL_2.assign_coords(time=jp_times_60)
#RSS=RSS.sel(time=slice('2020-02-17','2020-03-02'))

num=RSS.trajectory.values
num_1=RSS_1.trajectory.values
num_2=RSS_2.trajectory.values
print(num)
print(num_1)
print(num_2)

j=JPL.trajectory.values
j1=JPL_1.trajectory.values
j2=JPL_2.trajectory.values
print("JPL")
print(j)
print(j1)
print(j2)

# interp HYCOM data to saildrone dimensions
hysal = filled3.interp(lat=JPL.lat, longitude=JPL.lon, time=jp_times_61, method='nearest')
hysal_61 = hysal.salinity

# interp HYCOM data to saildrone dimensions
hysal_1 = filled3.interp(lat=JPL_1.lat, longitude=JPL_1.lon, time=jp_times_26, method='nearest')
hysal_26 = hysal_1.salinity

# interp HYCOM data to saildrone dimensions
hysal_2 = filled3.interp(lat=JPL_2.lat, longitude=JPL_2.lon, time=jp_times_60, method='nearest')
hysal_60 = hysal_2.salinity

#salinity data
sal_sbe37_61 = JPL['SAL_CTD_MEAN']
sal_jpl_61 = JPL.smap_SSS
sal_rss40_60 = RSS['smap_SSS_40km']
sal_rss_60 = RSS['smap_SSS']

sal_sbe37_26 = JPL_1['SAL_CTD_MEAN']
sal_jpl_26 = JPL_1['smap_SSS']
sal_rss40_61 = RSS_1['smap_SSS_40km']
sal_rss_61 = RSS_1['smap_SSS']

sal_sbe37_60 = JPL_2['SAL_CTD_MEAN']
sal_jpl_60 = JPL_2['smap_SSS']
sal_rss40_26 = RSS_2['smap_SSS_40km']
sal_rss_26 = RSS_2['smap_SSS']

d_fmt = DateFormatter("%m/%d")
legend_properties = {'weight':'semibold','size':'8'}
print(len(ss_times_26))
print(len(ss_times_61))

fig = plt.subplots(figsize=(28, 20))
ax = plt.subplot(3,1,1)
plt.plot(jp_times_26, sal_sbe37_26.values,color='r',label='SBE37')
plt.scatter(jp_times_26, sal_jpl_26.values,color='orange',label='JPL', s=16)
plt.scatter(ss_times_26, sal_rss40_26.values,color='g',label='RSS40', s=16)
plt.scatter(ss_times_26, sal_rss_26.values,color='c',label='RSS70', s=16)
plt.plot(jp_times_26, hysal_26.values, color='m', label='HYCOM')

ax.xaxis.set_major_formatter(d_fmt)
ax.set_ylim(33, 39)
ax.set_xticklabels(jp_times_26, fontdict=None, minor=False)
#ax.set_ylabel("Salinity (psu)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.legend(loc='center left', prop=legend_properties,  bbox_to_anchor=(1.0, 0.5))

plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
plt.text(0.05,0.05,'(a) SD1026', fontsize=15, color='k', style='normal', fontweight='semibold', transform=ax.transAxes)

ax1 = plt.subplot(3,1,2)
plt.plot(jp_times_60, sal_sbe37_60.values,color='r',label='SBE37')
plt.scatter(jp_times_60, sal_jpl_60.values,color='orange',label='JPL', s=16)
plt.scatter(ss_times_60, sal_rss40_60.values,color='g',label='RSS40', s=16)
plt.scatter(ss_times_60, sal_rss_60.values,color='c',label='RSS70', s=16)
plt.plot(jp_times_60, hysal_60.values, color='m', label='HYCOM')

ax1.xaxis.set_major_formatter(d_fmt)
ax1.set_ylim(33, 39)
#ax1.set_xlabel("Date (mm-dd)", fontsize=15, fontweight='semibold')
ax1.set_ylabel("Salinity (psu)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
ax1.tick_params(axis='y')
#plt.legend(loc='best', prop=legend_properties)
plt.legend(loc='center left', prop=legend_properties,  bbox_to_anchor=(1.0, 0.5))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
plt.text(0.05,0.05,'(b) SD1060', fontsize=15, color='k', style='normal', fontweight='semibold', transform=ax1.transAxes)

ax2 = plt.subplot(3,1,3)
plt.plot(jp_times_61, sal_sbe37_61.values,color='r',label='SBE37')
plt.scatter(jp_times_61, sal_jpl_61.values,color='orange',label='JPL', s=16)
plt.scatter(ss_times_61, sal_rss40_61.values,color='g',label='RSS40', s=16)
plt.scatter(ss_times_61, sal_rss_61.values,color='c',label='RSS70', s=16)
plt.plot(jp_times_61, hysal_61.values, color='m', label='HYCOM')

ax2.xaxis.set_major_formatter(d_fmt)
ax2.set_ylim(33, 39)
ax2.set_xlabel("Date (mm-dd)", fontsize=15, fontweight='semibold')
#ax2.set_ylabel("Salinity (PSU)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
ax2.tick_params(axis='y')
#plt.legend(loc='best', prop=legend_properties)
plt.legend(loc='center left', prop=legend_properties,  bbox_to_anchor=(1.0, 0.5))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
plt.text(0.05,0.05,'(c) SD1061', fontsize=15, color='k', style='normal', fontweight='semibold', transform=ax2.transAxes)

plt.savefig(PNG +'Salinityfour.png'
            , dpi=150, facecolor='w', edgecolor='w',
            orientation='lanscape', papertype=None, format='png',
              bbox_inches='tight', pad_inches=0.1)
plt.show()

exit()
