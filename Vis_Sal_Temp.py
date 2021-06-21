import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import cartopy.crs as ccrs                   # import projections
import cartopy.feature as cf                 # import features
from glob import glob
from matplotlib.dates import HourLocator, DateFormatter, DayLocator
from matplotlib.ticker import FormatStrFormatter
from math import sin, cos, sqrt, atan2, radians, pi


def merge(list1, list2):
    merged_list = [(list1[i], list2[i]) for i in range(0, len(list1))]
    return merged_list

def calculate_distance(positions):
    results = []
    for i in range(len(positions)):
        if i == 0 :
            loc1 = positions[i]
            loc2 = positions[i]
        else:
            loc1 = positions[i-1]
            loc2 = positions[i]

        lat1 = loc1[0]
        lng1 = loc1[1]

        lat2 = loc2[0]
        lng2 = loc2[1]

        degreesToRadians = (pi / 180)
        latrad1 = lat1 * degreesToRadians
        latrad2 = lat2 * degreesToRadians
        dlat = (lat2 - lat1) * degreesToRadians
        dlng = (lng2 - lng1) * degreesToRadians

        a = sin(dlat / 2) * sin(dlat / 2) + cos(latrad1) * \
        cos(latrad2) * sin(dlng / 2) * sin(dlng / 2)
        c = 2 * atan2(sqrt(a), sqrt(1 - a))
        r = 6371000

        results.append((r * c)/1000)
        
    return results #(sum(results) / 1000)

data_dir = '/home/alton/Github/paper_software/2020_ATOMIC_Salinity/data/'
adcp_dir = '/home/alton/Github/paper_software/2020_ATOMIC_Salinity/data/with_adcp/'
PNG = '/home/alton/Github/paper_software/2020_ATOMIC_Salinity/figures/'
# saildrone_filenames = glob(data_dir+'saildrone*.nc')
adcp_files = glob(adcp_dir+'combined*.nc')


sd = xr.open_dataset(adcp_files[2])
# sd = sd.isel(trajectory=0).swap_dims({'obs': 'time'})
sd['curr_spd'] = np.sqrt(sd.adcp_vel_east**2 + sd.adcp_vel_north**2)
sd['wspd']=np.sqrt(sd.UWND_MEAN**2+sd.VWND_MEAN**2)


# print(sd.latitude.values)
val = merge(sd.latitude.values, sd.longitude.values)
distan = np.cumsum(calculate_distance(val)) + 1574.29
times = pd.date_range("2020-01-31 00:00:00", freq="5T", periods=9216)
sd['distance'] = xr.DataArray(distan, coords=[times], dims=["time"])

#Salinity and Temperature Info
dist = sd['distance'].sel(time=slice('2020-02-15','2020-02-22'))
sal_rbr = sd['SAL_RBR_MEAN'].sel(time=slice('2020-02-15','2020-02-22'))
sal_sbe37 = sd['SAL_SBE37_MEAN'].sel(time=slice('2020-02-15','2020-02-22'))
air_temp = sd['TEMP_AIR_MEAN'].sel(time=slice('2020-02-15','2020-02-22'))
temp_sbe37 = sd['TEMP_SBE37_MEAN'].sel(time=slice('2020-02-15','2020-02-22'))
d_fmt = DateFormatter("%m/%d")   
legend_properties = {'weight':'semibold','size':'10'}

#Current and Wind Speed Info
cs = sd['curr_spd'].sel(time=slice('2020-02-15','2020-02-22'))
u_cs = sd['adcp_vel_east'].sel(time=slice('2020-02-15','2020-02-22'))
v_cs = sd['adcp_vel_north'].sel(time=slice('2020-02-15','2020-02-22'))
cs_lat = sd['latitude'].sel(time=slice('2020-02-15','2020-02-22'))
cs_lon = sd['longitude'].sel(time=slice('2020-02-15','2020-02-22'))

ws = sd['wspd'].sel(time=slice('2020-02-15','2020-02-22'))
u_ws = sd['UWND_MEAN'].sel(time=slice('2020-02-15','2020-02-22'))
v_ws = sd['VWND_MEAN'].sel(time=slice('2020-02-15','2020-02-22'))
ws_lat = sd['latitude'].sel(time=slice('2020-02-15','2020-02-22'))
ws_lon = sd['longitude'].sel(time=slice('2020-02-15','2020-02-22'))

fig = plt.subplots(figsize=(18,10))
ax = plt.subplot(2,2,1)
# plt.plot(sal_rbr.time.values, sal_rbr.values,color='b',label='RBR')
plt.plot(sal_sbe37.time.values, sal_sbe37.values,color='r',label='SBE37')
ax.xaxis.set_major_formatter(d_fmt)
#ax.set_xlabel("Time mm-dd (UTC)", fontsize=15, fontweight='semibold')
ax.set_ylabel("Salinity (psu)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.legend(loc='best', prop=legend_properties)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')

ax1 = plt.subplot(2,2,3)
plt.plot(air_temp.time.values, air_temp.values,color='r')
ax1.xaxis.set_major_formatter(d_fmt)
ax1.set_xlabel("Date (mm-dd)", fontsize=15, fontweight='semibold')
ax1.set_ylabel("Air Temperature", color='r', fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
ax1.tick_params(axis='y', labelcolor='r')
ax1.set_ylim(min(air_temp.values), max(temp_sbe37.values))
#ax1.set_ylim(26.5, max(temp_sbe37.values))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')

ax2 = ax1.twinx()
ax2.set_ylabel("Sea \n Surface \n Temperature", color='b', fontsize=15, fontweight='semibold')
plt.plot(temp_sbe37.time.values, temp_sbe37.values,color='b')
ax2.tick_params(axis='y', labelcolor='b')
#ax2.set_ylim(min(air_temp.values), max(temp_sbe37.values))
plt.ylim([26.5, max(temp_sbe37.values)])
plt.tick_params(axis='both', which='major', labelsize=15)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')


ax3 = ax1.twiny()
ax3.set_xticks(dist.values)
#plt.plot(dist.values,air_temp.values)
#newlabel = np.linspace(min(dist.values),max(dist.values),12)
newlabel = np.linspace(3600,4600,11)
ax3.set_xticks(newlabel)
ax3.set_xticklabels(newlabel)

ax3.xaxis.set_ticks_position('bottom')
ax3.xaxis.set_label_position('bottom')
ax3.spines['bottom'].set_position(('outward', 50))
ax3.set_xlabel('Distance Travelled (km)', fontsize=15, fontweight='semibold')
ax3.set_xlim(min(dist.values), max(dist.values))
plt.tick_params(axis='both', which='major', labelsize=15)
ax3.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.xticks(fontweight='semibold')

ax4 = plt.subplot(2,2,2)
wx = ws.time.values
wx = wx[::20]
wy = np.ones(len(wx)) * 8
dx, dy = u_ws[::30], v_ws[::30]
plt.quiver(wx, wy, dx, dy, pivot='tail', color='g',linewidths=0.005)
# plt.barbs(wx, wy, dx, dy, barbcolor=['b'])
# plt.quiver(ws.time.values, ws.values, u_ws.values, v_ws.values, 
#            pivot='tail', color='g',linewidths=0.005)
#plt.barbs(ws.time.values, ws.values*1.944,u_ws.values,v_ws.values, barbcolor=['b'])
#plt.plot(sal_sbe37.time.values, sal_sbe37.values,color='r',label='SBE37')
ax4.set_ylim(min(ws.values), max(ws.values))
ax4.xaxis.set_major_formatter(d_fmt)
ax4.set_ylabel("Wind Speed (m/s)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')



plt.savefig(PNG +'sal_temp_vs_datetime_distance_travelled_subplot_sst_min.png', 
            dpi=300, facecolor='w', edgecolor='w',
            orientation='lanscape', papertype=None, format='png',
            bbox_inches='tight', pad_inches=0.1)

plt.show()