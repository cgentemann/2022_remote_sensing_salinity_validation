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
adcp_files = glob(adcp_dir + 'combined*.nc')

d_fmt = DateFormatter("%m-%d")   
legend_properties = {'weight':'semibold','size':'10'}

print(adcp_files)
sd = xr.open_dataset(adcp_files[2])

#sd = sd.isel(trajectory=0).swap_dims({'obs': 'time'})
sd['curr_spd'] = np.sqrt(sd.adcp_vel_east**2 + sd.adcp_vel_north**2)
sd['wspd']=np.sqrt(sd.UWND_MEAN**2+sd.VWND_MEAN**2)
#curr = sd['curr_spd'].sel(time=slice('2020-02-15','2020-02-22'))
#sd['wspd']=np.sqrt(sd.UWND_MEAN**2+sd.VWND_MEAN**2)
#
#
val = merge(sd.latitude.values, sd.longitude.values)
distan = np.cumsum(calculate_distance(val)) + 1574.29
times = pd.date_range("2020-01-31 00:00:00", freq="5T", periods=9216)
sd['distance'] = xr.DataArray(distan, coords=[times], dims=["time"])
##
dist = sd['distance'].sel(time=slice('2020-02-15','2020-02-22'))
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
ax = plt.subplot(2,1,1)
wx = ws.time.values
wx = wx[::20]
wy = np.ones(len(wx)) * 0.5
dx, dy = u_ws[::30], v_ws[::30]
plt.quiver(wx, wy, dx, dy, pivot='tail', color='g',linewidths=0.005)
# plt.barbs(wx, wy, dx, dy, barbcolor=['b'])
# plt.quiver(ws.time.values, ws.values, u_ws.values, v_ws.values, 
#            pivot='tail', color='g',linewidths=0.005)
#plt.barbs(ws.time.values, ws.values*1.944,u_ws.values,v_ws.values, barbcolor=['b'])
#plt.plot(sal_sbe37.time.values, sal_sbe37.values,color='r',label='SBE37')
ax.xaxis.set_major_formatter(d_fmt)
ax.set_ylabel("Wind Speed (m/s)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')

ax1 = plt.subplot(2,1,2)
#plt.quiver(ws_lon.time.values, ws.values,u.values,v.values, pivot='tail', color='g',linewidths=1)
cx = cs.T[0].time.values
cx = cx[::20]
cy = np.ones(len(cx)) * 0.5
dx, dy = u_cs.T[0][::20], v_cs.T[0][::20]
plt.quiver(cx, cy, dx, dy, pivot='tail', color='g',linewidths=0.005)
# plt.quiver(cs.T[0].time.values, cs.T[0].values,u_cs.T[0].values,
#            v_cs.T[0].values, pivot='tail', color='g',linewidths=0.005)
#plt.barbs(cs.T[0].time.values, cs.T[0].values,u_cs.T[0].values,
#           v_cs.T[0].values, color='g',linewidths=0.005)
ax1.xaxis.set_major_formatter(d_fmt)
ax1.set_xlabel("Date (mm-dd)", fontsize=15, fontweight='semibold')
ax1.set_ylabel("Current Speed (m/s)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')

ax3 = ax1.twiny()
ax3.set_xticks(dist.values)
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

plt.savefig(PNG +'Wind_Curr_Vectors.png'
             , dpi=300, facecolor='w', edgecolor='w',
            orientation='lanscape', papertype=None, format='png',
              bbox_inches='tight', pad_inches=0.1)

plt.show()
