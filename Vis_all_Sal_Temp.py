import xarray as xr
import numpy as np
import pandas as pd
import datetime 
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import cartopy.crs as ccrs                   # import projections
import cartopy.feature as cf                 # import features
import cmocean.cm as cmo
from glob import glob
from matplotlib.dates import HourLocator, DateFormatter, DayLocator
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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
sat_dir = data_dir + 'sss_collocations_orbital_norepeat/'
adcp_dir = data_dir + 'with_adcp/'
PNG = '/home/alton/Github/paper_software/2020_ATOMIC_Salinity/figures/'

#IMPORT HYCOM DATA
#data directory for HYCOM data
hycom_dir1 = '/home/alton/Saildrone/cimh_saildrone/data/Sanola_Hycomm/'
hycom=xr.open_mfdataset(hycom_dir1+'*nc4',concat_dim='time').isel(depth=0)

files = glob(sat_dir+'saildrone*.nc')
adcp_files = glob(adcp_dir+'combined*.nc')


_, index = np.unique(hycom['time'], return_index=True)
hycom2=hycom.isel(time=index)

#Open Files
rss = xr.open_dataset(files[5], decode_times=False)
jpl = xr.open_dataset(files[4], decode_times=False)
sd = xr.open_dataset(adcp_files[2])
sd['curr_spd'] = np.sqrt(sd.adcp_vel_east**2 + sd.adcp_vel_north**2)
sd['wspd']=np.sqrt(sd.UWND_MEAN**2+sd.VWND_MEAN**2)
# print(files[4])
# print(files[5])
print(rss)
#Swap Dimensions
rss = rss.swap_dims({'ob': 'time'})
jpl = jpl.swap_dims({'ob': 'time'})

#Convert times from nano seconds to seconds 
ns=1e-9
rss.time.values = rss.time.values * ns
jpl.time.values = jpl.time.values * ns
# sv.time.values = sv.time.values * ns
# sw.time.values = sw.time.values * ns

ss_times = [ datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in rss.time.values ]
jp_times = [ datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in jpl.time.values ]

#change hycom coordinates to match with saildrone(0-359:-180-179)
hycom_lon=hycom2.assign_coords(longitude=(((hycom2.lon + 180) % 360) - 180))
hycom2=hycom_lon.swap_dims({'lon':'longitude'})

#remove nans from hycom data
filled=hycom2.chunk({'time':-1}).interpolate_na(dim="time",method="nearest",fill_value='extrapolate')
filled2=filled.interpolate_na(dim="lat", method="nearest",fill_value='extrapolate')
filled3=filled2.interpolate_na(dim="lon", method="nearest",fill_value='extrapolate')

val = merge(sd.latitude.values, sd.longitude.values)
distan = np.cumsum(calculate_distance(val)) + 1574.29
times = pd.date_range("2020-01-31 00:00:00", freq="5T", periods=9216)
sd['distance'] = xr.DataArray(distan, coords=[times], dims=["time"])

dist = sd['distance'].sel(time=slice('2020-02-15','2020-02-22'))
sal_rbr = sd['SAL_RBR_MEAN'].sel(time=slice('2020-02-15','2020-02-22'))
sal_sbe37 = sd['SAL_SBE37_MEAN'].sel(time=slice('2020-02-15','2020-02-22'))
air_temp = sd['TEMP_AIR_MEAN'].sel(time=slice('2020-02-15','2020-02-22'))
temp_sbe37 = sd['TEMP_SBE37_MEAN'].sel(time=slice('2020-02-15','2020-02-22'))

# interp HYCOM data to saildrone dimensions
hysal = filled3.interp(lat=sal_sbe37.latitude, longitude=sal_sbe37.longitude, time=sal_sbe37.time, method='nearest')
hysal2 = hysal.salinity

ws = sd['wspd'].sel(time=slice('2020-02-15','2020-02-22'))
u_ws = sd['UWND_MEAN'].sel(time=slice('2020-02-15','2020-02-22'))
v_ws = sd['VWND_MEAN'].sel(time=slice('2020-02-15','2020-02-22'))
ws_lat = sd['latitude'].sel(time=slice('2020-02-15','2020-02-22'))
ws_lon = sd['longitude'].sel(time=slice('2020-02-15','2020-02-22'))

cs = sd['curr_spd'].sel(time=slice('2020-02-15','2020-02-22'))
u_cs = sd['adcp_vel_east'].sel(time=slice('2020-02-15','2020-02-22'))
v_cs = sd['adcp_vel_north'].sel(time=slice('2020-02-15','2020-02-22'))
cs_lat = sd['latitude'].sel(time=slice('2020-02-15','2020-02-22'))
cs_lon = sd['longitude'].sel(time=slice('2020-02-15','2020-02-22'))

d_fmt = DateFormatter("%m-%d")   
legend_properties = {'weight':'semibold','size':'10'}
plt.rcParams['font.weight']='semibold'

fig = plt.subplots(figsize=(15,15))

# Salinity Plot
ax = plt.subplot(4,1,1)
# plt.plot(sal_rbr.time.values, sal_rbr.values,color='b',label='RBR')
plt.plot(sal_sbe37.time.values, sal_sbe37.values,color='r',label='SBE37')
plt.plot(ss_times[134:174], rss.smap_SSS.values[134:174],color='c',label='RSS70')
plt.plot(ss_times[134:174], rss.smap_SSS_40km.values[134:174],color='g',label='RSS40')
plt.plot(jp_times[96:127], jpl.smap_SSS.values[96:127],color='b',label='JPL')
plt.text(sal_sbe37.time.values[1], 34, 'a)', color='k', style='normal',fontsize='15')
plt.plot(hysal2.time.values, hysal2.values, color='m', label='HYCOM')
# ax.xaxis.set_major_formatter(d_fmt)
#ax.set_xlabel("Time mm-dd (UTC)", fontsize=15, fontweight='semibold')
ax.set_ylabel("Salinity (psu)", fontsize=15, fontweight='semibold')
# ax.set_ylim(min(air_temp.values), max(temp_sbe37.values))
plt.tick_params(axis='both', which='major', labelsize=15)
plt.legend(loc='lower right', prop=legend_properties)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')

# #SST and Air Temperature Plot
ax1 = plt.subplot(4,1,2)
plt.plot(air_temp.time.values, air_temp.values,color='r')
plt.text(air_temp.time.values[1], 25, 'b)', color='k', style='normal',fontsize='15')
ax1.xaxis.set_major_formatter(d_fmt)
# ax1.set_xlabel("Date (mm-dd)", fontsize=15, fontweight='semibold')
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
ax2.set_ylabel("Sea Surface Temperature", color='b', fontsize=15, fontweight='semibold')
plt.plot(temp_sbe37.time.values, temp_sbe37.values,color='b')
# plt.text(temp_sbe37.time.values[-1], 25, 'b)', color='k', style='normal',fontsize='15')
ax2.tick_params(axis='y', labelcolor='b')
#ax2.set_ylim(min(air_temp.values), max(temp_sbe37.values))
plt.ylim([26.5, max(temp_sbe37.values)])
plt.tick_params(axis='both', which='major', labelsize=15)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')

# ax3 = ax1.twiny()
# ax3.set_xticks(dist.values)
# #plt.plot(dist.values,air_temp.values)
# #newlabel = np.linspace(min(dist.values),max(dist.values),12)
# newlabel = np.linspace(3600,4600,11)
# ax3.set_xticks(newlabel)
# ax3.set_xticklabels(newlabel)

# ax3.xaxis.set_ticks_position('bottom')
# ax3.xaxis.set_label_position('bottom')
# ax3.spines['bottom'].set_position(('outward', 50))
# ax3.set_xlabel('Distance Travelled (km)', fontsize=15, fontweight='semibold')
# ax3.set_xlim(min(dist.values), max(dist.values))
# plt.tick_params(axis='both', which='major', labelsize=15)
# ax3.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
# plt.xticks(fontweight='semibold')

#Wind Speed Vectors
ax4 = plt.subplot(4,1,3)
wx = ws.time.values
wx = wx[::20]
wy = np.ones(len(wx)) * 8
dx, dy = u_ws[::30], v_ws[::30]
m = np.sqrt(np.power(dx, 2) + np.power(dy, 2))
q = plt.quiver(wx, wy, dx, dy, m/max(m), cmap=cmo.thermal, pivot='tail', width=0.0011)
plt.text(wx[1], 6, 'c)', color='k', style='normal',fontsize='15')
ax4.set_ylim(min(ws.values), max(ws.values))
ax4.xaxis.set_major_formatter(d_fmt)
ax4.set_ylabel("Wind Speed (m/s)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
cbaxes = inset_axes(ax4, width='30%', height='3%', loc='upper right')
plt.colorbar(cax=cbaxes, ticks=[0.5, 0.6, 0.7, 0.8, 0.9, 1.0], orientation='horizontal')

# al = [np.linspace(min(m/max(m)).values,max(m/max(m)).values,6)]
# print(al)

#Current Speed Vectors
ax5 = plt.subplot(4,1,4)
cx = cs.T[0].time.values
cx = cx[::20]
cy = np.ones(len(cx)) * 0.40
dx, dy = u_cs.T[0][::20], v_cs.T[0][::20]
m1 = np.sqrt(np.power(dx, 2) + np.power(dy, 2))
q = plt.quiver(cx, cy, dx, dy, m1/max(m1), cmap=cmo.thermal, pivot='tail', width=0.0011)
plt.text(cx[1], 0.2, 'd)', color='k', style='normal',fontsize='15')
ax5.set_ylim(min(cs.T[0].values), max(cs.T[0].values))
ax5.xaxis.set_major_formatter(d_fmt)
ax5.set_xlabel("Date (mm-dd)", fontsize=15, fontweight='semibold')
ax5.set_ylabel("Current Speed (m/s)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
cbaxes = inset_axes(ax5, width='30%', height='3%', loc='upper right')
plt.colorbar(cax=cbaxes, ticks=[0.0, 0.25, 0.50, 0.75, 1.0], orientation='horizontal')

ax6 = ax5.twiny()
ax6.set_xticks(dist.values)
newlabel = np.linspace(3600,4600,11)
ax6.set_xticks(newlabel)
ax6.set_xticklabels(newlabel)
ax6.xaxis.set_ticks_position('bottom')
ax6.xaxis.set_label_position('bottom')
ax6.spines['bottom'].set_position(('outward', 50))
ax6.set_xlabel('Distance Travelled (km)', fontsize=15, fontweight='semibold')
ax6.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax6.set_xlim(min(dist.values), max(dist.values))
plt.tick_params(axis='both', which='major', labelsize=15)
ax6.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.xticks(fontweight='semibold')


plt.savefig(PNG +' Salinity_Temperature_Time_Series.png'
            , dpi=150, facecolor='w', edgecolor='w',
            orientation='lanscape', papertype=None, format='png',
            bbox_inches='tight', pad_inches=0.1)

plt.show()

# print(rss)

# for i in range(len(jp_times)):
#     print(i, jp_times[i])

# rss['new_times'] = xr.DataArray(ss_times)

# print(rss.smap_SSS)

#Changes the time dimensions for each variable from nano seconds to datetime
# for i in list(rss.keys()):
#     rss[i] = xr.DataArray(rss[i].values, coords=[ss_times], dims=['time'])
# rss.assign_coords({"new_times": ss_times})
# rss.assign_coords(new_times = ss_times)
#This slice will not work because time is still reperesented in nanoseconds therefore I need to replace the old rss.time with ss_times
# rss_sal = rss['smap_SSS'].sel(time=slice(0,1))

# print(rss.smap_SSS.values)
# print(len(rss.time.values))
# new = list(rss.keys())
# print(new)
# print(rss)
# print(rss.coords)
# fig = plt.subplots(figsize=(18,10))
# ax = plt.subplot(2,1,1)
# # # plt.plot(sal_rbr.time.values, sal_rbr.values,color='b',label='RBR')
# plt.plot(ss_times[96:-1], rss.smap_SSS.values[96:-1],color='r',label='RSS')
# plt.plot(jp_times[134:-1], jpl.smap_SSS.values[134:-1],color='b',label='JPL')


# plt.show()
# print(rss_sal)

