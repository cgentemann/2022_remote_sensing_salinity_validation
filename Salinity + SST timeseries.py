# This script creates the timeseries for salinity and SST observed by the NASA Saildrone(SD1026) and satellite datasets for the entire ATOMIC cruise
import xarray as xr
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import HourLocator, DateFormatter, DayLocator
from matplotlib.ticker import FormatStrFormatter
import os
import datetime


data_dir = ''
PNG = 'Timeseries/'

# Opening data
saildrone_filenames = os.path.join(data_dir, 'saildrone-gen_5-atomic_eurec4a_2020-sd1026-20200117T000000-20200302T235959-1_minutes-v1.1595708344687.nc')
jpl_filenames = os.path.join(data_dir, 'saildrone-gen_5-atomic_eurec4a_2020-sd1026-20200117T000000-20200302T235959-1_minutes-v1.1589306725934jplv05.0_orbitalnorep.nc')
rss_filenames = os.path.join(data_dir, 'saildrone-gen_5-atomic_eurec4a_2020-sd1026-20200117T000000-20200302T235959-1_minutes-v1.1589306725934rssv04.0_orbitalnorep.nc')

sd = xr.open_dataset(saildrone_filenames)
ss = xr.open_dataset(jpl_filenames, decode_times=False)
st = xr.open_dataset(rss_filenames, decode_times=False)

# Converting to nanoseconds
ns=1e-9
ss.time.values = ss.time.values * ns
st.time.values = st.time.values * ns
ss_times = [ datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in ss.time.values ]
st_times = [ datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in st.time.values ]

# Switching dimensions
sd = sd.isel(trajectory=0).swap_dims({'obs': 'time'})
ss = ss.swap_dims({'ob': 'time'})
st = st.swap_dims({'ob': 'time'})

# Assigning salinity data
sal_sbe37 = sd['SAL_SBE37_MEAN']
sal_jpl_orbital=ss['smap_SSS']
sal_rss40_orbital=st['smap_SSS_40km']
sal_rss70_orbital=st['smap_SSS']

# Assigning SST data
seatemp_sbe37= sd['TEMP_SBE37_MEAN']
seatemp_jpl_orb = ss['TEMP_CTD_MEAN']
seatemp_rss_orb = st['TEMP_CTD_MEAN']

d_fmt = DateFormatter("%m/%d")
legend_properties = {'weight':'semibold','size':'10'}

fig = plt.subplots(figsize=(18, 10))
ax = plt.subplot(2,1,1)
plt.plot(sal_sbe37.time.values, sal_sbe37.values,color='r',label='SBE37')
plt.plot(ss_times, sal_jpl_orbital.values,color='b',label='JPL')
plt.plot(st_times, sal_rss40_orbital.values,color='g',label='RSS40')
plt.plot(st_times, sal_rss70_orbital.values,color='c',label='RSS70')
ax.xaxis.set_major_formatter(d_fmt)

ax.set_ylabel("Salinity (psu)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.legend(loc='best', prop=legend_properties)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
plt.text(0.05,0.05,'a)', fontsize=15, color='k', style='normal', fontweight='semibold', transform=ax.transAxes)

ax1 = plt.subplot(2,1,2)
plt.plot(seatemp_sbe37.time.values, seatemp_sbe37.values,color='r',label='SBE37')
plt.plot(ss_times, seatemp_jpl_orb.values,color='b',label='JPL')
plt.plot(st_times, seatemp_rss_orb.values,color='g',label='RSS')
ax1.xaxis.set_major_formatter(d_fmt)
ax1.set_xlabel("Date (mm-dd)", fontsize=15, fontweight='semibold')
ax1.set_ylabel("Sea Surface \n Temperature ($^\circ$C)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
ax1.tick_params(axis='y')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
plt.text(0.05,0.05,'b)', fontsize=15, color='k', style='normal', fontweight='semibold', transform=ax1.transAxes)
plt.savefig(PNG +'Salinity and SST timeseries .png'
             , dpi=150, facecolor='w', edgecolor='w',
            orientation='lanscape', papertype=None, format='png',
              bbox_inches='tight', pad_inches=0.1)
plt.show()

exit()
