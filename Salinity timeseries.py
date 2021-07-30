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

saildrone_filenames = os.path.join(data_dir, 'saildrone-gen_5-atomic_eurec4a_2020-sd1026-20200117T000000-20200302T235959-1_minutes-v1.1595708344687.nc')
jpl_filenames = os.path.join(data_dir, 'saildrone-gen_5-atomic_eurec4a_2020-sd1026-20200117T000000-20200302T235959-1_minutes-v1.1589306725934jplv05.0_orbitalnorep.nc')
rss_filenames = os.path.join(data_dir, 'saildrone-gen_5-atomic_eurec4a_2020-sd1026-20200117T000000-20200302T235959-1_minutes-v1.1589306725934rssv04.0_orbitalnorep.nc')

jpl_filename = os.path.join(data_dir,'saildrone-gen_5-atomic_eurec4a_2020-sd1026-20200117T000000-20200302T235959-1_minutes-v1.1589306725934_JPLv5.0_8dy_20210511norep_20210511.nc')
rss_filename = os.path.join(data_dir,'saildrone-gen_5-atomic_eurec4a_2020-sd1026-20200117T000000-20200302T235959-1_minutes-v1.1589306725934_RSSv4.0_8dy_20210511norep_20210511.nc')


sd = xr.open_dataset(saildrone_filenames)
ss = xr.open_dataset(jpl_filenames, decode_times=False)
st = xr.open_dataset(rss_filenames, decode_times=False)
sv = xr.open_dataset(jpl_filename, decode_times=False)
sw = xr.open_dataset(rss_filename, decode_times=False)

ns=1e-9
ss.time.values = ss.time.values * ns
st.time.values = st.time.values * ns
sv.time.values = sv.time.values * ns
sw.time.values = sw.time.values * ns

ss_times = [ datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in ss.time.values ]
st_times = [ datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in st.time.values ]
sv_times = [ datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in sv.time.values ]
sw_times = [ datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in sw.time.values ]

# ss_times = [ datetime.datetime(2020, 1, 17, 0, 48) + datetime.timedelta(seconds=s) for s in ss.time.values ]
# st_times = [ datetime.datetime(2020, 1, 17, 4) + datetime.timedelta(seconds=s) for s in st.time.values ]
# sv_times = [ datetime.datetime(2020, 1, 17, 0, 8) + datetime.timedelta(seconds=s) for s in sv.time.values ]
# sw_times = [ datetime.datetime(2020, 1, 17, 0, 8) + datetime.timedelta(seconds=s) for s in sw.time.values ]


sd = sd.isel(trajectory=0).swap_dims({'obs': 'time'})
ss = ss.swap_dims({'ob': 'time'})
st = st.swap_dims({'ob': 'time'})
sv = sv.swap_dims({'ob': 'time'})
sw = sw.swap_dims({'ob': 'time'})


sal_sbe37 = sd['SAL_SBE37_MEAN']
sal_jpl_orbital=ss['smap_SSS']
sal_rss40_orbital=st['smap_SSS_40km']
sal_rss70_orbital=st['smap_SSS']
sal_rss40=sw['sat_sss_smap_40km']
sal_rss70=sw['sat_sss_smap']
sal_jpl=sv['sat_smap_sss']

d_fmt = DateFormatter("%m/%d")
legend_properties = {'weight':'semibold','size':'10'}

fig = plt.subplots(figsize=(18, 10))
ax = plt.subplot(1,1,1)
plt.plot(sal_sbe37.time.values, sal_sbe37.values,color='r',label='SBE37')
#plt.plot(ss_times, sal_jpl_orbital.values,color='c',label='JPL orbital')
#plt.plot(st_times, sal_rss40_orbital.values,color='g',label='RSS40 orbital')
#plt.plot(st_times, sal_rss70_orbital.values,color='b',label='RSS70 orbital')
plt.plot(sv_times, sal_jpl.values,color='b',label='JPL 8 day no repeat')
plt.plot(sw_times, sal_rss40.values,color='g',label='RSS40 8 day no repeat')
plt.plot(sw_times, sal_rss70.values,color='c',label='RSS70 8 day no repeat')


ax.xaxis.set_major_formatter(d_fmt)
ax.set_xlabel("Time mm-dd (UTC)", fontsize=15, fontweight='semibold')
ax.set_ylabel("Salinity (psu)", fontsize=15, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.legend(loc='best', prop=legend_properties)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
plt.grid(True, lw=0.5, ls=':')
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')

plt.savefig(PNG +' Salinity time series with orbital no repeat data'
              , dpi=150, facecolor='w', edgecolor='w',
             orientation='lanscape', papertype=None, format='png',
               bbox_inches='tight', pad_inches=0.1)
plt.show()

exit()
