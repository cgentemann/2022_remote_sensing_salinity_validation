"""
Created Mar 29 2021, Author: Kashawn Hall
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cf
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmocean.cm as cmo
from glob import glob

'''Saildrone Data Pull'''
ds_1026 = xr.open_dataset('C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\saildrone-gen_5-atomic_eurec4a_2020-sd1026-20200117T000000-20200302T235959-1_minutes-v1.1589306725934.nc')  # Importing SD 1026
ds_1060 = xr.open_dataset('C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\saildrone-gen_5-atomic_eurec4a_2020-sd1060-20200117T000000-20200302T235959-1_minutes-v1.1589306886594.nc')  # Importing SD 1060
ds_1061 = xr.open_dataset('C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\saildrone-gen_5-atomic_eurec4a_2020-sd1061-20200117T000000-20200302T235959-1_minutes-v1.1589307121602.nc')  # Importing SD 1061
ds_1026 = ds_1026.isel(trajectory=0).swap_dims({'obs': 'time'})  # Switching dimensions from "obs" to "time"
ds_1060 = ds_1060.isel(trajectory=0).swap_dims({'obs': 'time'})  # Switching dimensions from "obs" to "time"
ds_1061 = ds_1061.isel(trajectory=0).swap_dims({'obs': 'time'})  # Switching dimensions from "obs" to "time"
ds_1026_dot = ds_1026.sel(time=slice('2020-02-18T00:00:00'))  # Slicing SD1026 time to showcase the SD in the middle of the fresh tongue

'''SMAP Data Pull'''
SMAP_JPL = xr.open_dataset("C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\RV_L'Atalante_data\\JPL_20200217_8DAYS_V5.0.nc")  # Importing JPL 8 day average
SMAP_RSS = xr.open_dataset("C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\RV_L'Atalante_data\\RSS_20200217_8day_running_v04.0.nc")  # Importing RSS 8 day average
# print(SMAP_JPL)

'''Extents'''
dx, dy = 3.05, 3.05
x1, x2 = ds_1026.longitude.min().data - dx, ds_1026.longitude.max().data + dx  # Setting X extents based on the max SD extents
y1, y2 = ds_1026.latitude.min().data - dy, ds_1026.latitude.max().data + dy  # Setting Y extents based on the max SD extents
# print('x1=', x1, 'x2=', x2, 'y1=', y1, 'y2=', y2)

'''Clipping the SMAP data to the extents to limit the color bar range'''
SMAP_JPL_clipped = SMAP_JPL.where((SMAP_JPL.latitude >= y1-0.125) & (SMAP_JPL.latitude <= y2+0.125) & (SMAP_JPL.longitude >= x1-0.125) & (SMAP_JPL.longitude <= x2+0.250), drop=True)
# print(SMAP_JPL_clipped)
sss_JPL = SMAP_JPL_clipped.smap_sss.values  # Setting the JPL Salinity variable
anc_JPL = SMAP_JPL_clipped.anc_sss.values  # Setting the JPL HYCOM Salinity variable
lat_JPL = SMAP_JPL_clipped.latitude.data  # Setting the JPL Latitude variable
lon_JPL = SMAP_JPL_clipped.longitude.data  # Setting the JPL Longitude variable

SMAP_RSS_clipped = SMAP_RSS.where((SMAP_RSS.lat >= y1-0.125) & (SMAP_RSS.lat <= y2+0.125) & (SMAP_RSS.lon >= 298.625-1.125) & (SMAP_RSS.lon <= 313.375+1.250), drop=True)
# print(SMAP_RSS_clipped)
sss_RSS = SMAP_RSS_clipped.sss_smap.values  # Setting the RSS Salinity variable
sss_RSS_40km = SMAP_RSS_clipped.sss_smap_40km.values  # Setting the RSS 40km Salinity variable
lat_RSS = SMAP_RSS_clipped.lat.data   # Setting the RSS Latitude variable
lon_RSS = SMAP_RSS_clipped.lon.data  # Setting the RSS Longitude variable

'''Plotting'''
cmap = cmo.haline  # Colormap choice
norm = mpl.colors.Normalize(vmin=None, vmax=None, clip=False)  # Normalizing colorbar between all subplots
feature = cf.NaturalEarthFeature(name='land', category='physical', scale='10m', edgecolor='#000000', facecolor='#FFFFFF')  # Import land features

# Plotting RSS SMAP Salinity
fig = plt.figure()
ax1 = plt.subplot(2, 2, 1, projection=ccrs.PlateCarree())
fd1 = ax1.contourf(lon_RSS, lat_RSS, sss_RSS, 60, cmap=cmap, norm=norm)
ax1.scatter(ds_1026.longitude, ds_1026.latitude, c=ds_1026.SAL_SBE37_MEAN, s=.15, transform=ccrs.PlateCarree(), cmap=cmap, label='Saildrone 1026', norm=norm)
ax1.scatter(ds_1060.longitude, ds_1060.latitude, c=ds_1060.SAL_SBE37_MEAN, s=.15, transform=ccrs.PlateCarree(), cmap=cmap, label='Saildrone 1060', norm=norm)
ax1.scatter(ds_1061.longitude, ds_1061.latitude, c=ds_1061.SAL_SBE37_MEAN, s=.15, transform=ccrs.PlateCarree(), cmap=cmap, label='Saildrone 1061', norm=norm)
dot = plt.plot(ds_1026_dot.longitude.data[-1], ds_1026_dot.latitude.data[-1],  markersize=5, marker='o', color='white')
ax1.set_xticks([-62, -58, -54, -50, -46], crs=ccrs.PlateCarree())
ax1.set_yticks([4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax1.xaxis.set_major_formatter(lon_formatter)
ax1.yaxis.set_major_formatter(lat_formatter)
ax1.text(-62, 6.75, 'South\nAmerica', fontsize=10)
ax1.text(-62, 5.05, 'A) RSS SMAP Salinity v4', fontsize=8)
ax1.text(-62, 4.55, '02-17-2020 8-day Average', fontsize=8)
ax1.coastlines(resolution='10m')

# Plotting RSS 40km SMAP Salinity
ax2 = plt.subplot(2, 2, 2, projection=ccrs.PlateCarree())
fd2 = ax2.contourf(lon_RSS, lat_RSS, sss_RSS_40km, 60, cmap=cmap, norm=norm)
ax2.scatter(ds_1026.longitude, ds_1026.latitude, c=ds_1026.SAL_SBE37_MEAN, s=.15, transform=ccrs.PlateCarree(), cmap=cmap, label='Saildrone 1026', norm=norm)
ax2.scatter(ds_1060.longitude, ds_1060.latitude, c=ds_1060.SAL_SBE37_MEAN, s=.15, transform=ccrs.PlateCarree(), cmap=cmap, label='Saildrone 1060', norm=norm)
ax2.scatter(ds_1061.longitude, ds_1061.latitude, c=ds_1061.SAL_SBE37_MEAN, s=.15, transform=ccrs.PlateCarree(), cmap=cmap, label='Saildrone 1061', norm=norm)
dot = plt.plot(ds_1026_dot.longitude.data[-1], ds_1026_dot.latitude.data[-1],  markersize=5, marker='o', color='white')
ax2.set_xticks([-62, -58, -54, -50, -46], crs=ccrs.PlateCarree())
ax2.set_yticks([4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax2.xaxis.set_major_formatter(lon_formatter)
ax2.yaxis.set_major_formatter(lat_formatter)
ax2.text(-62, 6.75, 'South\nAmerica', fontsize=10)
ax2.text(-62, 5.05, 'B) RSS SMAP Salinity v4-40km', fontsize=8)
ax2.text(-62, 4.55, '02-17-2020 8-day Average', fontsize=8)
ax2.coastlines(resolution='10m')

# Plotting JPL SMAP Salinity
ax3 = plt.subplot(2, 2, 3, projection=ccrs.PlateCarree())
fd3 = ax3.contourf(lon_JPL, lat_JPL, sss_JPL, 60, cmap=cmap, norm=norm)
ax3.scatter(ds_1026.longitude, ds_1026.latitude, c=ds_1026.SAL_SBE37_MEAN, s=.15, transform=ccrs.PlateCarree(), cmap=cmap, label='Saildrone 1026', norm=norm)
ax3.scatter(ds_1060.longitude, ds_1060.latitude, c=ds_1060.SAL_SBE37_MEAN, s=.15, transform=ccrs.PlateCarree(), cmap=cmap, label='Saildrone 1060', norm=norm)
ax3.scatter(ds_1061.longitude, ds_1061.latitude, c=ds_1061.SAL_SBE37_MEAN, s=.15, transform=ccrs.PlateCarree(), cmap=cmap, label='Saildrone 1061', norm=norm)
dot = plt.plot(ds_1026_dot.longitude.data[-1], ds_1026_dot.latitude.data[-1],  markersize=5, marker='o', color='white')
ax3.set_xticks([-62, -58, -54, -50, -46], crs=ccrs.PlateCarree())
ax3.set_yticks([4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax3.xaxis.set_major_formatter(lon_formatter)
ax3.yaxis.set_major_formatter(lat_formatter)
ax3.text(-62, 6.75, 'South\nAmerica', fontsize=10)
ax3.text(-62, 5.05, 'C) JPL SMAP Salinity v5', fontsize=8)
ax3.text(-62, 4.55, '02-17-2020 8-day Average', fontsize=8)
ax3.coastlines(resolution='10m')

# Plotting JPL HYCOM SMAP Salinity
ax4 = plt.subplot(2,2,4, projection=ccrs.PlateCarree())
fd4 = ax4.contourf(lon_JPL, lat_JPL, anc_JPL, 60, cmap=cmap, norm=norm)
im1026_4 = ax4.scatter(ds_1026.longitude, ds_1026.latitude, c=ds_1026.SAL_SBE37_MEAN, s=.15, transform=ccrs.PlateCarree(),cmap=cmap, label='Saildrone 1026', norm=norm)
dot = plt.plot(ds_1026_dot.longitude.data[-1], ds_1026_dot.latitude.data[-1],  markersize=5, marker='o', color='white')
ax4.set_xticks([-62, -58, -54, -50, -46], crs=ccrs.PlateCarree())
ax4.set_yticks([4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax4.xaxis.set_major_formatter(lon_formatter)
ax4.yaxis.set_major_formatter(lat_formatter)
ax4.text(-62, 6.75, 'South\nAmerica', fontsize=10)
ax4.text(-62, 5.05, 'D) JPL HYCOM Salinity v5', fontsize=8)
ax4.text(-62, 4.55, '02-17-2020 8-day Average', fontsize=8)
ax4.coastlines(resolution='10m')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.8, 0.11, 0.02, 0.77])  # Set Location of colorbar
cb = fig.colorbar(fd3, cax=cbar_ax, orientation="vertical")  # Create colorbar
cb.ax.set_ylabel('Salinity (PSU)', fontsize=10, fontweight='semibold')  # Set Label of colorbar
figManager = plt.get_current_fig_manager()
figManager.window.state('zoomed')
# plt.savefig('C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\figures\\Satellite_subplots.png', bbox_inches='tight')
plt.show()

