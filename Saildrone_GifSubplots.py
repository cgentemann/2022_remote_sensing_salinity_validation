"""
Created May 11 2021, Author: Kashawn Hall
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cf
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmocean.cm as cmo
from glob import glob

'''Saildrone and SMAP Data Pull'''
ds_1026 = xr.open_dataset('C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\saildrone-gen_5-atomic_eurec4a_2020-sd1026-20200117T000000-20200302T235959-1_minutes-v1.1589306725934.nc')  # Importing SD 1026
ds_1026 = ds_1026.isel(trajectory=0).swap_dims({'obs': 'time'})  # Switching dimensions from "obs" to "time"
ds_1026_dot = ds_1026.sel(time=slice('2020-02-18T00:00:00'))  # Slicing SD1026 time to showcase the SD in the middle of the fresh tongue

'''Extents'''
dx, dy = 14.05, 14.05
x1, x2 = ds_1026.longitude.min().data - dx, ds_1026.longitude.max().data + dx  # Setting X extents based on the max SD extents
y1, y2 = ds_1026.latitude.min().data - dy, ds_1026.latitude.max().data + dy  # Setting Y extents based on the max SD extents
# print('x1=', x1, 'x2=', x2, 'y1=', y1, 'y2=', y2)

'''Clipping the SMAP data to the extents to limit the color bar range'''
SMAP_JPL_01 = xr.open_dataset("C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\RV_L'Atalante_data\\SMAP_L3_SSS_20200201_8DAYS_V5.0.nc")  # Importing 2020/02/01 JPL 8 day average
SMAP_JPL_01_clipped = SMAP_JPL_01.where((SMAP_JPL_01.latitude >= y1-0.125) & (SMAP_JPL_01.latitude <= y2+2) & (SMAP_JPL_01.longitude >= x1-0.125) & (SMAP_JPL_01.longitude <= x2+0.250), drop=True)
sss_JPL_01 = SMAP_JPL_01_clipped.smap_sss.values  # Setting the JPL Salinity variable 2020/02/01
anc_JPL_01 = SMAP_JPL_01_clipped.anc_sss.values  # Setting the JPL HYCOM Salinity variable 2020/02/01
lat_JPL_01 = SMAP_JPL_01_clipped.latitude.data  # Setting the JPL Latitude variable 2020/02/01
lon_JPL_01 = SMAP_JPL_01_clipped.longitude.data  # Setting the JPL Longitude variable 2020/02/01

# Importing 2020/02/08 JPL 8 day average
SMAP_JPL_08 = xr.open_dataset("C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\RV_L'Atalante_data\\SMAP_L3_SSS_20200208_8DAYS_V5.0.nc")
SMAP_JPL_08_clipped = SMAP_JPL_08.where((SMAP_JPL_08.latitude >= y1-0.125) & (SMAP_JPL_08.latitude <= y2+2) & (SMAP_JPL_08.longitude >= x1-0.125) & (SMAP_JPL_08.longitude <= x2+0.250), drop=True)
sss_JPL_08 = SMAP_JPL_08_clipped.smap_sss.values  # Setting the JPL Salinity variable 2020/02/08
anc_JPL_08 = SMAP_JPL_08_clipped.anc_sss.values  # Setting the JPL HYCOM Salinity variable 2020/02/08
lat_JPL_08 = SMAP_JPL_08_clipped.latitude.data  # Setting the JPL Latitude variable 2020/02/08
lon_JPL_08 = SMAP_JPL_08_clipped.longitude.data  # Setting the JPL Longitude variable 2020/02/08

# Importing 2020/02/15 JPL 8 day average
SMAP_JPL_15 = xr.open_dataset("C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\RV_L'Atalante_data\\SMAP_L3_SSS_20200215_8DAYS_V5.0.nc")
SMAP_JPL_15_clipped = SMAP_JPL_15.where((SMAP_JPL_15.latitude >= y1-0.125) & (SMAP_JPL_15.latitude <= y2+2) & (SMAP_JPL_15.longitude >= x1-0.125) & (SMAP_JPL_15.longitude <= x2+0.250), drop=True)
sss_JPL_15 = SMAP_JPL_15_clipped.smap_sss.values  # Setting the JPL Salinity variable 2020/02/15
anc_JPL_15 = SMAP_JPL_15_clipped.anc_sss.values  # Setting the JPL HYCOM Salinity variable 2020/02/15
lat_JPL_15 = SMAP_JPL_15_clipped.latitude.data  # Setting the JPL Latitude variable 2020/02/15
lon_JPL_15 = SMAP_JPL_15_clipped.longitude.data  # Setting the JPL Longitude variable 2020/02/15

# Importing 2020/02/22 JPL 8 day average
SMAP_JPL_22 = xr.open_dataset("C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\RV_L'Atalante_data\\SMAP_L3_SSS_20200222_8DAYS_V5.0.nc")
SMAP_JPL_22_clipped = SMAP_JPL_22.where((SMAP_JPL_22.latitude >= y1-0.125) & (SMAP_JPL_22.latitude <= y2+2) & (SMAP_JPL_22.longitude >= x1-0.125) & (SMAP_JPL_22.longitude <= x2+0.250), drop=True)
sss_JPL_22 = SMAP_JPL_22_clipped.smap_sss.values  # Setting the JPL Salinity variable 2020/02/22
anc_JPL_22 = SMAP_JPL_22_clipped.anc_sss.values  # Setting the JPL HYCOM Salinity variable 2020/02/22
lat_JPL_22 = SMAP_JPL_22_clipped.latitude.data  # Setting the JPL Latitude variable 2020/02/22
lon_JPL_22 = SMAP_JPL_22_clipped.longitude.data  # Setting the JPL Longitude variable 2020/02/22

# Importing 2020/02/29 JPL 8 day average
SMAP_JPL_29 = xr.open_dataset("C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\RV_L'Atalante_data\\SMAP_L3_SSS_20200229_8DAYS_V5.0.nc")
SMAP_JPL_29_clipped = SMAP_JPL_29.where((SMAP_JPL_29.latitude >= y1-0.125) & (SMAP_JPL_29.latitude <= y2+2) & (SMAP_JPL_29.longitude >= x1-0.125) & (SMAP_JPL_29.longitude <= x2+0.250), drop=True)
sss_JPL_29 = SMAP_JPL_29_clipped.smap_sss.values  # Setting the JPL Salinity variable 2020/02/29
anc_JPL_29 = SMAP_JPL_29_clipped.anc_sss.values  # Setting the JPL HYCOM Salinity variable 2020/02/29
lat_JPL_29 = SMAP_JPL_29_clipped.latitude.data  # Setting the JPL Latitude variable 2020/02/29
lon_JPL_29 = SMAP_JPL_29_clipped.longitude.data  # Setting the JPL Longitude variable 2020/02/29

'''Plotting'''
cmap = cmo.haline  # Colormap choice

norm = mpl.colors.Normalize(vmin=31.85, vmax=max(sss_JPL_29.flatten()), clip=True)  # Normalizing colorbar between all subplots
feature = cf.NaturalEarthFeature(name='land', category='physical', scale='10m', edgecolor='#000000', facecolor='#FFFFFF')  # Import land features

# Plotting 2020/02/01 JPL SMAP Salinity
fig = plt.figure()
ax1 = plt.subplot(2, 5, 1, projection=ccrs.PlateCarree())
v = np.linspace(30.15, 38.28, 9, endpoint=True)
# print(v)
v = [round(num,2) for num in v]
fd1 = ax1.pcolormesh(lon_JPL_01, lat_JPL_01, sss_JPL_01, vmin=30.15, vmax=38.28, cmap=cmap)
ax1.set_xticks([-65, -55, -45], crs=ccrs.PlateCarree())
ax1.set_yticks([0, 10, 20], crs=ccrs.PlateCarree())
ax1.tick_params(axis='both', labelsize=7)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax1.xaxis.set_major_formatter(lon_formatter)
ax1.yaxis.set_major_formatter(lat_formatter)
ax1.xaxis.set_label_position('top')
ax1.set_xlabel('JPL Salinity v5\n02-01-2020 8-day Average')
ax1.text(-72, -5.05, 'South\nAmerica', fontsize=10)
ax1.coastlines(resolution='10m')

# Plotting 2020/02/08 JPL SMAP Salinity
ax2 = plt.subplot(2, 5, 2, projection=ccrs.PlateCarree())
ax2.add_feature(feature)
ax2.set_extent([x1, x2, y1, y2])
fd2 = ax2.pcolormesh(lon_JPL_08, lat_JPL_08, sss_JPL_08, vmin=30.15, vmax=38.28, cmap=cmap)
ax2.set_xticks([-65, -55, -45], crs=ccrs.PlateCarree())
ax2.set_yticks([0, 10, 20], crs=ccrs.PlateCarree())
ax2.tick_params(axis='both', labelsize=7)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax2.xaxis.set_major_formatter(lon_formatter)
ax2.yaxis.set_major_formatter(lat_formatter)
ax2.xaxis.set_label_position('top')
ax2.set_xlabel('JPL Salinity v5\n02-08-2020 8-day Average')
ax2.text(-72, -5.05, 'South\nAmerica', fontsize=10)
ax2.coastlines(resolution='10m')

# Plotting 2020/02/15 JPL SMAP Salinity
ax3 = plt.subplot(2, 5, 3, projection=ccrs.PlateCarree())
ax3.add_feature(feature)
ax3.set_extent([x1, x2, y1, y2])
fd3 = ax3.pcolormesh(lon_JPL_15, lat_JPL_15, sss_JPL_15, vmin=30.15, vmax=38.28, cmap=cmap)
ax3.set_xticks([-65, -55, -45], crs=ccrs.PlateCarree())
ax3.set_yticks([0, 10, 20], crs=ccrs.PlateCarree())
ax3.tick_params(axis='both', labelsize=7)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax3.xaxis.set_major_formatter(lon_formatter)
ax3.yaxis.set_major_formatter(lat_formatter)
ax3.xaxis.set_label_position('top')
ax3.set_xlabel('JPL Salinity v5\n02-15-2020 8-day Average')
ax3.text(-72, -5.05, 'South\nAmerica', fontsize=10)
ax3.coastlines(resolution='10m')

# Plotting 2020/02/22 JPL SMAP Salinity
ax4 = plt.subplot(2, 5, 4, projection=ccrs.PlateCarree())
ax4.add_feature(feature)
ax4.set_extent([x1, x2, y1, y2])
fd4 = ax4.pcolormesh(lon_JPL_22, lat_JPL_22, sss_JPL_22, vmin=30.15, vmax=38.28, cmap=cmap)
ax4.set_xticks([-65, -55, -45], crs=ccrs.PlateCarree())
ax4.set_yticks([0, 10, 20], crs=ccrs.PlateCarree())
ax4.tick_params(axis='both', labelsize=7)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax4.xaxis.set_major_formatter(lon_formatter)
ax4.yaxis.set_major_formatter(lat_formatter)
ax4.xaxis.set_label_position('top')
ax4.set_xlabel('JPL Salinity v5\n02-22-2020 8-day Average')
ax4.text(-72, -5.05, 'South\nAmerica', fontsize=10)
ax4.coastlines(resolution='10m')

# Plotting 2020/02/29 JPL SMAP Salinity
ax5 = plt.subplot(2, 5, 5, projection=ccrs.PlateCarree())
ax5.add_feature(feature)
ax5.set_extent([x1, x2, y1, y2])
fd5 = ax5.pcolormesh(lon_JPL_29, lat_JPL_29, sss_JPL_29, vmin=30.15, vmax=38.28, cmap=cmap, norm=norm)
ax5.set_xticks([-65, -55, -45], crs=ccrs.PlateCarree())
ax5.set_yticks([0, 10, 20], crs=ccrs.PlateCarree())
ax5.tick_params(axis='both', labelsize=7)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax5.xaxis.set_major_formatter(lon_formatter)
ax5.yaxis.set_major_formatter(lat_formatter)
ax5.xaxis.set_label_position('top')
ax5.set_xlabel('JPL Salinity v5\n02-29-2020 8-day Average')
ax5.text(-72, -5.05, 'South\nAmerica', fontsize=10)
ax5.coastlines(resolution='10m')

# Plotting 2020/02/01 JPL HYCOM SMAP Salinity
ax6 = plt.subplot(2, 5, 6, projection=ccrs.PlateCarree())
ax6.add_feature(feature)
ax6.set_extent([x1, x2, y1, y2])
fd6 = ax6.pcolormesh(lon_JPL_01, lat_JPL_01, anc_JPL_01, vmin=30.15, vmax=38.28, cmap=cmap)
ax6.set_xticks([-65, -55, -45], crs=ccrs.PlateCarree())
ax6.set_yticks([0, 10, 20], crs=ccrs.PlateCarree())
ax6.tick_params(axis='both', labelsize=7)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax6.xaxis.set_major_formatter(lon_formatter)
ax6.yaxis.set_major_formatter(lat_formatter)
ax6.xaxis.set_label_position('top')
ax6.set_xlabel('JPL HYCOM Salinity v5\n02-01-2020 8-day Average')
ax6.text(-72, -5.05, 'South\nAmerica', fontsize=10)
ax6.coastlines(resolution='10m')

# Plotting 2020/02/08 JPL HYCOM SMAP Salinity
ax7 = plt.subplot(2, 5, 7, projection=ccrs.PlateCarree())
ax7.add_feature(feature)
ax7.set_extent([x1, x2, y1, y2])
fd7 = ax7.pcolormesh(lon_JPL_08, lat_JPL_08, anc_JPL_08, vmin=30.15, vmax=38.28, cmap=cmap)
ax7.set_xticks([-65, -55, -45], crs=ccrs.PlateCarree())
ax7.set_yticks([0, 10, 20], crs=ccrs.PlateCarree())
ax7.tick_params(axis='both', labelsize=7)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax7.xaxis.set_major_formatter(lon_formatter)
ax7.yaxis.set_major_formatter(lat_formatter)
ax7.xaxis.set_label_position('top')
ax7.set_xlabel('JPL HYCOM Salinity v5\n02-08-2020 8-day Average')
ax7.text(-72, -5.05, 'South\nAmerica', fontsize=10)
ax7.coastlines(resolution='10m')

# Plotting 2020/02/15 JPL HYCOM SMAP Salinity
ax8 = plt.subplot(2, 5, 8, projection=ccrs.PlateCarree())
ax8.add_feature(feature)
ax8.set_extent([x1, x2, y1, y2])
fd8 = ax8.pcolormesh(lon_JPL_15, lat_JPL_15, anc_JPL_15, vmin=30.15, vmax=38.28, cmap=cmap)
ax8.set_xticks([-65, -55, -45], crs=ccrs.PlateCarree())
ax8.set_yticks([0, 10, 20], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax8.xaxis.set_major_formatter(lon_formatter)
ax8.yaxis.set_major_formatter(lat_formatter)
ax8.tick_params(axis='both', labelsize=7)
ax8.xaxis.set_label_position('top')
ax8.set_xlabel('JPL HYCOM Salinity v5\n02-15-2020 8-day Average')
ax8.text(-72, -5.05, 'South\nAmerica', fontsize=10)
ax8.coastlines(resolution='10m')

# Plotting 2020/02/22 JPL HYCOM SMAP Salinity
ax9 = plt.subplot(2, 5, 9, projection=ccrs.PlateCarree())
ax9.add_feature(feature)
ax9.set_extent([x1, x2, y1, y2])
fd9 = ax9.pcolormesh(lon_JPL_22, lat_JPL_22, anc_JPL_22, vmin=30.15, vmax=38.28, cmap=cmap)
ax9.set_xticks([-65, -55, -45], crs=ccrs.PlateCarree())
ax9.set_yticks([0, 10, 20], crs=ccrs.PlateCarree())
ax9.tick_params(axis='both', labelsize=7)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax9.xaxis.set_major_formatter(lon_formatter)
ax9.yaxis.set_major_formatter(lat_formatter)
ax9.xaxis.set_label_position('top')
ax9.set_xlabel('JPL HYCOM Salinity v5\n02-22-2020 8-day Average')
ax9.text(-72, -5.05, 'South\nAmerica', fontsize=10)
ax9.coastlines(resolution='10m')

# Plotting 2020/02/29 JPL HYCOM SMAP Salinity
ax10 = plt.subplot(2, 5, 10, projection=ccrs.PlateCarree())
ax10.add_feature(feature)
ax10.set_extent([x1, x2, y1, y2])
fd10 = ax10.pcolormesh(lon_JPL_29, lat_JPL_29, anc_JPL_29, vmin=30.15, vmax=38.28, cmap=cmap, norm=norm)
ax10.set_xticks([-65, -55, -45], crs=ccrs.PlateCarree())
ax10.set_yticks([0, 10, 20], crs=ccrs.PlateCarree())
ax10.tick_params(axis='both', labelsize=7)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax10.xaxis.set_major_formatter(lon_formatter)
ax10.yaxis.set_major_formatter(lat_formatter)
ax10.xaxis.set_label_position('top')
ax10.set_xlabel('JPL HYCOM Salinity v5\n02-29-2020 8-day Average')
ax10.text(-72, -5.05, 'South\nAmerica', fontsize=10)
ax10.coastlines(resolution='10m')

fig.subplots_adjust(bottom=0.31)
cbar_ax = fig.add_axes([0.122, 0.275, 0.78125, 0.015])  # Set Location of colorbar
cb = fig.colorbar(fd1, cax=cbar_ax, ticks=v, orientation="horizontal")  # Create colorbar
cb.ax.set_xlabel('Salinity (PSU)', fontsize=10, fontweight='semibold')  # Set Label of colorbar
figManager = plt.get_current_fig_manager()
figManager.window.state('zoomed')
# plt.savefig('C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\figures\\GIFsubplotset.png', bbox_inches='tight')
plt.show()
