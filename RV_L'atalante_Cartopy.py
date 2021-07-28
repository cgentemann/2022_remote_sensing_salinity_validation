"""
Created April 10, 2021, Author: Kashawn Hall
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
import datetime as dt

'''Research Vessel Data Pull'''
RV_latalante = xr.open_dataset("C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\RV_L'Atalante_data\\Latalante_TSG_20200206.nc")  # Importing RV TSG data 2020/02/06
RV_latalante1 = xr.open_dataset("C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\RV_L'Atalante_data\\Latalante_TSG_20200207.nc")  # Importing RV TSG data 2020/02/07
RV_latalante2 = xr.open_dataset("C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\RV_L'Atalante_data\\Latalante_TSG_20200208.nc")  # Importing RV TSG data 2020/02/08
rv_psal = RV_latalante.PSAL.values  # Setting the RV Salinity variable 2020/02/06
rv_lat = RV_latalante.LATITUDE.values  # Setting the RV Latitude variable 2020/02/06
rv_lon = RV_latalante.LONGITUDE.values  # Setting the RV Longitude variable 2020/02/06
rv1_psal = RV_latalante1.PSAL.values  # Setting the RV Salinity variable 2020/02/07
rv1_lat = RV_latalante1.LATITUDE.values  # Setting the RV Latitude variable 2020/02/07
rv1_lon = RV_latalante1.LONGITUDE.values  # Setting the RV Longitude variable 2020/02/07
rv2_psal = RV_latalante2.PSAL.values  # Setting the RV Salinity variable 2020/02/08
rv2_lat = RV_latalante2.LATITUDE.values  # Setting the RV Latitude variable 2020/02/08
rv2_lon = RV_latalante2.LONGITUDE.values  # Setting the RV Longitude variable 2020/02/08

'''CTD Data Pull'''
RV_CTD = xr.open_dataset("C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\RV_L'Atalante_data\\Latalante_CTD_20200207.nc")  # Importing RV CTD data 2020/02/07
RV_CTD1 = xr.open_dataset("C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\RV_L'Atalante_data\\Latalante_CTD_20200208.nc")  # Importing RV CTD data 2020/02/08

'''SMAP Data Pull'''
SMAP_JPL = xr.open_dataset("C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\RV_L'Atalante_data\\SMAP_L3_SSS_20200207_8DAYS_V5.0.nc")  # Importing 2020/02/07 JPL 8 day average
SMAP_JPL2 = xr.open_dataset("C:\\Users\\khall\\Documents\\GitHub\\paper_software\\2020_ATOMIC_Salinity\\data\\RV_L'Atalante_data\\SMAP_L3_SSS_20200208_8DAYS_V5.0.nc")  # Importing 2020/02/08 JPL 8 day average
# print(SMAP_JPL)

'''Extents'''
x1, x2 = -62.3922848, -45.5022528  # Extents taken from Saildrone_SMAP_Subplots.py
y1, y2 = 4.4158584, 15.101028  # Extents taken from Saildrone_SMAP_Subplots.py

'''Clipping the SMAP data to the extents to limit the color bar range'''
SMAP_JPL_clipped = SMAP_JPL.where((SMAP_JPL.latitude >= y1-0.125) & (SMAP_JPL.latitude <= y2+0.125) & (SMAP_JPL.longitude >= x1-0.125) & (SMAP_JPL.longitude <= x2+0.250), drop=True)
SMAP_JPL2_clipped = SMAP_JPL2.where((SMAP_JPL.latitude >= y1-0.125) & (SMAP_JPL.latitude <= y2+0.125) & (SMAP_JPL.longitude >= x1-0.125) & (SMAP_JPL.longitude <= x2+0.250), drop=True)
# print(SMAP_JPL_clipped)
sss_JPL = SMAP_JPL_clipped.smap_sss.values
anc_JPL = SMAP_JPL_clipped.anc_sss.values
lat_JPL = SMAP_JPL_clipped.latitude.data
lon_JPL = SMAP_JPL_clipped.longitude.data
sss_JPL2 = SMAP_JPL2_clipped.smap_sss.values
anc_JPL2 = SMAP_JPL2_clipped.anc_sss.values
lat_JPL2 = SMAP_JPL2_clipped.latitude.data
lon_JPL2 = SMAP_JPL2_clipped.longitude.data

'''Plotting'''
cmap = cmo.haline
norm = mpl.colors.Normalize(vmin=None, vmax=None, clip=False)
feature = cf.NaturalEarthFeature(name='land', category='physical', scale='10m', edgecolor='#000000', facecolor='#FFFFFF')

grid = plt.GridSpec(2, 3, wspace=0.2, hspace=0.1)  # Creating the layout of the figure
fig = plt.figure()

# Creating Top Left CTD Profile Plot
ax1 = plt.subplot(grid[0, 0])
ax1.plot(RV_CTD.PSAL[0][:], RV_CTD.PRES.data[0][:], label='01:01:59')
ax1.plot(RV_CTD.PSAL[1][:], RV_CTD.PRES.data[1][:], label='02:28:52')
ax1.plot(RV_CTD.PSAL[2][:], RV_CTD.PRES.data[2][:], label='05:00:06')
ax1.plot(RV_CTD.PSAL[3][:], RV_CTD.PRES.data[3][:], label='06:28:25')
ax1.plot(RV_CTD.PSAL[4][:], RV_CTD.PRES.data[4][:], label='08:48:47')
ax1.plot(RV_CTD.PSAL[5][:], RV_CTD.PRES.data[5][:], label='10:21:22')
ax1.plot(RV_CTD.PSAL[6][:], RV_CTD.PRES.data[6][:], label='12:15:12')
ax1.plot(RV_CTD.PSAL[7][:], RV_CTD.PRES.data[7][:], label='13:41:10')
ax1.set_ylim(-5, 500)
ax1.xaxis.set_label_position('top')
ax1.set_xlabel('Salinity (PSU)\n02-07-2020')
ax1.set_ylabel('Ocean Depth (Meters)')
ax1.xaxis.set_ticks_position('top')
ax1.grid(which='major', linestyle='-', linewidth='0.2', color='black')
ax1.minorticks_on()
ax1.grid(which='minor', linestyle=':', linewidth='0.1', color='black')
ax1.invert_yaxis()
ax1.text(34.55, 30.0, 'a)', fontsize=11)
ax1.legend(fontsize=9.2)

# Creating Top Right CTD Profile Plot
ax2 = plt.subplot(grid[0, 1])
ax2.plot(RV_CTD1.PSAL[0][:], RV_CTD1.PRES.data[0][:], label='04:05:29')
ax2.plot(RV_CTD1.PSAL[1][:], RV_CTD1.PRES.data[1][:], label='04:34:54')
ax2.plot(RV_CTD1.PSAL[2][:], RV_CTD1.PRES.data[2][:], label='06:39:23')
ax2.plot(RV_CTD1.PSAL[3][:], RV_CTD1.PRES.data[3][:], label='07:59:26')
ax2.plot(RV_CTD1.PSAL[4][:], RV_CTD1.PRES.data[4][:], label='10:33:08')
ax2.plot(RV_CTD1.PSAL[5][:], RV_CTD1.PRES.data[5][:], label='11:03:18')
ax2.plot(RV_CTD1.PSAL[6][:], RV_CTD1.PRES.data[6][:], label='13:08:13')
ax2.plot(RV_CTD1.PSAL[7][:], RV_CTD1.PRES.data[7][:], label='13:39:50')
ax2.set_ylim(-5, 500)
ax2.xaxis.set_label_position('top')
ax2.set_xlabel('Salinity (PSU)\n02-08-2020')
ax2.set_ylabel('Ocean Depth (Meters)')
ax2.xaxis.set_ticks_position('top')
ax2.grid(which='major', linestyle='-', linewidth='0.2', color='black')
ax2.minorticks_on()
ax2.grid(which='minor', linestyle=':', linewidth='0.1', color='black')
ax2.invert_yaxis()
ax2.text(34.58, 30.0, 'b)', fontsize=11)
ax2.legend(fontsize=9.2)

# Creating Bottom Left TSG/CTD Plot
ax3 = plt.subplot(grid[1, :2])
ax3.plot(RV_CTD.TIME.data[0], RV_CTD.PSAL[0][0], 'kX')  # Plotting Upper most CTD salinity reading of the cast 2020/02/07
ax3.plot(RV_CTD.TIME.data[1], RV_CTD.PSAL[1][0], 'kX')  # Plotting Upper most CTD salinity reading of the cast 2020/02/07
ax3.plot(RV_CTD.TIME.data[2], RV_CTD.PSAL[2][0], 'kX')  # Plotting Upper most CTD salinity reading of the cast 2020/02/07
ax3.plot(RV_CTD.TIME.data[3], RV_CTD.PSAL[3][0], 'kX')  # Plotting Upper most CTD salinity reading of the cast 2020/02/07
ax3.plot(RV_CTD.TIME.data[4], RV_CTD.PSAL[4][0], 'kX')  # Plotting Upper most CTD salinity reading of the cast 2020/02/07
ax3.plot(RV_CTD.TIME.data[5], RV_CTD.PSAL[5][0], 'kX')  # Plotting Upper most CTD salinity reading of the cast 2020/02/07
ax3.plot(RV_CTD.TIME.data[6], RV_CTD.PSAL[6][0], 'kX')  # Plotting Upper most CTD salinity reading of the cast 2020/02/07
ax3.plot(RV_CTD.TIME.data[7], RV_CTD.PSAL[7][0], 'kX')  # Plotting Upper most CTD salinity reading of the cast 2020/02/07
ax3.plot(RV_latalante1.TIME.values, rv1_psal, c='b', label='Thermosaliograph Salinity (3.5 Depth)')  # Plotting TSG vs DateTime 2020/02/07
ax3.plot(RV_latalante2.TIME.values, rv2_psal, c='b')  # Plotting TSG vs DateTime 2020/02/08
ax3.plot(RV_CTD1.TIME.data[0], RV_CTD1.PSAL[0][0], 'kX', label="CTD Salinty (approx. 5.0m Depth)")  # Plotting Upper most CTD salinity reading of the cast 2020/02/08
ax3.plot(RV_CTD1.TIME.data[1], RV_CTD1.PSAL[1][0], 'kX')  # Plotting Upper most CTD salinity reading of the cast 2020/02/08
ax3.plot(RV_CTD1.TIME.data[2], RV_CTD1.PSAL[2][0], 'kX')  # Plotting Upper most CTD salinity reading of the cast 2020/02/08
# ax3.plot(RV_CTD1.TIME.data[3], RV_CTD1.PSAL[3][0], 'kX')  #Erroneous data point
ax3.plot(RV_CTD1.TIME.data[4], RV_CTD1.PSAL[4][0], 'kX')  # Plotting Upper most CTD salinity reading of the cast 2020/02/08
ax3.plot(RV_CTD1.TIME.data[5], RV_CTD1.PSAL[5][0], 'kX')  # Plotting Upper most CTD salinity reading of the cast 2020/02/08
ax3.plot(RV_CTD1.TIME.data[6], RV_CTD1.PSAL[6][0], 'kX')  # Plotting Upper most CTD salinity reading of the cast 2020/02/08
ax3.plot(RV_CTD1.TIME.data[7], RV_CTD1.PSAL[7][0], 'kX')  # Plotting Upper most CTD salinity reading of the cast 2020/02/08
ax3.set_xlabel('Date time (mm-dd hr)')
ax3.set_ylabel('Salinity (PSU)')
ax3.grid(which='major', linestyle='-', linewidth='0.2', color='black')
ax3.minorticks_on()
ax3.grid(which='minor', linestyle=':', linewidth='0.1', color='black')
ax3.text(RV_latalante.TIME.data[615], 34.65, 'd)', fontsize=11)
ax3.legend()

# Creating Bottom right Cartopy Plot 2020/02/08
ax4 = plt.subplot(grid[1, 2], projection=ccrs.PlateCarree())
ax4.add_feature(feature)
ax4.set_extent([x1, x2, y1, y2])
fd1 = ax4.contourf(lon_JPL2, lat_JPL2, sss_JPL2, 60, cmap=cmap, norm=norm)
im_rv = ax4.scatter(rv_lon, rv_lat, c=rv_psal, s=.15, transform=ccrs.PlateCarree(), cmap=cmo.haline_r, norm=norm, label="RV L'atalante")
im_rv1 = ax4.scatter(rv1_lon, rv1_lat, c=rv1_psal, s=.15, transform=ccrs.PlateCarree(), cmap=cmo.haline_r, norm=norm)
im_rv2 = ax4.scatter(rv2_lon, rv2_lat, c=rv2_psal, s=.15, transform=ccrs.PlateCarree(), cmap=cmo.haline_r, norm=norm)
ctd2 = plt.plot(RV_CTD1.LONGITUDE, RV_CTD1.LATITUDE,  markersize=2.5, marker='o', color='red', label="CTD Scan Locations")
ax4.set_xticks([-62, -58, -54, -50, -46], crs=ccrs.PlateCarree())
ax4.set_yticks([4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax4.xaxis.set_major_formatter(lon_formatter)
ax4.yaxis.set_major_formatter(lat_formatter)
ax4.legend(loc='upper right')
ax4.xaxis.set_label_position('top')
ax4.set_xlabel('JPL Salinity v5\n02-08-2020 8-day Average')
ax4.text(-62, 5.75, 'e)', fontsize=11)
ax4.text(-62, 4.65, 'South America', fontsize=11)
ax4.coastlines(resolution='10m')

# Creating Top right Cartopy Plot 2020/02/07
ax5 = plt.subplot(grid[:1, 2], projection=ccrs.PlateCarree())
fd1 = ax5.contourf(lon_JPL, lat_JPL, sss_JPL, 60, cmap=cmap, norm=norm)
im_rv = ax5.scatter(rv_lon, rv_lat, c=rv_psal, s=.15, transform=ccrs.PlateCarree(), cmap=cmo.haline_r, norm=norm, label="RV L'atalante")
im_rv1 = ax5.scatter(rv1_lon, rv1_lat, c=rv1_psal, s=.15, transform=ccrs.PlateCarree(), cmap=cmo.haline_r, norm=norm)
im_rv2 = ax5.scatter(rv2_lon, rv2_lat, c=rv2_psal, s=.15, transform=ccrs.PlateCarree(), cmap=cmo.haline_r, norm=norm)
ctd1 = plt.plot(RV_CTD.LONGITUDE, RV_CTD.LATITUDE,  markersize=2.5, marker='o', color='red', label="CTD Scan Locations")
ax5.set_xticks([-62, -58, -54, -50, -46], crs=ccrs.PlateCarree())
ax5.set_yticks([4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax5.xaxis.set_major_formatter(lon_formatter)
ax5.yaxis.set_major_formatter(lat_formatter)
ax5.legend(loc='upper right')
ax5.xaxis.set_label_position('top')
ax5.set_xlabel('JPL Salinity v5\n02-07-2020 8-day Average')
ax5.text(-62, 5.75, 'c)', fontsize=11)
ax5.text(-62, 4.65, 'South America', fontsize=11)
ax5.coastlines(resolution='10m')

# fig.subplots_adjust(right=0.1)
cbar_ax = fig.add_axes([0.935, 0.1, 0.01, 0.8])
cb = fig.colorbar(fd1, cax=cbar_ax, orientation="vertical")
cb.ax.set_ylabel('Salinity (PSU)', fontsize=10, fontweight='semibold')
figManager = plt.get_current_fig_manager()
figManager.window.state('zoomed')
# plt.savefig('C:\\Users\\khall\\Desktop\\PythonMPIM\\OTECandSWAC\\EUREC4A\\Saildrone\\Figures\\RV_Subplot.png', bbox_inches='tight')
plt.show()
