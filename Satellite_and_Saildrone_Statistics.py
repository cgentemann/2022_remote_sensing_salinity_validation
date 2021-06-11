#!/usr/bin/env python
##This script does calculations on saildrone data and satellite data
import xarray as xr
import xskillscore as xs
from glob import glob

#data directory for 8 day nearest satellite data
data_dir = './../paper_software/2020_ATOMIC_Salinity/data/sss_collocations_8day_nearest/'

#data directory for saildrone and satellite orbital data
data_sat = './../paper_software/2020_ATOMIC_Salinity/data/sss_collocations_orbital/'

##data directory for raw saildrone data
data_sail = './../paper_software/2020_ATOMIC_Salinity/data/'

#create xarray dataset with saildrone filenames
saildrone_filenames = glob(data_sail+'saildrone*.nc')
ds_usv=xr.open_mfdataset(saildrone_filenames,concat_dim='trajectory').rename({'latitude':'lat','longitude':'lon'})
files = [x for x in glob(data_sail+'saildrone*.nc')]

#Comparing RBR and SBE37 Saildrone Instruments
for x in files:

    sail = xr.open_dataset(x).isel(trajectory=0).swap_dims({'obs': 'time'})
    num=sail.trajectory.values

    ##Data stats
    #Difference in Standard Deviation
    std_rbr= sail["SAL_RBR_MEAN"].std(dim="time", skipna=None)
    std_sbe37= sail["SAL_SBE37_MEAN"].std(dim="time",skipna=None)
    sail_diff=std_sbe37-std_rbr

    #RMSE
    rmse_sal=xs.rmse(sail.SAL_SBE37_MEAN,sail.SAL_RBR_MEAN, dim="time",skipna="True")

    #Mean Difference
    mean_sbe=sail.SAL_SBE37_MEAN.mean(dim='time')
    mean_rbr=sail.SAL_RBR_MEAN.mean(dim='time')
    sail_mean=mean_sbe-mean_rbr

    print("SAILDRONE " +str(num))
    print("STD Difference SBE37 - RBR = " + str(sail_diff.values))
    print("RMSE = " + str(rmse_sal.values))
    print("Mean Difference SBE37 - RBR = " +str(sail_mean.values))

#Comparing 8 day and orbital satellite data
##Netcdf Satellite files into xarray
SMAP_8RSS=xr.open_mfdataset(data_dir+'*RSS*.nc',concat_dim="trajectory",combine="nested")
SMAPRSS=xr.open_mfdataset(data_sat+'*RSS*.nc',concat_dim="trajectory",combine="nested")

SMAP_8JPL=xr.open_mfdataset(data_dir+'*JPL*.nc',concat_dim="trajectory",combine="nested")
SMAPJPL=xr.open_mfdataset(data_sat+'*JPL*.nc',concat_dim="trajectory",combine="nested")

#Mean Difference
rss_mean = SMAPRSS.smap_SSS.mean(dim='trajectory')
rss_mean2 = rss_mean.mean(dim='time')

rss8_mean= SMAP_8RSS.sat_sss_smap.mean(dim='trajectory')
rss8_mean2 = rss8_mean.mean(dim='time')
rssdiff = rss8_mean2 - rss_mean2

jpl_mean = SMAPJPL.smap_SSS.mean(dim='trajectory')
jpl_mean2 = jpl_mean.mean(dim='time')

jpl8_mean = SMAP_8JPL.sat_smap_sss.mean(dim='trajectory')
jpl8_mean2 = jpl8_mean.mean(dim='time')
jpldiff = jpl8_mean2 - jpl_mean2

print("Mean Difference RSS 8 day - RSS Orbital")
print(rssdiff.values)
print("Mean Difference JPL 8 day - JPL Orbital")
print(jpldiff.values)

##Standard Deviation Difference
rss_std = rss_mean.std(dim='time')
rss8_std = rss8_mean.std(dim='time')

jpl_std = jpl_mean.std(dim='time')
jpl8_std = jpl8_mean.std(dim='time')

std_rss= rss8_std - rss_std
std_jpl= jpl8_std - jpl_std

print("STD Difference RSS 8 day - RSS Orbital")
print(std_rss.values)

print("STD Difference JPL 8 day - JPL Orbital")
print(std_jpl.values)

##RMSE
rmse_rss=xs.rmse(rss8_mean, rss_mean, dim="time", skipna="True")
rmse_jpl=xs.rmse(jpl8_mean, jpl_mean, dim="time", skipna="True")

print("RMSE RSS 8 day and RSS Orbital")
print(rmse_rss.values)

print("RMSE JPL 8 day and JPL Orbital")
print(rmse_jpl.values)


