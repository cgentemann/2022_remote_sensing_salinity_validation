#!/usr/bin/env python
##This script does calculations on hycom and saildrone data
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import datetime
import math
from scipy import stats
import xskillscore as xs
import sys
import csv

#data directory for saildrone and satellite orbital data
data_sat = './../paper_software/2020_ATOMIC_Salinity/data/sss_collocations_orbital/'

#data directory for HYCOM data
data_dir1 = './../paper_software/2020_ATOMIC_Salinity/data/'
files = glob(data_dir1+'*nc4')
hycom=xr.open_mfdataset(data_dir1+'*nc4',concat_dim='time').isel(depth=0)

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

#Collocated Saildrone and Satellite data
SMAP = [x for x in glob(data_sat+'*.nc')]

for x in SMAP:
    usv_2 = xr.open_dataset(x)
    usv_2.close()
    num = usv_2.trajectory.values
    sail = usv_2

    #interp HYCOM data to saildrone dimensions
    hysal = filled3.interp(lat=sail.lat, longitude=sail.lon, time=sail.time, method='nearest')
    hysal2 = hysal.salinity

    #Comparing RSS SMAP data with saildrone data
    if np.sum(usv_2.smap_name.str.contains('RSS')):

        # RMSE
        rmse_RSS_sbe = xs.rmse(sail.SAL_CTD_MEAN, sail.smap_SSS, dim="time", skipna="True")
        rmse_RSS_rbr = xs.rmse(sail.SAL_RBR_MEAN, sail.smap_SSS, dim="time", skipna="True")

        # Mean Diference
        mean_sbe = sail.SAL_CTD_MEAN.mean(dim='time')
        mean_rbr = sail.SAL_RBR_MEAN.mean(dim='time')
        mean_RSS = sail.smap_SSS.mean(dim='time')

        # Difference in Standard Deviation
        std_rbr = sail["SAL_RBR_MEAN"].std(dim="time", skipna=None)
        std_sbe37 = sail["SAL_CTD_MEAN"].std(dim="time", skipna=None)
        std_RSS = sail["smap_SSS"].std(dim="time", skipna=None)
        std_RSS_sbe = (std_sbe37.values - std_RSS.values)
        std_RSS_rbr = (std_rbr.values - std_RSS.values)

        with open('RSS_HYCOM_' + str(num) + '.csv', 'w') as csvfile1:
            fieldnames = ['stat', 'value']
            writer = csv.DictWriter(csvfile1, fieldnames)
            writer.writeheader()
            writer.writerow({'stat': 'RSS and Saildrone', 'value': str(num)})
            writer.writerow({'stat': 'Saildrone Number', 'value': str(num)})
            writer.writerow({'stat': 'Standard Deviation RBR - RSS', 'value': str(std_RSS_rbr)})
            writer.writerow({'stat': 'Standard Deviation SBE37 - RSS', 'value': str(std_RSS_sbe)})
            writer.writerow({'stat': 'RMSE SBE37 and RSS', 'value': str(rmse_RSS_sbe.values)})
            writer.writerow({'stat': 'RMSE RBR and RSS', 'value': str(rmse_RSS_rbr.values)})
            writer.writerow({'stat': 'Mean Difference RBR - RSS', 'value': str(mean_rbr.values - mean_RSS.values)})
            writer.writerow({'stat': 'Mean Difference SBE37 - RSS', 'value': str(mean_sbe.values - mean_RSS.values)})


            #Comparing HYCOM to saildrone and RSS data
            ##Saildrone
            hycom_std = hysal2.std(dim='time')
            sbe_std = (std_sbe37.values - hycom_std.values)
            rbr_std = (std_rbr.values - hycom_std.values)

            # RMSE
            rmse_sal_sbe = xs.rmse(sail.SAL_CTD_MEAN, hysal2, dim="time", skipna="True")
            rmse_sal_rbr = xs.rmse(sail.SAL_RBR_MEAN, hysal2, dim="time", skipna="True")

            # Mean Diference
            mean_sbe = sail.SAL_CTD_MEAN.mean(dim='time')
            mean_rbr = sail.SAL_RBR_MEAN.mean(dim='time')
            hycom_mean = hysal2.mean(dim='time')

            writer.writerow({'stat': 'HYCOM and Saildrone', 'value': str(num)})
            writer.writerow({'stat': 'Saildrone Number', 'value': str(num)})
            writer.writerow({'stat': 'Standard Deviation RBR - HYCOM', 'value': str(rbr_std)})
            writer.writerow({'stat': 'Standard Deviation SBE37 - HYCOM', 'value': str(sbe_std)})
            writer.writerow({'stat': 'RMSE SBE37 and HYCOM', 'value': str(rmse_sal_sbe.values)})
            writer.writerow({'stat': 'RMSE RBR and HYCOM', 'value': str(rmse_sal_rbr.values)})
            writer.writerow({'stat': 'Mean Difference RBR - HYCOM', 'value': str(mean_rbr.values - hycom_mean.values)})
            writer.writerow({'stat': 'Mean Difference SBE37 - HYCOM', 'value': str(mean_sbe.values - hycom_mean.values)})


            ##Comparing RSS and HYCOM
            # STD
            hycom_std = hysal2.std(dim='time')
            std_RSS = sail["smap_SSS"].std(dim="time", skipna=None)
            RSS_std = (std_RSS.values - hycom_std.values)

            # RMSE
            rmse_hycom_rss = xs.rmse(sail.smap_SSS, hysal2, dim="time", skipna="True")

            # Mean Difference
            hycom_mean = hysal2.mean(dim='time')
            mean_RSS = sail.smap_SSS.mean(dim='time')

            writer.writerow({'stat': 'HYCOM and RSS'})
            writer.writerow({'stat': 'Saildrone Number', 'value': str(num)})
            writer.writerow({'stat': 'Standard Deviation RSS - HYCOM', 'value': str(RSS_std)})
            writer.writerow({'stat': 'RMSE HYCOM and RSS', 'value': str(rmse_hycom_rss.values)})
            writer.writerow({'stat': 'Mean Difference RSS - HYCOM', 'value': str(mean_RSS.values - hycom_mean.values)})

    else:
        # Comparing JPL SMAP data with saildrone data
        # RMSE
        print(x)
        rmse_JPL_sbe = xs.rmse(sail.SAL_CTD_MEAN, sail.smap_SSS, dim="time", skipna="True")
        rmse_JPL_rbr = xs.rmse(sail.SAL_RBR_MEAN, sail.smap_SSS, dim="time", skipna="True")

        # Mean Difference
        mean_sbe = sail.SAL_CTD_MEAN.mean(dim='time')
        mean_rbr = sail.SAL_RBR_MEAN.mean(dim='time')
        mean_JPL = sail.smap_SSS.mean(dim='time')

        # Difference in Standard Deviation
        std_rbr = sail["SAL_RBR_MEAN"].std(dim="time", skipna=None)
        std_sbe37 = sail["SAL_CTD_MEAN"].std(dim="time", skipna=None)
        std_JPL = sail["smap_SSS"].std(dim="time", skipna=None)
        std_JPL_sbe = (std_sbe37.values - std_JPL.values)
        std_JPL_rbr = (std_rbr.values - std_JPL.values)

        with open('JPL_Stats' + str(num) + '.csv', 'w') as csvfile1:
            fieldnames = ['stat', 'value']
            writer = csv.DictWriter(csvfile1, fieldnames)
            writer.writeheader()

            writer.writerow({'stat': 'JPL and Saildrone', 'value': str(num)})
            writer.writerow({'stat': 'Saildrone Number', 'value': str(num)})
            writer.writerow({'stat': 'Standard Deviation RBR - JPL', 'value': str(std_JPL_rbr)})
            writer.writerow({'stat': 'Standard Deviation SBE37 - JPL', 'value': str(std_JPL_sbe)})
            writer.writerow({'stat': 'RMSE SBE37 and JPL', 'value': str(rmse_JPL_sbe.values)})
            writer.writerow({'stat': 'RMSE RBR and JPL', 'value': str( rmse_JPL_rbr.values)})
            writer.writerow({'stat': 'Mean Difference RBR - JPL', 'value': str(mean_rbr.values - mean_JPL.values)})
            writer.writerow({'stat': 'Mean Difference SBE37 - JPL', 'value': str(mean_sbe.values - mean_JPL.values)})

            #Comparing JPL Data to HYCOM data
            hycom_std = hysal2.std(dim='time')
            JPL_std = (std_JPL.values - hycom_std.values)

            # RMSE
            rmse_hycom_jpl = xs.rmse(sail.smap_SSS, hysal2, dim="time", skipna="True")

            # Mean Difference
            hycom_mean = hysal2.mean(dim='time')
            mean_JPL = sail.smap_SSS.mean(dim='time')
            JPL_hy_diff = mean_JPL.values - hycom_mean.values

            writer.writerow({'stat':'HYCOM and JPL'})
            writer.writerow({'stat': 'Standard Deviation JPL-HYCOM', 'value': str(std_JPL_rbr)})
            writer.writerow({'stat': 'RMSE HYCOM and JPL', 'value': str(rmse_hycom_jpl.values)})
            writer.writerow({'stat': 'Mean Difference JPL - HYCOM', 'value': str(JPL_hy_diff)})



