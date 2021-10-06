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
data_sat = 'C:/Users/intern-1/Documents/GitHub/paper_software/2020_ATOMIC_Salinity/data/sss_collocations_orbital_norepeat/'
data_sat8 = 'C:/Users/intern-1/Documents/GitHub/paper_software/2020_ATOMIC_Salinity/data/sss_collocations_8day_nearest_norepeat/'

#data directory for HYCOM data
data_dir1 = 'C:/Users/intern-1/Documents/hycom_files/'
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
SMAP = [x for x in glob(data_sat+'*1026*.nc')]

#Uncomment to do comparison with 8 day average data
SMAP8 = [x for x in glob(data_sat8+'*1026*.nc')]
print(SMAP8)
for x in SMAP8:
    usv_2 = xr.open_dataset(x, decode_times=False)
    usv_2.close()
    num = usv_2.trajectory.values
    sail = usv_2
    sail = sail.swap_dims({'ob': 'time'})

    # Convert times from nano seconds to seconds
    ns = 1e-9
    sat_time = sail.time.values * ns
    test = sat_time.astype(np.float)

    ss_times = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test]
    sail.assign_coords(time=ss_times)

    # interp HYCOM data to saildrone dimensions
    hysal = filled3.interp(lat=sail.lat, longitude=sail.lon, time=ss_times, method='nearest')
    hysal2 = hysal.salinity

    # mean_sbe_RSS = RSS.SAL_CTD_MEAN
    # mean_RSS_40km = RSS.smap_SSS_40km

    # Comparing RSS SMAP data with saildrone data
    if "RSS" in x:
        print(x)
        #hysal2 = sail.sat_anc_sss
        SBE = sail.SAL_CTD_MEAN
        RBR = sail.SAL_RBR_MEAN
        RSS = sail.sat_sss_smap
        RSS_40 = sail.sat_sss_smap_40km

        # RMSE
        rmse_RSS_sbe = xs.rmse(SBE, RSS, dim="time", skipna="True")
        rmse_RSS_rbr = xs.rmse(RBR, RSS, dim="time", skipna="True")
        rmse_RSS_40 = xs.rmse(SBE, RSS_40, dim="time", skipna="True")

        # Mean Difference
        RBR_RSS=(RBR - RSS).mean(dim='time').values
        SBE_RSS=(SBE - RSS).mean(dim='time').values
        SBE_RSS40=(SBE - RSS_40).mean(dim='time').values

        # Difference in Standard Deviation
        RBR_RSS_std=(RBR - RSS).std(dim='time').values
        SBE_RSS_std=(SBE - RSS).std(dim='time').values
        SBE_RSS40_std=(SBE - RSS_40).std(dim='time').values

        with open('RSS_8DAY_Stats'
                  '' + str(num) + '.csv', 'w') as csvfile1:
            fieldnames = ['stat', 'value']
            writer = csv.DictWriter(csvfile1, fieldnames)
            writer.writeheader()
            writer.writerow({'stat': 'RSS and Saildrone', 'value': str(num)})
            writer.writerow({'stat': 'Saildrone Number', 'value': str(num)})
            writer.writerow({'stat': 'Standard Deviation RBR - RSS', 'value': str(RBR_RSS_std)})
            writer.writerow({'stat': 'Standard Deviation SBE37 - RSS', 'value': str(SBE_RSS_std)})
            writer.writerow({'stat': 'Standard Deviation SBE37 - RSS 40 km', 'value': str(SBE_RSS40_std)})

            writer.writerow({'stat': 'RMSE SBE37 and RSS', 'value': str(rmse_RSS_sbe.values)})
            writer.writerow({'stat': 'RMSE RBR and RSS', 'value': str(rmse_RSS_rbr.values)})
            writer.writerow({'stat': 'RMSE SBE37 and RSS 40 km', 'value': str(rmse_RSS_40.values)})

            writer.writerow({'stat': 'Mean Difference RBR - RSS', 'value': str(RBR_RSS)})
            writer.writerow({'stat': 'Mean Difference SBE37 - RSS', 'value': str(SBE_RSS)})
            writer.writerow({'stat': 'Mean Difference SBE37 - RSS 40km', 'value': str(SBE_RSS40)})

            # Comparing HYCOM to saildrone and RSS data
            ##Saildrone
            hycom_std_sbe = (SBE-hysal2).std(dim='time').values
            hycom_std_rbr = (RBR - hysal2).std(dim='time').values

            # RMSE
            rmse_sal_sbe = xs.rmse(SBE, hysal2, dim="time", skipna="True")
            rmse_sal_rbr = xs.rmse(RBR, hysal2, dim="time", skipna="True")

            # Mean Diference
            RBR_HYC = (RBR - hysal2).mean(dim='time').values
            SBE_HYC = (SBE - hysal2).mean(dim='time').values

            writer.writerow({'stat': 'HYCOM and Saildrone', 'value': str(num)})
            writer.writerow({'stat': 'Saildrone Number', 'value': str(num)})
            writer.writerow({'stat': 'Standard Deviation RBR - HYCOM', 'value': str(hycom_std_rbr)})
            writer.writerow({'stat': 'Standard Deviation SBE37 - HYCOM', 'value': str(hycom_std_sbe)})
            writer.writerow({'stat': 'RMSE SBE37 and HYCOM', 'value': str(rmse_sal_sbe.values)})
            writer.writerow({'stat': 'RMSE RBR and HYCOM', 'value': str(rmse_sal_rbr.values)})
            writer.writerow({'stat': 'Mean Difference RBR - HYCOM', 'value': str(RBR_HYC)})
            writer.writerow({'stat': 'Mean Difference SBE37 - HYCOM', 'value': str(SBE_HYC)})

            ##Comparing RSS and HYCOM
            # STD
            RSS_std = (RSS - hysal2).std(dim='time').values

            # RMSE
            rmse_hycom_rss = xs.rmse(RSS, hysal2, dim="time", skipna="True")

            # Mean Difference
            RSS_HYC = (RSS - hysal2).mean(dim='time').values

            writer.writerow({'stat': 'HYCOM and RSS'})
            writer.writerow({'stat': 'Saildrone Number', 'value': str(num)})
            writer.writerow({'stat': 'Standard Deviation RSS - HYCOM', 'value': str(RSS_std)})
            writer.writerow({'stat': 'RMSE HYCOM and RSS', 'value': str(rmse_hycom_rss.values)})
            writer.writerow({'stat': 'Mean Difference RSS - HYCOM', 'value': str(RSS_HYC)})

    elif "JPL" in x:

        SBE = sail.SAL_CTD_MEAN
        RBR = sail.SAL_RBR_MEAN
        JPL = sail.sat_smap_sss
        # Comparing JPL SMAP data with saildrone data
        # RMSE
        rmse_JPL_sbe = xs.rmse(SBE, JPL, dim="time", skipna="True")
        rmse_JPL_rbr = xs.rmse(RBR, JPL, dim="time", skipna="True")

        # Mean Difference
        RBR_JPL=(RBR - JPL).mean(dim='time').values
        SBE_JPL=(SBE - JPL).mean(dim='time').values

        # Difference in Standard Deviation
        RBR_JPL_std=(RBR - JPL).std(dim='time').values
        SBE_JPL_std=(SBE - JPL).std(dim='time').values

        with open('JPL_8day_Stats' + str(num) + '.csv', 'w') as csvfile1:
            fieldnames = ['stat', 'value']
            writer = csv.DictWriter(csvfile1, fieldnames)
            writer.writeheader()

            writer.writerow({'stat': 'JPL and Saildrone', 'value': str(num)})
            writer.writerow({'stat': 'Saildrone Number', 'value': str(num)})
            writer.writerow({'stat': 'Standard Deviation RBR - JPL', 'value': str(RBR_JPL_std)})
            writer.writerow({'stat': 'Standard Deviation SBE37 - JPL', 'value': str(SBE_JPL_std)})
            writer.writerow({'stat': 'RMSE SBE37 and JPL', 'value': str(rmse_JPL_sbe.values)})
            writer.writerow({'stat': 'RMSE RBR and JPL', 'value': str(rmse_JPL_rbr.values)})
            writer.writerow({'stat': 'Mean Difference RBR - JPL', 'value': str(RBR_JPL)})
            writer.writerow({'stat': 'Mean Difference SBE37 - JPL', 'value': str(SBE_JPL)})

            # Comparing JPL Data to HYCOM data
            #STD
            JPL_HYC_std = (JPL - hysal2).std(dim='time').values

            # RMSE
            rmse_hycom_jpl = xs.rmse(sail.sat_smap_sss, hysal2, dim="time", skipna="True")

            # Mean Difference
            JPL_HYC= (JPL - hysal2).mean(dim='time').values

            writer.writerow({'stat': 'HYCOM and JPL'})
            writer.writerow({'stat': 'Standard Deviation JPL-HYCOM', 'value': str(JPL_HYC_std)})
            writer.writerow({'stat': 'RMSE HYCOM and JPL', 'value': str(rmse_hycom_jpl.values)})
            writer.writerow({'stat': 'Mean Difference JPL - HYCOM', 'value': str(JPL_HYC)})

for x in SMAP:
    print("SMAP")
    usv_2 = xr.open_dataset(x, decode_times=False)
    usv_2.close()
    num = usv_2.trajectory.values
    sail = usv_2
    sail = sail.swap_dims({'ob': 'time'})

    # Convert times from nano seconds to seconds
    ns = 1e-9
    sat_time = sail.time.values * ns
    test = sat_time.astype(np.float)

    ss_times = [datetime.datetime(2020, 1, 17, 0, 0) + datetime.timedelta(seconds=s) for s in test]
    sail.assign_coords(time=ss_times)

    # interp HYCOM data to saildrone dimensions
    hysal = filled3.interp(lat=sail.lat, longitude=sail.lon, time=ss_times, method='nearest')
    hysal2 = hysal.salinity

    #Comparing RSS SMAP data with saildrone data
    if "rss" in x:
        print(x)
        SBE = sail.SAL_CTD_MEAN
        RBR = sail.SAL_RBR_MEAN
        RSS = sail.smap_SSS
        RSS_40 = sail.smap_SSS_40km

        # RMSE
        rmse_RSS_sbe = xs.rmse(SBE, RSS, dim="time", skipna="True")
        rmse_RSS_rbr = xs.rmse(RBR, RSS, dim="time", skipna="True")
        rmse_RSS_40 = xs.rmse(SBE, RSS_40, dim="time", skipna="True")

        # Mean Difference
        RBR_RSS = (RBR - RSS).mean(dim='time').values
        SBE_RSS = (SBE - RSS).mean(dim='time').values
        SBE_RSS40 = (SBE - RSS_40).mean(dim='time').values

        # Difference in Standard Deviation
        RBR_RSS_std = (RBR - RSS).std(dim='time').values
        SBE_RSS_std = (SBE - RSS).std(dim='time').values
        SBE_RSS40_std = (SBE - RSS_40).std(dim='time').values

        with open('RSS_orb_Stats'
                  '' + str(num) + '.csv', 'w') as csvfile1:
            fieldnames = ['stat', 'value']
            writer = csv.DictWriter(csvfile1, fieldnames)
            writer.writeheader()
            writer.writerow({'stat': 'RSS and Saildrone', 'value': str(num)})
            writer.writerow({'stat': 'Saildrone Number', 'value': str(num)})
            writer.writerow({'stat': 'Standard Deviation RBR - RSS', 'value': str(RBR_RSS_std)})
            writer.writerow({'stat': 'Standard Deviation SBE37 - RSS', 'value': str(SBE_RSS_std)})
            writer.writerow({'stat': 'Standard Deviation SBE37 - RSS 40km', 'value': str(SBE_RSS40_std)})

            writer.writerow({'stat': 'RMSE SBE37 and RSS', 'value': str(rmse_RSS_sbe.values)})
            writer.writerow({'stat': 'RMSE RBR and RSS', 'value': str(rmse_RSS_rbr.values)})
            writer.writerow({'stat': 'RMSE SBE37 and RSS 40km', 'value': str(rmse_RSS_40.values)})

            writer.writerow({'stat': 'Mean Difference RBR - RSS', 'value': str(RBR_RSS)})
            writer.writerow({'stat': 'Mean Difference SBE37 - RSS', 'value': str(SBE_RSS)})
            writer.writerow({'stat': 'Mean Difference SBE37 - RSS 40km', 'value': str(SBE_RSS40)})

            # Comparing HYCOM to saildrone and RSS data
            ##Saildrone
            hycom_std_sbe = (SBE-hysal2).std(dim='time').values
            hycom_std_rbr = (RBR - hysal2).std(dim='time').values

            # RMSE
            rmse_sal_sbe = xs.rmse(SBE, hysal2, dim="time", skipna="True")
            rmse_sal_rbr = xs.rmse(RBR, hysal2, dim="time", skipna="True")

            # Mean Diference
            RBR_HYC = (RBR - hysal2).mean(dim='time').values
            SBE_HYC = (SBE - hysal2).mean(dim='time').values

            writer.writerow({'stat': 'HYCOM and Saildrone', 'value': str(num)})
            writer.writerow({'stat': 'Saildrone Number', 'value': str(num)})
            writer.writerow({'stat': 'Standard Deviation RBR - HYCOM', 'value': str(hycom_std_rbr)})
            writer.writerow({'stat': 'Standard Deviation SBE37 - HYCOM', 'value': str(hycom_std_sbe)})
            writer.writerow({'stat': 'RMSE SBE37 and HYCOM', 'value': str(rmse_sal_sbe.values)})
            writer.writerow({'stat': 'RMSE RBR and HYCOM', 'value': str(rmse_sal_rbr.values)})
            writer.writerow({'stat': 'Mean Difference RBR - HYCOM', 'value': str(RBR_HYC)})
            writer.writerow({'stat': 'Mean Difference SBE37 - HYCOM', 'value': str(SBE_HYC)})

            ##Comparing RSS and HYCOM
            # STD
            RSS_std = (RSS - hysal2).std(dim='time').values

            # RMSE
            rmse_hycom_rss = xs.rmse(RSS, hysal2, dim="time", skipna="True")

            # Mean Difference
            RSS_HYC = (RSS - hysal2).mean(dim='time').values

            writer.writerow({'stat': 'HYCOM and RSS'})
            writer.writerow({'stat': 'Saildrone Number', 'value': str(num)})
            writer.writerow({'stat': 'Standard Deviation RSS - HYCOM', 'value': str(RSS_std)})
            writer.writerow({'stat': 'RMSE HYCOM and RSS', 'value': str(rmse_hycom_rss.values)})
            writer.writerow({'stat': 'Mean Difference RSS - HYCOM', 'value': str(RSS_HYC)})

    else:
        # Comparing JPL SMAP data with saildrone data
        # RMSE
        print(x)
        SBE = sail.SAL_CTD_MEAN
        RBR = sail.SAL_RBR_MEAN
        JPL = sail.smap_SSS
        # Comparing JPL SMAP data with saildrone data
        # RMSE
        rmse_JPL_sbe = xs.rmse(SBE, JPL, dim="time", skipna="True")
        rmse_JPL_rbr = xs.rmse(RBR, JPL, dim="time", skipna="True")

        # Mean Difference
        RBR_JPL=(RBR - JPL).mean(dim='time').values
        SBE_JPL=(SBE - JPL).mean(dim='time').values

        # Difference in Standard Deviation
        RBR_JPL_std=(RBR - JPL).std(dim='time').values
        SBE_JPL_std=(SBE - JPL).std(dim='time').values

        with open('JPL_orb_Stats' + str(num) + '.csv', 'w') as csvfile1:
            fieldnames = ['stat', 'value']
            writer = csv.DictWriter(csvfile1, fieldnames)
            writer.writeheader()

            writer.writerow({'stat': 'JPL and Saildrone', 'value': str(num)})
            writer.writerow({'stat': 'Saildrone Number', 'value': str(num)})
            writer.writerow({'stat': 'Standard Deviation RBR - JPL', 'value': str(RBR_JPL_std)})
            writer.writerow({'stat': 'Standard Deviation SBE37 - JPL', 'value': str(SBE_JPL_std)})
            writer.writerow({'stat': 'RMSE SBE37 and JPL', 'value': str(rmse_JPL_sbe.values)})
            writer.writerow({'stat': 'RMSE RBR and JPL', 'value': str( rmse_JPL_rbr.values)})
            writer.writerow({'stat': 'Mean Difference RBR - JPL', 'value': str(RBR_JPL)})
            writer.writerow({'stat': 'Mean Difference SBE37 - JPL', 'value': str(SBE_JPL)})

            #Comparing JPL Data to HYCOM data
            HYC_JPL_std=(JPL - hysal2).std(dim='time').values

            # RMSE
            rmse_hycom_jpl = xs.rmse(sail.smap_SSS, hysal2, dim="time", skipna="True")

            # Mean Difference
            JPL_HYC = (JPL - hysal2).mean(dim='time').values

            writer.writerow({'stat':'HYCOM and JPL'})
            writer.writerow({'stat': 'Standard Deviation JPL-HYCOM', 'value': str(HYC_JPL_std)})
            writer.writerow({'stat': 'RMSE HYCOM and JPL', 'value': str(rmse_hycom_jpl.values)})
            writer.writerow({'stat': 'Mean Difference JPL - HYCOM', 'value': str(JPL_HYC)})

