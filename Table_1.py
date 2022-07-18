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



