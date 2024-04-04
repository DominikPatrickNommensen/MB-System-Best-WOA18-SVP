#!/usr/bin/env python3
# coding: utf-8

from pathlib import Path
import subprocess
import xarray as xr
import numpy as np
import pandas as pd

class SwathFileInfo():
    '''Should only get mbinfo and return the files needed for the swathfile of the project'''
    
    def __init__(self, period='both', resolution='both', svpfolder=None, 
                 swathfile='datalist.mb-1', nearest=False):
        
        self.svpfolder = Path(svpfolder)
        self.nearest = nearest
        self.base_file = 'woa18_{0}_svp{1}_{2}.nc'
        self.swathfile = swathfile
        self.resolution = resolution
        self.period = period
         
    def mbinfo(self, end=False):
        
        info = subprocess.run(["mbinfo", "-I {}".format(self.swathfile)],
                              capture_output=True, text=True).stdout

        self.latmin_swathfile = float(info.split(
                                            "Minimum Latitude:")[1].split()[0])
        self.latmax_swathfile = float(info.split(
                                            "Maximum Latitude:")[1].split()[0])
        self.lonmin_swathfile = float(info.split(
                                            "Minimum Longitude:")[1].split()[0])
        self.lonmax_swathfile = float(info.split(
                                            "Maximum Longitude:")[1].split()[0])
    
        self.maximum_depth = float(info.split("Maximum Depth:")[1].split()[0])
             
        
        if not end:    
            self.month = info.split("Start of Data:")[1].split()[1]
        else:
            self.month = info.split("End of Data:")[1].split()[1]
   
        if self.month in ["01", "02", "03"]:
            self.season = "13"
        elif self.month in ["04", "05", "06"]:
            self.season = "14"
        elif self.month in ["07", "08", "09"]:
            self.season = "15"
        else:
            self.season = "16"
    
        self.annual = "00"
        
        self.time = [self.season, self.month, self.annual]
   
        
    def extract_svp(self):
        
        subprocess.run(["mbsvplist", "-I {}".format(self.swathfile), "-O", "-V"])
        
    def get_filenames(self):
        
        filenames = []
        for resolution in self.resolution:
            for period in self.period:
                for time in self.time:
                    filenames.append(self.svpfolder / self.base_file.format(period, time, 
                                                           resolution))
                    
        return filenames
    
    def convert_lon_to_360(self, x):
        return x % 360
    
    def calculate_mean_profiles(self, file):
        
        ds = xr.open_dataset(file)
        if any(ds.lon > 180):
            self.lonmin_swathfile = self.convert_lon_to_360(self.lonmin_swathfile)
            self.lonmax_swathfile = self.convert_lon_to_360(self.lonmax_swathfile)
            
        if self.nearest:
            lat_mean = np.mean([self.latmin_swathfile, self.latmax_swathfile])        
            lon_mean = np.mean([self.lonmin_swathfile, self.lonmax_swathfile])
            
            ds = xr.open_dataset(file).sel(lat=lat_mean, lon=lon_mean, 
                                           method='nearest')
            self.latmin = ds.lat.data
            self.latmax = self.latmin
            self.lonmin = ds.lon.data
            self.lonmax = self.lonmin
            
        else:
            ds = xr.open_dataset(file)
            ffill_bounds = ds.sel(lat=self.latmin_swathfile,
                                  lon=self.lonmin_swathfile, method='ffill')
            bfill_bounds = ds.sel(lat=self.latmax_swathfile,
                                  lon=self.lonmax_swathfile, method='bfill')
    
            self.latmin = ffill_bounds.lat.data
            self.lonmin = ffill_bounds.lon.data
            self.latmax = bfill_bounds.lat.data
            self.lonmax = bfill_bounds.lon.data
            
            ds = ds.sel(lat=slice(self.latmin, self.latmax),
                        lon=slice(self.lonmin, self.lonmax))
            
            ds = ds.mean(dim=['lat', 'lon'], skipna=True)
        
        df = ds.to_dataframe()
        
        #the append part
        if (file.stem[-5:-3] == self.season) and (self.maximum_depth > 1500):
            self.season_to_append = df.loc[1550:]
        elif (file.stem[-5:-3] == self.month) and (self.maximum_depth > 1500):
            df = pd.concat([df, self.season_to_append])
        df.dropna(inplace=True)
        df.loc["12000"] = 1670.864
        df.index = df.index.astype("int")
        
        return df