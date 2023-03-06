#!/usr/bin/env python3
# coding: utf-8

#import xarray as xr
#from pathlib import Path
import numpy as np

#TODO:
#capture the user input if any passed, think about it?
#make some constrains to the user input oustside the class like only 01 and 04 etc.

class AccessOpendap18():
    
    base_url = ('https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/{0}/{1}'
                '/{2}/woa18_{1}_{3}{4}_{5}.nc')
    
    def __init__(self, period=None, resolution=None, time=None):
        
        '''Constructer to assign variables and setup desired parameters.'''
        
        self.url_dict = {"variable_long" : ["temperature", "salinity"],
                    "period" : ["decav", "A5B7"],
                    "time" : ["00", "01", "02", "03", "04", "05", "06", "07", "08",
                              "09", "10", "11", "12", "13", "14", "15", "16"],
                    "resolution" : ["01", "04"]}
        self._set_arguments(period, resolution, time)
        self._fill_url()
 
    def _set_arguments(self, period, resolution, time):
        
        '''Method to check the user input and update some attributes.'''

        self.period = ("period", period)
        self.resolution = ("resolution", resolution)
        self.time = ("time", time) 
        userargs = np.array([self.period, self.resolution, self.time],
                            dtype=object)
        userargs_mask = [userarg[-1]  is not None for userarg in userargs]
        if any(userargs_mask):
            self._update_dict(userargs[userargs_mask]) 
        
    def _update_dict(self, userarg_list):
        
        '''Method to update the url_dict based on user specified parameters.'''
        
        for name, value in userarg_list:
            if type(value) == list:
                self.url_dict[name] = value
            else:
                self.url_dict[name] = [value]
    
    def _fill_url(self):
        
        #TODO: replace with itertools?
        '''Method to fill the base_url with the fitting arguments from the url_dict.'''

        self.url_list = []
        for time in self.url_dict["time"]:
            for period in self.url_dict["period"]:
                for resolution_int in self.url_dict["resolution"]:
                    if resolution_int == "01":
                        resolution_float = "1.00"
                    elif resolution_int == "04":
                        resolution_float = "0.25"
                    for parameter in self.url_dict["variable_long"]:
                        parameter_short = parameter[0]
                        
                        url = self.base_url.format(parameter, period,
                                             resolution_float, 
                                             parameter_short, time,
                                             resolution_int)
                        self.url_list.append(url)
        self.url_list = np.reshape(self.url_list, (-1,2))                           
         
    
    
    