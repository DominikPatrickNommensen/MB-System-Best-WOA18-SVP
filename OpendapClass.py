#!/usr/bin/env python3
# coding: utf-8

#import xarray as xr
#from pathlib import Path
import numpy as np

class AccessOpendap18():
    
    base_url = ('https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/{0}/{1}'
                '/{2}/woa18_{1}_{3}{4}_{5}.nc')
    
    def __init__(self, period, resolution, time):
        '''Constructer to assign variables and setup desired parameters.'''
        
        self.url_dict = {"variable_long" : ["temperature", "salinity"],
                    "period" : period,
                    "time" : time,
                    "resolution" : resolution}
        self._fill_url()
     
    def _update_dict(self, userarg_list):
        '''Method to update the url_dict based on user specified parameters.'''
        
        for name, value in userarg_list:
            if type(value) == list:
                self.url_dict[name] = value
            else:
                self.url_dict[name] = [value]
    
    def _fill_url(self):
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
         
    
    
    