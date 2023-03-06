#!/usr/bin/env python3
# coding: utf-8

import xarray as xr
import numpy as np
import gc
import subprocess
from pathlib import Path
import pandas as pd
from io import StringIO
from numpy.polynomial.polynomial import polyfit
import matplotlib.pyplot as plt

class CalculateSvp():
    '''Class that handles the calculation of pressure and different sound 
    velocities using the World Ocean Atlas 2018 temperature and salinity files.
    '''
    c00 =  1402.388  
    c01 =  5.03830   
    c02 = -5.81090e-2
    c03 =  3.3432e-4 
    c04 = -1.47797e-6
    c05 =  3.1419e-9 
    c10 =  0.153563  
    c11 =  6.8999e-4 
    c12 = -8.1829e-6 
    c13 =  1.3632e-7 
    c14 = -6.1260e-10
    c20 =  3.1260e-5 
    c21 = -1.7111e-6 
    c22 =  2.5986e-8 
    c23 = -2.5353e-10
    c24 =  1.0415e-12
    c30 = -9.7729e-9 
    c31 =  3.8513e-10
    c32 = -2.3654e-12
    
    a00 =  1.389     
    a01 = -1.262e-2  
    a02 =  7.166e-5  
    a03 =  2.008e-6  
    a04 = -3.21e-8   
    a10 =  9.4742e-5
    a11 = -1.2583e-5 
    a12 = -6.4928e-8
    a13 =  1.0515e-8
    a14 = -2.0142e-10
    a20 = -3.9064e-7 
    a21 =  9.1061e-9
    a22 = -1.6009e-10
    a23 =  7.994e-12  
    a30 =  1.100e-10 
    a31 =  6.651e-12 
    a32 = -3.391e-13 
    
    b00 = -1.922e-2  
    b01 = -4.42e-5   
    b10 =  7.3637e-5 
    b11 =  1.7950e-7 
    
    d00 =  1.727e-3  
    d10 = -7.9836e-6
    
    
    def __init__(self, lon_lat, url_list, corrective_term="NONE"):
        '''Sets the initial parameters.
        
        Parameters
        ----------
        lon_lat : string or list of 4 strings
            Longitude und latitude extent 4 float 
            list with min/max lon/lat.
            
        url_list : nested list of string pairs
            Nested list of URL World Ocean Atlas temperature and salinty pairs.
            
        corrective_term : string, optional
            Corrective term to use for pressure calculation. 
            The default is "NONE".

        Returns
        -------
        None.
        '''

        self.temp_file, self.sali_file = url_list
        self.lolon, self.uplon, self.lolat, self.uplat = lon_lat

        self.corrective_term = corrective_term
        self._term_dict = {"NONE" : self.d2p_without_corrective_terms,
                          "COM" : self.d2p_60N_40S_correction,
                          "NEA" : self.d2p_North_Eastern_Atlantic_correction,
                          "CPAW" : self.d2p_CPAW_correction,
                          "MED" : self.d2p_Mediterranean_Sea_correction,
                          "JAP" : self.d2p_Sea_of_Japan_correction,
                          "SULU" : self.d2p_Sulu_Sea_correction,
                          "HAMA" : self.d2p_Halmahera_Basin_correction,
                          "CELE" : self.d2p_Celebes_Basin_Weber_Sea_correction,
                          "BLACK" : self.d2p_Black_Sea_correction,
                          "BALTIC" : self.d2p_Baltic_Sea_correction 
                          }
        
        try:
            self._pressure_formula = self._term_dict[corrective_term]
        except KeyError:
            print("No corrective term will be applied. A corrective term was "
                  "given that does not exist. Must be one of {}".format(
                  list(self._term_dict.keys())))
            self._pressure_formula = self._term_dict["NONE"]
            
    def __h_Z45(self, Z):
        """Returns the Leroy and Parthiot 1998 depth to pressure term (9)""" 
        return (1.00818 * 10**-2 * Z + 2.465 * 10**-8 * Z**2 -1.25 * 10**-13 
                * Z**3 + 2.8 * 10**-19 * Z**4)
    
    def __g_Phi(self, L):
        """Returns the Leroy and Parthiot 1998 depth to pressure term (11)"""
        return (9.7803 * (1 + 5.3 * 10**-3 * (np.sin(L*np.pi/180))**2))
    
    def __k_ZPhi(self, Z, L):
        """"Returns the Leroy and Parthiot 1998 depth to pressure term (10)"""
        return ((self.__g_Phi(L) - 2 * 10**-5 * Z)/(9.80612 - 2 * 10**-5 * Z))
    
    def d2p_without_corrective_terms(self, Z, L):
        """Returns the result of the Leroy and Parthiot 1998 depth to pressure 
        formula without corrective term
        """
        return (self.__h_Z45(Z) * self.__k_ZPhi(Z, L) *10) 
    
    def d2p_60N_40S_correction(self, Z, L):
        """Returns the result of the Leroy and Parthiot 1998 depth to pressure 
        formula with 60°N to 40°S corrective term for common ocean
        """
        return ((self.__h_Z45(Z) * self.__k_ZPhi(Z, L) - (1.0 * 10**-2 * Z 
                                                         / (Z + 100) + 6.2 
                                                         * 10**-6 * Z))*10)
    
    def d2p_North_Eastern_Atlantic_correction(self, Z, L):
        """Returns the result of the Leroy and Parthiot 1998 depth to pressure 
        formula with corrective term for the Northestern_Atlantic 
        (30-35° latitude)
        """
        return ((self.__h_Z45(Z) * self.__k_ZPhi(Z, L) - (8 * 10**-3 * Z 
                                                         / (Z + 200) + 4 
                                                         * 10**-6 * Z))*10)

    def d2p_CPAW_correction(self, Z, L):
        """Returns the result of the Leroy and Parthiot 1998 depth to pressure 
        formula with corrective term for the circumpolar antartic water
        """
        return ((self.__h_Z45(Z) * self.__k_ZPhi(Z, L) - (8 * 10**-3 * Z 
                                                         / (Z + 1000) + 1.6 
                                                         * 10**-6 * Z))*10)

    def d2p_Mediterranean_Sea_correction(self, Z, L):
        """Returns the result of the Leroy and Parthiot 1998 depth to pressure 
        formula with corrective term for the Mediterranean Sea
        """
        return ((self.__h_Z45(Z) * self.__k_ZPhi(Z, L) 
                 - (-8.5 * 10**-6 * Z + 1.4 * 10**-9 * Z**2)) * 10) 

    def d2p_Sea_of_Japan_correction(self, Z, L):
        """Returns the result of the Leroy and Parthiot 1998 depth to pressure 
        formula with corrective term for the Sea of Japan alternatively the 
        common ocean corrective term can be used as well
        """
        return ((self.__h_Z45(Z) * self.__k_ZPhi(Z, L) 
                 - (7.8 * 10**-6 * Z)) * 10) 

    def d2p_Sulu_Sea_correction(self, Z, L):
        """Returns the result of the Leroy and Parthiot 1998 depth to pressure 
        formula with corrective term for the Sulu Sea (8° latitude)
        """
        return ((self.__h_Z45(Z) * self.__k_ZPhi(Z, L) 
                 - (1 * 10**-2 * Z / (Z + 100) + 1.6 * 10**-5 * Z + 1 * 10**-9 
                    * Z**2)) * 10) 

    def d2p_Halmahera_Basin_correction(self, Z, L):
        """Returns the result of the Leroy and Parthiot 1998 depth to pressure 
        formula with corrective term for the Hamahera Basin (0° latitude)
        """
        return ((self.__h_Z45(Z) * self.__k_ZPhi(Z, L) 
                 - (8 * 10**-3 * Z / (Z + 50) + 1.3 * 10**-5 * Z)) * 10) 

    def d2p_Celebes_Basin_Weber_Sea_correction(self, Z, L):
        """Returns the result of the Leroy and Parthiot 1998 depth to pressure 
        formula with corrective term for the Celebes Basin (4°) and 
        the Weber Deep (6°)
        """
        return ((self.__h_Z45(Z) * self.__k_ZPhi(Z, L) 
                 - (1.2 * 10**-2 * Z / (Z + 100) + 7 * 10**-6 * Z)) * 10) 
    
    def d2p_Black_Sea_correction(self, Z, L):
        """Returns the result of the Leroy and Parthiot 1998 depth to pressure 
        formula with corrective term for the Black Sea (43°)
        """
        return ((self.__h_Z45(Z) * self.__k_ZPhi(Z, L) 
                 - (1.13 * 10**-4 * Z)) * 10) 
    
    def d2p_Baltic_Sea_correction(self, Z, L):
        """Returns the result of the Leroy and Parthiot 1998 depth to pressure 
        formula with corrective term for the Baltic Sea (60°)
        """
        return ((self.__h_Z45(Z) * self.__k_ZPhi(Z, L) 
                 - (1.8 * 10**-4 * Z)) * 10)

    def __chen_C_term(self, T, P):
        """Returns Chen and Millero 1977 C term"""
        return (self.c00 + self.c01 * T + self.c02 * T**2 + self.c03 * T**3 
                + self.c04 * T**4 + self.c05 * T**5 + (self.c10 + self.c11 * T 
                + self.c12 * T**2 + self.c13 * T**3 + self.c14 * T**4) * P 
                + (self.c20 + self.c21 * T + self.c22 * T**2 + self.c23 * T**3 
                + self.c24 * T**4) * P**2 + (self.c30 + self.c31 * T 
                + self.c32 * T**2) * P**3)

    def __chen_A_term(self, T, P):
        """Returns Chen and Millero 1977 A term"""
        return (self.a00 + self.a01 * T + self.a02 * T**2 + self.a03 * T**3 
                + self.a04 * T**4 + (self.a10 + self.a11 * T + self.a12 * T**2 
                + self.a13 * T**3 + self.a14 * T**4) * P + (self.a20 
                + self.a21 * T + self.a22 * T**2 + self.a23 * T**3) * P**2 
                + (self.a30 + self.a31 * T + self.a32 * T**2) * P**3)

    def __chen_B_term(self, T, P):
        """Returns Chen and Millero 1977 B term"""
        return (self.b00 + self.b01 * T + (self.b10 + self.b11*T) * P)

    def __chen_D_term(self, T, P):
        """Returns Chen and Millero 1977 D term"""
        return (self.d00 + self.d10*P)

    def __chen_millero_terms(self, T, S, P):
        """Returns the complete result of the Chen Millero formula 1977 based on 
        above terms
        """
        return (self.__chen_C_term(T, P) + self.__chen_A_term(T, P) * S 
              + self.__chen_B_term(T, P) * S**(3/2) + self.__chen_D_term(T, P) 
              * S**2)

    def convert_lon_to_360(self, x):
        return x % 360
    
    def mimic_circular_axis(self, ds):
        ds = ds.assign_coords(lon=(self.convert_lon_to_360(ds.lon)))
        ds = ds.roll(lon=int(len(ds['lon']) / 2), roll_coords=True)
        
        return ds 
    
    def crop_and_combine(self):
        '''Opens the temperature and salinity netcdf files from the WOA18 
        and creates xarray datasets to crop them to the latitude and longitude 
        extent. Both datasets are merged in the ds_merged dataset.
        '''
        self.ds_temp = xr.open_dataset(self.temp_file, decode_times=False).t_an
        self.ds_sali = xr.open_dataset(self.sali_file, decode_times=False).s_an
        
        if self.uplon > 180:
            self.ds_merged = (xr.merge([self.ds_temp, self.ds_sali])
                            .isel(time=0)
                            .reset_coords('time')
                            .drop_vars("time"))
            
            self.ds_merged = self.mimic_circular_axis(self.ds_merged)
            self.ds_merged = self.ds_merged.sel(lon=slice(self.lolon, self.uplon),
                                             lat=slice(self.lolat, self.uplat))
        else:
            self.ds_temp = self.ds_temp.sel(lon=slice(self.lolon, self.uplon),
                                             lat=slice(self.lolat, self.uplat))
            self.ds_sali = self.ds_sali.sel(lon=slice(self.lolon, self.uplon),
                                             lat=slice(self.lolat, self.uplat))
        
            self.ds_merged = (xr.merge([self.ds_temp, self.ds_sali])
                            .isel(time=0)
                            .reset_coords('time')
                            .drop_vars("time"))
            
        del self.ds_temp
        del self.ds_sali
        gc.collect()
        
    def calculate_pressure(self):
        '''Calculates the pressure from depth and latitude based on the 
        pressure formula specified and adds it to ds_merged.
        '''
        self.ds_merged["p_an"] = (self._pressure_formula(self.ds_merged.depth,
                                                         self.ds_merged.lat)
                                                         .expand_dims(
                                                         dim={"lon": 
                                                         self.ds_merged.lon},
                                                         axis=-1))
        
        self.ds_merged.p_an.attrs["standard_name"] = 'sea_water_pressure'
        self.ds_merged.p_an.attrs["long_name"] = ('Pressure calculated from '
                                                 'temperature and salinity '
                                                 'with corrective term: '
                                                 '{}'
                                                 .format(self.corrective_term))
        self.ds_merged.p_an.attrs["units"] = 'pressure_MPa'
    
  
    def __leroy_terms(self, t, s, z, lat):
        '''Calculates the sound velocity from temperature, salinity, depth and 
        latitude using the Leroy formula.
        '''
        return (1402.5 + 5 * t - 5.44 * 10**-2 * t**2 + 2.1 * 10**-4 
                * t**3 + 1.33 * s - 1.23 * 10**-2 * s * t + 8.7 * 10**-5
                * s * t**2 + 1.56 * 10**-2 * z + 2.55 * 10**-7 * z**2 
                - 7.3 * 10**-12 * z**3 + 1.2 * 10**-6 * z * (lat - 45)
                - 9.5 * 10**-13 * t * z**3 + 3 * 10**-7 * t**2 * z
                + 1.43 * 10**-5 * s * z)
    

    def calculate_sound_velocity(self):
        '''Calls the Chen_Millero and the Leroy formula using the ds_merged 
        parameters, adds them to the ds_merged and removes temperature, salinity 
        and pressure from the ds_merged.
        '''
        self.ds_merged["v_an_chen"] = self.__chen_millero_terms(self.ds_merged.t_an, 
                                                           self.ds_merged.s_an, 
                                                           self.ds_merged.p_an)
        
        self.ds_merged.v_an_chen.attrs["standard_name"] = 'sea_water_velocity_chen'
        self.ds_merged.v_an_chen.attrs["long_name"] = ('Sound velocity calculated '
                                                 'from temperature, salinity '
                                                 'and pressure (corrective '
                                                 'term {}) using the Chen '
                                                 'Millero formula 1977'
                                                 .format(self.corrective_term))

        self.ds_merged["v_an_leroy"] = self.__leroy_terms(self.ds_merged.t_an, 
                                                           self.ds_merged.s_an, 
                                                           self.ds_merged.depth,
                                                           self.ds_merged.lat)
        
        self.ds_merged.v_an_leroy.attrs["standard_name"] = 'sea_water_velocity_leroy'
        self.ds_merged.v_an_leroy.attrs["long_name"] = ('Sound velocity '
                                                 'calculated from temperature, '
                                                 'salinity, depth and latitude '
                                                 'using the Leroy formula.')
        
        self.ds_merged = self.ds_merged.drop_vars(['t_an', 's_an', 'p_an'])




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
            

class ApplySVP():
    
    def __init__(self, swathfile, png=True):
        
        self.png = png
        self.metrics = {}
        self.swathfile = swathfile
        self.processed_swathfile = "p.".join(self.swathfile.split("."))
        
    def apply_svp(self, svpfile):
        
        subprocess.run(["mbset", "-I {}".format(self.swathfile),
                        "-PSVPFILE:{}".format(svpfile)])
        subprocess.run(["mbprocess", "-I {}".format(self.swathfile)])
        current_ping = subprocess.run(["mblist", "-I {}".format(self.processed_swathfile),
                                       "-ON#XYZ", "-MA", "-G,"],
                                      capture_output=True, text=True).stdout
        
        
        
        df_current_ping = pd.read_csv(StringIO(current_ping),
                    names=["Ping number", "Beam number", "Longitude",
                           "Latitude", "Depth"])
        
        df_current_ping = self.adjust_pings(df_current_ping)

        x, y, ystd, yline = self.average_residual(df_current_ping)
        
        rmse = self.calculate_rmse(y, yline)
        mse = self.calculate_mse(y, yline)
        me = self.calculate_me(y, yline)
        self.metrics[svpfile.name] = rmse

        
        if self.png:
            pngfolder = Path.cwd() / 'svppngs'
            if not pngfolder.is_dir():
                pngfolder.mkdir()
            self.plot_metrics_single(svpfile, x, y, ystd, yline, rmse, mse, me, pngfolder)
            
    def svp_statistics(self, directory):
        #print RMSE info
        metrics_sorted = {k: v for k, v in sorted(self.metrics.items(), key=lambda item: item[1])}
        best_svp_file = Path(directory) / (list(metrics_sorted.keys())[0])
        
        #apply the best svp at the end (again)
        print("\nApplying the best SVP: {}".format(best_svp_file.name))
        subprocess.run(["mbset", "-I {}".format(self.swathfile), "-PSVPFILE:{}".format(best_svp_file)])
        subprocess.run(["mbprocess", "-I" + self.swathfile])
        
        print("\n\nSVP ranking:\n")
        for key, value in metrics_sorted.items():
            print("Rmse for {}: {}".format(key, value))
        print("WARNING: All parameter files of the specified swathfiles/datalist"
              " are set to use {} for the PSVPFILE parameter! Use mbset to set"
              " no or other SVP file.".format(best_svp_file.name))
        
    #define metrics functions
    def calculate_rmse(self, reference, data):
        """Calculates the root mean square error of data against zero line."""
        return np.sqrt(np.square(np.subtract(reference, data)).mean())

    def calculate_mse(self, reference, data):
        """Calculates the mean square error of data against zero line."""
        return np.abs(np.subtract(reference, data)).mean()

    def calculate_me(self, reference, data):
        """Calculates the mean error of data against zero line."""
        return np.subtract(reference, data).mean()
    
    def average_residual(self, df):

        df = df.sort_values(['Beam number'], ascending=[True])
        
        #get the ybeammean and ybeamstd with group by beam number
        ybeammean = np.array(df.groupby(["Beam number"])["Yfits"].mean())
        ybeamstd = np.array(df.groupby(["Beam number"])["Yfits"].std())
        
        x= df["Beam number"].unique()
        yline = np.zeros(len(x))
        
        return x, ybeammean, ybeamstd, yline
    
    def adjust_pings(self, df):
        
        number_of_pings = len(df["Ping number"].unique())+1
        
        for i, j in zip(np.arange(1, number_of_pings), df["Ping number"].unique()):
            df.loc[df["Ping number"] == j, "Ping number"] = i
            
        df.set_index("Ping number", inplace=True)
            
        difftotalunshaped = []
        for i in range(1, number_of_pings):
            x= np.array(df.loc[i, "Beam number"])
            y= np.array(df.loc[i, "Depth"])
            
            #get the yfit values
            b, m = polyfit(x, y, 1)
            yfit = b + m * x
            
            #get the residuals of y minus yfit
            diffy = yfit - y
            difftotalunshaped.append(diffy)
        
        difftotal = []
        for i in difftotalunshaped:
            for j in i:
                difftotal.append(j)
            
        df["Yfits"] = difftotal
        
        return df 
    
    def plot_metrics_single(self, filename, x, y, ystd, yline, rmse, mse, me, pngfolder):
        
        fig = plt.figure(figsize=(12,9))
        
        ax = fig.add_subplot()
        ax.set_title("Used swathfile: " + self.swathfile + '\n' + filename.name,
                     pad=12, fontsize=15)
                
        ax.errorbar(x, y, yerr=ystd, ecolor="grey", color='black')
                    
        ax.hlines(0, 0, x.max(), color="k")
        ax.set_xlim([0,x.max()])
        ax.set_xticks([0, int(x.max()*0.25), int(x.max()/2), int(x.max()*0.75), x.max()])
        ax.set_xlabel("Beam number", fontsize=12)
        ax.set_ylim([-15,15])
        ax.set_ylabel("Difference [$m$]", fontsize=12)
        ax.invert_yaxis()
        ax.text(0.5,0.02, "ME: " + str(round(me,3)) + "\n MSE: " + str(round(mse,3)) + "\n RMSE: " + 
                str(round(rmse,3)), ha="center", va="bottom", transform=ax.transAxes)
        
        fig.tight_layout()    
        fig.savefig(pngfolder / "{}.png".format(filename.name), bbox_inches="tight", dpi= 300)
        plt.close(fig)