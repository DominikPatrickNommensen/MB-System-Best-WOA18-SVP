#!/usr/bin/env python3
# coding: utf-8

import xarray as xr
import numpy as np
import gc

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