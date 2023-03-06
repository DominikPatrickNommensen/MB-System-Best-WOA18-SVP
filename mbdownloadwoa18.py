#!/usr/bin/env python3
# coding: utf-8

import argparse
from OpendapClass import AccessOpendap18
from CalculateSvpClass import CalculateSvp
from pathlib import Path
import re
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# import tracemalloc
import gc
import time


def progress_bar(progress, total):
    percent = 100 * (progress / float(total))
    bar ='█' * int(percent) + '-' * (100 - int(percent))
    print(f"\r|{bar}| {percent:.2f}%", end="\r")


def overview_map(area, folder):
    if area[1] > 180:
        fig, ax = plt.subplots(figsize=(12,9),
                               subplot_kw=dict(projection=ccrs.PlateCarree(
                                   central_longitude=180)))
        ax.set_extent([0, 359.9, -90, 90])
    else:
        fig, ax = plt.subplots(figsize=(12,9),
                               subplot_kw=dict(projection=ccrs.PlateCarree(
                                   central_longitude=0)))
        ax.set_extent([-180, 180, -90, 90])
        
    width = area[1] - area[0]
    if width == 360:
        width = 359.9
    height = area[3] - area[2]
    ax.add_feature(cfeature.LAND, zorder=10, color="grey")
    ax.coastlines(zorder=11)
    ax.add_patch(mpatches.Rectangle(xy=[area[0], area[2]],
                                    width=width, height=height,
                                        facecolor='blue',
                                        alpha=0.2,
                                        transform=ccrs.PlateCarree()))
    
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    
    fig.savefig(folder / 'cropped_area.png', bbox_inches='tight')

def allowed_area(area):
    '''Validate the area specified by the user to gurantee plausible longitudinal 
    and latitudinal extent.
    
    Parameters
    ----------
    area : list of 4 string elements or string
        List of 4 string elements list or string

    Raises
    ------
    argparse
        Error message if passed argument is not valid.

    Returns
    -------
    list of floats
        After successful validation returns either a list of 4 float elements 
        representing the latitudinal and longitudinal extent of the desired area
        as passed to the function or the whole globe (cut through the pacific). 
    '''
    
    msg = ("%r is not a valid argument. Please use 'globe' "
           "or 'min lon, max lon, min lat, max lat' with -180 <= lon <= 180 "
           ", -90 <= lat <= 90 and min lon/lat < max lon/lat. "
           "If you provide a positive lon min and a negative lon max then it's "
           "assumed you want to extract an area around the -180 to 180 zone, "
           "the dataset will be then cut at 0 degree."% area)
    
    if area == 'globe':
        return [-180, 180, -90, 90]
    
    elif (len(area) == 4):
        if ((float(area[0]) >= -180) and 
            (float(area[0]) <=  180) and 
            (float(area[1]) <=  180) and 
            (float(area[1]) >= -180) and
            (float(area[2]) >=  -90) and 
            (float(area[3]) <=   90) and
            (float(area[2]) < float(area[3]))):
            
            print('Valid longitude and latitude provided.')
            if ((float(area[0]) >= float(area[1])) and
                (float(area[0]) >= 0) and 
                (float(area[1]) <= 0)):
                
                print('Longitude were passed in decreasing order. This will '
                      'shift the dataset by 180 degree allowing to consisitently '
                      'slice areas extending around 180 degree e.g. Pacific.')
                
                area[1] = str(float(area[1]) % 360)
                if float(area[1]) == 0:
                    area[1] = 360
                
                return [float(x) for x in area]
            elif (float(area[0]) < float(area[1])):
                print('Longitude passed in increasing order. No shifting of the '
                      'dataset will be done (cut through pacific).')
                return [float(x) for x in area]
            else:
                raise argparse.ArgumentTypeError(msg)
        else:
            raise argparse.ArgumentTypeError(msg)
    else:
        raise argparse.ArgumentTypeError(msg)


if __name__=='__main__':
    start_time = time.time()
    parser = argparse.ArgumentParser(formatter_class=
                                     argparse.ArgumentDefaultsHelpFormatter,
                                     description=('Access the World Ocean Atlas'
                                     '2018 Temperature and Salinity netcdf '
                                     'files, calculate the sound velocity '
                                     'profiles and store them.'))
      
    parser.add_argument('-A', '--area', required=False,
                        default='globe', nargs=4, 
                        metavar=('min_lon', 'max_lon', 'min_lat', 'max_lat'),
                        help='Specify which longitude and latitude extent '
                        'should be used from the netcdf files.\n The user can specify '
                        'minima and maxima longitudes and latitudes [-180 to 180] '
                        'like: [min lon, max lon, min lat, max lat]. This '
                        'option might be useful to extract a subset e.g. a '
                        'basin used for multiple missions. If specifying a ' 
                        'min lon > max lon it is assumed the user wants to '
                        'download the data with no cut through the pacific e.g. '
                        '160 -120 would result in cropping the dataset to only '
                        'get the area around 160 to 240 (360 degree notation). '
                        'If the user wants to get the entire world but with a '
                        'cut at 0 degree instead through the Pacific they should '
                        'use 0 for both min and max longitude.\n NOTE: Larger '
                        'areas will take longer times. Also not cutting through '
                        'the Pacific takes longer due to "rolling" of the dataset.\n '
                        'Omitting this argument does not crop the original extent of the '
                        'WOA temperature and salinity files and use the '
                        'whole globe with a cut through the Pacific.\n')

    parser.add_argument('-O', '--outputfolder', type=str, required=True,
                        metavar="SVP folderpath",
                        help='Output folder path to where the calculated '
                        'sound velocity netcdf files should be stored.')
    
    parser.add_argument('-P', '--period', type=str, required=False,
                        metavar="Averaged period", nargs='+', 
                        choices=["decav", "A5B7"],
                        help='The user can choose between averaged decades '
                        '[decav] and the 2005 to 2017 average [A5B7] periods '
                        'If not using this argument both periods will be '
                        'used.')
    
    parser.add_argument('-R', '--resolution', type=str, required=False, 
                        metavar="Grid resolution", nargs='+',
                        choices=["01", "04"],
                        help='The user can choose between the coarser 1° [01] '
                        'or the finer 0.25° [04] grid resolution. If not '
                        'using this argument both grid resolutions will be '
                        'used.')
    
    parser.add_argument('-T', '--time', type=str, required=False, 
                        metavar="Time period", nargs='+',
                        choices=["00", "01", "02", "03", "04", "05", "06", 
                                 "07", "08", "09", "10", "11", "12", "13", 
                                 "14", "15", "16"],
                        help='The user can specify which months '
                        '[01-12], which seasons [winter=13, spring=14, '
                        'summer=15, autumn=16] and/or the annual [00] should '
                        'be used. If not using this argument all times will '
                        'be used.')
    
    args = parser.parse_args()
    args.area = allowed_area(args.area)
    ao = AccessOpendap18(period=args.period, resolution=args.resolution,
                          time=args.time)

    progress_bar(0, len(ao.url_list))
    for i, pair in enumerate(ao.url_list):
        svp = CalculateSvp(args.area, pair)
        svp.crop_and_combine()
        svp.calculate_pressure()
        svp.calculate_sound_velocity()
        
        current_file = Path(pair[0])
        svp_filename = "_svp".join(re.split('_s|_t', current_file.name))
        outputfolder = Path(args.outputfolder) / svp_filename
        svp.ds_merged.to_netcdf(outputfolder)
        del svp
        gc.collect()
        progress_bar(i + 1, len(ao.url_list))
        
    overview_map(args.area, Path(args.outputfolder))
    print(time.time() - start_time)