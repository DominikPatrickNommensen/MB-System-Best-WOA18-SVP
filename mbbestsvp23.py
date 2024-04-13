#!/usr/bin/env python3
# coding: utf-8

import argparse
from SwathfileInfoClass23 import SwathFileInfo 
from ApplySvpClass import ApplySVP
from pathlib import Path
import datetime
import subprocess

program_name = 'mbbestsvp'
mb_version = 'VERSION'

if __name__=='__main__':

    date = datetime.datetime.now()
    #create a parser which takes some arguments as specified below.
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=('''Calculates the World Ocean 
                                                  Atlas 2023 SVPs for a specified 
                                                  flat part of a multibeam survey 
                                                  based on the sound velocity 
                                                  netcdf files created by "mbdownloadwoa18", 
                                                  applies them on the swathfile and 
                                                  evaluates its capabilities of 
                                                  flattening the seafloor. The 
                                                  WOA18 and SVPs stored in the 
                                                  swathfile are ranked and the 
                                                  best SVP is applied to the swathfile.
                                                  The best SVP could than be used 
                                                  for a sound speed correction of 
                                                  the whole survey. In addition, 
                                                  a folder with pngs showing the 
                                                  flattening of the seafloor by each 
                                                  SVP is also created.'''))
      
    parser.add_argument('-N', '--nearest',  required=False, action='store_true',
                        help='''The user can use this flag to use the nearest 
                        grid point of the center of the provided swathfile area 
                        to extract SVPs. If the user does not use this flag 
                        the surrounding 4 grid points are used to extract mean 
                        SVPs that bound the area.''')
    
    parser.add_argument('-P', '--period',type=str, required=False,
                        metavar="Averaged period", nargs='+', 
                        choices=["decav", "decav91C0", "decav81B0", "decav71A0", "B5C2", "A5B4", "95A4", "8594", "7584", "6574", "5564"],
                        default=["decav", "B5C2"],
                        help='''The user can choose which period(s) to use
                        If not using this argument decav (averaged decades)
                        and B5C2 (2015 to 2022) will be 
                        used.''')
    
    parser.add_argument('-R', '--resolution', type=str, required=False, 
                        metavar="Grid resolution", nargs='+',
                        choices=["01", "04"],
                        default=["01", "04"],
                        help='''The user can choose between the coarser 1° [01] 
                        or the finer 0.25° [04] grid resolution. If not 
                        using this argument both grid resolutions will be 
                        used.''')

    parser.add_argument('-I', '--swathfile', type=str, required=False,
                        default="datalist.mb-1", metavar="Swathfile/Datalist",
                        help='''The user can specify a swathfile or a datalist 
                        from a realitvely flat area of a survey to evaluate 
                        the WOA23 SVPS based on their capabilities to "flatten" 
                        the seafloor.''')
    
    parser.add_argument('-S', '--svpfolder', type=str, required=True,
                        metavar="SVP folder",
                        help='''The user needs to specify the path where the 
                        calculated (and cropped) SVP netcdf files are stored, 
                        that were created using "mbdownloadwoa23".''')
    
    parser.add_argument('-O', '--outputfolder', type=str, required=False,
                        default='./svpprofiles', metavar="SVP profile folder",
                        help='''The user can specify the name of a new or existing 
                        folder where all the individual SVP profiles will be 
                        stored in. If not a new "svpprofiles" folder will be 
                        created in the current working directory.''')
                        
    parser.add_argument('-E', '--enddate', action='store_true', 
                        required=False,
                        help='''If using this flag the end date of the swathfile 
                        will be used instead of the start date. This is only 
                        important if your transect happens to start at the end 
                        of one month into another month so another month 
                        might be equally interesting to check.''')
    
    args = parser.parse_args()
    sf = SwathFileInfo(period=args.period, resolution=args.resolution, 
                       swathfile=args.swathfile, svpfolder=args.svpfolder,
                       nearest=args.nearest)
    sf.mbinfo(end=args.enddate)
    files = sf.get_filenames()
    
    swathfile_lolon = sf.lonmin_swathfile
    swathfile_uplon = sf.lonmax_swathfile
    swathfile_lolat = sf.latmin_swathfile
    swathfile_uplat = sf.latmax_swathfile
    
    profile_folder = Path(args.outputfolder)
    if not profile_folder.is_dir():
        profile_folder.mkdir()
    
    for file in files:
        df = sf.calculate_mean_profiles(file)
        if args.nearest:
            method = '## WOA23 Longitude/Latitude method used: Nearest\n'
            extent = ('## WOA23 SVP longitude point: {}\n'.format(sf.lonmin) +
                      '## WOA23 SVP Latitude point: {}\n'.format(sf.latmin))
        else:
            method = '## WOA23 Longitude/Latitude method used: Mean area\n' 
            extent = ('## WOA23 SVP min/max longitudes: {}/{}\n'.format(sf.lonmin, sf.lonmax) +
                      '## WOA23 SVP min/max Latitudes: {}/{}\n'.format(sf.latmin, sf.latmax))
            
        leroy_file = profile_folder / "{}_leroy.svp".format(file.stem)
        with leroy_file.open(mode='w') as f:
            f.write("## World Ocean Atlas 2023 (WOA23) SVP\n" +
                    "## WOA23 grid resolution: {}\n".format(file.stem[-2:]) +
                    "## WOA23 averaged decades: {}\n".format(file.stem[6:10]) +
                    "## WOA23 time period: {}\n".format(file.stem[-5:-3]) +
                    "## Water sound velocity formula used: Leroy\n" +
                    "## Output by Program {}\n".format(program_name) +
                    "## Swath File: {}\n".format(args.swathfile) +
                    "## Swath File min/max longitude and latitude: {} {} {} {}\n".format(swathfile_lolon,
                                                                                         swathfile_uplon,
                                                                                         swathfile_lolat,
                                                                                         swathfile_uplat) +
                    method +
                    extent +
                    "## Start time: {}\n".format(date.strftime("%Y/%m/%d %H:%M:%S.%f")) +
                    "## Water Sound Velocity Profile (SVP)\n" +
                    "## SVP Count: {}\n".format(len(df))) 
            f.write(df["v_an_leroy"].to_string(header=False, index=True))
            
        chen_file = profile_folder / "{}_chen.svp".format(file.stem)
        with chen_file.open(mode='w') as f:
            f.write("## World Ocean Atlas 2023 (WOA23) SVP\n" +
                    "## WOA23 grid resolution: {}\n".format(file.stem[-2:]) +
                    "## WOA23 averaged decades: {}\n".format(file.stem[6:10]) +
                    "## WOA23 time period: {}\n".format(file.stem[-5:-3]) +
                    "## Water sound velocity formula used: Chen\n" +
                    "## Output by Program {}\n".format(program_name) +
                    "## Swath File: {}\n".format(args.swathfile) +
                    "## Swath File min/max longitude and latitude: {} {} {} {}\n".format(swathfile_lolon,
                                                                                         swathfile_uplon,
                                                                                         swathfile_lolat,
                                                                                         swathfile_uplat) +
                    method +
                    extent +
                    "## Start time: {}\n".format(date.strftime("%Y/%m/%d %H:%M:%S.%f")) +
                    "## Water Sound Velocity Profile (SVP)\n" +
                    "## SVP Count: {}\n".format(len(df)))
            f.write(df["v_an_chen"].to_string(header=False, index=True))
    
    # sf.extract_svp()
    # internal_svp = subprocess.run(["mbsvplist", "-I {}".format(args.swathfile), "-V"],
    #                           capture_output=True, text=True).stdout


    # internal_svp_list = internal_svp.split('## MB-SVP')
    # internal_svp_file = profile_folder / (Path(internal_svp.split('## Swath File: ')[-1].split('\n')[0]).stem + '.svp')
    
    # internal_file = profile_folder / 
    # with Path(internal_svp_file).open(mode='w') as f:
    #     f.write(internal_svp)
    # for internal in Path.cwd().glob('*svp'):
    #     internal.rename(profile_folder / internal.name)
    
    subprocess.run(["mbsvplist", "-I {}".format(args.swathfile), "-O", "-V"])
    for internal in Path(args.swathfile).parents[0].glob('*svp'):
        internal.rename(profile_folder / internal.name)
    
    #apply svp part
    asvp = ApplySVP(args.swathfile)
    for svp in profile_folder.glob('*.svp'):
        # print(svp)
        asvp.apply_svp(svp)
    asvp.svp_statistics(profile_folder)
        
        
        
        
        