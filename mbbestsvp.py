#!/usr/bin/env python3
# coding: utf-8

import argparse
from CalculateSvpClass import SwathFileInfo, ApplySVP
from pathlib import Path
import datetime

program_name = 'mbbestsvp'

mb_version = 'VERSION'

if __name__=='__main__':

    date = datetime.datetime.now()
    #create a parser which takes some arguments as specified below.
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=(''))
      
    parser.add_argument('-N', '--nearest',  required=False, action='store_true',
                        help='')
    
    parser.add_argument('-P', '--period', type=str, required=False,
                        default="both", metavar="",
                        help='')
    
    parser.add_argument('-R', '--resolution', type=str, required=False,
                        default="both", metavar="", 
                        help='')

    parser.add_argument('-I', '--swathfile', type=str, required=False,
                        default="datalist.mb-1", metavar="",
                        help='')
    
    parser.add_argument('-S', '--svpfolder', type=str, required=True,
                        metavar="SVP folder",
                        help='')
    
    parser.add_argument('-O', '--outputfolder', type=str, required=False,
                        default='./svpprofiles', metavar="SVP profile folder",
                        help='')
    
    args = parser.parse_args()
    sf = SwathFileInfo(period=args.period, resolution=args.resolution, 
                       swathfile=args.swathfile, svpfolder=args.svpfolder,
                       nearest=args.nearest)
    sf.mbinfo()
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
            method = '## WOA18 Longitude/Latitude method used: Nearest\n'
            extent = ('## WOA18 SVP longitude point: {}\n'.format(sf.lonmin) +
                      '## WOA18 SVP Latitude point: {}\n'.format(sf.latmin))
        else:
            method = '## WOA18 Longitude/Latitude method used: Mean area\n' 
            extent = ('## WOA18 SVP min/max longitudes: {}/{}\n'.format(sf.lonmin, sf.lonmax) +
                      '## WOA18 SVP min/max Latitudes: {}/{}\n'.format(sf.latmin, sf.latmax))
            
        leroy_file = profile_folder / "{}_leroy.svp".format(file.stem)
        with leroy_file.open(mode='w') as f:
            f.write("## World Ocean Atlas 2018 (WOA18) SVP\n" +
                    "## WOA18 grid resolution: {}\n".format(file.stem[-2:]) +
                    "## WOA18 averaged decades: {}\n".format(file.stem[6:10]) +
                    "## WOA18 time period: {}\n".format(file.stem[-5:-3]) +
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
            f.write("## World Ocean Atlas 2018 (WOA18) SVP\n" +
                    "## WOA18 grid resolution: {}\n".format(file.stem[-2:]) +
                    "## WOA18 averaged decades: {}\n".format(file.stem[6:10]) +
                    "## WOA18 time period: {}\n".format(file.stem[-5:-3]) +
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
    
    sf.extract_svp()
    for internal in Path.cwd().glob('*svp'):
        internal.rename(profile_folder / internal.name)
        
    #apply svp part
    asvp = ApplySVP(args.swathfile)
    for svp in profile_folder.glob('*.svp'):
        asvp.apply_svp(svp)
    asvp.svp_statistics(args.svpfolder)
        
        
        
        
        