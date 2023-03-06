#!/usr/bin/env python3
# coding: utf-8

import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from io import StringIO
from numpy.polynomial.polynomial import polyfit
import matplotlib.pyplot as plt

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