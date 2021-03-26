# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 15:50:46 2019

@author: Thomas
"""

import numpy as np
import matplotlib.pyplot as plt
import re
import scipy.signal as ss
import scipy.optimize as so
import math
import os.path


def lorentz(x, a, s, m, offset):
    return a/(1+(x-m)**2/s**2)+offset


def poly_lorentz(x, *args):
    N = int ( np.size(args)/4 )

    poly = np.zeros((x.size, N))
    
    for i in range(N):
        poly[:,i] = lorentz(x, *args[ 4*i:4*i+4 ] )
    
    return poly.sum(axis=1)

def get_no_peaks(signal, n):
    all_peaks = ss.find_peaks(signal)[0]
    proms = ss.peak_prominences(signal, all_peaks)[0]
    
    # Taking the n most prominent peaks
    index = np.argsort(proms)
    all_peaks = all_peaks[index]
    
    return all_peaks[-n:]

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def magnitude(x):
    return int(math.floor(math.log10(x)))

def load_log(filename):
    with open(filename) as f:
        s = f.read()
    return s.splitlines()

class Experiment(object):
    def __init__(self, filename):
        with open(filename) as f:
            data = f.read()
        
        # Assigning each line to an array element
        data = data.splitlines()
        
        self.title = re.split("\s{2,}", data[0])[0]
        
        print("Loaded experiment: ", self.title)
        
        # Have to convert it back to string and split by whitespace
        info_str = ''.join(data[1])
        info = np.asarray( info_str.split(), dtype=float)
        
        # Loading calibrations from file
        self.channels = info[0]
        self.fold_ch = 511
        
        self.cal = 1
        self.zero_ch = 0
        
        counts_str = ''.join(data[2:])
        self.counts = np.asarray( counts_str.split(), dtype=float )
    
#        self.relative()
        self.fold()
        
        
        # Loading log file if it exists
        logname = filename[:-3]+'log'
        if os.path.isfile(logname):
            info = load_log(logname)
            
            # Time elapsed in hours
            time_str = info[2]
            self.time = time_str.split(':')[-1]
            
            # Mean Temperature
            temp_str = info[3]
            self.temp = temp_str.split(':')[-1] 
            
            # Min and max temperature
            max_str = info[4]
            min_str = info[5]
            
            max_temp = float( max_str.split(':')[-1] )
            min_temp = float( min_str.split(':')[-1] )
            print("The difference in maximum and minimum temperature is %.3f K" %(max_temp-min_temp))
            
        else:
            self.time = "No log file found"
            self.temp = "No log file found"
    
#    def relative(self):
#        baseline = (self.counts[:5].mean()+self.counts[-5:].mean() ) /2
#        self.counts_rel = self.counts / baseline
    
    def fold(self):               
        remain = self.channels-self.fold_ch
        
        ch_skip = 5
        if np.abs(remain) > ch_skip:
            start_ch = remain
        else:
            start_ch = ch_skip
        
        N = np.array( (self.fold_ch - start_ch)/2, dtype=int)
        
        self.counts_fold = np.zeros(N)
        for i in range(N):
            ch = start_ch + i
            
            neg = int(self.fold_ch-ch-2)
            count = self.counts[ch] + self.counts[neg]
            self.counts_fold[i] = count
        
        add = np.sum( self.counts_fold[:5] ) + np.sum( self.counts_fold[-5:])
        baseline = add / (2*5)

        
        self.counts_rel = self.counts_fold / baseline

        self.velocity = np.zeros(N)
        for i in range(N):
            v = (i+1+5-self.zero_ch)*self.cal
            self.velocity[i] = v
    
    def plot(self, N):
        popt = self.fit_data(N, return_popt=True)
        x = np.linspace(self.velocity[0], self.velocity[-1], 1000)
        y = poly_lorentz(x, *popt)
        
        plt.figure(figsize=(12,9))

        plt.plot(self.velocity, self.counts_rel, 'ro-', lw=1, ms=2)
        plt.plot(x,y, color='blue', lw=1)
        plt.ylabel("Relative number of counts")

        plt.title(self.title)
        
        plt.xlabel("Velocity [mm/s]")
        plt.tick_params('both', direction='inout', top=True, right=True)
                
        
        plt.grid(axis='y')
        plt.show()


    def fit_data(self, N, return_pos = False, return_popt = False):

        p = np.sort( get_no_peaks(1/self.counts_rel, N) )

            
        p0 = np.zeros(4*N)
        for i in range(N):
            p0[4*i] =  -self.counts_rel[ p[i] ]/2
            p0[4*i+1] = .1
            p0[4*i+2] = self.velocity[ p[i] ]   
            p0[4*i+3] = 1
        
        popt, pcov = so.curve_fit(poly_lorentz, self.velocity, self.counts_rel.ravel(), p0)
        
               
        if return_pos:
            pos = np.zeros(N)
            for i in range(N):
                pos[i] = popt[4*i+2]
            return pos
        
        elif return_popt:
            return popt
        
        
        # Extracting widths from popt
        widths = np.zeros(N)
        for i in range(N):
            widths[i] = popt[4*i+1]*2        

        return widths
    
    def optimal_fold(self, N):
        
        fold_ch = [510, 511, 512]
        tmp = np.zeros( np.size(fold_ch) )
        
        for i, f in enumerate(fold_ch):
            self.fold_ch=f
            self.fold()
    
            widths = self.fit_data(N)
            tmp[i] = widths.mean()
#            print(tmp[i])
        
        indx = np.argmin(tmp)
        self.fold_ch = fold_ch[indx]  
        print("Most optimal fold channel is: %i" %self.fold_ch) 
    
    def calibration(self, N, weight):
        self.fold()
        
        if N%2 != 0:
            raise Exception("There has to be an equal number of peaks!")
        P = int(N/2)
        if P > 3:
            raise Exception("There cannot be more than 3 peaks for Fe-57")
        

        
        line_pos = self.fit_data(N, return_pos=True)
        print("Peaks are positioned at channel:")        
        print(line_pos)
        

        
        direct = np.zeros(P)
        zero_direct = np.zeros(P)
        # From Gutlich
        values = np.array([10.657, 6.167, 1.677])
        
        for i in range(P):
            direct[i] = values[i]/(line_pos[N-i-1] - line_pos[i])
            zero_direct[i] = (line_pos[N-i-1] + line_pos[i])/2
            
        self.cal = np.dot(direct, weight[:P])
        self.zero_ch = np.dot(zero_direct, weight[:P])
        
#        if N == 6:            
#            #From Cathrine
##            values = np.array([10.648, 6.166, 1.68])
##            weight = np.array([0.5, 0.3, 0.2])
#            
#            for i in range(direct.size):
#                direct[i] = values[i]/(line_pos[N-i-1] - line_pos[i])
#                zero_direct[i] = (line_pos[N-i-1] + line_pos[i])/2
#            
#            self.cal = np.dot(direct, weight)
#            self.zero_ch = np.dot(zero_direct, weight)
#            
#        elif N == 4:
#            #From Cathrine
##            values = np.array([6.152, 1.6794])
##            weight = np.array([0.6, 0.4])
#            
#            for i in range(direct.size):
#                direct[i] = values[i]/(line_pos[N-i-1] - line_pos[i])
#                zero_direct[i] = (line_pos[N-i-1] + line_pos[i])/2
#                
#            self.cal = np.dot(direct, weight[:P])
#            self.zero_ch = np.dot(zero_direct, weight[:P])
#    
#        else:
#            raise Exception("Calibration values have not been added for %i peaks" %N)
        
        print("Calibration constant: %.4f" %self.cal)
        print("Zero channel:         %.4f" %self.zero_ch)       
