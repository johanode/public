# -*- coding: utf-8 -*-
"""
Created on Thu May 25 11:39:36 2023

@author: johode
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import scipy

# Define path to measurement folder - MODIFY 
pread = r'C:\Data'

# Select measurement folders
folders = ['01NOV22B1','24MAY23B', '01NOV22A2','24MAY23A']

#%% Subfunctions
def read_calsheet(folder):
    # Read projectname and date from matlab file calsheet.m
    filename = os.path.join(pread, folder, 'calsheet.m')
    
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    # Extract project name and project date meting
    projectname = lines[0].split('=')[1].split(';')[0].replace("'","").strip()
    projectdate = lines[1].split('=')[1].split(';')[0].replace("'","").strip()
    
    data = {
            'projectname' : projectname,
            'projectdate' : projectdate
            }

    return data

def load_frf(folder, idx, freq=None, freqlim=[100, 5000], oct3=True):
    # Load data from computed FRF stored in matlab .mat files
    
    # Read freqeuncies
    leading_zeros = 4
    if freq is None:
        meas_number = str(1).zfill(leading_zeros)
        data = scipy.io.loadmat(os.path.join(pread, folder, 'tfrSVAN'+meas_number+'.mat'))
        if oct3:
            freq = data['okt3']
        else:
            freq = data['G']
    freq = np.squeeze(freq)        
    
    # Define frequency range
    if freqlim is not None:
        i_frange = (freq >= freqlim[0]) & (freq <= freqlim[1])
    else:
        i_frange = np.isfinite(freq)
    
    freq = freq[i_frange].astype(int)
    Nfreq = len(freq)
    
    # Load tf data
    TF = {}        
    for direction in idx:  
        Npos = len(idx[direction])
        TF[direction] = np.zeros((Nfreq, Npos))
        
        # Adapt to Python indexing, i.e. channel 1 is index 0
        ch = channel[direction]-1
        
        for k in range(Npos):    
            #data = scipy.io.loadmat(os.path.join(pread, d[idx[direction][k]]))
            meas_number = str(idx[direction][k]).zfill(leading_zeros)
            data = scipy.io.loadmat(os.path.join(pread, folder, 'tfrSVAN'+meas_number+'.mat'))
            
            if oct3:
                TF[direction][:, k] = data['tfrokt3'][i_frange, ch]
            else:
                TF[direction][:,k] = data['tfr'][i_frange,ch]    
    
    return TF, freq


def load_results(folder):
    # Load initial measurement results computed by matlab script
    T1 = pd.read_csv(os.path.join(pread, folder, 'TDRprocessed1.txt'), header=None, delimiter=',')    
    T2 = pd.read_csv(os.path.join(pread, folder, 'TDRprocessed2.txt'), header=None, delimiter=',')
    
    TDR = {
        'horizontal' : T1.iloc[:,1].values,
        'vertical' : T2.iloc[:,1].values
        }

    return TDR    
    
def load_limits():
    # Limits according to IS0 3095
    limits = {
        'freq' : [250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000],
        'vertical' : [2, 2, 6, 6, 6, 2.19, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8],
        'horizontal' : [2.04, 1.38, 0.94, 0.64, 0.43, 0.29, 0.2, 0.2, 0.32, 0.5, 0.5, 0.5, 0.5, 0.5]
        }
    return limits

def compute_tdr(A, dx):
    # Compute decay rates according to EN 15461
    DR = np.zeros(A.shape[0])
        
    for i in range(A.shape[0]):
        A0 = np.abs(A[i, 0]) ** 2
        An = np.abs(A[i,:]) ** 2
        DR[i] = 4.343 / np.sum(An / A0 * dx)
            
    return DR


def rms(y):
    return np.sqrt(np.mean(y**2))

#%% Main
# Track design
sleeper_distance = 0.65

# Measurement positions in accordance with EN 15461 
sleeper_pos = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5,
              3, 3.5, 4, 5, 6, 7, 8, 10, 12, 16, 20, 24, 30, 36, 42, 48, 54, 66])

# Measurement channel setup
channel = {
     'horizontal' : 2,
     'vertical' : 3
     }

# Loop through folders
TDR = {}
FRF = {}
for folder in folders: 
    print(folder)
    
    if folder.startswith('01NOV22B'):
        # Correction position B 1/11 2022
        sleeper_pos[18:] = sleeper_pos[18:] - 1
        print('change sleepers')
   
    # Define measurement order
    if folder.startswith('01NOV22'):
        idx = {
             'horizontal' : np.arange(1, 58+1, 2),
             'vertical' : np.arange(2, 58+1, 2)
             }
    else:
        idx = {
             'horizontal' : np.arange(2, 58+1, 2),
             'vertical' : np.arange(1, 58+1, 2)
             }
    
    # Define positions in m
    x = sleeper_pos * sleeper_distance
    dx = np.concatenate(([0], x[1:] - x[:-1]))
    Npos = len(sleeper_pos)
     
    # Load 1/3 octave frfs with center frequencies
    TFoct3, fc = load_frf(folder, idx, oct3=True, freqlim=[100, 5000])
    FRF[folder] = TFoct3
    Nfc = len(fc)
    
    # Compute decay rates for both directions
    DR = {direction : compute_tdr(TFoct3[direction], dx) for direction in TFoct3}
    TDR[folder] = DR
    
    # Load S&M matlab results
    DR['mat'] = load_results(folder)
    TDR[folder]['mat'] = DR['mat']
     
    
#%% Plot results
# Plotting results, will create a new figure for each measurement

# Uncomment and mofidy if you want to change folders to plot 
# folders = ['01NOV22B1','24MAY23B', '01NOV22A2','24MAY23A']

# Get line colors for every position using colormap jet
colors = plt.cm.jet(np.linspace(0, 1, Npos))[::-1]

# Get TDR limits
limits = load_limits()    


for folder in folders:    
    fig, ax = plt.subplots(3, 2)
    
    # Read name and create figure title
    cal = read_calsheet(folder)
    fig.suptitle(f"{cal['projectname']} ({cal['projectdate']})")
    
    # Load decay rates and transfer functions (FRFs)
    TFoct3 = FRF[folder]
    DR = TDR[folder]
    
    # Loop both directions
    for j, direction in enumerate(['horizontal','vertical']):
        ax1 = ax[0,j]
        ax2 = ax[1,j]
        ax3 = ax[2,j]
        
        # Plot transfer functions (FRF)
        for k in range(Npos):
            ax1.semilogx(fc, 10*np.log10((TFoct3[direction][:,k]*1e9)**2), color=colors[k])
        
        # Set axex properteis
        ax1.set_ylim([85,200])
        ax1.set_xlim([100, 6000])
        ax1.grid(True, which='both', linestyle='dotted')
        ax1.set_title('Frequency response function')
        ax1.set_xlabel(r'Frequency $(Hz)$')
        ax1.set_ylabel(r'$dB~rel~1e^{-9}~m/s^2/N$')
    
        
        # Spectrogram (contour) plot    
        # Compute levels for the contour plot
        TFoct3dB = 10*np.log10(TFoct3[direction]**2)+180
        min_level = np.ceil(np.min(TFoct3dB) / 10) * 10
        max_level = np.ceil(np.max(TFoct3dB))//10*10
        levels = np.arange(min_level,max_level,2)
        
        # Plot contours of transfer functions
        contour = ax2.contourf(TFoct3dB.T, levels=levels, extend='both', cmap='jet', origin='lower')
        
        # Set axes
        ax2.set_ylabel(r'distance $(m)$')
        ax2.set_xlabel(r'1/3 octave center frequency $(Hz)$')
        ax2.set_xticks(range(1, Nfc, 2), fc[1::2])
        ax2.set_yticks(range(0,Npos,2),np.round(x[0::2],4))
        
        # Colorbar
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cbar = plt.colorbar(contour, cax=cax, extendrect=True, ticks=np.arange(min_level,max_level,10))
    
        # Plot track decay rate 
        ax3.semilogy(DR[direction], 'b-o', label='TDR')
        
        # Uncomment this to show TDR results from S&M Matlab scripts
        # ax3.semilogy(DR['mat'][direction], 'r--*', label='TDR S&M')
        
        # Plot limits
        nLimit = [list(fc).index(f) for f in limits['freq']]
        ax3.semilogy(nLimit, limits[direction], 'k-', label=f"ref. {direction}")
        
        # Set axes properties
        ax3.set_ylabel(r'$dB/m$')
        ax3.set_xlabel(r'1/3 octave center frequency $(Hz)$')
        ax3.set_title('Track Decay Rate')
        ax3.set_xlim([0, Nfc-1])
        ax3.set_xticks(range(1, Nfc, 2), fc[1::2])
        ax3.set_ylim([0.1, 100])
        ax3.grid(True, which='both', linestyle='dotted')
        ax3.legend(loc='upper right')
        
    #Adjust distance between subplots
    plt.subplots_adjust(hspace=0.5)
    

#%% Plot comparison of TDRs for different measurement
# Uncomment and mofidy if you want to change folders to plot
# folders = ['01NOV22B1', '01NOV22A2', '24MAY23B', '24MAY23A']

fig, ax = plt.subplots(2, 1)

# Define linecolors
linecolors = plt.cm.tab20([1,0, 7,6, 3,2 ,5,4, 9,8, 11,10, 13,12, 15,14, 17,16])

for j, direction in enumerate(['horizontal','vertical']):
    ax1 = ax[j]
    for i, folder in enumerate(folders): 
        cal = read_calsheet(folder)
        label = f"{cal['projectname']} {cal['projectdate']}"
        
        # Define style
        if 'ref' in label.lower():
            linestyle = '-o'
            fillstyle = 'full'
        else:
            linestyle = '-o'
            fillstyle = 'full'
        
        ax1.semilogy(TDR[folder][direction], linestyle, fillstyle=fillstyle, label=label, color=linecolors[i])
        
        # Uncomment this if you want to plot result from Matlab S&M TDR programme
        #ax1.semilogy(TDR[folder]['mat'][direction], '--*', label=f"{label} S&M", color=linecolors[i+len(folders)])
    
    #Plot Limits
    ax1.semilogy(nLimit, limits[direction], 'k-', label=f"ref. {direction}")
    
    # Set axis properties
    ax1.set_ylabel(r'$dB/m$')
    ax1.set_xlabel(r'1/3 octave center frequency $(Hz)$')
    ax1.set_xlim([0, Nfc-1])
    ax1.set_xticks(range(1, Nfc, 2), fc[1::2])
    ax1.set_ylim([0.1, 100])
    ax1.grid(True, which='both', linestyle='dotted')
    ax1.legend(loc='upper right')
    ax1.set_title('Track Decay Rate')

plt.show()
