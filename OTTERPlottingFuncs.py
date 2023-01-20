#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Plotting code supporting the OTTERPy modeling code by Clarke Delisle and
    Brian J. Yanites, supporting the manuscript 'Rethinking variability in bedrock rivers:
    variability of hillslope sediment supply modulates bedrock incision during floods'
    to be submitted to JGR Earth Surface
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.special as sc

eta_colormap = 'YlOrRd'
kv_colormap = 'YlGnBu_r'
BedrockColor = "gray"
WidthColor = "forestgreen"
ErosionColor = "firebrick"
SedimentColor = "burlywood"
QwColor = "powderblue"
SlopeColor = "plum"
ElevationColor = "black"

""" READ THIS

To test these plotting functions for the two datasets provided on the github,
use the following lines:"""

#This will plot a model run for 1.5Myr at 0.1mm/yr with kv = 1 and eta = 1
uplift = 0.0001 # Rock uplift rate for a single model
kv = 1 # discharge variabilty parameter
eta = 1 # sediment variablity parameter
totaltime = 1500000 # total model time
dt = 1/24 # time step

#Uncomment this to plot a model run for 1.5Myr at 0.1mm/yr with kv = 1 and eta = 1
# =============================================================================
# uplift = 0.005 # Rock uplift rate for a single model
# kv = 1 # discharge variabilty parameter
# eta = 1 # sediment variablity parameter
# totaltime = 1500000 # total model time
# dt = 1/24 # time step
# =============================================================================

""" The above variables should be edited to plot user-generated model data"""

#What node should we plot? for any function that plots data only from one node
timeseries_location = 30 
    
nsaves = int(10000) # How many timesteps were saved in the model run?

modelname = 'uplift' + str(int(uplift*10000)) + '_kv' +str(int(kv*10)) \
     + '_eta' +str(int(eta*10)) + '_time' + str(totaltime) + ".pkl"

modelname_eq = 'uplift' + str(int(uplift*10000)) + '_kv' +str(int(kv*10)) \
     + '_eta' +str(int(eta*10)) + '_time' + str(totaltime) + "Eq.pkl"

def OTTERoverviewplotting(modelname,upliftrate,location,totaltime,firsttime):
    
    """ This functions plots a time series of channel slope, elevation, width
    sediment cover, and bedrock incision for the entire model time series plotted
    at a frequency determined by the number of saves specified in the OTTERPy code"""
    
    df = pd.read_pickle(modelname)  
    uplift = upliftrate
    savestep = totaltime/nsaves
    time = np.arange(0,totaltime+savestep,savestep)
    
    Width = np.asarray(df.Width)
    Width = np.asarray(Width[0,])
    Slope = np.asarray(df.Slope)
    Slope = np.asarray(Slope[0,])
    Z = np.asarray(df.Z)
    Z = np.asarray(Z[0,])
    Sed = np.asarray(df.Sed_Depth)
    Sed = np.asarray(Sed[0,])
    Qw = np.asarray(df.Qw)
    Qw = np.asarray(Qw[0,])
       
    Zdiff = -uplift * dt - np.diff(Z[location,:]) 
    Zdiff[Zdiff<0] = 0
    Zdiff = abs(Zdiff)
    Zdiff = Zdiff/(totaltime/nsaves)
    
    fig, (ax,ax2) = plt.subplots(2,1,figsize = (12,6))
    ax.set_title("U = " + str(uplift*1000) + "mm/yr, k = " + str(kv) + ", eta = " + str(eta))
    fig.subplots_adjust(right=0.75)
    twin1 = ax.twinx()
    twin2 = ax.twinx()
    twin2.spines.right.set_position(("axes", 1.1))
    p1 = ax.scatter(time[firsttime:-2], Qw[location,firsttime:-2], 10,color = "powderblue",edgecolors="cornflowerblue", label="Water Discharge (cms)")
    p2, = twin1.plot(time[firsttime:-2], Sed[location,firsttime:-2], color ="burlywood", label="Sediment Cover (m)")
    twin1.fill_between(time[firsttime:-2], Sed[location,firsttime:-2],alpha = 0.4, color="tan")
    p3, = twin2.plot(time[firsttime:-2], Zdiff[firsttime:-1]*1000/dt, linewidth = 0.5, color ="crimson", label="River Incision Rate (mm/yr)")
    twin2.axhline(y=uplift, xmin=min(time), xmax=max(time), c="black", linewidth=2, zorder=0)
    ax.set_ylabel("River Discharge (cms)")
    twin1.set_ylabel("Sediment Cover (m)")
    twin2.set_ylabel("River Incision Rate (mm/yr)")
    ax.yaxis.label.set_color("cornflowerblue")
    twin1.yaxis.label.set_color(p2.get_color())
    twin2.yaxis.label.set_color(p3.get_color())
    tkw = dict(size=4, width=1.5)
    ax.tick_params(axis='y', color="cornflowerblue", **tkw)
    twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    ax.tick_params(axis='x', **tkw)
    twin2.legend(handles=[p1, p2, p3], loc='upper right')
    
    twin21 = ax2.twinx()
    twin22 = ax2.twinx()
    twin22.spines.right.set_position(("axes", 1.1))
    p21, = ax2.plot(time[firsttime:-2], Slope[location,firsttime:-2],linewidth=2, color = "plum", label="Channel Slope")
    p22, = twin21.plot(time[firsttime:-1], Width[location,firsttime:-1], linewidth=2,color ="forestgreen", label="Channel Width (m)")
    p23, = twin22.plot(time[firsttime:-1], Z[location,firsttime:-1], linewidth=2.5, color ="black", label="Elevation (m)")
    ax2.set_xlabel("Time (years)")
    ax2.set_ylabel("Channel Slope")
    twin21.set_ylabel("Channel Width (m)")
    twin22.set_ylabel("Elevation (m)")
    ax2.yaxis.label.set_color(p21.get_color())
    twin21.yaxis.label.set_color(p22.get_color())
    twin22.yaxis.label.set_color(p23.get_color())
    tkw = dict(size=4, width=1.5)
    ax2.tick_params(axis='y', colors=p21.get_color(), **tkw)
    twin21.tick_params(axis='y', colors=p22.get_color(), **tkw)
    twin22.tick_params(axis='y', colors=p23.get_color(), **tkw)
    ax2.tick_params(axis='x', **tkw)
    twin22.legend(handles=[p21, p22, p23],loc='upper right')
    
    return 

def OTTEREqplotting(modelname,upliftrate,location,lasttime, kv, eta):
    
    """ This functions plots a short time series of channel slope, elevation, width
    sediment cover, and bedrock incision for a short time series plotted every time
    step once the model reaches equilibrium. """
    
    df = pd.read_pickle(modelname_eq)  
    uplift = upliftrate
    time = np.arange(0,12000)
    Width = np.asarray(df.Width)
    Width = np.asarray(Width[0,])
    Slope = np.asarray(df.Slope)
    Slope = np.asarray(Slope[0,])
    Z = np.asarray(df.Z)
    Z = np.asarray(Z[0,])
    Sed = np.asarray(df.Sed_Depth)
    Sed = np.asarray(Sed[0,])
    Qw = np.asarray(df.Qw)
    Qw = np.asarray(Qw[0,])
        
    Zdiff = -uplift * dt - np.diff(Z[location,:]) 
    Zdiff[Zdiff<0] = 0
    Zdiff = abs(Zdiff)
        
    fig, (ax,ax2) = plt.subplots(2,1,figsize = (12,6))
    fig.subplots_adjust(right=0.8)
    ax.ticklabel_format(useOffset=False)
    twin1 = ax.twinx()
    twin2 = ax.twinx()
    twin1.ticklabel_format(useOffset=False)
    twin2.ticklabel_format(useOffset=False)
    twin2.spines.right.set_position(("axes", 1.1))
    
    p1 = ax.scatter(time[0:lasttime]*dt, Qw[location,0:lasttime], 10,color = QwColor, edgecolors="cornflowerblue", label="Water Discharge (cms)")
    p2, = twin1.plot(time[0:lasttime]*dt, Sed[location,0:lasttime], color = SedimentColor, label="Sediment Cover (m)")
    twin1.fill_between(time[0:lasttime]*dt, Sed[location,0:lasttime],alpha = 0.4, color=SedimentColor)
    p3, = twin2.plot(time[0:lasttime]*dt, Zdiff[0:lasttime]*1000/dt, linewidth = 0.5, color =ErosionColor, label="River Incision Rate (mm/yr)")
   
    ax.set_ylabel("River Discharge (cms)")
    ax.set_title("U = " + str(uplift*1000) + "mm/yr, k = " + str(kv) + ", eta = " + str(eta))
    twin1.set_ylabel("Sediment Cover (m)")
    twin2.set_ylabel("River Incision Rate (mm/yr)")
    ax.yaxis.label.set_color("cornflowerblue")
    twin1.yaxis.label.set_color(p2.get_color())
    twin2.yaxis.label.set_color(p3.get_color())
    tkw = dict(size=4, width=1.5)
    ax.tick_params(axis='y', color="cornflowerblue", **tkw)
    twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    ax.tick_params(axis='x', **tkw)
    twin2.legend(handles=[p1, p2, p3], loc='upper right')
            
    twin21 = ax2.twinx()
    twin22 = ax2.twinx()
    ax2.ticklabel_format(useOffset=False)
    twin21.ticklabel_format(useOffset=False)
    twin22.ticklabel_format(useOffset=False)
    twin22.spines.right.set_position(("axes", 1.1))
    
    p21, = ax2.plot(time[0:lasttime]*dt, Slope[location,0:lasttime],linewidth=2, color = SlopeColor, label="Channel Slope")
    p22, = twin21.plot(time[0:lasttime]*dt, Width[location,0:lasttime], linewidth=2,color = WidthColor, label="Channel Width (m)")
    p23, = twin22.plot(time[0:lasttime]*dt, Z[location,0:lasttime], linewidth=2.5, color = ElevationColor, label="Elevation (m)")
    
    ax2.set_xlabel("Time (years)")
    ax2.set_ylabel("Channel Slope")
    twin21.set_ylabel("Channel Width (m)")
    twin22.set_ylabel("Elevation (m)")
    ax2.yaxis.label.set_color(p21.get_color())
    twin21.yaxis.label.set_color(p22.get_color())
    twin22.yaxis.label.set_color(p23.get_color())
    tkw = dict(size=4, width=1.5)
    ax2.tick_params(axis='y', colors=p21.get_color(), **tkw)
    twin21.tick_params(axis='y', colors=p22.get_color(), **tkw)
    twin22.tick_params(axis='y', colors=p23.get_color(), **tkw)
    ax2.tick_params(axis='x', **tkw)
    twin22.legend(handles=[p21, p22, p23],loc='upper right').set_zorder(5)
    
    return 
    
def OTTERsingletimeplotting(modelname,totaltime,singletime):
    
    """ This function plots relationships of channel width, slope, and 
    sediment cover vs distance downstream for a single time in the model """
    
    df = pd.read_pickle(modelname)  
    x = np.arange(1,100000,2000)
    
    Width = np.asarray(df.Width)
    Width = np.asarray(Width[0,])
    Slope = np.asarray(df.Slope)
    Slope = np.asarray(Slope[0,])
    Z = np.asarray(df.Z)
    Z = np.asarray(Z[0,])
    Sed = np.asarray(df.Sed_Depth)
    Sed = np.asarray(Sed[0,])
    Qw = np.asarray(df.Qw)
    Qw = np.asarray(Qw[0,])
    
    fig, ax  = plt.subplots(figsize = (12,4))
    fig.subplots_adjust(right=0.75)
    twin1 = ax.twinx()
    twin2 = ax.twinx()
    twin2.spines.right.set_position(("axes", 1.1))
    p1, = ax.plot(x[2:-2], Sed[2:-2,singletime],linewidth=2, color = "tan", label="Sediment Thickness (m)")
    p2, = twin1.plot(x[2:-2], Width[2:-2,singletime], linewidth=2,color ="forestgreen", label="Channel Width (m)")
    p3, = twin2.plot(x[2:-2], Z[2:-2,singletime], linewidth=2, color ="black", label="Bedrock Elevation (m)")
    ax.set_xlabel("Distance (m)")
    ax.set_ylabel("Channel Elevation [Bedrock + Sediment] (m)")
    twin1.set_ylabel("Channel Width (m)")
    twin2.set_ylabel("Bedrock Elevation (m)")
    ax.set_ylim(0,np.max(Sed[:,singletime]))
    twin2.set_ylim(0,np.max(Z[:,singletime]+Sed[:,singletime]))
    ax.yaxis.label.set_color(p1.get_color())
    twin1.yaxis.label.set_color(p2.get_color())
    twin2.yaxis.label.set_color(p3.get_color())
    tkw = dict(size=4, width=1.5)
    ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
    twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    ax.tick_params(axis='x', **tkw)
    ax.legend(handles=[p1, p2, p3],loc='upper center')
    
    return 

# Calling the functions above

OTTERoverviewplotting(modelname,uplift,timeseries_location,totaltime,1000)

OTTEREqplotting(modelname_eq,uplift,timeseries_location,4000,kv,eta)

OTTERsingletimeplotting(modelname,totaltime,9000)

