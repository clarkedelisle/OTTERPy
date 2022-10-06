#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 12:25:41 2022

@author: clarkedelisle
"""

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.special as sc
from matplotlib.lines import Line2D
from scipy import stats
from scipy.optimize import fsolve
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)

eta_colormap = 'YlOrRd'
kv_colormap = 'YlGnBu_r'
BedrockColor = "gray"
WidthColor = "forestgreen"
ErosionColor = "firebrick"
SedimentColor = "burlywood"
QwColor = "powderblue"
SlopeColor = "plum"
ElevationColor = "black"

uplift = 0.005 # Rock uplift rate for a single model
kv = 3 # discharge variabilty parameter
eta = 1 # sediment variablity parameter
totaltime = 2000000 # total model time
dt = 1/24 # time step

#What node should we plot? for any function that plots data only from one node
timeseries_location = 30 
    
nsaves = int(10000) # How many timesteps were saved in the model run?

modelname = 'uplift' + str(int(uplift*10000)) + '_kv' +str(int(kv*10)) \
     + '_eta' +str(int(eta*10)) + '_time' + str(totaltime) + ".pkl"

modelname_eq = 'uplift' + str(int(uplift*10000)) + '_kv' +str(int(kv*10)) \
     + '_eta' +str(int(eta*10)) + '_time' + str(totaltime) + "Eq.pkl"

def OTTERoverviewplotting(modelname,upliftrate,location,totaltime,firsttime):
    
    """ This functions plots a short time series of channel slope, elevation, width
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

def OTTERSadler(time):
    
    """ This functions calculates and plots coefficients of variations for distributions
    of incision rates measured over different time windows. We use this to investigate 
    the strength of the Sadler effect in bedrock rivers with different model 
    parameters """
    
    k_arr = [3, 2.5, 2, 1.5, 1, 0.6, 0.3]
    eta_arr = [2.0,1.75,1.5, 1.25,1.0,0.75,0.5,0.25]
    U_arr = [0.0001]

    timestep = time/nsaves
    EqLength = 5000

    fig, ax1 = plt.subplots(figsize = (12,5))
   
    ax2 = plt.axes([0,0,1,1])
    ip = InsetPosition(ax1, [0.45,0.45,0.45,0.45])
    ax2.set_axes_locator(ip)
   
    pointsize = 75
    
    counter = -1
    
    column_names = ["ModelName","uplift","kv","eta","std1","std2",\
                    "std3","std4","std5","std6","std7","std8","std9","std10"]
    
    df1 = pd.DataFrame(columns = column_names)
            
    for i in range(0,(len(U_arr))):
        uplift = U_arr[i]
        
        for j in range(0,(len(k_arr))):
            kv = k_arr[j]

            for k in range(0,(len(eta_arr))):
                eta = eta_arr[k]
                counter = int(counter + 1)
                
                ModelName = 'uplift' + str(int(uplift*10000)) + '_kv' +str(int(kv*10)) \
                    + '_eta' +str(int(eta*10)) + '_time' + str(time)
                    
                if ModelName == 'uplift1_kv25_eta17_time2000000':
                    continue
                
                if ModelName == 'uplift1_kv25_eta15_time2000000':
                    continue
                
                if ModelName == 'uplift1_kv25_eta2_time2000000':
                    continue
                
                if ModelName == 'uplift1_kv10_eta2_time2000000':
                    continue
                
                if ModelName == 'uplift1_kv3_eta12_time2000000':
                    continue
                    
                df = pd.read_pickle(str(ModelName) + '.pkl')
                Z = np.asarray(df.Z)
                Z = np.asarray(Z[0,])
                Z_Eq = Z[30,(-EqLength-2):-2]
                
                nsamples = 70
     
                store1 = np.zeros(nsamples)
                store2 = np.zeros(nsamples)
                store3 = np.zeros(nsamples)
                store4 = np.zeros(nsamples)
                store5 = np.zeros(nsamples)
                store6 = np.zeros(nsamples)
                store7 = np.zeros(nsamples)
                store8 = np.zeros(nsamples)
                store9 = np.zeros(nsamples)
                store10 = np.zeros(nsamples)
                
                d1 = int(200 // timestep)
                d2 = int(600 // timestep)
                d3 = int(1000 // timestep)
                d4 = int(1500 // timestep)
                d5 = int(2000 // timestep)
                d6 = int(2500 // timestep)
                d7 = int(3000 // timestep)
                d8 = int(3500 // timestep)
                d9 = int(4000 // timestep)
                d10 = int(4500 // timestep)

                
                for a in range(0,(len(store1)-1)):
                    if a == 0: 
                        store1[a] = (Z_Eq[a+d1] - Z_Eq[a])/(timestep*d1) + uplift
                    else:
                        store1[a] = (Z_Eq[(a+1)*d1] - Z_Eq[a*d1])/(timestep*d1) + uplift 

                for b in range(0,(len(store2)-1)):
                    if b == 0: 
                        store2[b] = (Z_Eq[b+d2] - Z_Eq[b])/(timestep*d2) + uplift
                    else:
                        store2[b] = (Z_Eq[(b+1)*d2] - Z_Eq[b*d2])/(timestep*d2) + uplift 
                        
                for c in range(0,(len(store3)-1)):
                    if c == 0: 
                        store3[c] = (Z_Eq[c+d3] - Z_Eq[c])/(timestep*d3) + uplift
                    else:
                        store3[c] = (Z_Eq[(c+1)*d3] - Z_Eq[c*d3])/(timestep*d3) + uplift 
                
                for d in range(0,(len(store4)-1)):
                    if d == 0: 
                        store4[d] = (Z_Eq[d+d4] - Z_Eq[d])/(timestep*d4) + uplift
                    else:
                        store4[d] = (Z_Eq[(d+1)*d4] - Z_Eq[d*d4])/(timestep*d4) + uplift 
                        
                for e in range(0,(len(store5)-1)):
                    if e == 0: 
                        store5[e] = (Z_Eq[e+d5] - Z_Eq[e])/(timestep*d5) + uplift
                    else:
                        store5[e] = (Z_Eq[(e+1)*d5] - Z_Eq[e*d5])/(timestep*d5) + uplift 
                
                for f in range(0,(len(store6)-1)):
                    if f == 0: 
                        store6[f] = (Z_Eq[f+d6] - Z_Eq[f])/(timestep*d6) + uplift
                    else:
                        store6[f] = (Z_Eq[(f+1)*d6] - Z_Eq[f*d6])/(timestep*d6) + uplift 
                        
                for g in range(0,(len(store7)-1)):
                    if b == 0: 
                        store7[g] = (Z_Eq[g+d7] - Z_Eq[g])/(timestep*d7) + uplift
                    else:
                        store7[g] = (Z_Eq[(g+1)*d7] - Z_Eq[g*d7])/(timestep*d7) + uplift 
                        
                for h in range(0,(len(store8)-1)):
                    if h == 0: 
                        store8[h] = (Z_Eq[h+d8] - Z_Eq[h])/(timestep*d8) + uplift
                    else:
                        store8[h] = (Z_Eq[(h+1)*d8] - Z_Eq[h*8])/(timestep*d8) + uplift 
                        
                for l in range(0,(len(store9)-1)):
                    if l == 0: 
                        store9[l] = (Z_Eq[l+d9] - Z_Eq[l])/(timestep*d9) + uplift
                    else:
                        store9[l] = (Z_Eq[(l+1)*d9] - Z_Eq[l*d9])/(timestep*d9) + uplift 
                        
                for m in range(0,(len(store10)-1)):
                    if m == 0: 
                        store10[m] = (Z_Eq[m+d10] - Z_Eq[m])/(timestep*d10) + uplift
                    else:
                        store10[m] = (Z_Eq[(m+1)*d10] - Z_Eq[m*d10])/(timestep*d10) + uplift 
                        
                      
                std1 = np.std(store1) / uplift
                std2 = np.std(store2) / uplift   
                std3 = np.std(store3) / uplift   
                std4 = np.std(store4) / uplift   
                std5 = np.std(store5) / uplift   
                std6 = np.std(store6) / uplift   
                std7 = np.std(store7) / uplift  
                std8 = np.std(store8) / uplift  
                std9 = np.std(store9) / uplift  
                std10 = np.std(store10) / uplift  

                df1 = df1.append({"ModelName":ModelName,"uplift":uplift,"kv":kv,"eta":eta,"std1":std1,\
                                "std2":std2,"std3":std3,"std4":std4,"std5":std5,"std6":std6,"std7":std7,\
                                "std8":std8,"std9":std9,"std10":std10\
                                 },ignore_index = True)
                               
                plot1 = np.zeros(len(df1.std1))+(d1*timestep)
                plot2 = np.zeros(len(df1.std2))+(d2*timestep)
                plot3 = np.zeros(len(df1.std3))+(d3*timestep)
                plot4 = np.zeros(len(df1.std4))+(d4*timestep)
                plot5 = np.zeros(len(df1.std5))+(d5*timestep)
                plot6 = np.zeros(len(df1.std6))+(d6*timestep)
                plot7 = np.zeros(len(df1.std7))+(d7*timestep)
                plot8 = np.zeros(len(df1.std8))+(d8*timestep)
                plot9 = np.zeros(len(df1.std9))+(d9*timestep)
                plot10 = np.zeros(len(df1.std10))+(d10*timestep)
                
                colored_by = df1.kv
                colormap = kv_colormap
                
                ax1.set_axisbelow(True)
                ax1.grid(True, which="both",color='lightgray', linestyle='--')
                ax1.scatter(plot1,df1.std1,s=pointsize,c=colored_by,cmap=colormap,edgecolor='gray')
                ax1.scatter(plot2,df1.std2,s=pointsize,c=colored_by,cmap=colormap,edgecolor='gray')
                ax1.scatter(plot3,df1.std3,s=pointsize,c=colored_by,cmap=colormap,edgecolor='gray')
                ax1.scatter(plot4,df1.std4,s=pointsize,c=colored_by,cmap=colormap,edgecolor='gray')
                ax1.scatter(plot5,df1.std5,s=pointsize,c=colored_by,cmap=colormap,edgecolor='gray')
                ax1.scatter(plot6,df1.std6,s=pointsize,c=colored_by,cmap=colormap,edgecolor='gray')
                ax1.scatter(plot7,df1.std7,s=pointsize,c=colored_by,cmap=colormap,edgecolor='gray')
                ax1.scatter(plot8,df1.std8,s=pointsize,c=colored_by,cmap=colormap,edgecolor='gray')
                ax1.scatter(plot9,df1.std9,s=pointsize,c=colored_by,cmap=colormap,edgecolor='gray')
                ax1.scatter(plot10,df1.std10,s=pointsize,c=colored_by,cmap=colormap,edgecolor='gray')
                ax1.set_xlabel('Measurement Window (Years)')
                ax1.set_ylabel('Coefficient of Variation in Incision Rate')
                ax1.set_title('Ur = 5mm/yr')
                ax1.set_ylim([0.1,0.5])
               
                ax2.grid(True, which="both",color='lightgray', linestyle='--', zorder=0)
                ax2.scatter(df1.eta, df1.std1, s=pointsize, c=colored_by,cmap=colormap,edgecolor='gray',zorder = 1000)
                ax2.set_xlabel('Sediment Variability eta')
                ax2.set_ylabel('CV')
                ax2.set_title('t=200years')
                ax2.set_ylim([0.1,0.5])
                                
    if colormap == kv_colormap:
        fig3, ax3 = plt.subplots(figsize=(6, 1))
        fig3.subplots_adjust(bottom=0.5)
        fig3.colorbar(mpl.cm.ScalarMappable(cmap=kv_colormap),
            cax=ax3, orientation='horizontal', label='Discharge Variability k')
          
    if colormap == eta_colormap:
        fig4, ax4 = plt.subplots(figsize=(6, 1))
        fig4.subplots_adjust(bottom=0.5)
        fig4.colorbar(mpl.cm.ScalarMappable(cmap=eta_colormap),
            cax=ax4, orientation='horizontal', label='Sediment Variability eta')
                     
    return

def OTTERTrends(time):
    
    """
    OTTERTrends reads in data from an entire suite of model runs with varying
    k, eta, and uplift and plots trends in width, slope, and sediment cover
    at a one upstream location and one downstream location. Inputs are the 
    k, eta, and U arrays defined in this script globally, and the time over 
    which all models were run in years. 
    """
    
    k_arr = [0.3, 0.6, 1, 1.5, 2, 2.5, 3]
    eta_arr = [0.25,0.5,0.75,1.0,1.25,1.5,1.75,2]
    U_arr = [0.005]
    
    column_names = ["ModelName","uplift","kv","eta","S_US","S_DS","W_US",\
                    "W_DS","SedThicknessUS","SedThicknessDS","QwMax"]
        
    df1 = pd.DataFrame(columns = column_names)
        
    counter = -1
    
    fig1,ax1 = plt.subplots(1,2,figsize = (12,5))
    fig2,ax2 = plt.subplots(1,2,figsize = (12,5))
    
    for i in range(0,(len(U_arr))):
        uplift = U_arr[i]
        
        for j in range(0,(len(k_arr))):
            kv = k_arr[j]

            for k in range(0,(len(eta_arr))):
                eta = eta_arr[k]
                counter = int(counter + 1)
                
                ModelName = 'uplift' + str(int(uplift*10000)) + '_kv' +str(int(kv*10)) \
                    + '_eta' +str(int(eta*10)) + '_time' + str(time)
                
                df = pd.read_pickle(str(ModelName) + '.pkl')
                Width = np.asarray(df.Width)
                Width = np.asarray(Width[0,])
                Slope = np.asarray(df.Slope)
                Slope = np.asarray(Slope[0,])
                Sed = np.asarray(df.Sed_Depth)
                Sed = np.asarray(Sed[0,])
                Z = np.asarray(df.Z)
                Z = np.asarray(Z[0,])
                Qw = np.asarray(df.Qw)
                Qw = np.asarray(Qw[0,])
                
                SlopeUS = (np.mean(Slope[5,int(0.8*nsaves):nsaves]))
                SlopeDS = (np.mean(Slope[45,int(0.8*nsaves):nsaves]))

                WidthUS = (np.mean(Width[5,int(0.8*nsaves):nsaves]))
                WidthDS = (np.mean(Width[45,int(0.8*nsaves):nsaves]))
                
                SedUS = (np.mean(Sed[5,int(0.8*nsaves):nsaves]))
                SedDS = (np.mean(Sed[45,int(0.8*nsaves):nsaves]))
                
                QwMax = np.max(Qw)
                
                df1 = df1.append({"ModelName":ModelName,"uplift":uplift,"kv":kv,"eta":eta,\
                                "S_US":SlopeUS,"S_DS":SlopeDS,"W_US":WidthUS,"W_DS":WidthDS,\
                                "SedThicknessUS":SedUS,"SedThicknessDS":SedDS,"QwMax":QwMax},ignore_index = True)

    ax1[0].set_axisbelow(True)
    ax1[0].grid(color='lightgray', linestyle='dashed')
    ax1[0].scatter(df1.eta,df1.S_DS,s=(abs(df1.uplift*30000+10)),c=df1.kv,cmap=kv_colormap,edgecolor='gray')
    ax1[0].set_xlabel("Sediment Variability Exponent")
    ax1[0].set_ylabel("Equilibrium Channel Slope")
    ax1[0].set_yscale('log')
    
    ax1[1].set_axisbelow(True)
    ax1[1].grid(color='lightgray', linestyle='dashed')
    ax1[1].scatter(df1.eta,df1.W_DS,s=(abs(df1.uplift*30000+10)),c=df1.kv,cmap=kv_colormap,edgecolor='gray')
    ax1[1].set_xlabel("Sediment Variability Exponent")
    ax1[1].set_ylabel("Equilibrium Channel Width (m)")
    
    fig1.suptitle('10km Upstream of Channel Mouth')
        
    ax2[0].set_axisbelow(True)
    ax2[0].grid(color='lightgray', linestyle='dashed')
    ax2[0].scatter(df1.eta,df1.S_US,s=(abs(df1.uplift*30000+10)),c=df1.kv,cmap=kv_colormap,edgecolor='gray')
    ax2[0].set_xlabel("Sediment Variability Exponent")
    ax2[0].set_ylabel("Equilibrium Channel Slope Upstream")
    ax2[0].set_yscale('log')
   
    ax2[1].set_axisbelow(True)
    ax2[1].grid(color='lightgray', linestyle='dashed')
    ax2[1].scatter(df1.eta,df1.W_US,s=(abs(df1.uplift*30000+10)),c=df1.kv,cmap=kv_colormap,edgecolor='gray')
    ax2[1].set_xlabel("Sediment Variability Exponent")
    ax2[1].set_ylabel("Equilibrium Channel Width Upstream")
    
    fig2.suptitle('10km Downstream of Headwaters')
    
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)
    fig.colorbar(mpl.cm.ScalarMappable(cmap=kv_colormap),
             cax=ax, orientation='horizontal', label='Discharge Variability k')
   
    return

def OTTERHiatuses(time):
    
    """ This function plots distributions of hiatuses in bedrock incision
    for models run with different endmember model parameters"""
     
    k_arr = [3.0, 0.3]
    eta_arr = [0.5, 2.0]
    U_arr = [0.0001,0.001, 0.005]

    colors = ['lightsteelblue','dodgerblue','navy']
    
    fig, axes = plt.subplots(figsize=(10,8),ncols=2, nrows=2)
    
    for i in range(0,(len(U_arr))):
        uplift = U_arr[i]
        
        for j in range(0,(len(k_arr))):
            kv = k_arr[j]

            for k in range(0,(len(eta_arr))):
                eta = eta_arr[k]
                
                ModelName = 'uplift' + str(int(uplift*10000)) + '_kv' +str(int(kv*10)) \
                    + '_eta' +str(int(eta*10)) + '_time' + str(time)
                    
                df = pd.read_pickle(str(ModelName) + '.pkl')
                F = np.asarray(df.F)
                F = np.asarray(F[0,])
                F = F[20,:]
                
                Hiatuses = np.load(ModelName + 'Hiatus.npy')
                Hiatuses = Hiatuses * dt
                Hiatuses = Hiatuses[Hiatuses>(0)]
                Hiatuses = Hiatuses[Hiatuses<(250000)]
                Hiatuses = Hiatuses[-50000:]

                
                custom_lines = [Line2D([0], [0], color='lightsteelblue', lw=4),
                                Line2D([0], [0], color='dodgerblue', lw=4),
                                Line2D([0], [0], color='navy', lw=4)]

                axes[0,0].legend(custom_lines, ['U = 0.1mm/yr', 'U = 5mm/yr'])
            
                sns.kdeplot(data=Hiatuses,ax=axes[j,k],bw_adjust=3,color=colors[i])
                axes[j,k].set_xlabel('Incision Hiatius (Years)')
                axes[j,k].set_ylabel('CDF')
                axes[j,k].set_xlim([0.03,50])
                               
                plt.show()

    return 

def OTTERGif(ModelName,time):
    
    """
    OTTERGif plots and saves a .gif for the river evolution of one specific model 
    run. Inputs are the Model Name as defined above, and the time over which the
    model was run in years. 
    """
    
    fig, (ax1, ax3, ax5) = plt.subplots(3,1)
    ax2 = ax1.twinx()
    ax4 = ax3.twinx()
    ax6 = ax5.twinx()
    
    dsave = 1
    dx = 2000
    x = np.arange(5000,105000,dx)
    
    df = pd.read_pickle(ModelName)
    Width = np.asarray(df.Width)
    Width = np.asarray(Width[0,])
    Slope = np.asarray(df.Slope)
    Slope = np.asarray(Slope[0,])
    Z = np.asarray(df.Z)
    Z = np.asarray(Z[0,])
    Sed = np.asarray(df.Sed_Depth)
    Sed = np.asarray(Sed[0,])
    E = np.asarray(df.E)
    E = np.asarray(E[0,])
    
    Width = Width[:,-1000]
    Slope = Slope[:,-1000]
    Z = Z[:,-1000]
    Sed = Sed[:,-1000]
    E = E[:,-1000]
    
    def animate(i):
        timevar = i*dsave + dsave
        
        ax1.cla()
        ax1.plot(x[3:len(x)-2]/1000,Z[3:len(x)-2,i]/1000 + Sed[3:len(x)-2,i]/1000,lw=3,c=SedimentColor) 
        ax1.plot(x[3:len(x)-2]/1000,Z[3:len(x)-2,i]/1000,lw=3,c=BedrockColor)
        ax1.set_ylabel('Z (km)')
        ax1.set_title('Time = %i' %timevar)
        
        ax2.cla()
        ax2.plot(x[1:len(x)-2]/1000,Width[1:len(x)-2,i],lw=3,c=WidthColor)
        ax2.axvline(x=timeseries_location*2,lw=2,c = 'k')
        ax2.set_ylabel('Channel Width (m)',c=WidthColor)
        
        ax3.cla()
        ax3.plot(x[1:len(x)-2]/1000,-1*E[1:len(x)-2,i],lw=3,c=ErosionColor)
        ax3.set_ylabel('Bedrock Incision (m)',c=ErosionColor)
        ax3.set_xlabel('X (km)')
        
        ax4.cla()
        ax4.plot(x[1:len(x)-2]/1000,Sed[1:len(x)-2,i],lw=3,c=SedimentColor)
        ax4.set_ylabel('Sediment Cover (m)',c=SedimentColor)
        ax4.set_xlabel('X (km)')
          
        ax5.cla()
        ax5.plot((np.linspace(1,i*dsave,i)),-1*E[timeseries_location,0:i],lw=0.75,c=ErosionColor)
        ax5.scatter(i*dsave,-1*E[timeseries_location,i],c='red')
        ax5.set_ylabel('Bedrock Incision (m)',c=ErosionColor)
        ax5.set_xlim(0,time)
        
        ax6.cla()
        ax6.plot((np.linspace(1,i*dsave,i)),Sed[timeseries_location,0:i],lw=0.75,c=SedimentColor)
        ax6.scatter(i*dsave,Sed[timeseries_location,i],c='tan')
        ax6.set_ylabel('Sediment Cover (m)',c=SedimentColor)
        ax6.set_xlim(0,time)
    
    anim = animation.FuncAnimation(fig, animate, frames = 15000, interval = 1, blit = False)
    anim.save(str(ModelName) + '.gif', writer='PillowWriter', fps=8)
    
    return 

def OTTERQw_E(modelname,upliftrate,location,lasttime, kv, eta):
    
    """ This function plots individual and cumulative bedrock erosion and sediment
    aggradation/removal for floods discharge events across the distribution of Qw """
    
    df = pd.read_pickle(modelname_eq)  
    uplift = upliftrate
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
    F = np.asarray(df.F)
    F = np.asarray(F[0,])
        
    Zdiff = -uplift * dt - np.diff(Z[location,:]) 
    Zdiff[Zdiff<0] = 0
    Zdiff = abs(Zdiff) * 1000
    
    Seddiff = np.diff(Sed[location,:]) 
    
    pointsize = 60
    numbins = 30
    
    fig1, (ax1, ax2) = plt.subplots(1,2,figsize = (15,4))
    fig1.tight_layout(pad=3)
    fig1.subplots_adjust(wspace=0.3)
    
    ax1.set_axisbelow(True)
    ax1.grid(color='gray', linestyle='dashed')
    sc1 = ax1.scatter(Qw[location,0:lasttime],Zdiff[0:lasttime],c=Sed[location,0:lasttime],cmap = 'viridis',edgecolor='black',s=pointsize)    
    ax1.set_xlabel('Qw (cms)')
    ax1.set_ylabel('Bedrock Erosion (mm)')
    ax1.set_xscale('log')
    cbar = fig1.colorbar(sc1, ax=ax1,orientation = 'horizontal')
    cbar.set_label('Sediment Cover (m)')
    ax1.set_title("U = " + str(uplift*1000) + "mm/yr, k = " + str(kv) + ", eta = " + str(eta))
    
    ax2.set_axisbelow(True)
    ax2.grid(color='gray', linestyle='dashed')
    sc2 = ax2.scatter(Qw[location,0:lasttime],Seddiff[0:lasttime],c=Zdiff[0:lasttime],cmap = 'magma',edgecolor='black',s=pointsize)    
    ax2.set_xlabel('Qw (cms)')
    ax2.set_ylabel('Change in Sediment Thickness (m)')
    ax2.set_xscale('log')
    ax2.set_title("U = " + str(uplift*1000) + "mm/yr, k = " + str(kv) + ", eta = " + str(eta))
    cbar = fig1.colorbar(sc2, ax=ax2, orientation = 'horizontal')
    cbar.set_label('Bedrock Erosion (mm)')
    
    fig2, (ax3, ax4) = plt.subplots(1,2,figsize = (15,4))
    fig2.tight_layout(pad=3)
    fig2.subplots_adjust(wspace=0.3)
    
    Ebin_sums, Ebin_edges, Ebinnumber = stats.binned_statistic(np.log10(Qw[location,0:lasttime]), Zdiff[0:lasttime], statistic = 'sum', bins = numbins)
    Ebin_means, Ebin_edges, Ebinnumber = stats.binned_statistic(np.log10(Qw[location,0:lasttime]), Zdiff[0:lasttime], statistic = 'mean', bins = numbins)
    Ebin_counts, Ebin_edges, Ebinnumber = stats.binned_statistic(np.log10(Qw[location,0:lasttime]), Zdiff[0:lasttime], statistic = 'count', bins = numbins)
    center = (Ebin_edges[:-1] + Ebin_edges[1:]) / 2
    
    twin3 = ax3.twinx()
    ax3.set_axisbelow(True)
    ax3.grid(color='gray', linestyle='dashed')
    twin3.set_ylim([0,(np.max(Ebin_counts))*1.5])
    twin3.bar(10**center, Ebin_counts,width=0.15*np.array(10**center), color = 'gray', edgecolor='k')
    twin3.set_ylabel('Discharge Event Counts')
    sc3 = ax3.scatter(10**center, Ebin_sums,s = pointsize, c=Ebin_means,cmap = 'magma', edgecolor = 'k')
    ax3.set_xlabel('Qw (cms)')
    ax3.set_ylabel('Cumulative Erosion (mm)')
    ax3.set_title("U = " + str(uplift*1000) + "mm/yr, k = " + str(kv) + ", eta = " + str(eta))
    ax3.set_xscale('log')
    cbar = fig2.colorbar(sc3, ax=ax3, orientation = 'horizontal')
    cbar.set_label('Mean Bedrock Erosion (mm)')
    
    Sedbin_sums, Sedbin_edges, Sedbinnumber = stats.binned_statistic(np.log10(Qw[location,0:lasttime]), Seddiff[0:lasttime], statistic = 'sum', bins = numbins)
    Sedbin_means, Sedbin_edges, Sedbinnumber = stats.binned_statistic(np.log10(Qw[location,0:lasttime]), Seddiff[0:lasttime], statistic = 'mean', bins = numbins)

    center = (Sedbin_edges[:-1] + Sedbin_edges[1:]) / 2
    
    ax4.set_axisbelow(True)
    PhaseChange1 = (10**center[np.argmax(Sedbin_sums)] - 10**center[np.argmin(Sedbin_sums)]) * 0.3
    PhaseChange2 = (10**center[np.argmax(Sedbin_sums)] - 10**center[np.argmin(Sedbin_sums)]) * 0.7
    ax4.grid(color='gray', linestyle='dashed', zorder = 0)
    ax4.axvspan(PhaseChange1, PhaseChange2, facecolor = 'lightgray', alpha=0.5, zorder = 1)
    ax4.axhline(y=0,linewidth = 2, color = 'k', zorder = 2)
    sc4 = ax4.scatter(10**center, Sedbin_sums,s = pointsize, c=Sedbin_means,cmap = 'viridis', edgecolor = 'k', zorder = 5)
    ax4.set_xlabel('Qw (cms)')
    ax4.set_xscale('log')
    ax4.set_ylabel('Cumulative Change in Sed Cover (m)')
    cbar = fig2.colorbar(sc4, ax=ax4, orientation = 'horizontal')
    cbar.set_label('Mean Change in Sed Cover (m)')
    ax4.set_title("U = " + str(uplift*1000) + "mm/yr, k = " + str(kv) + ", eta = " + str(eta))
    
    return
           
def OTTERSedSupplyVsTrans(location):
    
    """ This function plots the tradeoff between sediment supply and tranpsort for discharge
    events across the distribution of Qw"""
    
    x = np.arange(5000, 105000, 2000, dtype='float')
    Hc = 1
    He = 1.8598
    Ah = Hc * ((x) ** He)
    dA = np.zeros(len(x))
    dA[1:100] = np.diff(Ah)
    kv = 0.3
    eta_sedvar = 2
    pw = 1000
    ps = 2650
    n = 0.04
    W = 8
    S = 0.04
    Htest = 1
    D = 0.005
    Tc = 0.045
    kQ = 10 ** (-8)
    eQ = 1
    Qwo = kQ * (Ah ** eQ)
    
    tarray = np.arange(1, 10000, 1, dtype='float')

    def calc_E(a,n,simpleorcomplex,eta_sedvar):
 
        prb = 1 / n
    
        if simpleorcomplex == 1:
            sum = 0
            for i in range(0, n):
                sum += (a[i] * prb)
                
        elif simpleorcomplex == 2:
            sum = 0
            for i in range(0, n):
                sum += (prb * (a[i] ** eta_sedvar))
         
        return float(sum)
    
    def QwVar(kv, eta_sedvar):
        Qstar = np.arange(0.01, 10, 0.01)
        pdf = (kv ** (kv + 1) / sc.gamma(kv + 1)) * np.exp (-( kv / Qstar))\
         * Qstar ** (-(2 + kv)) * 0.01
        pdf[np.isnan(pdf)] = 0
        pdf = pdf / np.sum(pdf)
        tsteps = len(tarray)
        QwScale = np.random.choice(Qstar, size=tsteps, p=pdf)
        EQw = calc_E(QwScale,len(QwScale),1,eta_sedvar)
        QwScale = QwScale / EQw
        return QwScale
 
    QwScale = QwVar(kv,eta_sedvar)
         
    def QsVar(kv,eta_sedvar):
        EQs = calc_E(QwScale,len(QwScale),2,eta_sedvar)
        ETA_sednorm = 1/EQs
        QsScale = (QwScale ** eta_sedvar) * ETA_sednorm
        return QsScale, ETA_sednorm
             
    QsScale,ETA_sednorm = QsVar(3,eta_sedvar)
    
    Qw = Qwo[location] * QwScale
    Qs = QsScale*dA[location] * 0.001 * 2650

    H = np.zeros(len(Qw))

    for i in range(len(Qw)):
        
        Qwi = (Qw[i])
        
        def Hw(x):
            return 1/n * ((W * x)/(2 * x + W))**(2/3) * ((S) ** 0.5) * W * x - Qwi
        
        H[i] = fsolve(Hw, Htest)
        
    Tb = pw*9.81*H*S
    
    yr2sec = 3.14 * (10 ** 7)
    dt = 1/24
        
    tauNonD = Tb/((ps-pw)*9.81*D) # non dimensional shear stress
    tauNonD[tauNonD<0.000001] = 0.000001
    tau_diff = tauNonD - Tc
    
    qssf = 3.97*((tau_diff)**(1.5)) # Wong and Parker, 2006 re-analysis of MPM
    qsf = np.sqrt((ps-pw)*9.81*D/pw)*D*qssf # sediment transport per unit width
    qsf[qsf < 0] = 0
    Qtf = ps * qsf * W
    Qtf = Qtf*yr2sec*dt
        
    zipped_lists = zip(Qw, Qs,Qtf)
    sorted_pairs = sorted(zipped_lists)

    tuples = zip(*sorted_pairs)
    Qw,Qs,Qtf = [ list(tuple) for tuple in  tuples]
    
    fig, ax2 = plt.subplots()
    ax1 = ax2.twinx()
    ax2.hist(Qw/np.mean(Qw),bins = 30, zorder = 0 , color = 'gray', edgecolor='k', alpha = 0.2)
    ax1.plot(Qw/np.mean(Qw),Qs/np.mean(Qs),c = 'tan', linewidth = 4, label = 'Sediment Supply' , zorder = 100000)
    ax1.plot(Qw/np.mean(Qw),Qtf/np.mean(Qs),c = 'black', linewidth = 4, label = 'Sediment Transport Capacity')
    ax1.set_xlabel('Norm Water Discharge')
    ax1.set_ylabel('Norm Sediment Supply / Transport')
    ax1.set_ylim([0,75])
    ax2.set_yscale('log')
    ax2.set_ylabel('histogram counts')
    
    Qw1 = Qw/np.mean(Qw)
    Qw1 = Qw1[Qs/np.mean(Qs) < Qtf/np.mean(Qs)]
    Qs1 = Qs/np.mean(Qs)
    Qs1 = Qs1[Qs/np.mean(Qs) < Qtf/np.mean(Qs)]
    Qtf1 = Qtf/np.mean(Qs)
    Qtf1 = Qtf1[Qs/np.mean(Qs) < Qtf/np.mean(Qs)]
    
    Qw2 = Qw/np.mean(Qw)
    Qw2 = Qw2[Qs/np.mean(Qs) > Qtf/np.mean(Qs)]
    Qs2 = Qs/np.mean(Qs)
    Qs2 = Qs2[Qs/np.mean(Qs) > Qtf/np.mean(Qs)]
    Qtf2 = Qtf/np.mean(Qs)
    Qtf2 = Qtf2[Qs/np.mean(Qs) > Qtf/np.mean(Qs)]

    ax1.fill_between(Qw1,Qs1,Qtf1, color='firebrick', alpha=.3)
    ax1.fill_between(Qw2,Qs2,Qtf2, color='b', alpha=.3)
    
    ax1.legend()
    
    return

#OTTERoverviewplotting(modelname,uplift,timeseries_location,totaltime,0)

#OTTERoverviewplotting(modelname,uplift,timeseries_location,totaltime,1000)

#OTTEREqplotting(modelname_eq,uplift,timeseries_location,4000,kv,eta)

#OTTERsingletimeplotting(modelname,totaltime,9000)

#OTTERSadler(totaltime)

#OTTERTrends(totaltime)

#OTTERHiatuses(totaltime)

#OTTERGif(modelname_eq,1000)

#OTTERQw_E(modelname_eq,uplift,timeseries_location,8400,kv,eta)

#OTTERSedSupplyVsTrans(30)

