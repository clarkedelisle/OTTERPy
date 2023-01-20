"""
    This code support the manuscript titled 'Rethinking variability in bedrock rivers:
    variability of hillslope sediment supply modulates bedrock incision during floods'
    by Clarke Delisle and Brian J. Yanites. 
"""

import sys
import numpy as np
import scipy.special as sc
from scipy.optimize import fsolve
import pandas as pd

kv = float(sys.argv[1])
eta_sedvar = float(sys.argv[2])
uprate = float(sys.argv[3])
beta = float(sys.argv[4])
onoffsedvar = float(sys.argv[5])
time = float(sys.argv[6])
dt = 1/float(sys.argv[7])
dx = float(sys.argv[8])
    
def otterpystochastic(kv,eta_sedvar,uprate,beta,onoffsedvar,time,dt,dx):
        
    def DomainStart(dx, time, dt):
        
        """Initialize space and time domains"""
        
        x = np.arange(5000, 105000, dx, dtype='float')
        time = int(time)
        tarray = np.arange(1, time, dt, dtype='float')
        
        return x, time, tarray
    
    x, time, tarray = DomainStart(dx, time, dt)
    
    def DrainageAreaStart(x):
        
        """Estimate drainage area from Hacks Law. This will be used to
        dimensionalize nondimensional water discharge values."""
        
        Hc = 1
        He = 1.8598
        Ah = Hc * ((x) ** He)
        dA = np.zeros(len(x))
        dA[1:100] = np.diff(Ah)
        
        return dA, Ah
    
    dA, Ah = DrainageAreaStart(x)
    
    def Constants():
        
        """define important constants"""
        
        g = 9.81
        rho_w = 1000
        rho_s = 2650
        kf = -2 * 10 ** (-4)
        kw = 0.03
        lmbda = 0.2
        n = 0.04
        tau_c = 0.04
        yr2sec = 3.14 * (10 ** 7)
        
        return g, rho_w, rho_s, kf, kw, lmbda, n, tau_c, yr2sec
    
    g, rho_w, rho_s, kf, kw, lmbda, n, tau_c, yr2sec = Constants()
    
    def UpliftStart(x, uprate):
        
        """Defining the initial uplift field, the multiplicative increase in uplift
        and the timing for that increase, if desired"""
        
        uplift = np.zeros(len(x), dtype='float') + uprate
        
        return uplift
    
    uplift = UpliftStart(x,uprate)
    
    def QwStart(Ah):
        
        """Long-term average Qw setup"""
        
        kQ = 10 ** (-8)
        eQ = 1
        Qw = kQ * (Ah ** eQ)
        Qwo = Qw
        
        return Qw, Qwo
        
    Qw, Qwo = QwStart(Ah)
    
    def WidthStart(Qw):
       
        """Channel width setup"""
        kWid = 5
        eW = 0.3
        W = kWid * (Qw ** eW)
        
        return W
    
    W = WidthStart(Qw)
    
    def DStart(x):
        
        """Grain size setup - sternberg's law"""
        a = 0.00005
        Do = 0.3
        D = Do * (np.exp(-a*x))
        D[D <= 0.001] = 0.001
        depth_thresh = D * 3
        
        return D, depth_thresh
    
    D, depth_thresh = DStart(x)
    
    def BedrockSetup(uprate):
        
        """Initial bedrock topo setup - I generally start models with an oversteep
        channel profile. It tends to be faster to erode bedrock down to an
        equilibrium profile rather than building topography to equilibrium
        through rock uplift"""
        
        z = (5400, 5031, 4636, 4312, 4039, 3801, 3585, 3386, 3203, 3034, \
        2877, 2729.68901125, 2591.20636389, 2398.20347555, 2337.38758255, 2278.15752649, \
        2055.63552662, 2003.32119626, 1853.28290795, 1805.4614463,  1758.66371036,\
        1580.75558004, 1538.46519751, 1416.37139168, 1377.1541044,  1338.63632907,\
        1191.25407471, 1156.01386186, 1053.72297386, 1020.71322264,  988.22288788,\
        863.28899696,  833.24447839, 745.76217526,  717.44807561,  689.53528337,\
        581.67300318,  555.61712391, 479.5095483,  454.80576469,  430.42067691,\
        335.93360109,  313.04603787, 246.01856164,  224.2025779,   202.64508769,\
        118.99316647,   98.69778219, 58.71886933,  0) 
        z = np.asarray(z)
        
        if uprate < 0.0002:
            zmultiplier = 0.35
        if ((uprate > 0.0002) and (uprate < 0.003)):
            zmultiplier = 0.7
        if uprate > 0.003:
            zmultiplier = 1.3
        
        z = z * zmultiplier
        
        topo = z
        
        return z, topo
    
    z, topo = BedrockSetup(uprate)
    
    sed_depth = np.zeros(len(x), dtype='float')
    Qs_down = np.zeros(len(x), dtype='float')
    slope = np.zeros(len(x), dtype='float')
    eroded_sedsup = np.zeros(len(x), dtype='float') + uplift
    initial_sedsup = np.zeros(len(x), dtype='float') + uplift
    dz_b = np.zeros(len(x), dtype='float')
    dz_b_t = np.zeros((int(time/dt), int(len(x))))
    tau_b = np.zeros(len(x), dtype='float')
    H = np.zeros(len(x), dtype='float') + 0.5
    F = np.zeros(len(x), dtype='float') + 1
    Qt = np.zeros(len(x), dtype='float')
    dz_s = np.zeros(len(x), dtype='float')
    F3 = np.zeros((3,int(len(x))))
    Qtf3 = np.zeros((3,int(len(x))))
    tau3 = np.zeros((3,int(len(x))))
    taubQSQT3 = np.zeros((3,int(len(x))))  
    dz_b_store = dz_b
    dz_s_store = dz_s
    
    nsaves = 1000
    savecounter = 0
    savesteps = time/nsaves
    WidthStore = np.zeros((len(x), nsaves+1), dtype='float')
    SlopeStore = np.zeros((len(x), nsaves+1), dtype='float')
    ZStore = np.zeros((len(x), nsaves+1), dtype='float')
    SedStore = np.zeros((len(x), nsaves+1), dtype='float')
    EStore = np.zeros((len(x), nsaves+1), dtype='float')
    TauStore = np.zeros((len(x), nsaves+1), dtype='float')
    QwStore = np.zeros((len(x),nsaves+1), dtype='float')
    FStore = np.zeros((len(x),nsaves+1), dtype='float')
    WidthStoreEq = np.zeros((len(x),12000), dtype='float')
    SlopeStoreEq = np.zeros((len(x),12000), dtype='float')
    ZStoreEq = np.zeros((len(x),12000), dtype='float')
    SedStoreEq = np.zeros((len(x),12000), dtype='float') 
    EStoreEq = np.zeros((len(x),12000), dtype='float')
    TauStoreEq = np.zeros((len(x),12000), dtype='float')
    QwStoreEq = np.zeros((len(x),12000), dtype='float')
    FStoreEq = np.zeros((len(x),12000), dtype='float')
    
    column_names = ["Width","Slope","Z","Sed_Depth","E","Tau","Qw","F"]

    df = pd.DataFrame(columns = column_names)
    df_Eq = pd.DataFrame(columns = column_names)
            
    def variableQwQs(kv,eta_sedvar):
        
        """
        variable QwQs is the function that will generate distributions of 
        nondimensional water and sediment discharge values, based on the user 
        input parameters. 
        """
    
        def calc_E(a,n,simpleorcomplex,eta_sedvar):
            
            """ calculate the expectation value of a distribution"""
     
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
            
            """ Generate non-dimensional Qw Values"""
            
            Qstar = np.arange(0.01, 10, 0.01)
            pdf = (kv ** (kv + 1) / sc.gamma(kv + 1)) * np.exp (-( kv / Qstar))\
             * Qstar ** (-(2 + kv)) * 0.01 # defining the PDF
            pdf[np.isnan(pdf)] = 0
            pdf = pdf / np.sum(pdf)
            tsteps = len(tarray)
            QwScale = np.random.choice(Qstar, size=tsteps, p=pdf) #drawing from dist
            EQw = calc_E(QwScale,len(QwScale),1,eta_sedvar)
            QwScale = QwScale / EQw
            
            return QwScale
     
        QwScale = QwVar(kv,eta_sedvar)
        
        def QsVar(kv,eta_sedvar):
            
            """Generate non-dimensional Qs Values"""
            
            EQs = calc_E(QwScale,len(QwScale),2,eta_sedvar) # normalization constant
            ETA_sednorm = 1/EQs
            QsScale = (QwScale ** eta_sedvar) * ETA_sednorm # values of Qs*
            return QsScale, ETA_sednorm
             
        QsScale,ETA_sednorm = QsVar(kv,eta_sedvar)
         
        return QwScale,QsScale,ETA_sednorm
 
    def width_test(W, Qw, S, D, QSS):
        
        """ Function that calculates flow depths and bedload sediment
        transport capacities """
      
        # test wide channel assumption to get new H:
        Htest = (n * Qw/W) ** 0.6 * (S) ** (0.3)
    
        #Manning's
        def Hw(x):
            return 1/n * ((W * x)/(2 * x + W))**(2/3) * ((S) ** 0.5) * W * x - Qw
        
        H = fsolve(Hw, Htest) # solve for flow depth
        H[np.isnan(H)] = 0.01
        H[-1] = H[-2]
        tau = rho_w * g * S * (W*H/(2*H+W)) # shear stress

        tauNonD = tau/((rho_s-rho_w)*g*D) # non dimensional shear stress
        tauNonD[tauNonD<0.000001] = 0.000001
        tau_diff = tauNonD - tau_c # excess shear stress
        tau_diff[tau_diff<0] = 0
        
        qssf = 3.97*((tau_diff)**(1.5)) # Wong and Parker, 2006 re-analysis of MPM
        qsf = np.sqrt((rho_s-rho_w)*g*D/rho_w)*D*qssf # sediment transport per unit width
        
        qsf[qsf < 0] = 0
    
        Qtf = rho_s * qsf * W # mass transport capacity kg/s
        
        Qtf[np.isnan(Qtf)] = 0
        
        F = np.zeros(len(x), dtype='float') + 1
        F[Qtf == 0] = 0
        F[Qtf != 0] = 1 - QSS[Qtf != 0]/Qtf[Qtf != 0] # bedrock fraction
        F[F<0] = 0
     
        taubQSQT = tau*F
            
        return F, tau, Qtf, taubQSQT, H
    
    def calculate_new_width(Wo, tau, Qtf, taubQSQT, QSS, Widening_rate, Narrowing_rate):
        
        """ Function that test three scenarios of channel width """
    
        maxtau = np.argmax(tau, axis=0)
        maxtaubQSQT = np.argmax(taubQSQT, axis=0)
        Qtf0 = Qtf[0,:]
        W = np.zeros(len(x))
        QSS = np.nan_to_num(QSS)
        maxtau = np.nan_to_num(maxtau)
        Qtf0 = np.nan_to_num(Qtf0)
        maxtaubQSQT = np.nan_to_num(maxtaubQSQT)
            
        mask1 = np.logical_and(QSS <= 0, maxtau == 2)
        mask2 = np.logical_and(QSS <= 0, maxtau == 1)
        mask3 = np.logical_and(QSS <= 0, maxtau == 0)
        W[mask1] = Wo[mask1] + (dt*Narrowing_rate[0,mask1])
        W[mask2] = Wo[mask2] + (dt*Widening_rate[mask2])
        W[mask3] = Wo[mask3]
        
        mask4 = np.logical_and(Qtf0 >= QSS, maxtaubQSQT == 2)
        mask5 = np.logical_and(Qtf0 >= QSS, maxtaubQSQT == 1)
        mask6 = np.logical_and(Qtf0 >= QSS, maxtaubQSQT == 0)
        W[mask4] = Wo[mask4] + (dt*Narrowing_rate[0,mask4])
        W[mask5] = Wo[mask5] + (dt*Widening_rate[mask5])
        W[mask6] = Wo[mask6]
        
        mask7 = np.logical_and(Qtf0 < QSS, maxtaubQSQT == 2)
        mask8 = np.logical_and(Qtf0 < QSS, maxtaubQSQT == 1)
        mask9 = np.logical_and(Qtf0 < QSS, maxtaubQSQT == 0)
        W[mask7] = Wo[mask7] + (dt*Narrowing_rate[0,mask7])
        W[mask8] = Wo[mask8] + (dt*Widening_rate[mask8])
        W[mask9] = Wo[mask9]
        
        return W
    
    def calc_width(W, H, Qw, slope, EE, DS, QSS, D, kw, dt):
        
        """ Function that calculates the new channel width """
          
        S = -1*slope  # make slope positive
    
        Narrowing_rate_rock = (W / (H/(dt*EE)))
        Narrowing_rate_sed = (W / (H/(dt*DS)))
        Narrowing_rate = np.minimum([Narrowing_rate_rock], [Narrowing_rate_sed])
    
        # make sure narrowing rate is negative
        Narrowing_rate = abs(Narrowing_rate) * -1
            
        # limit narrowing to 1% of current width:
        MinNarrowing = np.array(-0.01*W,ndmin=2)
        Narrowing_rate[Narrowing_rate < MinNarrowing] = MinNarrowing[Narrowing_rate < MinNarrowing]
    
        Widening_rate = (rho_w * g * S * (W*H/(2*H+W))) * kw
    
        # push width to wider and lower to test optimization
        Wo = W
        Wbig = W * 1.05
        Wsmall = W * 0.95
       
        F3[0, :], tau3[0, :], Qtf3[0, :], taubQSQT3[0, :], H = width_test(Wo, Qw, S, D, QSS)
        F3[1, :], tau3[1, :], Qtf3[1, :], taubQSQT3[1, :], H = width_test(Wbig, Qw, S, D, QSS)
        F3[2, :], tau3[2, :], Qtf3[2, :], taubQSQT3[2, :], H = width_test(Wsmall, Qw, S, D, QSS)
    
        #find the new width
        W = calculate_new_width(Wo, tau3, Qtf3, taubQSQT3, QSS, Widening_rate, Narrowing_rate)
        F = F3
        tau = tau3
        Qtf = Qtf3
        taubQSQT = taubQSQT3
        
        return F, W, H, tau, Qtf, taubQSQT
    
    # simulate all values of Qw*, Qs* before starting time loop
    QwScale, QsScale, ETA_sednorm = variableQwQs(kv, eta_sedvar)
    
    for t in range(1, (len(tarray)) - 1):
    
        # determine water discharge
        Qw = Qwo * QwScale[t]
                
        # apply rock uplift to present channel elevation to get new elevation
        z = z + (uplift * dt)
        z[-1] = 0
        sed_depth[-1] = 0
        
        # topography is bedrock elevation plus sediment depth
        topo = z + sed_depth
        
        # calculate channel slope
        slope[0:(len(x)-1)] = (np.diff(topo) / np.diff(x))
    
        # determine sediment input from hillslopes
        
        if onoffsedvar == 0:
            Qs_n = (rho_s * beta * dA * eroded_sedsup * dt) + Qs_down
        else:
           # Qs_n = (rho_s * beta * dA * eroded_sedsup * dt * QsScale[t]) + Qs_down
           Qs_n = (rho_s * beta * dA * eroded_sedsup * dt * QsScale[t]) + Qs_down
            
        EE = (dz_s_store / dt)
        QSS = (dz_b_store / dt)
        DS = (Qs_n / (yr2sec*dt))
        
        # calculate a new channel width via the optimization algorithm
        F, W, H, tau, Qtf, taubQSQT = calc_width(W, H, Qw, slope, QSS, EE, DS, D, kw, dt)  
        Qtf0 = np.amax(Qtf, 0)
        F0 = np.amax(F, 0)
        Qtfreal = np.isreal(Qtf0)
        Qt[Qtfreal == 1] = Qtf0[Qtfreal == 1]
        Qt[Qtfreal == 0] = 0
        
        # pass sediment downstream
        dz_s = (1 - (1 - lmbda))*((dt*Qt*yr2sec)-Qs_n)/(W*rho_s*dx)
        Qs_pass = dt*Qt*yr2sec
        mask = np.logical_and(dz_s > sed_depth, sed_depth > 0)
        dz_s[mask] = sed_depth[mask]
        Qs_down = (rho_s*dz_s*W*dx)+(Qs_pass)
        Qs_down[Qs_down < 0] = 0
        
        # calculate bedrock erosion
        sed_depth = sed_depth-(dz_s)
        sed_depth[sed_depth < 0] = 0
        tau_b = rho_w*g*(-slope)*(W*H/(W+(2*H)))
        F0[F0<0] = 0
        dz_b[sed_depth < depth_thresh] = F0[sed_depth < depth_thresh]*kf*tau_b[sed_depth < depth_thresh] * dt
        dz_b[dz_b > 0] = 0
        dz_b[sed_depth >= depth_thresh] = 0
        dz_b[-1] = dz_b[-2]
        dz_b[0] = dz_b[1]    
        dz_b_t[t, :] = dz_b

        # calculate new topography
        z[sed_depth < depth_thresh] = z[sed_depth < depth_thresh] + dz_b[sed_depth < depth_thresh]
        topo = z + sed_depth
        
        dz_s_store = dz_s
        dz_b_store = dz_b
                 
        # Kind of messy, but this is just saving the important variables which 
        # we want to output and examine. 
        if t == 1:
            WidthStore[:, 0] = W
            SlopeStore[:, 0] = slope * -1
            ZStore[:, 0] = z
            SedStore[:, 0] = sed_depth
            EStore[:, 0] = dz_b
            TauStore[:,0] = tau_b
            QwStore[:,0] = Qw
            FStore[:,0] = F0
        elif (t*dt) % savesteps == 0:
            saveindex = int((t*dt)/savesteps)
            WidthStore[:, saveindex] = W
            SlopeStore[:, saveindex] = slope * -1
            ZStore[:, saveindex] = z
            SedStore[:, saveindex] = sed_depth
            EStore[:, saveindex] = dz_b
            TauStore[:,saveindex] = tau_b
            QwStore[:,saveindex] = Qw
            # I usually leave this in just if I'm testing model parameters, it will print these variables
            # for one node (here I use node 40) through time in the terminal
            print('t = ' + str(t*dt) + ', Qw = ' + str(Qw[40]) + ', H = ' + str(H[40])+ ', Z = ' + str(z[40]) + ', Sed = ' + str(sed_depth[40]))
            FStore[:,saveindex] = F0
        elif t == (len(tarray)) - 2:
            WidthStore[:, -2] = W
            SlopeStore[:, -2] = slope * -1
            ZStore[:, -2] = z
            SedStore[:, -2] = sed_depth
            EStore[:, -2] = dz_b
            TauStore[:,-2] = tau_b
            QwStore[:,-2] = Qw
            FStore[:,-2] = F0
            
        # Saving variables at high - resolution for last 2k years
        if t > (time/dt-12000):
            WidthStoreEq[:, savecounter] = W
            SlopeStoreEq[:, savecounter] = slope * -1
            ZStoreEq[:, savecounter] = z
            SedStoreEq[:, savecounter] = sed_depth
            EStoreEq[:, savecounter] = dz_b
            TauStoreEq[:,savecounter] = tau_b
            QwStoreEq[:,savecounter] = Qw
            FStoreEq[:,savecounter] = F0
            savecounter = savecounter + 1
                        
        # reset some variables to zero before the loop repeats    
        eroded_sedsup = initial_sedsup
        Qs_down = 0
        dz_s[0:-1] = 0
        dz_b[0:-1] = 0
   
    # Once the time loop has completed, we add the important variables to the
    # dataframe and save it with pickle. 
    ModelName = 'uplift' + str(int(uprate*10000)) + '_kv' +str(int(kv*10))\
        + '_eta' +str(int(eta_sedvar*10)) + '_time' + str(time)
        
    # save long-timescale variables
    df = df.append({"Width":WidthStore,"Slope":SlopeStore,"Z":ZStore,"Sed_Depth":SedStore,\
                    "E":EStore,"Tau":TauStore,"Qw":QwStore,"F":FStore},ignore_index = True)

    # save models for 5k years at equilibrium
    df_Eq = df_Eq.append({"Width":WidthStoreEq,"Slope":SlopeStoreEq,"Z":ZStoreEq,\
                        "Sed_Depth":SedStoreEq,"E":EStoreEq,"Tau":TauStoreEq,"Qw":QwStoreEq,\
                            "F":FStoreEq},ignore_index = True)
    
    # save data to a .pkl file
    df.to_pickle(ModelName + ".pkl") 
    df_Eq.to_pickle(ModelName + "Eq.pkl")
    
    return
  
# Here, at the end of the code, we finally call the huge function that we just
# built. The inputs for this function (kv,eta_sedvar,uprate,beta,onoffsedvar,
# time,dt,dx) need to be input by the model user in the command line. 

otterpystochastic(kv,eta_sedvar,uprate,beta,onoffsedvar,time,dt,dx)

# Here is an example of how I would go about actually running this model. Make
# sure that this file is in a directory that you are okay saving new files to, 
# open a terminal window on your Mac and enter
#
#     python OTTERPy.py 0.3 0.5 0.0001 0.3 1 1000 24 2000
#
# this will initiate a model where 
# k = 0.3 (this controls water discharge variability)
# eta = 0.5 (this controls sediment supply variability)
# uplift = 0.1mm/yr 
# beta = 0.3 (30% of sediment is bedload, the rest suspended)
# onoffsedvar = 1 (this turns on dynamic sediment supply, 0 would turn it off)
# time = 1000 years
# dt = 1/24 years
# dx = 2000m
