# OTTERPy
This is the repository for model code, sample data, and plotting scripts for the OTTER model for Python. 

The file OTTERPy.py contains the actual model code

The file OTTERPlottingFuncs.py contains a few scripts for plotting model results.

The four sample data files in the repository are results from two long (1.5million model year) runs at different
rock uplift rates which took several days to run on the IU Bloomington Quartz supercomputing cluster. 
This data should be viewed and plotted using the OTTERPlottingFuncs.py script. See that code for information on 
how to select which model to plot.
    
    Glossary of variables - this is pretty close to exhaustive, there are a
    probably a couple of throwaway variables not listed here.
    
    a:                Sternberg's law multiplier
    Ah:               drainage area from Hack's Law (m^2)
    beta:             fraction of sediment transported as bedload from 0 to 1
    D:                grain size at each node (m)
    dA:               change in drainage area from node to node (m^2)
    depth_thresh:     sediment depth at which bedrock erosion can occur
    Do:               initial grain size
    Do:               grain size at headwaters (m)
    dQs:              change in sediment supply from node to node
    dt:               time step (years)
    dx:               space step (m)
    dz_b:             bedrock erosion at each node
    dz_s:             change in sed depth
    dz_b_t:           bedrock erosion saved through time
    dz_b_store:       bedrock erosion stored
    dz_s_store:       sed depth change stored
    eQ:               exponent for scaling Qw from Ah
    eroded_sedsup:    sediment created by bedrock erosion
    eW:               exponent for scaling W from Qw
    F:                fractional bedrock exposure
    g:                gravitational acceleration (M/s^2)
    H:                water depth
    Hc:               Hack's Law constant
    He:               Hack's law Exponent
    initial_sedsup:   Sediment supply in first 10k yr (equal to uplift)
    kf:               bedrock erodibility
    kQ:               constant for scaling Qw from Ah
    kv:               parameter controlling discharge variability
                        - Crave & Davy, 2001. Lower kv is more variable
    kw:               lateral bedrock erosion constant
    kWid:             constant for scaling W from Qw
    last_z:           elevation at channel mouth
    lmbda:            sediment porosity
    n:                manning's roughness coefficient
    onlyonce:         changes from 0 to 1 to make sure the uplift pulse only happens one time
    onoff:            turns sediment on and off entirely
    onoffsedvar:      turns stochastic sediment on and off
    Qs:               sediment supply along the length of the river
    Qs_down:          downstream sediment flux at each node
    Qt:               sediment transport capacity
    Qw:               water discharge
    Qwo:              initial water discharge
    QwScale:          array of scalar multipliers for Qw, from the chosen distribution
    R:                effective density of submerged sediment
    rho_w:            density of water (kg/m^3)
    rho_s:            density of sediment (kg/m^3)
    savesteps:        at what time should we save outputs?
    sed_depth:        depth of seidment at each node (m)
    SedStore:         storing snapshots of sediment thickness through time
    slope:            river gradient (negative usually)
    SlopeStore:       storing snapshots of channel slope through time
    tarray:           array of model time
    tau_b:            basal shear stress
    tau_c:            critical shields criterion for initiation of sediment motion
    time:             model duration (years)
    topo:             channel elevation (bedrock + sediment)
    upinc:            size up uplift pulse, as a multiplier of initial uplift rate
    upinctime:        when should the pulse happen? (years)
    uplift:           uplift field across the entire domain
    uprate:           initial uplift rate, if using constant uplift across domain
    W:                channel width
    WidthStore:       storing snapshots of channel width through time
    Wo:               initial channel width
    yr2sec:           number of seconds in a year
    z:                bedrock channel elevation
    ZStore:           storing snapshots of bedrock elevation through time
