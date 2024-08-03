# BARRIER INPUT PARAMETERS ##############
############################################
B = 0.002  # Mainland Slope [m/m]
Dt = 10.  # Toe depth [m]
We = 350.  # Equilibrium (or critical) width [m]
Wo = 350.
He = 2.25  # Equilibrium (or critical) height [m]
Ho = 2.25
Ae = 0.015  # Equilibrium shoreface slope
zdot = 0.0115
Qow_max = 50.  # Maximum overwash flux [m^3/m/yr]
K = 5050  # Shoreface Flux constant [m^2/year]

# MARSH INPUT PARAMETERS ################
#############################################
bm1c = 275.  # Critical marsh width [m]
rhos = 1000.  # marsh and mudlfat density [kg/m^3]
rhom = 1000.  # marsh and mudlfat density [kg/m^3]
P = 12.5/(24.*365.)  # tidal period [yr]
ws = 0.275e-3*(365.*24.*3600.)  # settling velocity [m/yr]
lamda = 0.0001  # 0.0001 sediment erodibility parameter
dist = 10  # parameter for marsh erosion [m]
rng = 1.75  # tidal range, double the tidal amplitude [m]
Dmax = .7167*rng-.0483  # maximum depth that allows vegetation growth [m]
Dmin = 0.  # minimum vegetation growth depth [m]
tcr = 0.125  # critical bed shear stress [Pa or N/m^2]
wind = 7.5  # wind speed [m/s]
Ka = 2.  # marsh progradation coeff.
Ke = 0.16  # 0.15 bank erosion coeff.
Co = 115e-3  # (25–300 mg/l)  # external (open ocean) sed. conc. [g/l =kg/m3]

# SALT MARSH GROWTH #####################
############################################
Bpeak = 2.500  # peak biomass [kg/m^2]
por = 1000./2650.  # organic matter porosity = 0.37
chiref = 0.15  # fraction of refractory accumulated carbon
po = 1000.  # dens of organic matter [kg/m^3]

# MARSH-LAGOON CONFIGURATION ##########
#############################################
bfo = 5500.  # lagoon width between marshes
bm1o = 525.  # initial backbarrier marsh width
bm2o = 300.  # initial interior (mainland) marsh width
dfo = 2.  # Initial lagoon depth w/ respect to MHW
