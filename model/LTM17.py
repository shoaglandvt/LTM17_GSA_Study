""" LTM17 Python Code Used in GSA Study."""

import pickle
from numpy import *
from model.modelsubfunctions.model_subfunctions import *

import copy
import os


def model(num_years, time_step, runarray):
    ##################################################################
    # ASSIGN VARIABLES #####
    ##################################################################

    from barrier_systems.GSAstudy import B, Dt, We, He, Wo, Ho, Ae, zdot, Qow_max, K, bm1c, rhos, rhom, P, ws, lamda, dist, rng, Dmax, Dmin, tcr, wind, Ka, Ke, Co, Bpeak, por, chiref, po, bfo, bm1o, bm2o, dfo

    print_interm = False
    print_warnings = False

    B = runarray[0]
    Dt = runarray[1]
    We = runarray[2]
    He = runarray[3]
    Ae = runarray[4]
    bm1c = runarray[5]
    Qow_max = runarray[6]
    zdot = runarray[7]
    K = runarray[8]
    wind = runarray[9]
    rng = runarray[10]
    ws = runarray[11]
    tcr = runarray[12]
    Co = runarray[13]
    Bpeak = runarray[14]

    ##################################################################
    # INITIAL CONDITIONS #####
    ##################################################################
    Vd_max = We*He  # Maximum deficit volume (m^2/year)
    amp = rng/2  # tidal amplitude, equal to 0.5*tidal range, m
    Dmax = .7167*rng-.0483  # max backbarrier depth that allows veg. growth [m]

    A = 1.0*Ae  # Shoreface slope, barrier initially in equlibrium
    W = Wo  # Barrier width, initially in equlibrium
    H = Ho  # Barrier height, initially in equlibrium
    Ao = A  # initial slope

    xt = 0.  # toe position [m]
    xtv = 0-Dt  # vertical toe position [m]
    xs = Dt/A  # shoreline position, {depth of toe / slope} [m]
    xb = xs+W  # backbarrier position {shoreline position + barrier width} [m]
    xto = xt  # initial toe position
    xso = xs  # initial shoreline position
    xbo = xb  # initial backbarrier position

    Z = 0.  # 1.0*Dt  # sea level [m] {based off of toe depth}
    Zo = 0.  # 1.0*Dt  # initial sea level [m]

    BL = bfo+bm1o+bm2o  # total lagoon width (lagoon + marshes) [m]
    xm1o = xbo+bm1o  # backbarrier marsh position [m]
    xm2o = xbo+bm1o+bfo  # interior marsh position [m]
    xm3o = xm2o+bm2o  # mainland position [m]
    xLo = xm3o+(Zo+amp)/B
    dmo = Dmax/2  # initial marsh depth with respect to MHW (optimal depth for peak biomass)
    bf = bfo
    bm1 = bm1o  # backbarrier marsh width
    bm2 = bm2o  # interior (mainland) marsh width
    xm1 = xm1o  # backbarrier marsh position [m]
    xm2 = xm2o  # interior marsh position [m]
    df = dfo  # lagoon (or mudflat) depth [m]
    dm1 = dmo  # backbarrier marsh depth [m]
    dm2 = dmo  # interior marsh depth [m]
    dm = dmo  # backbarrier marsh depth [m]
    xm3 = xm3o  # position where interior marsh meets mainland slope [m]
    xL = xLo

    ##################################################################
    # COMPUTATIONAL PARAMETERS #####
    ##################################################################
    Tmax = num_years  # simulation time [yrs]
    dt = time_step  # time step [yrs]
    n = int(Tmax/dt)  # calculate number of time steps

    # initilize arrays
    t = linspace(dt, Tmax, n)
    t = insert(t, 0, 0)
    lent = len(t)

    Vol_array = zeros(lent)
    Transgression = zeros(lent)
    Qow_max_array = zeros(lent)
    wind_array = zeros(lent)
    Height = zeros(lent)  # subaerial barrier height [m]
    Width = zeros(lent)  # subaerial barrier width [m]
    dHeight = zeros(lent)
    dWidth = zeros(lent)
    QSF = zeros(lent)
    OverwashFlux = zeros(lent)  # overwash flux [m^3/m/yr]
    XB = zeros(lent)  # backbarrier position [m]
    XT = zeros(lent)  # toe position [m]
    XTV = zeros(lent)  # vertical toe position [m] starts at -Dt
    XS = zeros(lent)  # shoreline position [m]
    SHrate = zeros(lent)  # shoreline position change rate [m/yr]
    XTrate = zeros(lent)  # toe position change rate [m/yr]
    XBrate = zeros(lent)  # backbarrier position change rate [m/yr]

    XM1 = zeros(lent)  # backbarrier marsh position [m]
    XM2 = zeros(lent)  # interior marsh position [m]
    DF = zeros(lent)  # lagoon depth [m]
    dDF = zeros(lent)
    DM = zeros(lent)  # marsh depth [m]

    DM1 = zeros(lent)  # backbarrier marsh depth [m]
    DM2 = zeros(lent)  # interior marsh depth [m]

    BM1 = zeros(lent)  # backbarrier marsh width [m]
    dBM1 = zeros(lent)
    BM2 = zeros(lent)  # interior marsh width [m]
    XL = zeros(lent)  # mainland position [m]
    BF = zeros(lent)  # lagoon width [m]
    dBF = zeros(lent)
    ZZ = zeros(lent)  # sea level [m]
    QOWBM = zeros(lent)  # overwash flux into marsh [m^3/m/yr]
    QOWBL = zeros(lent)  # overwash flux into lagoon [m^3/m/yr]
    QOWB = zeros(lent)  # total overwash flux [m^3/m/yr]
    QOWH = zeros(lent)  # total overwash flux [m^3/m/yr]
    WO = zeros(lent)  # initial barrier width

    XBDOT = zeros(lent)  # backbarrier position change rate [m/yr]
    DMDOT = zeros(lent)  # lagoon depth change rate [m/yr]
    XM3 = zeros(lent)
    CR = zeros(lent)

    ZDOT_A = zeros(lent)

    # time markers
    tdrown_W = 0
    tdrown_H = 0
    t_Lfill = 0
    t_vmc = 0
    tbbm = 0
    tintm = 0
    breakreturn = 0
    codereturn = 0

    # assign initial non-zero values to arrays
    Height[0] = Ho
    Width[0] = Wo
    XB[0] = xb
    XT[0] = xt
    XTV[0] = xtv
    XS[0] = xs
    XM1[0] = xm1
    XM2[0] = xm2
    XM3[0] = xm3
    DF[0] = df
    DM[0] = dmo
    DM1[0] = dm1
    DM2[0] = dm2
    BM1[0] = bm1
    BM2[0] = bm2
    XL[0] = xL
    BF[0] = bf
    ZZ[0] = Z
    Vol_array[0] = Ho*Wo

    ##################################################################
    # MAIN CODE #####
    ##################################################################
    for i in range(1, n+1):

        if print_interm is True:
            print('####### Time Step ', i, '#########')

        if breakreturn == 0:

            A = Dt/(xs-xt)  # calc new shoreface slope value
            Z = Z+zdot*dt  # calc new sea level according to SLR rate zdot

        ##################################################################
        # DEFICIT VOLUME #####
        ##################################################################

            Phi = min(1.0, bm1/bm1c)
            EDm = max(0, dm1-amp)  # depth of backbarrier marsh
            EDl = max(0, df-amp)  # depth of lagoon
            DD = Phi*EDm+EDl*(1.0-Phi)  # backbarrier depth for vol. considerations

            # deficit volume calculation
            Vd_H = max((He-H)*W, 0.0)  # top of barrier deficit
            Vd_B = max((We-W)*(H+DD), 0.0)  # back of barrier deficit
            Vd = Vd_H+Vd_B  # total deficit volume

        ##################################################################
        # OVERWASH FORMULATION #####
        ##################################################################
            if Vd < Vd_max:
                Qow_H = Qow_max*Vd_H/Vd_max
                Qow_B = Qow_max*Vd_B/Vd_max
            else:
                Qow_H = Qow_max*Vd_H/Vd
                Qow_B = Qow_max*Vd_B/Vd

            Qow = Qow_H+Qow_B  # calculate total overwash
            Qow_B_m = Qow_B*Phi  # distribute part of barrier OW to marsh
            Qow_B_l = Qow_B*(1.0-Phi)  # distribute part of barrier OW to lagoon

        ##################################################################
        # BARRIER EQUATIONS #####
        ##################################################################
            Qsf = K*(Ae-A)
            Hdot = Qow_H/W-zdot  # barrier height change rate [m/yr]
            xbdot = Qow_B_m/(H+EDm)  # backbarrier position change rate [m/yr]
            xsdot = 2.0*Qow/(Dt+2.0*H)-4.0*Qsf*(H+Dt) / \
                (2.0*H+Dt)**2.0  # shoreline pos. change [m/yr]
            xtdot = 4. * Qsf * ((Dt + H)/(Dt*(2.0*H+Dt))) + 2. * zdot / A  # toe pos. change [m/yr]

            H = H+Hdot*dt  # calculate new barrier height [m]
            if H <= 0.0:  # if height is less than 0, mark timestep as height drowned
                H = 0.0
                if tdrown_H == 0:
                    tdrown_H = i*dt
                    if print_warnings is True:
                        print('WARNING: Barrier Height drowned @ t =',
                              i*dt, 'yrs.')
                    breakreturn = 1

            xb = xb+xbdot*dt  # new backbarrier position [m]
            xs = xs+xsdot*dt  # new shoreline position [m]
            xt = xt+xtdot*dt  # new toe position [m]

            A = Dt/(xs-xt)  # new shoreface slope [m/m]
            W = xb-xs  # new barrier width [m]
            if W <= 0.0:  # if width is less than 0, mark timestep as width drowned
                W = 0.0
                if tdrown_W == 0:
                    tdrown_W = i*dt
                    if print_warnings is True:
                        print('WARNING: Barrier Width drowned @ t =',
                              i*dt, 'yrs.')
                    breakreturn = 2

            if print_interm is True:
                print('- - - - OW Calculations - - - -')
                print('Qow: ', Qow)
                print('Qow_H: ', Qow_H)
                print('Qow_B: ', Qow_B)
                print('Qow_B_m: ', Qow_B_m)
                print('Qow_B_l: ', Qow_B_l)
                print('EDm: ', EDm)

        ##################################################################
        # VARIABLE STORAGE #####
        ##################################################################
            Height[i] = H
            Width[i] = W
            dHeight[i] = Ho-H
            dWidth[i] = Wo-W
            XB[i] = xb
            XT[i] = xt
            XTV[i] = XTV[i-1] + zdot*dt
            XS[i] = xs
            Transgression[i] = XS[i]-XS[0]
            ZZ[i] = Z
            QSF[i] = Qsf
            OverwashFlux[i] = Qow
            QOWBM[i] = Qow_B_m
            QOWBL[i] = Qow_B_l
            QOWB[i] = Qow_B
            QOWH[i] = Qow_H
            XBrate[i] = xbdot
            XTrate[i] = xtdot
            ZDOT_A[i] = zdot
            Qow_max_array[i] = Qow_max
            wind_array[i] = wind
            Vol_array[i] = H*W

        ##################################################################
        # BACKBARRIER EQUATIONS #####
        ##################################################################
            AA = .25*(Dmax-Dmin)**2.0
            BMax = Bpeak*(Dmax-dm)*(dm-Dmin)/AA  # salt marsh biomass [kg/m^2]
            if (BMax <= 1e-3):
                BMax = 0.0

            Bfrac = (BMax/Bpeak)
            nuGp = .0138
            AMC = (180.)*BMax*(nuGp)  # accumulated carbon per year
            Rref = AMC*chiref
            FFm = (1./por)*(Rref/po)

            Df = (df+(df-min(rng, df))) / 2.  # lagoon depth - tidally avg'd [m]
            Dm = (dm+(dm-min(rng, dm)))/2.  # marsh depth - tidally avg'd [m]

            tw = wavetau(bf, wind, Df)  # calculate bed shear stress in lagoon
            if Dm > 1e-4:
                twm = wavetauBmod(bf, wind, Dm, Bfrac)
            else:
                twm = 0

            tau = max((tw-tcr)/tcr, 0.0)*lamda
            taum = max((twm-tcr)/tcr, 0.0)*lamda
            Cr = rhos*tau/(1.+tau)  # reference sed. concentration in mudflat [kg/m^3]
            Cm = rhos*taum/(1.+taum)  # reference sed. concentration on marsh [kg/m^3]
            hb = dm+(df-dm)*(1.-exp(-dist*0.1/df))  # scarp height according to profile

            Fc = (Cr-Co)*min(rng, df)/P/rhom

            if df > dm:
                WP = waveTRNS(amp, wind, bf, hb)
                E = (Ke*WP/(hb-dm)-Ka*Cr*ws/rhom)
                Fm = (Cr-Cm)*min(rng, dm)/P/rhom
            else:
                E = 0
                Fm = 0

            if print_interm is True:
                print('- - - - Calculations - - - -')
                print('Bfrac: ', Bfrac)
                print('tw: ', tw)
                print('tau: ', tau)
                print('taum: ', taum)
                print('Cr: ', Cr)
                print('Cm: ', Cm)
                print('hb: ', hb)
                print('WP: ', WP)
                print('Fc: ', Fc)
                print('E: ', E)
                print('Fm: ', Fm)

            dmdot = -Fm-FFm+zdot  # backbarrier marsh depth rate of change [m/yr]
            xm1dot = -E+Qow_B_l/(df-dm)
            xm2dot = E
            xm3dot = (zdot-dmdot)/B
            dfdot = -E*(df-dm1)/bf-E*(df-dm2)/bf+Fm*(bm1+bm2)/bf + \
                Fc+zdot

            xm1 = xm1+xm1dot*dt  # backbarrier marsh position [m]
            xm2 = xm2+xm2dot*dt  # interior marsh position [m]
            xm3 = xm3+xm3dot*dt  # interior marsh/mainland intersect position [m]
            xL = xLo+(Z-Zo)/B  # mainland shoreline position [m]

            dm = dm+dmdot*dt
            df = df+dfdot*dt

            if print_interm is True:
                print('xbdot: ', xbdot)
                print('xm1dot: ', xm1dot)

            # WARNINGS
            if xm1 <= xb:
                xm1 = xb
                if tbbm == 0:
                    tbbm = i*dt
                    if print_warnings is True:
                        print('WARNING: BB Marsh to Zero @ time',
                              i*dt, 'yrs')

            if xm3 <= xm2:
                xm2 = xm3
                if tintm == 0:
                    tintm = i*dt
                    if print_warnings is True:
                        print('WARNING: Interior Marsh to Zero @ time',
                              i*dt, 'yrs')

            if dm > Dmax:
                if t_vmc == 0:
                    t_vmc = i*dt
                    if print_warnings is True:
                        print('WARNING: Vertical Marsh Collapse @ time',
                              i*dt, 'yrs')
                    breakreturn = 3

            if df <= dm:
                if t_Lfill == 0:
                    if print_warnings is True:
                        print('WARNING: Lagoon filled @ time',
                              i*dt, 'yrs')
                    t_Lfill = i*dt
                    codereturn = 4
                xm1 = xb
                xm2 = xm3

            dm1 = dm  # backbarrier marsh depth beomes marsh depth
            dm2 = dm  # interior marsh depth becomes marsh depth

            bm1 = xm1-xb  # calculate backbarrier marsh width
            bm2 = xm3-xm2  # calculate interior marsh width
            bf = xm2-xm1  # calculate lagoon width

            if print_interm is True:
                print('- - - - Final Values - - - -')
                print('H: ', H)
                print('W: ', W)
                print('xb: ', xb)
                print('xs: ', xs)
                print('xt: ', xt)
                print('xm1: ', xm1)
                print('xm2: ', xm2)
                print('xm3: ', xm3)
                print('df: ', df)
                print('dm1: ', dm1)
                print('dm2: ', dm2)
                print('bm1: ', bm1)
                print('bm2: ', bm2)
                print('xL: ', xL)
                print('bf: ', bf)
                print('Z: ', Z)
                print('zdot: ', zdot)
                print('Cr: ', Cr)
                print('Dt: ', Dt)
                print('A: ', A)

        ##################################################################
        # VARIABLE STORAGE #####
        ##################################################################

            XM1[i] = xm1  # backbarrier marsh position [m]
            XM2[i] = xm2  # interior marsh position [m]
            XM3[i] = xm3  # interior marsh/mainland slope intersection [m]
            DF[i] = df  # lagoon depth [m]
            dDF[i] = dfo - df
            DM[i] = dm  # marsh depth [m]
            DM1[i] = dm1  # backbarrier marsh depth [m]
            DM2[i] = dm2  # interior marsh depth [m]
            BF[i] = bf  # lagoon width [m]
            dBF[i] = bfo - bf
            BM1[i] = bm1  # backbarrier marsh width [m]
            dBM1[i] = bm1o - bm1
            BM2[i] = bm2  # interior marsh width [m]
            XL[i] = xL  # mainland shoreline position [m]
            CR[i] = Cr  # reference sed. concentration in mudflat [kg/m^3]
            DMDOT[i] = dmdot  # backbarrier marsh depth rate of change [m/yr]
        else:
            Height[i] = Height[i-1]
            Width[i] = Width[i-1]
            dHeight[i] = dHeight[i-1]
            dWidth[i] = dWidth[i-1]
            XB[i] = XB[i-1]
            XT[i] = XT[i-1]
            XTV[i] = XTV[i-1]
            XS[i] = XS[i-1]
            Transgression[i] = Transgression[i-1]
            ZZ[i] = ZZ[i-1] + zdot*dt
            QSF[i] = 0
            OverwashFlux[i] = 0
            QOWBM[i] = 0
            QOWBL[i] = 0
            QOWB[i] = 0
            QOWH[i] = 0
            XBrate[i] = 0
            XTrate[i] = 0
            ZDOT_A[i] = zdot
            Qow_max_array[i] = 0
            Vol_array[i] = 0

            XM1[i] = XM1[i-1]  # backbarrier marsh position [m]
            XM2[i] = XM2[i-1]  # interior marsh position [m]
            XM3[i] = XM3[i-1]  # interior marsh/mainland slope intersection [m]
            DF[i] = DF[i-1]  # lagoon depth [m]
            dDF[i] = dDF[i-1]
            DM[i] = DM[i-1]  # marsh depth [m]
            DM1[i] = DM1[i-1]  # backbarrier marsh depth [m]
            DM2[i] = DM2[i-1]  # interior marsh depth [m]
            BF[i] = BF[i-1]  # lagoon width [m]
            dBF[i] = dBF[i-1]
            BM1[i] = BM1[i-1]  # backbarrier marsh width [m]
            dBM1[i] = dBM1[i-1]
            BM2[i] = BM2[i-1]  # interior marsh width [m]
            XL[i] = XL[i-1]  # mainland shoreline position [m]
            CR[i] = CR[i-1]  # reference sed. concentration in mudflat [kg/m^3]
            DMDOT[i] = 0  # backbarrier marsh depth rate of change [m/yr]

    return [XM1, XM2, XM3, DF, DM, Height, XT, XS, XB, breakreturn, codereturn]
