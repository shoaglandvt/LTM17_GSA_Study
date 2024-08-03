import numpy as np


def wavetau(fetch, wind, Df):
    # calculates bed shear stress (tw) based off of fetch (bf for this model),
    # wind speed (wind), and tidally avg'd lagoon depth (Df)
    # (note: uses linear wave theory)

    [Hs, Tp] = YeV(fetch, wind, Df)  # calculate significant wave height and wave period
    if Hs == 0:
        tw = 0
    else:
        kk = wavek(1/Tp, Df)  # compute wave number
        Um = (np.pi*Hs/Tp/np.sinh(kk*Df))
        aw = Tp*Um/np.pi
        ko = 0.001
        fw = 0.4*(aw/ko)**(-0.75)
        tw = 1/2*1020*fw*Um**2
    return np.round(tw, 4)


def YeV(fetch, wind, D):
    # this function computes significant wave height (Hs) and wave period (Tp)
    # based off of fetch (bf in this model), wind speed, and tidally avg'd
    # lagoon depth (Df); uses semi-empirical equations from Young and
    # Verhagen (1996)
    if (D < 0.01) or (fetch < 25):
        Hs = 0
        Tp = 0
    else:
        g = 9.81  # gravitational acceleration [m/s^2]
        delta = D*g/wind**2
        chi = fetch*g/wind**2
        epsilon = 3.64*10**(-3)*(np.tanh(0.493*delta**0.75)*np.tanh(3.13*10 **
                                                                    (-3)*chi**0.57/np.tanh(0.493*delta**0.75)))**1.74
        ni = 0.133*(np.tanh(0.331*delta**1.01)*np.tanh(5.215*10**(-4)
                                                       * chi**0.73/np.tanh(0.331*delta**1.01)))**(-0.37)
        Hs = 4*np.sqrt(wind**4*epsilon/g**2)
        Tp = wind/ni/g
    return [Hs, Tp]


def wavek(F, H):
    # computes wave number (K) based on wave frequency (F) and depth (H)

    # % function K=wavek(F,H);
    # % where K is wavenumber (rad/m)
    # %       F is frequency (Hz)
    # %       H is depth (m)
    # %
    # % Copyright (C) 2001, Lee Gordon, NortekUSA LLC
    # % https://github.com/csdms-contrib/wetland3p/blob/master/wavek.m
    # % http://bbl.ancl.hawaii.edu/wsvn/filedetails.php?repname=BBL+Code+Repository&path=%2Fbbl%2Ftrunk%2Fsrc%2Fmatlab%2Fkilonalu%2Fprocessing%2FDirSpec0%2FWAVEK.M&rev=1088

    g = 9.80171  # this value of gravity seems off??? why not 9.80665?

#     % This routine use an approximate equation, then sharpens the result with
#     % one interpolation. The result is good to around 1 part in 10^-6.

#     % The equation came out of a textbook, but I have long since forgotton
#     % which one. If you know, please tell me! lgordon@nortekusa.com

    e1 = 4*np.pi**2*F**2*H/g  # f4 = omega**2 * h1/g
    e2 = 1+0.6666666*e1 + 0.355555555*e1**2 + 0.1608465608*e1**3 + \
        0.0632098765*e1**4 + 0.0217540484*e1**5 + 0.0065407983*e1**6
    e3 = +e1**2 + e1/e2
    K1 = np.sqrt(e3)/H

#     %compute error as basis for interpolation

    o1 = np.sqrt(g*K1*np.tanh(K1*H))
    e1 = o1**2*H/g
    e2 = 1+0.6666666*e1 + 0.355555555*e1**2 + 0.1608465608*e1**3 + \
        0.0632098765*e1**4 + 0.0217540484*e1**5 + 0.0065407983*e1**6
    e3 = +e1**2 + e1/e2
    K2 = np.sqrt(e3)/H

#     %interpolate
    K = 2*K1-K2
    return np.round(K, 4)


def wavetauBmod(fetch, wind, Dm, Bfrac):
    if (Bfrac == 0):  # if there is no biomass on the marsh

        # calculate bed shear stress on the marsh using marsh depth (Dm)
        [Hs, Tp] = YeV(fetch, wind, Dm)
        if Hs == 0:
            tw = 0
        else:
            kk = wavek(1/Tp, Dm)
            Um = (np.pi*Hs/Tp/np.sinh(kk*Dm))
            aw = Tp*Um/(2*np.pi)  # changed to 2*pi instead of pi to match LTM original code
            ko = 0.001
            fw = 0.4*(aw/ko)**(-0.75)
            tw = 1/2*1020*fw*Um**2

    else:  # if there is biomass
        tw = 0
    return tw


def waveTRNS(amp, wind, fetch, hb):
    depth = hb  # scarp height
    fac = min(1, depth/(2*amp))
    D = (depth+(depth-fac*2*amp))/2
    [Hs, Tp] = YeV(fetch, wind, D)
    if Hs == 0:
        WP = 0
    else:
        kk = wavek(1/Tp, D)
        cg = 2*np.pi/kk/Tp*0.5*(1+2*kk*D/(np.sinh(2*kk*D)))
        WP = cg*9800/16*abs(Hs)**2
    return WP
