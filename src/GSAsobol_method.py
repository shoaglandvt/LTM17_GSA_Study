""" THIS SCRIPT RUNS THE SOBOL METHOD FOR VARIABLES LISTED
IN SVARDICT BELOW. THE VERSION OF THE LTM MODEL CODE  IS
IMPORTED BELOW.  THE SOBOL INDEX CALCS ARE
MADE FOR ALL SIMULATIONS THAT DO NOT EXPERIENCE DROWNING
OF THE BARRIER WIDTH, BARRIER HEIGHT, OR MARSH. LAGOON FILLING
SIMULATIONS ARE INCLUDED IN THE SOBOL CALCS."""

from functools import partial
import argparse
import pickle
from numpy import *
from tools.misc_functions import *
from concurrent.futures import ProcessPoolExecutor
import time
import copy
import os

# This import statement imports the LTM model
from model.LTM17 import *

#########################################################
### START OF SCRIPT #####################################
#########################################################
#########################################################


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "numyrs", help="The number of years in the simulation.", type=int)
    parser.add_argument(
        "ts", help="The simulation timestep.", type=float)
    parser.add_argument(
        "samples", help="Number of samples.", type=int)
    parser.add_argument(
        "save_text", help="What to name the output .pkl file", type=str)

    args = parser.parse_args()


#########################################################
### USER INPUT #########################################
#########################################################
#########################################################
    print('Welcome! Preparing to run Sobol Method...')
    print(' ')
    statevarextractlist = ['XM1', 'XM2', 'XM3', 'DF', 'DM',
                           'Height', 'XT', 'XS', 'XB', 'breakreturn', 'codereturn']
    rvarextractlist = ['DFr', 'DMr', 'BFr', 'BM1r', 'BM2r', 'Heightr', 'Widthr', 'Transgressionr', 'Volr', 'dHeightr', 'dWidthr',
                       'dBM1r', 'dBFr', 'dDFr', 'breakreturn', 'codereturn']
    rvarcalclist = ['DFr', 'DMr', 'BFr', 'BM1r', 'BM2r', 'Heightr', 'Widthr', 'Transgressionr', 'Volr', 'dHeightr', 'dWidthr',
                    'dBM1r', 'dBFr', 'dDFr']

    extractlist = statevarextractlist

    svardict = {}
    svardict['B'] = np.array([0.0001, 0.005])
    svardict['Dt'] = np.array([5, 15])
    svardict['We'] = np.array([100, 600])
    svardict['He'] = np.array([0.5, 4.])
    svardict['Ae'] = np.array([0.005, 0.025])
    svardict['bm1c'] = np.array([50, 500])
    svardict['Qow_max'] = np.array([1, 100])
    svardict['zdot'] = np.array([0.003, 0.02])
    svardict['K'] = np.array([100, 10000])
    svardict['wind'] = np.array([5., 10.])
    svardict['rng'] = np.array([0.7, 2.8])
    svardict['ws'] = np.array([0.05e-3*(365.*24.*3600.), 0.5e-3*(365.*24.*3600.)])
    svardict['tcr'] = np.array([0.05, 0.2])
    svardict['Co'] = np.array([30e-3, 200e-3])
    svardict['Bpeak'] = np.array([1.5, 3.5])

    svarlist = list(svardict.keys())
    k = len(svardict)  # number of variables
    N = args.samples  # number of samples

#########################################################
### MODEL SETUP #######################################
#########################################################
#########################################################
    np.random.seed(14)
    ai = np.random.random([N, k])
    bi = np.random.random([N, k])
    A = np.zeros([N, k])
    B = np.zeros([N, k])

    for i in range(N):
        for j in range(k):
            A[i, j] = svardict[svarlist[j]][0] + ai[i, j] * \
                (svardict[svarlist[j]][1]-svardict[svarlist[j]][0])
            B[i, j] = svardict[svarlist[j]][0] + bi[i, j] * \
                (svardict[svarlist[j]][1]-svardict[svarlist[j]][0])

#########################################################
### RUN A/B MATRICES #######################################
#########################################################
#########################################################
    tic = time.perf_counter()
    print('Running Simulations on A and B Matrices...')

    Aout = []
    count = 1
    with ProcessPoolExecutor() as executor:
        func = partial(model, args.numyrs, args.ts)
        results = executor.map(func, A)
        for result in results:
            if (count == N/4) or (count == N/2) or (count == N*3/4) or (count == N):
                print('..... run %i of %i complete on A' % (count, N))
            count += 1
            Aout.append(result)

    Bout = []
    count = 1
    with ProcessPoolExecutor() as executor:
        func = partial(model, args.numyrs, args.ts)
        results = executor.map(func, B)
        for result in results:
            if (count == N/4) or (count == N/2) or (count == N*3/4) or (count == N):
                print('..... run %i of %i complete on B' % (count, N))
            count += 1
            Bout.append(result)

    Adict = {}
    Bdict = {}
    for i in range(len(extractlist)):
        dum1 = []
        dum2 = []
        for j in range(N):
            dum1.append(Aout[j][i])
            dum2.append(Bout[j][i])
        Adict[extractlist[i]] = dum1
        Bdict[extractlist[i]] = dum2

    toc = time.perf_counter()
    print('Done. (%s)' % displaytimeshort(tic, toc))


#########################################################
### RUN C MATRIX #########################
#########################################################
#########################################################
    print('Running Simulations on C Matrices...')
    tic = time.perf_counter()

    Cdict = {}
    for fi in range(k):
        C = np.zeros([N, k])
        for i in range(k):
            if i == fi:
                C[:, i] = A[:, i]
            else:
                C[:, i] = B[:, i]

        Cout = []

        count = 1
        with ProcessPoolExecutor() as executor:
            func = partial(model, args.numyrs, args.ts)
            results = executor.map(func, C)
            for result in results:
                if (count == np.round(N/4)) or (count == np.round(N/2)) or (count == np.round(N*3/4)) or (count == N):
                    print('..... run %i of %i complete on C' % (count, N))
                count += 1
                Cout.append(result)

        Csubdict = {}
        for i in range(len(extractlist)):
            dum1 = []
            for j in range(N):
                dum1.append(Cout[j][i])
            Csubdict[extractlist[i]] = dum1

        Cdict[svarlist[fi]] = Csubdict
        print('..... Variable %i of %i complete' % (fi+1, k))

    toc = time.perf_counter()
    print('Done. (%s)' % displaytimeshort(tic, toc))

########################################################
## SAVE MATRICES ###########################
########################################################
########################################################
    print('Saving Matrices...', end=" ")
    tic = time.perf_counter()

    filename = open("../output/ms2/Adict%s.pkl" % args.save_text, "wb")
    pickle.dump(Adict, filename)
    filename.close()

    filename = open("../output/ms2/Bdict%s.pkl" % args.save_text, "wb")
    pickle.dump(Bdict, filename)
    filename.close()

    filename = open("../output/ms2/Cdict%s.pkl" % args.save_text, "wb")
    pickle.dump(Cdict, filename)
    filename.close()

    toc = time.perf_counter()
    print('Done. (%s)' % displaytimeshort(tic, toc))

    print('End Program.')
