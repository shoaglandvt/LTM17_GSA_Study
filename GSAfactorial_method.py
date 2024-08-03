from functools import partial
import argparse
import pickle
from numpy import *
from tools.misc_functions import *
from concurrent.futures import ProcessPoolExecutor
import time
import os
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
        "save_text", help="What to name the output .pkl file", type=str)

    args = parser.parse_args()

#########################################################
#########################################################
#########################################################
#########################################################

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

    k = len(svardict)  # number of variables
    svarlist = list(svardict.keys())
    numlevels = 2
    numvars = k
    numsims = numlevels**numvars
    numEEs = int(numsims/2)

    print('Welcome! Preparing to run %i simulations.' % numsims)
    tic = time.perf_counter()
    testarray = np.zeros([numsims, numvars])
    print('Creating index array...', end=" ")
    for i in range(numvars):
        curvar = numvars-i
        varind = numlevels**(numvars-curvar)

        varcounter = 0
        for j in range(varind, numsims, 2*varind):
            for k in range(varind):
                testarray[j+k, curvar-1] = 1
                varcounter += 1

            if varcounter == numsims/2:
                break
    toc = time.perf_counter()
    print('Done. (%s)' % displaytimeshort(tic, toc))
    print('Assigning values to index array...', end=" ")
    tic = time.perf_counter()
    runarray = np.zeros([numsims, numvars])
    for i in range(numsims):
        for j in range(numvars):
            runarray[i, j] = svardict[svarlist[j]][int(testarray[i, j])]
    toc = time.perf_counter()
    print('Done. (%s)' % displaytimeshort(tic, toc))

#########################################################
#########################################################
#########################################################
#########################################################
    tic = time.perf_counter()
    print('Running Simulations:')
    output = []
    simnum = 0
    simnumlist = [int(numsims*1/10), int(numsims*2/10), int(numsims*3/10), int(numsims*4/10),
                  int(numsims*5/10), int(numsims*6/10), int(numsims*7/10), int(numsims*8/10),
                  int(numsims*9/10), int(numsims)]
    #
    with ProcessPoolExecutor() as executor:
        func = partial(model, args.numyrs, args.ts)
        results = executor.map(func, runarray)
        for result in results:
            output.append(result)
            simnum += 1
            if (simnum in simnumlist):
                print('.......Completed', simnum, 'of', numsims, 'simulations.')
    toc = time.perf_counter()
    print('Simulations Completed (%s)' % displaytimeshort(tic, toc))

    tic = time.perf_counter()
    print('Postprocessing and saving data...', end=" ")
    output = np.array(output, dtype=object)

    rvarlist = ['DFr', 'DM1r', 'BFr', 'BM1r', 'BM2r', 'Heightr', 'Widthr',
                'Transgressionr', 'Vol_arrayr', 'dHeightr', 'dWidthr', 'dBM1r', 'dBFr', 'dDFr']
    out = {}
    for i in range(len(rvarlist)):
        out[rvarlist[i]] = list(output[:, i])

    filename = open("../output/%s.pkl" % args.save_text, "wb")
    pickle.dump(out, filename)
    filename.close()

    toc = time.perf_counter()
    print('Done. (%s)' % displaytimeshort(tic, toc))

    #########################################################
    # CALCULATE EFFECTS ###########################
    #########################################################
    #########################################################

    outkeylist = list(out.keys())
    for var in outkeylist:
        out[var] = np.array(out[var]).T.tolist()

    tic = time.perf_counter()
    print('Calculating Effects ...')

    eed = {}
    unique_list = []
    for rvar in outkeylist:
        eed[rvar] = {}
        temp1 = []
        for i in range(numsims):
            temp1.append(out[rvar][i])
        unique_list.append(temp1)

    tempresults = []
    with ProcessPoolExecutor() as executor:
        func = partial(ee_calc, numEEs, numvars, numsims, svarlist)
        results = executor.map(func, unique_list)
        for result in results:
            tempresults.append(result)

    for i in range(len(outkeylist)):
        eed[outkeylist[i]] = tempresults[i]

    del tempresults

    toc = time.perf_counter()
    print('Done. (%s)' % displaytimeshort(tic, toc))

    print('Saving EE data...', end=" ")
    tic = time.perf_counter()

    filename = open("../output/%s.pkl" % args.save_text, "wb")
    pickle.dump(eed, filename)
    filename.close()

    toc = time.perf_counter()
    print('Done. (%s)' % displaytimeshort(tic, toc))

    print('End program.')
