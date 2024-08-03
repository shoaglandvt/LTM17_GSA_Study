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
    parser.add_argument("numtraj", help="Number of trajectories.", type=int)
    parser.add_argument("numlevels", help="Number of levels (p).", type=int)
    parser.add_argument("numEEs", help="Number of EEs/variable (sets Delta value).", type=int)
    parser.add_argument(
        "save_text", help="What to name the output .pkl file", type=str)
    args = parser.parse_args()

    numyrs = args.numyrs
    ts = args.ts


#########################################################
### USER INPUT #########################################
#########################################################
#########################################################
    print('Welcome! Preparing to run Morris Method...')
    print(' ')

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
    p = args.numlevels  # number of levels
    x = args.numEEs  # number of EEs per variable (sets the delta value)
    r = args.numtraj  # number of trajectories
    rlist = np.linspace(0, r-1, r)

    #########################################################
    ### MODEL SETUP #######################################
    #########################################################
    #########################################################

    levels = np.zeros(p)
    for i in range(1, p):
        levels[i] = levels[i-1]+(1/(p-1))

    '''Note x must be between 2 and (p-1)'''
    if (x < 1) or (x > (p-1)):
        print('Error: x must be between 2 and (p-1)')
        sys.exit()

    numEEs = x*p**(k-1)
    numsims = 2*numEEs
    del_base = 1/(p-1)
    delta = del_base*(p-x)

    print('Number of Variables: ', k)
    print('Number of Levels: ', p)
    print('Delta: ', delta)
    print('Number of EEs in Parameter Space: ', numEEs)
    print('Number of Simulations in Parameter Space: ', numsims)
    print(' ')
    print('Running Trajectories...')
    tic = time.perf_counter()

#########################################################
### DEFINE ELEM EFFECTS DICT #########################
#########################################################
#########################################################
    eed = {}

    model_results_list = ['DFr', 'DM1r', 'BFr', 'BM1r',
                          'BM2r', 'Heightr', 'Widthr', 'Transgressionr', 'Volr', 'dHeightr', 'dWidthr']
    model_results_list_nums = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    for traj in range(r):
        eed[traj] = {}
        for rvar in model_results_list:
            eed[traj][rvar] = {}
            for i in range(k):
                eed[traj][rvar][svarlist[i]] = []

#########################################################
### RUN TRAJECTORIES #########################
#########################################################
#########################################################
    def create_runarray(k, levels, delta, svardict, svarlist):

        #########################################################
        ### CREATE B* MATRIX #########################
        #########################################################

        # randomize starting level (value) for each variable
        xstar = np.zeros([1, k])
        for i in range(k):
            xstar[0, i] = np.random.choice(levels)

        # lower triangular identity matrix
        B = np.zeros([k+1, k])
        for i in range(1, k+1):
            B[i, 0:i] = 1

        # ones matrix
        J = np.ones([k+1, k])

        # random diagonal matrix of +1 or -1
        Dstar = np.zeros([k, k])
        for i in range(k):
            dv = np.random.random([1])
            if dv <= 0.5:
                Dstar[i, i] = -1
            else:
                Dstar[i, i] = 1

        # zeros matrix - each row has single 1 in it, randomly assigned non repeating
        perm = np.random.permutation(k)
        Pstar = np.zeros([k, k])
        for i in range(k):
            Pstar[i, perm[i]] = 1

        Omega = delta/2*(np.matmul((2*B-J), Dstar)+J)

        J1 = np.zeros([k+1, 1])
        for i in range(k+1):
            J1[i, 0] = J[i, 0]

        J1X = np.matmul(J1, xstar)

        J1X_Omega = J1X + Omega

        for i in J1X_Omega:
            for j in range(len(i)):
                if i[j] > 1:
                    for k in range(len(J1X_Omega)):
                        J1X_Omega[k, j] -= delta
                elif i[j] < 0:
                    for k in range(len(J1X_Omega)):
                        J1X_Omega[k, j] += delta

        Bstar = np.matmul(J1X_Omega, Pstar)

        varorder = []
        for i in Pstar:
            for j in range(k):
                if i[j] == 1:
                    varorder.append(j)

        #########################################################
        ### ASSIGN MODEL VALUES TO B* MATRIX ################
        #########################################################

        runarray = np.zeros([k+1, k])
        for i in range(k+1):
            for j in range(k):
                runarray[i, j] = svardict[svarlist[j]][0] + \
                    (svardict[svarlist[j]][1]-svardict[svarlist[j]][0])*Bstar[i, j]

        return runarray, varorder

    def run_trajectory(numyrs, ts, subrunarray):
        numruns = len(subrunarray)
        trajoutput = []
        for i in range(numruns):
            trajoutput.append(model(numyrs, ts, subrunarray[i]))

        return trajoutput

    #########################################################
    ### CREATE RUNARRAY MASTER LIST ####################
    #########################################################
    #########################################################
    runarray_masterlist = []
    varorder_masterlist = []
    for i in range(r):
        runarray_i, varorder_i = create_runarray(k, levels, delta, svardict, svarlist)
        runarray_masterlist.append(runarray_i)
        varorder_masterlist.append(varorder_i)

    #########################################################
    ### RUN SIMULATIONS ###################################
    #########################################################
    #########################################################
    output = []
    count = 1
    with ProcessPoolExecutor() as executor:
        func = partial(run_trajectory, numyrs, ts)
        results = executor.map(func, runarray_masterlist)
        for result in results:
            if (count == r/4) or (count == r/2) or (count == r*3/4) or (count == r):
                print('..... trajectory %i of %i complete' % (count, r))
            count += 1
            output.append(result)

    #########################################################
    ### CALCULATE ELEM EFFECTS ##########################
    #########################################################
    #########################################################

    for traj in range(r):
        for i in range(len(model_results_list)):
            for j in range(k):
                numerator = output[traj][1+j][i] - output[traj][j][i]
                eed[traj][model_results_list[i]][svarlist[varorder_masterlist[traj][j]]
                                                 ].append(numerator)
                # divide by delta value (denominator) in post-processing

    filename = open("../output/eed_%s.pkl" % args.save_text, "wb")
    pickle.dump(eed, filename)
    filename.close()

    filename = open("../output/output_%s.pkl" % args.save_text, "wb")
    pickle.dump(output, filename)
    filename.close()

    toc = time.perf_counter()
    print('End program. (%s)' % displaytimeshort(tic, toc))
