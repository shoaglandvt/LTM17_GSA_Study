import numpy as np


def slr_fill(slr, t, dt):
    tslr = np.linspace(1, 200, 200)
    tslr = np.insert(tslr, 0, 0)
    slradj = np.interp(t, tslr, slr)
    slr_rate = np.diff(slradj)
    slr_rate = np.insert(slr_rate, 0, 0)
    slr_rate = slr_rate/dt

    return slradj, slr_rate


def displaytime(tic, toc):
    if toc-tic > 60:
        print('Code Executed in ', (toc-tic)//60, 'minute(s) and',
              round((toc-tic)-(toc-tic)//60*60, 2), 'seconds')
    else:
        print('Code Executed in ', round(toc-tic, 2), ' seconds')


def displaytimeshort(tic, toc):
    if toc-tic > 60:
        x = '%f minute(s) %f seconds' % ((toc-tic)//60, round((toc-tic)-(toc-tic)//60*60, 2))
    else:
        x = '%f seconds' % round(toc-tic, 2)

    return x


def displayplottime(tic, toc):
    if toc-tic > 60:
        print('Plotting Executed in ', (toc-tic)//60, 'minute(s) and',
              round((toc-tic)-(toc-tic)//60*60, 2), 'seconds')
    else:
        print('Plotting Executed in ', round(toc-tic, 2), ' seconds')


def ee_calc(numEEs, varcount, numsims, svarlist, yvs):
    rdict = {}

    for i in range(len(svarlist)):
        rvar = svarlist[i]
        rarray = np.zeros(numEEs)
        vcicounter = i+1
        vci = varcount - vcicounter
        vcindex = 2**vci

        basearray = []
        EEcount = 0
        i = 0

        while i < numsims:
            for j in range(vcindex):
                # basearray.append(i)
                rarray[EEcount] = np.abs(yvs[i+vcindex] - yvs[i])
                EEcount += 1
                i += 1

            i += vcindex

        if (rvar == 'Widthr') or (rvar == 'Height'):
            rarray = np.nan_to_num(rarray)

        rdict[rvar] = rarray

    return rdict
