"""
This is for benchmark testing to see how the program did with the final
classifications.
"""
from astropy.io import fits
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import progressBar as pb
matplotlib.rc('text', usetex=True)
matplotlib.rc('font',size=24)
import sys
import glob

def bench_combine():
    test=False
    benchfiles = np.sort(glob.glob('./bench_files/benchmark*.dat'))
    b_array = np.loadtxt(benchfiles[0],delimiter=',',skiprows=1,dtype=bytes).astype(str)
    if test:
        nf = len(benchfiles) - 1
    else:
        nf = len(benchfiles)
    for i in range(1,nf):
        bench_temp = np.loadtxt(benchfiles[i],delimiter=',',skiprows=1,dtype=bytes).astype(str)
        b_array = np.append(b_array,bench_temp,axis=0)
    bst = np.zeros(len(b_array),dtype=[('CTIME','f8'),('PMF','U15'),
                                    ('AUTO_CLASS','U6'),('REAL_CLASS','U6')])

    bst['CTIME'] = b_array[:,0].astype('f8')
    bst['PMF'] = b_array[:,1]
    bst['AUTO_CLASS'] = b_array[:,2]
    bst['REAL_CLASS'] = b_array[:,3]

    return bst

def bench_describe():
    bstruc = bench_combine()
    num_rec = len(bstruc)

    wlost = np.where((bstruc['AUTO_CLASS']=='NotQSO')&(bstruc['REAL_CLASS']=='QSO'))[0]
    wcontam = np.where((bstruc['AUTO_CLASS']=='QSO')&(bstruc['REAL_CLASS']=='NotQSO'))[0]
    wqsoC = np.where((bstruc['AUTO_CLASS']=='QSO')&(bstruc['REAL_CLASS']=='QSO'))[0]
    wnqC = np.where((bstruc['AUTO_CLASS']=='NotQSO')&(bstruc['REAL_CLASS']=='NotQSO'))[0]
    wqsoT = np.where(bstruc['REAL_CLASS']=='QSO')[0]
    wnqT = np.where(bstruc['REAL_CLASS']=='NotQSO')[0]
    wmeq = np.where(bstruc['AUTO_CLASS']=='QSO')[0]
    wmeNQ = np.where(bstruc['AUTO_CLASS']=='NotQSO')[0]

    num_correct = len(wqsoC) + len(wnqC)
    completeness = (len(wqsoC) / len(wqsoT)) * 100
    contaminants = (len(wcontam) / (len(wcontam) + len(wqsoC))) * 100
    pct_correct = (num_correct / num_rec) * 100
    pct_lost = (len(wlost) / len(wqsoT)) * 100
    print('RQSO: '+str(len(wqsoT)))
    print('AQSO:' + str(len(wmeq)))
    print('RNQ: '+str(len(wnqT)))
    print('ANQ: '+str(len(wmeNQ)))
    print(len(wnqC))

    qso_times = bstruc['CTIME'][wqsoT]
    qsoC_times = bstruc['CTIME'][wqsoC]
    nq_times = bstruc['CTIME'][wnqT]
    nqC_times = bstruc['CTIME'][wnqC]
    lost_times = bstruc['CTIME'][wlost]
    contam_times = bstruc['CTIME'][wcontam]
    qsoA_times = bstruc['CTIME'][wmeq]
    nqA_times = bstruc['CTIME'][wmeNQ]

    qCmean,qCmed,qCstd = np.mean(qso_times),np.median(qso_times),np.std(qso_times)
    nqCmean,nqCmed,nqCstd = np.mean(nq_times),np.median(nq_times),np.std(nq_times)

    qbench_str = 'Mean: {0:6.3f} s | Median: {1:6.3f} s | StdDev: {2:6.3f} s'.format(qCmean,qCmed,qCstd)
    nqbench_str = 'Mean: {0:6.3f} s | Median: {1:6.3f} s | StdDev: {2:6.3f} s'.format(nqCmean,nqCmed,nqCstd)
    cc_str = 'Completeness: {0:6.2f}%  | Contam: {1:6.2f}%  | Correct: {2:6.2f}%'.format(completeness,contaminants,pct_correct)
    print('\n')
    print('  QSOs  ' + qbench_str)
    print('NotQSOs ' + nqbench_str)
    print(cc_str)
    print('Percent Lost: {0:6.2f}%'.format(pct_lost))
    print('\n')
    print('Number of Objects: {}'.format(num_rec))
    print('Number of QSOs: {}'.format(len(wqsoT)))
    print('Number of Contaminants: {}'.format(len(wcontam)))
    print('Number of QSOs IDed: {}'.format(len(wmeq)))
    print('\n')


    fig1,ax1 = plt.subplots(nrows=2,ncols=2,sharex='col',sharey='row',figsize=(15,12))
    ax1[0,0].hist(qso_times,bins=100,histtype='stepfilled',color='black',label='R-QSOs',alpha=0.7)
    ax1[0,0].hist(nq_times,bins=100,histtype='stepfilled',color='red',label='R-Non QSOs',linewidth=1.2,alpha=0.6)
    ax1[0,1].hist(qsoC_times,bins=100,histtype='stepfilled',color='black',label='AR-QSOs',alpha=0.7)
    ax1[0,1].hist(nqC_times,bins=100,histtype='stepfilled',color='red',label='AR-Non QSOs',linewidth=1.2,alpha=0.6)
    ax1[1,0].hist(qsoA_times,bins=100,histtype='stepfilled',color='black',label='A-QSOs',alpha=0.7)
    ax1[1,0].hist(nqA_times,bins=100,histtype='stepfilled',color='red',label='A-Non QSOs',linewidth=1.2,alpha=0.6)
    ax1[1,1].hist(lost_times,bins=20,histtype='stepfilled',color='black',label='Lost QSOs',alpha=0.7)
    ax1[1,1].hist(contam_times,bins=5,histtype='stepfilled',color='red',label='Contaminants',linewidth=1.2,alpha=0.6)
    ax1[1,0].set_xlabel('Computation Time (s)')
    ax1[1,1].set_xlabel('Computation Time (s)')
    ax1[0,0].set_ylabel('Number of Objects')
    ax1[1,0].set_ylabel('Number of Objects')
    for i in range(2):
        for j in range(2):
            ax1[i,j].tick_params(axis='both',direction='in')
            ax1[i,j].tick_params(axis='both',which='minor',direction='in')
            ax1[i,j].tick_params(top=True,right=True)
            ax1[i,j].tick_params(which='minor',top=True,right=True)
            ax1[i,j].legend(loc='upper right')
            if i == 1:
                ax1[i,j].xaxis.set_major_locator(ticker.MultipleLocator(2))
                ax1[i,j].xaxis.set_minor_locator(ticker.MultipleLocator(1))
            if j == 0:
                ax1[i,j].yaxis.set_major_locator(ticker.MultipleLocator(250))
                ax1[i,j].yaxis.set_minor_locator(ticker.MultipleLocator(50))
    fig1.tight_layout()
    fig1.subplots_adjust(wspace=0,hspace=0)

    fig2,ax2 = plt.subplots(figsize=(5,5))
    ax2.hist(qsoC_times,bins=100,histtype='stepfilled',color='black',label='AR-QSOs',alpha=0.7)
    ax2.hist(lost_times,bins=100,histtype='stepfilled',color='red',label='Lost QSOs',alpha=0.7)
    ax2.set_xlabel('Computation Time (s)')
    ax2.set_ylabel('Number of Objects')
    ax2.tick_params(axis='both',direction='in')
    ax2.tick_params(axis='both',which='minor',direction='in')
    ax2.tick_params(top=True,right=True)
    ax2.tick_params(which='minor',top=True,right=True)
    ax2.legend(loc='upper right')
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(2))
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(250))
    ax2.yaxis.set_minor_locator(ticker.MultipleLocator(50))
    fig2.tight_layout()

    #plt.show()
    #fig1.savefig('4xhist2.png')
    #fig2.savefig('singlehist.png')

bench_describe()
