# Sk.py
#

import readcenters
import shutil
import csv
import os
import pylab as pyl
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np

from subprocess import Popen, PIPE
from scipy.optimize import curve_fit

# parameters for plotting figures
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['font.size'] = 8.0
mpl.rcParams['figure.subplot.top'] = 0.8
mpl.rcParams['figure.subplot.wspace'] = 0.25
mpl.rcParams['figure.subplot.wspace'] = 0.25
mpl.rcParams['text.usetex'] = True


maxx = 1
maxy = 1

code_Sk         = './bin/calc_Sk.out'
out_Sk           = 'Sk_bin.txt'

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def CalcSk(curr_id):
    
    fout_Sk    = './dat/' + curr_id + '/Sk/' + curr_id + '_Sk_'
    
    ensure_dir(fout_Sk)
    
    [n_s, n_c, r_c, coords] = readcenters.read(curr_id)

    N = 0
    for i in range(0,n_s):
        N += n_c[i]
    
    r_beg = 3
    
    total_in = [maxx, maxy, N]
    
    # calc Sk for each species
    
    tot_x = []
    tot_y = []
    
    for i in range(0,n_s):

        proc = Popen(code_Sk, stdin = PIPE)
        
        inputs = [maxx, maxy, n_c[i]]
        
        x = []
        y = []
        
        for j in range(0,n_c[i]):
            tmpx = coords[i][j][0]
            tmpy = coords[i][j][1]
            x.append(tmpx)
            y.append(tmpy)
            inputs.extend( [ tmpx, tmpy ] )
            total_in.extend( [ tmpx, tmpy ] )
        
        tot_x.extend(x)
        tot_y.extend(y)
        
        proc.communicate('\n '.join(str(input) for input in inputs) )
        proc.wait()
        
        shutil.copyfile(out_Sk, fout_Sk + str(i) + '.txt')

    proc = Popen(code_Sk, stdin = PIPE)
    proc.communicate('\n '.join(str(input) for input in total_in) )
    proc.wait()
    
    shutil.move(out_Sk, fout_Sk + 'T' + '.txt')

def FitSk(k,Sk):
    
    num_pts = len(k)
    
    cutoff_Sk = .7
    cutoff_k = 0
    cutoff_ind = 0
    
    for i in range(0,num_pts):
        if Sk[i] > cutoff_Sk:
            cutoff_ind = i
            cutoff_k = k[i]
            break

    def func(x,a,b):
        return a * x**b
    
    popt, pcov = curve_fit(func, k[0:cutoff_ind], Sk[0:cutoff_ind])

    fit_k = np.linspace(0, cutoff_k, 50)
    fit_Sk = func(fit_k, popt[0], popt[1])
    
    return [popt, pcov, fit_k, fit_Sk]

def PlotSkFit(popt, pcov, fit_k, fit_Sk, axes):

    a = popt[0]
    b = popt[1]
    sig_a = pcov[0][0]
    sig_b = pcov[1][1]
    
    stats = '$\displaystyle S(k) \sim a k^b$\n\
            $a = $ {:0.2f} , $\sigma_a = ${:0.2f}\n\
            $b = $ {:0.2f} , $\sigma_b = ${:0.2f}'.format(a,sig_a,b,sig_b)

    plt.text(0.1, 0.8, stats, horizontalalignment='left',
             verticalalignment='center',
             transform=axes.transAxes)

    plt.loglog(fit_k, fit_Sk)


def PlotSk(curr_id):
    
    [n_s, n_c, r_c, coords] = readcenters.read(curr_id)
    
    print 'plotting Sk...'
    
    plt.figure(figsize=(12.0,4.0))
    
    for i in range(0,n_s):
        
        axes_indiv = plt.subplot(1, n_s + 1, i+1)
        axes_indiv.grid()
        
        plt.xlabel('k')
        plt.ylabel('S(k)')
        plt.title('species {:d}'.format(i))
        #plt.tick_params(axis='y', which='major', labelsize=8)
        #plt.tick_params(axis='x', which='major', labelsize=8)
        #plt.tick_params(axis='both', which='minor', labelsize=4)
        
        fname = './dat/' + curr_id + '/Sk/' + curr_id + '_Sk_' + str(i) + '.txt'
        
        lines = list(csv.reader( open(fname, 'rb') , delimiter = '\t' ))
        k = []
        Sk = []

        for i in range(len(lines)):
            el = lines[i][0].split(' ')
            k.append( float(el[0]) )
            Sk.append( float(el[1]) )

        [popt, pcov, fit_k, fit_Sk] = FitSk(k,Sk)
        PlotSkFit(popt, pcov, fit_k, fit_Sk, axes_indiv)
        plt.loglog(k,Sk)
    
    fname = './dat/' + curr_id + '/Sk/' + curr_id + '_Sk_T.txt'
    
    axes_combined = plt.subplot(1, n_s + 1, n_s+1)
    axes_combined.grid()
    plt.xlabel('k')
    plt.ylabel('S(k)')
    plt.title('all species')
    
    lines = list(csv.reader( open(fname, 'rb') , delimiter = '\t' ))
    k = []
    Sk = []
    for i in range(len(lines)):
        el = lines[i][0].split(' ')
        k.append( float(el[0]) )
        Sk.append( float(el[1]) )

    [popt, pcov, fit_k, fit_Sk] = FitSk(k,Sk)
    PlotSkFit(popt, pcov, fit_k, fit_Sk, axes_combined)
    plt.loglog(k,Sk)
    
    fname = curr_id + '_Sk.eps'
    fout = './dat/' + curr_id + '/Plots/' + fname
    
    ensure_dir(fout)
    
    plt.suptitle('structure factor, S(k)', fontsize=12)
    plt.savefig(fout, bbox_inches=0, dpi = 300)

