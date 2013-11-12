# birds.py
#

import readcenters
import csv
import sys
import os
from lib import tm
import pylab as pyl
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np

# parameters for plotting figures
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['font.size'] = 8.0
mpl.rcParams['figure.subplot.top'] = 0.8
mpl.rcParams['figure.subplot.wspace'] = 0.25
mpl.rcParams['figure.subplot.wspace'] = 0.25

# makes sure folder for file [f] exists
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

# draws circles on plots
def my_circle_scatter(axes, x_array, y_array, radius=0.5, **kwargs):
    for x, y in zip(x_array, y_array):
        circle = plt.Circle((x,y), radius=radius, **kwargs)
        axes.add_patch(circle)
    return True

# returns distance between pos1 and pos2
def GetDist(pos1, pos2):
    return ((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)**0.5

# plots nearest neighbor distribution, with statistics
def GetNN(curr_id, n_s, n_c, r_c, coords):

    print 'getting nearest neighbor information...'
    
    plt.figure(figsize=(12.0,4.0))

    # do species first
    for i in range(0,n_s):
        NN_sp = []
        for j in range(0, n_c[i]):
            curr_nbr = []
            pos1 = coords[i][j]
            for l in range(0 , n_c[i]):
                if (l == j): continue
                pos2 = coords[i][l]
                dist = GetDist(pos1,pos2)
                curr_nbr.append(dist)
            tmp_a = np.array(curr_nbr)
            min = np.min(tmp_a)
            NN_sp.append(min)

        axes_indiv = plt.subplot(1, n_s + 1, i+1)
        axes_indiv.grid()
        plt.title('species {:d}'.format(i))
        plt.hist(NN_sp,normed=1, bins=20)
        
        tmp_a = np.array(NN_sp)
        N = len(NN_sp)
        mean = np.mean(tmp_a)
        med = np.median(tmp_a)
        std = np.std(tmp_a)
        stats = 'N = {:d}\n\
                mean = {:0.4f}\n\
                med = {:0.4f}\n\
                std = {:0.4f}'.format(N,mean,med,std)
        plt.text(0.8, 0.8, stats, horizontalalignment='right',
                            verticalalignment='center',
                            transform=axes_indiv.transAxes)

    NN_all = []
    for i in range(0,n_s):
        for k in range(0, n_c[i]):
            curr_nbr = []
            pos1 = coords[i][k]
            for j in range(0,n_s):
                for l in range(0, n_c[j]):
                    if (i==j and k==l): continue
                    pos2 = coords[j][l]
                    dist = GetDist(pos1,pos2)
                    curr_nbr.append(dist)
            tmp_a = np.array(curr_nbr)
            min = np.min(tmp_a)
            NN_all.append(min)
    
    axes_combined = plt.subplot(1, n_s + 1, n_s+1)
    axes_combined.grid()
    plt.title('all species')
    plt.hist(NN_all, normed=1, bins=20)

    tmp_a = np.array(NN_all)
    N = len(NN_all)
    mean = np.mean(tmp_a)
    med = np.median(tmp_a)
    std = np.std(tmp_a)
    stats = 'N = {:d}\n\
    mean = {:0.4f}\n\
    med = {:0.4f}\n\
    std = {:0.4f}'.format(N,mean,med,std)
    plt.text(0.8, 0.8, stats, horizontalalignment='right',
             verticalalignment='center',
             transform=axes_combined.transAxes)

    fname = curr_id + '_NN.eps'
    fout = './dat/' + curr_id + '/Plots/' + fname

    ensure_dir(fout)

    plt.suptitle('nearest neighbor distribution', fontsize=12)
    plt.savefig(fout, bbox_inches=0, dpi = 300)


# plots centers
def PlotCenters(curr_id, n_s, n_c, r_c, coords):
    
    print 'plotting centers...'

    plt.figure(figsize=(12.0,4.0))

    axes_combined = plt.subplot(1, n_s + 1, n_s+1)
    axes_combined.grid()
    plt.axis('scaled')
    plt.axis([0,1,0,1])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('all species')


    jet = cm = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=0, vmax=n_s)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    for i in range(0,n_s):
        
        rad = r_c[i]
        
        axes_indiv = plt.subplot(1, n_s + 1, i+1)
        axes_indiv.grid()

        colorVal = scalarMap.to_rgba(i)
        
        for j in range(n_c[i]):
            x = coords[i][j][0]
            y = coords[i][j][1]
            my_circle_scatter(axes_indiv, [x], [y], radius=rad, alpha=0.5, color=colorVal)
            my_circle_scatter(axes_combined, [x], [y], radius=rad, alpha=0.5, color=colorVal)
        
        plt.axis('scaled')
        plt.axis([0,1,0,1])
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('species {:d}, rad {:0.4f}'.format(i,rad))

    fname = curr_id + '_centers.eps'
    fout = './dat/' + curr_id + '/Plots/' + fname

    ensure_dir(fout)

    plt.suptitle('centers', fontsize=12)
    plt.savefig(fout, bbox_inches=0, dpi = 300)

def PlotG2():
    print 'plotting g2...'

def PlotSk(curr_id, n_s):
    
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

        data = list(csv.reader( open(fname, 'rb') , delimiter = '\t' ))
        k = []
        Sk = []
        for i in range(len(data)):
            el = data[i][0].split(' ')
            k.append( float(el[0]) )
            Sk.append( float(el[1]) )
        plt.loglog(k,Sk)
        
    fname = './dat/' + curr_id + '/Sk/' + curr_id + '_Sk_T.txt'

    axes_combined = plt.subplot(1, n_s + 1, n_s+1)
    axes_combined.grid()
    plt.xlabel('k')
    plt.ylabel('S(k)')
    plt.title('all species')
    
    data = list(csv.reader( open(fname, 'rb') , delimiter = '\t' ))
    k = []
    Sk = []
    for i in range(len(data)):
        el = data[i][0].split(' ')
        k.append( float(el[0]) )
        Sk.append( float(el[1]) )
    plt.loglog(k,Sk)

    fname = curr_id + '_Sk.eps'
    fout = './dat/' + curr_id + '/Plots/' + fname

    ensure_dir(fout)

    plt.suptitle('structure factor, S(k)', fontsize=12)
    plt.savefig(fout, bbox_inches=0, dpi = 300)


def PlotBands(curr_id, res, n_s):
    
    print 'plotting bands...'

    rad = tm.MakeRad(res, n_s)

    fname_out = './dat/' + curr_id + '/Download/' + curr_id + '_out'
    fname_band = './dat/' + curr_id + '/Plots/Bands/' + curr_id + '_bands'
    
    ensure_dir(fname_band)

    for r in rad:

        curr_fout = fname_out
        curr_fband = fname_band
        
        for i in range(0, n_s):
            curr_fout += '_r{:0.4f}'.format(float(r[i]))
            curr_fband += '_r{:0.4f}'.format(float(r[i]))
        curr_fout += '.OUT'
        curr_fband += '.png'
        print 'plotting bands for ' + curr_fout

        tmp_bands = 'tmp_bands.del'
        tmp_bands2 = 'tmp.del'
        cmd0 = 'cat ' + curr_fout + ' | grep -i \"Range\" > ' + tmp_bands2
        cmd1 = 'awk \'{print $4, $10}\' tmp.del > ' + tmp_bands
        os.system(cmd0)
        os.system(cmd1)

        fig = plt.figure(figsize=(10,10))

        ax = fig.gca()
        ax.set_aspect(5)
        ax.xaxis.set_ticks_position('none')
        plt.tick_params(axis='y', which='major', labelsize=8)
        plt.tick_params(axis='x', which='major', labelsize=0)
        plt.tick_params(axis='both', which='minor', labelsize=0)
        plt.setp(ax.get_xticklabels(), visible=False)

        with open(tmp_bands,'r') as fin_bands:
            bands = [band.split() for band in fin_bands]
            numbands_f = len(bands)
            print 'num bands...' + str(numbands_f)
            
            prv_lo = 0.0
            prv_hi = 0.0
            
            for band in bands:
                cur_lo = float(band[0])
                cur_hi = float(band[1])
                
                #print '(lo, hi) = ({:0.4f}, {:0.4f})'.format(cur_lo, cur_hi)
                
                if(cur_lo > prv_hi):
                    plt.vlines(0, prv_hi, cur_lo, color='r', linewidths = 4)
                
                prv_lo = cur_lo
                prv_hi = cur_hi

        cmd0 = 'rm ' + tmp_bands
        cmd1 = 'rm ' + tmp_bands2
        
        os.system(cmd0)
        os.system(cmd1)

        plt.axis([-0.25,0.25,0,1])
        #        plt.text(0,5.05,fname_out,size='8',horizontalalignment='center')
        plt.title(r)
        plt.savefig(curr_fband, format='png', dpi=300, bbox_inches='tight')

def MakePlots(curr_id):

    [n_s, n_c, r_c, coords] = readcenters.read(curr_id)

    N = 0
    for i in range(0,n_s):
        N += n_c[i]

#    GetNN(curr_id, n_s, n_c, r_c, coords)
#    PlotCenters(curr_id, n_s, n_c, r_c, coords)
#    PlotSk(curr_id, n_s)
#    PlotBands(curr_id, 1, n_s)
#    PlotBands(curr_id, 2, n_s)