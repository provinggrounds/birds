# birds.py
#

#import pylab as plt
import readcenters
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def my_circle_scatter(axes, x_array, y_array, radius=0.5, **kwargs):
    for x, y in zip(x_array, y_array):
        circle = plt.Circle((x,y), radius=radius, **kwargs)
        axes.add_patch(circle)
    return True

# gets statistics of nearest neighbors
# plots distribution
# also prints boxlength, number, and statistics (avg, stdev)
# also outputs to file the average and stdev

def GetDist(pos1, pos2):
    return ((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)**0.5

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

    plt.suptitle('nearest neighbor distribution')
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

    plt.suptitle('centers')
    plt.savefig(fout, bbox_inches=0, dpi = 300)

def PlotG2():
    print 'plotting g2...'

def PlotSk():
    print 'plotting Sk...'

def PlotBands():
    print 'plotting bands...'

def MakePlots(curr_id):

    [n_s, n_c, r_c, coords] = readcenters.read(curr_id)

    N = 0
    for i in range(0,n_s):
        N += n_c[i]

    GetNN(curr_id, n_s, n_c, r_c, coords)
    PlotCenters(curr_id, n_s, n_c, r_c, coords)