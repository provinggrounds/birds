# birds.py
#

import readcenters
import csv
import sys
import os
from PIL import Image
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

def min(a,b):
    if (a < b):
        return a
    else:
        return b

# returns distance between pos1 and pos2
#     (taking into account PBC, assuming boxlength 1)
def GetDist(pos1, pos2):
    
    boxlen = 1
    
    [x_1, y_1] = pos1
    [x_2, y_2] = pos2
    
    x_distsq = min( (x_1 - x_2)**2 , boxlen**2 - 2 * boxlen * (x_1 - x_2) + (x_1 - x_2)**2 )
    y_distsq = min( (y_1 - y_2)**2 , boxlen**2 - 2 * boxlen * (y_1 - y_2) + (y_1 - y_2)**2 )
    
    return (x_distsq + y_distsq)**0.5



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

        lines = list(csv.reader( open(fname, 'rb') , delimiter = '\t' ))
        k = []
        Sk = []
        for i in range(len(lines)):
            el = lines[i][0].split(' ')
            k.append( float(el[0]) )
            Sk.append( float(el[1]) )
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
    plt.loglog(k,Sk)

    fname = curr_id + '_Sk.eps'
    fout = './dat/' + curr_id + '/Plots/' + fname

    ensure_dir(fout)

    plt.suptitle('structure factor, S(k)', fontsize=12)
    plt.savefig(fout, bbox_inches=0, dpi = 900)

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

def PlotSomeBands(curr_id, res, n_s, subbands):
    
    print 'plotting bands...'
    
    rad = tm.MakeRad(res, n_s)
    
    fname_out = './dat/' + curr_id + '/Download/' + curr_id + '_out'
    fname_band = './dat/' + curr_id + '/Plots/Bands/' + curr_id + '_subbands'
    
    ensure_dir(fname_band)
    
    min = 0.1
    max = 0.6
    dpi_in = 400

    jet = cm = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=0, vmax=len(subbands))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    
    for r in rad:
        
        curr_fout = fname_out
        curr_fband = fname_band
        
        for i in range(0, n_s):
            curr_fout += '_r{:0.4f}'.format(float(r[i]))
            curr_fband += '_r{:0.4f}'.format(float(r[i]))
        curr_fout += '.OUT'
        curr_fband = curr_fband + '_dpi' + str(dpi_in) + '.png'
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
        ax.yaxis.set_ticks_position('none')
        plt.xticks([])
        plt.yticks([0.2, 0.4, 0.6])
        plt.tick_params(axis='y', which='major', labelsize=2)
        plt.tick_params(axis='x', which='major', labelsize=0)
        plt.tick_params(axis='both', which='minor', labelsize=0)
        plt.setp(ax.get_xticklabels(), visible=False)
        
        with open(tmp_bands,'r') as fin_bands:
            bands = [b.split() for b in fin_bands]
            
            for i in range(len(subbands)):
                
                colorVal = scalarMap.to_rgba(i)

                cur_hi = float(bands[int(subbands[i-1])][1])
                nxt_lo = float(bands[int(subbands[i-1])+1][0])
                
                if(nxt_lo > cur_hi):
                    plt.vlines(0, cur_hi, nxt_lo, color=colorVal, linewidths = 100)
                        
        cmd0 = 'rm ' + tmp_bands
        cmd1 = 'rm ' + tmp_bands2
        
        os.system(cmd0)
        os.system(cmd1)
        
        plt.axis([-0.1,0.1,min,max])
        ax.yaxis.labelpad = -1

        ax.text(0.5,0.97,r,
                horizontalalignment='center',
                transform=ax.transAxes)
        
        #        plt.text(0,5.05,fname_out,size='8',horizontalalignment='center')
        plt.savefig(curr_fband, format='png', dpi=dpi_in, bbox_inches='tight')

def CalcBands(curr_id, res, n_s):
    
    print 'calculating bands...'

    rad = tm.MakeRad(res, n_s)

    fname_out = './dat/' + curr_id + '/Download/' + curr_id + '_out'
    fname_band = './dat/' + curr_id + '/' + curr_id + '_bands.csv'

    ensure_dir(fname_band)

    initiate = 0
    
    if not os.path.exists(fname_band):
        lines = [[]]
        for i in range(n_s):
            lines[0].append('r' + str(i))
        numrows = 0
        numcols = 0
        numbands_i = 0
        initiate = 1

    for r in rad:

        cur_fout = fname_out
        for i in range(0, n_s):
            cur_fout += '_r{:0.4f}'.format(float(r[i]))
        cur_fout += '.OUT'
        print 'opening file ' + cur_fout

        print 'adding bands for',
        print r

        if not initiate:
            lines = list(csv.reader( open(fname_band, 'rb') , delimiter = '\t' ))

            lines = [line[0].split(',') for line in lines]
            numrows = len(lines)
            numcols = len(lines[0])
            numbands_i = ((numcols-2)+1)/3

        cur_row = numrows

        print('(rows,cols,bands) = {:d}, {:d}, {:d}'.format(numrows,numcols,numbands_i))

        print 'checking if already in file...'

        changerad = 0
        for row in range(1,numrows):
            checknum = 0

            print '(cur, new) = ',
            print [ float(i) for i in lines[row][0:n_s] ],
            print r

            for c in range(n_s):
                cur_rad = float('{:0.4f}'.format(float(lines[row][c])))
                my_rad = float('{:0.4f}'.format(float(r[c])))


                if(cur_rad != my_rad):
                    break
                checknum = checknum + 1
                if(checknum == n_s):
                    print('ALREADY IN FILE')
                    changerad = 1
                    break
            if(changerad): break
        if(changerad): continue

        tmp_bands = 'tmp_bands.del'
        tmp_bands2 = 'tmp.del'
        cmd0 = 'cat ' + cur_fout + ' | grep -i \"Range\" > ' + tmp_bands2
        cmd1 = 'awk \'{print $4, $10}\' tmp.del > ' + tmp_bands
        os.system(cmd0)
        os.system(cmd1)

        # write all band lines to cur_line
        cur_line = []

        with open(tmp_bands,'r') as fin_bands:
            bands = [band.split() for band in fin_bands]
            numbands_f = len(bands)
            
            print 'num bands...' + str(numbands_f)
        
            # first write all radii
            for r_sp in r:
                cur_line.append(r_sp)
            
            prv_lo = 0.0
            prv_hi = 0.0
    
            for band in bands:
                cur_lo = float(band[0])
                cur_hi = float(band[1])
                if(cur_lo > prv_hi):
                    gap = cur_lo - prv_hi
                    midfreq = (cur_lo + prv_hi)/2
                    normgap = gap/midfreq
                else:
                    normgap = 0

                if(cur_lo>0):
                    cur_line.append(normgap)
                cur_line.append(cur_lo)
                cur_line.append(cur_hi)

                prv_lo = cur_lo
                prv_hi = cur_hi
                
        cmd0 = 'rm ' + tmp_bands2
        cmd1 = 'rm ' + tmp_bands

        os.system(cmd0)
        os.system(cmd1)

        lines.append(cur_line)
        
        # add band numbers to topline if exceed existing
        if(initiate):
            print 'initating bands'
            lines[0].append('1L')
            lines[0].append('1H')
            for i in range(2,numbands_f+1):
                col1 = '{:d}to{:d}'.format(i-1,i)
                col2 = '{:d}L'.format(i)
                col3 = '{:d}H'.format(i)
                lines[0].extend([col1,col2,col3])
            initiate = 0
        elif(numbands_f > numbands_i):
            print 'adding new bands'
            for i in range(numbands_i+1,numbands_f+1):
                col1 = '{:d}to{:d}'.format(i-1,i)
                col2 = '{:d}L'.format(i)
                col3 = '{:d}H'.format(i)
                lines[0].extend([col1,col2,col3])

        with open(fname_band, 'w') as f:
            for line in lines:
                for el in line:
                    f.write('%s,' % el)
                f.write('\n')
            f.close()
    cur_row = numrows + 1

def CalcBandsMin(curr_id, min, n_s):
    
    print 'outputting all gaps for ' + curr_id + ' greater than {:0.4f}'.format(float(min))

    
    fname_band = './dat/' + curr_id + '/' + curr_id + '_bands.csv'
    lines = list(csv.reader( open(fname_band, 'rb') , delimiter = '\t' ))
    lines = [line[0].split(',') for line in lines]
    numrows = len(lines)
    numcols = len(lines[0])
    
    y_ind = []

    for r in range(1, numrows):
        for c in range(n_s+3-1, numcols, 3):
            if c>=len(lines[r]):
                break
            if lines[r][c] == "":
                break
            cur_num = float(lines[r][c])
            if cur_num > float(min):
                print 'comparing {:0.4f} to {:0.4f}'.format(cur_num, float(min))
                y_ind.append(c)
    y_ind = list(set(y_ind))
    y_ind.sort()

    out_name = './dat/' + curr_id + '/' + curr_id + '_bands_min{:0.2f}'.format(float(min)) + '.csv'
    
    with open(out_name,'w+') as fout:
        
        for x in range(numrows):
            print lines[x][0],
            print lines[x][1],
            fout.write(lines[x][0])
            fout.write(',')
            fout.write(lines[x][1])
            fout.write(',')
            print '|',
            for y in y_ind:
                if y>=len(lines[x]):
                    break
                print lines[x][y-1],
                print lines[x][y+1],
                print lines[x][y],
                fout.write(lines[x][y-1])
                fout.write(',')
                fout.write(lines[x][y+1])
                fout.write(',')
                fout.write(lines[x][y])
                fout.write(',')
                print '|',
            fout.write('\n')
            print '\n'

        fout.close()

def GridPlots(curr_id, res, n_s):

    min = 0
    max = 0.5
    inc = 0.05
    
    dpi = 400
    
    rad = [i * inc for i in range(0, int(max/inc)+1)]
    
    fname_band = './dat/' + curr_id + '/Plots/Bands/' + curr_id + '_subbands'
    
    max_x = 0
    max_y = 0
    
    for i in rad:
        for j in rad:
            curr_fband = fname_band + '_r{:0.4f}_r{:0.4f}'.format(i,j)
            curr_fband = curr_fband + '_dpi' + str(dpi) + '.png'
            cur_x = Image.open(curr_fband).size[0]
            cur_y = Image.open(curr_fband).size[1]
            if cur_x > max_x: max_x = cur_x
            if cur_y > max_y: max_y = cur_y
    print 'max (x,y) = ({:d},{:d})'.format(max_x,max_y)

    blank_img = Image.new("RGB", (len(rad) * max_x, len(rad) * max_y), "white")
    
    for i in range(len(rad)):
        for j in range(len(rad)):
            cur_x = i * max_x
            cur_y = j * max_y
            curr_fband = fname_band + '_r{:0.4f}_r{:0.4f}'.format(rad[i],rad[j])
            curr_fband = curr_fband + '_dpi' + str(dpi) + '.png'
            print 'pasting ' + curr_fband
            blank_img.paste(Image.open(curr_fband), (cur_x, cur_y))

    out_img = './dat/' + curr_id + '/Plots/' + curr_id + 'grid_dpi' + str(dpi) + '.png'

    ensure_dir(out_img)

    blank_img.save(out_img)


def MakePlots(curr_id, opt):

    [n_s, n_c, r_c, coords] = readcenters.read(curr_id)

    N = 0
    for i in range(0,n_s):
        N += n_c[i]

#if opt==1:
        #GetNN(curr_id, n_s, n_c, r_c, coords)
        #PlotCenters(curr_id, n_s, n_c, r_c, coords)
        #PlotSk(curr_id, n_s)
        #       elif opt==2:
        subbands = [97, 98, 99, 100, 101, 102, 198, 199, 200, 201, 202, 297, 298, 299, 300, 301, 302]
#    PlotBands(curr_id, 1, n_s)
#PlotSomeBands(curr_id, 1, n_s, subbands)
#    PlotSomeBands(curr_id, 2, n_s, subbands)
    GridPlots(curr_id, 2, n_s)
#    PlotBands(curr_id, 2, n_s)
#    CalcBands(curr_id, 2, n_s)
#    CalcBandsMin(curr_id, 0.03, n_s)
#    CalcBandsMin(curr_id, 0.05, n_s)