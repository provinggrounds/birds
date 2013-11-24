# delanalyze.py
#
# returns 1D plot

import readcenters
import csv
import sys
import os

from PIL import Image
from lib import tm, delaunaymap, NN
from readcenters import ensure_dir

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

# should modify to include different res
def MakeRad(res):
    
    min = 0.0
    max = 0.5
    inc = 0.05
    
    rad = []
    cur = min
    while cur <= max:
        rad.append(cur)
        cur += inc
    
    return rad

def CalcBands(curr_id, res):
    
    print 'calculating bands...'

    rad = MakeRad(1)
    
    fname_out = './dat/' + curr_id + '/Del/Download/' + curr_id + '_out'
    
    fname_band = './dat/' + curr_id + '/Del/' + curr_id + '_del_bands.csv'

    ensure_dir(fname_band)

    initiate = 0
    
    if not os.path.exists(fname_band):
        lines = [[]]
        lines[0].append('rad')
        numrows = 0
        numcols = 0
        numbands_i = 0
        initiate = 1

    for r in rad:

        cur_fout = fname_out + '_r{:0.4f}'.format(float(r)) + '.OUT'
        
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

            print '(cur, new) = ',

            cur_rad = float('{:0.4f}'.format(float(lines[row][0])))
            my_rad = float('{:0.4f}'.format(float(r)))
            
            print cur_rad,
            print ',',
            print my_rad

            if(cur_rad == my_rad):
                print('ALREADY IN FILE')
                changerad = 1
                break

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
            cur_line.append(r)
            
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

def CalcBandsMin(curr_id, min):
    
    print 'outputting all gaps for ' + curr_id + ' greater than {:0.4f}'.format(float(min))

    fname_band = './dat/' + curr_id + '/Del/' + curr_id + '_del_bands.csv'
    
    lines = list(csv.reader( open(fname_band, 'rb') , delimiter = '\t' ))
    lines = [line[0].split(',') for line in lines]
    numrows = len(lines)
    numcols = len(lines[0])
    
    y_ind = []

    for r in range(1, numrows):
        for c in range(3, numcols, 3):
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

    out_name = './dat/' + curr_id + '/Del/' + curr_id + '_del_bands_min{:0.2f}.csv'.format(float(min))

    with open(out_name,'w+') as fout:
        
        for x in range(numrows):
            print lines[x][0],
            print lines[x][1],
            fout.write(lines[x][0])
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

def PlotSomeBands(curr_id, res, subbands):
    
    print 'plotting bands...'
    
    rad = MakeRad(res)
    
    fname_out = './dat/' + curr_id + '/Del/Download/' + curr_id + '_out'
    
    fname_band = './dat/' + curr_id + '/Del/Plots/Bands/' + curr_id + '_del_subbands'
    
    ensure_dir(fname_band)
    
    min = 0.15
    max = 0.60
    dpi_in = 400
    
    jet = cm = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=0, vmax=len(subbands))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    
    
    for r in rad:
        
        curr_fout = fname_out + '_r{:0.4f}.OUT'.format(float(r))
        curr_fband = fname_band + '_r{:0.4f}'.format(float(r)) + '_dpi' + str(dpi_in) + '.png'
        
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
        plt.yticks([0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6])
        plt.tick_params(axis='y', which='major', labelsize=2)
        plt.tick_params(axis='x', which='major', labelsize=0)
        plt.tick_params(axis='both', which='minor', labelsize=0)
        plt.setp(ax.get_xticklabels(), visible=False)
        
        with open(tmp_bands,'r') as fin_bands:
            bands = [b.split() for b in fin_bands]
            
            for i in range(len(subbands)):
                
                colorVal = scalarMap.to_rgba(i)
                
                print int(subbands[i-1])+1,

                cur_hi = float(bands[int(subbands[i-1])][1])
                nxt_lo = float(bands[int(subbands[i-1])+1][0])
                print cur_hi,
                print ',',
                print nxt_lo,
                if(nxt_lo > cur_hi):
                    print 'gap!'
                    plt.vlines(0, cur_hi, nxt_lo, color=colorVal, linewidths = 100)
                print '\n'
        
        cmd0 = 'rm ' + tmp_bands
        cmd1 = 'rm ' + tmp_bands2
        
        #os.system(cmd0)
        #os.system(cmd1)
        
        plt.axis([-0.1,0.1,min,max])
        ax.yaxis.labelpad = -1
        
        ax.text(0.5,0.97,r,
                horizontalalignment='center',
                transform=ax.transAxes)
                
                #        plt.text(0,5.05,fname_out,size='8',horizontalalignment='center')
        plt.savefig(curr_fband, format='png', dpi=dpi_in, bbox_inches='tight')

def GridPlots(curr_id, res):

    rad = MakeRad(res)

    fname_band = './dat/' + curr_id + '/Del/Plots/Bands/' + curr_id + '_del_subbands'
    
    max_x = 0
    max_y = 0
    dpi_in = 400

    for i in rad:
        curr_fband = fname_band + '_r{:0.4f}'.format(i)
        curr_fband = curr_fband + '_dpi' + str(dpi_in) + '.png'
        cur_x = Image.open(curr_fband).size[0]
        cur_y = Image.open(curr_fband).size[1]
        if cur_x > max_x: max_x = cur_x
        if cur_y > max_y: max_y = cur_y
    print 'max (x,y) = ({:d},{:d})'.format(max_x,max_y)

    blank_img = Image.new("RGB", (len(rad) * max_x, max_y), "white")
    
    for i in range(len(rad)):
        cur_x = i * max_x
        curr_fband = fname_band + '_r{:0.4f}'.format(rad[i])
        curr_fband = curr_fband + '_dpi' + str(dpi_in) + '.png'
        print 'pasting ' + curr_fband
        blank_img.paste(Image.open(curr_fband), (cur_x, 0))

    out_img = './dat/' + curr_id + '/Del/Plots/' + curr_id + '_del_grid_dpi' + str(dpi_in) + '.png'

    ensure_dir(out_img)

    blank_img.save(out_img)

def MakePlots(curr_id):

#CalcBands(curr_id, 1)
#CalcBandsMin(curr_id, 0.03)
#CalcBandsMin(curr_id, 0.05)
#subbands = range(580,591)
    #PlotSomeBands(curr_id, 1, subbands)
    GridPlots(curr_id, 1)
