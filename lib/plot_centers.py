#!/Users/ccl/Library/Enthought/Canopy_64bit/User/bin/python

import pylab as plt
import shutil
import csv
import sys

from subprocess import Popen, PIPE
from numpy import *
from matplotlib import rc, rcParams

# formatting plots
rc('font',**{'size':4.0,'family':'Times New Roman','serif':['Computer Modern']})
#rc('figure',**{'figsize':[12,10]})

curr_id             = 'twospeciesB' #unique id to identify pattern and associated Sk

# for packing code
num_species         = 2
num_stagegrowth     = 1000
num_stagerelax      = 30
#num_cells           = [250, 120, 95, 90, 60]
#radii               = [0.012, 0.0095, 0.0085, 0.009, 0.008]
#num_cells           = [50,100,800]
#radii               = [0.01,0.0001,0.00005]
num_cells           = [200, 200]
radii               = [0.01, 0.01]
growthrate          = 1
transmod            = 0.0001
start_quench        = 0.58
start_config        = 0 # generate (0) or read (1)
num_MC              = 200
len_cell            = 0.08

# for structure factor code
maxx = 1
maxy = 1

kmax = 50
Nmax = 75
mesh = 100

###########

class SaveFile():
    def __init__(self, folder, curr_id):
        self.params = './' + folder + '/' + curr_id + '_params.txt'
        self.centers= './' + folder + '/' + curr_id + '_centers.txt'
        self.Sk     = './' + folder + '/' + curr_id + '_Sk_'
        self.Skfig  = './' + folder + '/' + curr_id + '_fig_Sk_'
        self.cenfig = './' + folder + '/' + curr_id + '_fig_centers_'

filenames = SaveFile(folder, curr_id)

code_packing        = './bin/pack.out'
code_Sk             = './bin/calc_Sk.out'

out_packing         = 'centers.txt'
out_Sk              = 'Sk_bin.txt'

args = ([num_species, num_stagegrowth, num_stagerelax] +
        num_cells + radii +
        [growthrate, transmod, start_quench, start_config, num_MC, len_cell])

######################
# SAVE PARAMS
fo = open(filenames.params, 'w')
for var in args:
    fo.write(str(var)+'\n')
fo.close()

######################
# RUN code_packing 
# format of centers.txt
#    N_tot [1] N_c (x,y) //##### [2] 
if(enable_packing):
    with open(filenames.params,'r') as inf:
        proc = Popen([code_packing], stdin = inf)
        proc.wait()
    shutil.copyfile(out_packing, filenames.centers)

######################
# RUN code_Sk
if(enable_Sk):

    data = list(csv.reader(open(out_packing,'rb'),delimiter='\t'))
    
    n_s =  int(float(data[0][0]))
    N = int(float(data[1][0]))
    r_beg = 3
    total_in = [maxx, maxy, N]
    
    f = plt.figure(1)
    
    # calc Sk for each indiv species

    tot_x = []
    tot_y = []

    for i in range(0,n_s):
        n_c = int(float(data[r_beg][0]))
        print n_c
        print r_beg
        proc = Popen(code_Sk, stdin = PIPE)
        inputs = [maxx, maxy, n_c]
        x = []
        y = []

        for j in range(3,n_c+3):
            tmpx = data[r_beg+j][0]
            tmpy = data[r_beg+j][1]
            x.append(tmpx)
            y.append(tmpy)
            print '(' + str(x) + ',' + str(y) + ')' 
            inputs.extend( [ tmpx, tmpy ] )
            total_in.extend( [ tmpx, tmpy ] )

        tot_x.extend(x)
        tot_y.extend(y)

        proc.communicate('\n '.join(str(input) for input in inputs) )
        proc.wait()
        k, Sk = loadtxt(out_Sk, unpack=True, usecols=[0,1])

        plt.subplot(2,n_s+1,i+1,aspect='equal')
        plt.scatter(x,y,s=5)
        plt.title('Species '+str(i))
        plt.xlabel('x')
        plt.ylabel('y')

        rng = [0,1]
        plt.xlim(rng)
        plt.ylim(rng)

        plt.subplot(2,n_s+1,i+1+1+n_s, aspect = 'equal')
        plt.plot(k,Sk)
        plt.xlabel('k')
        plt.ylabel('S(k)')

        plt.xlim([0,3])
        plt.ylim([0,3])

        shutil.copyfile(out_Sk, filenames.Sk + str(i) + '.txt')
        print "inputs were: "+ ', '.join(str(input) for input in inputs)
        r_beg = r_beg + n_c + 6 #the 3 includes one dummy line and two empty lines
    
    print "total inputs were: "+ ', '.join(str(input) for input in total_in)

    # calc Sk for all species together
    proc = Popen(code_Sk, stdin = PIPE)
    proc.communicate('\n '.join(str(input) for input in total_in) )
    proc.wait()

    k, Sk = loadtxt(out_Sk, unpack=True, usecols=[0,1])

    plt.subplot(2,n_s+1,n_s+1,aspect = 'equal')
    plt.scatter(tot_x,tot_y,s=5)
    plt.title('all species')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(rng)
    plt.ylim(rng)
   
    plt.subplot(2,n_s+1,2*(n_s+1),aspect = 'equal')
    plt.plot(k,Sk)
    plt.xlabel('k')
    plt.ylabel('S(k)')

    plt.xlim([0,3])
    plt.ylim([0,3])

    plt.tight_layout(pad=0.5, h_pad=None, w_pad = None, rect = [0,0,1,1])
    plt.savefig(filenames.Skfig+'all.eps', bbox_inches=0)
    
    shutil.copyfile(out_Sk, filenames.Sk + 'tot' + '.txt')
