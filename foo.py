# foo.py

import pylab as pl
import shutil
import csv
import sys

from subprocess import Popen, PIPE
from numpy import *

enable_packing      = 0
enable_Sk           = 1

curr_id             = 'paper' #unique id to identify pattern and associated Sk
folder              = 'data'    #name of folder in which data is stored

# for packing code
num_species         = 5
num_stagegrowth     = 100
num_stagerelax      = 5
num_cells           = [250, 120, 95, 90, 60]
radii               = [0.012, 0.0095, 0.0085, 0.009, 0.008]
growthrate          = 1
transmod            = 0.0001
start_quench        = 0.58
start_config        = 0 # generate (0) or read (1)
num_MC              = 100
len_cell            = 0.08

# for structure factor code
maxx = 1
maxy = 1

kmax = 80
Nmax = 75
mesh = 100

###########

class SaveFile():
    def __init__(self, folder, curr_id):
        self.params = './' + folder + '/' + curr_id + '_params.txt'
        self.centers= './' + folder + '/' + curr_id + '_centers.txt'
        self.Sk     = './' + folder + '/' + curr_id + '_Sk.txt'
        self.Skfig  = './' + folder + '/' + curr_id + '_fig_Sk.png'
        self.cenfig = './' + folder + '/' + curr_id + '_fig_centers.png'


filenames = SaveFile(folder, curr_id)

code_packing        = './pack.out'
code_Sk             = './calc_Sk.out'

out_packing         = 'centers.txt'
out_Sk              = 'Sk_bin.txt'

args = ([num_species, num_stagegrowth, num_stagerelax] +
        num_cells + radii +
        [growthrate, transmod, start_quench, start_config, num_MC, len_cell])

# save params

fo = open(filenames.params, 'w')
for var in args:
    fo.write(str(var)+'\n')
fo.close()

# RUN code_packing 
# format of centers.txt
#    N_tot [1] N_c (x,y) //##### [2] 
if(enable_packing):
    with open(filenames.params,'r') as inf:
        proc = Popen([code_packing], stdin = inf)
        proc.wait()
    shutil.copyfile(out_packing, filenames.centers)


# RUN code_Sk
if(enable_Sk):
    data = list(csv.reader(open(out_packing,'rb'),delimiter='\t'))
    n_s =  int(float(data[0][0]))
    N = int(float(data[1][0]))
    r_beg = 3
    for i in range(0,n_s):
        n_c = int(float(data[r_beg][0]))
        print n_c
        print r_beg
        proc = Popen(code_Sk, stdin = PIPE)
        inputs = [maxx, maxy, n_c]
        
        for j in range(3,n_c+3):
            x = data[r_beg+j][0]
            y = data[r_beg+j][1]
            print '(' + str(x) + ',' + str(y) + ')' 
            inputs.extend( [ x, y ] )
        proc.communicate('\n '.join(str(input) for input in inputs) )
        
        print "inputs were: "+ ', '.join(str(input) for input in inputs)

        r_beg = r_beg + n_c + 6 #the 3 includes one dummy line and two empty lines

# PLOT 

x,y = loadtxt(out_Sk, unpack=True, usecols=[0,1])

plot1 = pl.plot(x,y,'r')
pl.title('S(k) vs k')
pl.xlabel('k')
pl.ylabel('S(k)')
#pl.show()

# SAVE DATA

pl.savefig(filenames.Skfig, bbox_inches=0)
