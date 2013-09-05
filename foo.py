# foo.py


import pylab as pl
import shutil
import csv
import sys
import matplotlib as mpl

from subprocess import Popen, PIPE
from numpy import *

enable_packing      = 0
enable_Sk           = 1

plot_indiv          = 0

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
        self.Sk     = './' + folder + '/' + curr_id + '_Sk_'
        self.Skfig  = './' + folder + '/' + curr_id + '_fig_Sk_'
        self.cenfig = './' + folder + '/' + curr_id + '_fig_centers_'

# plotSk(load_file, save_file, [spec_id])
def plotSk(*arg):
    if len(arg) != 2 and len(arg) != 3: 
        sys.exit(0)
    load_file = arg[0]
    save_file = arg[1]
    if len(arg) == 3:
        save_file += arg[2]
    else:
        save_file += 'tot'
    save_file += '.png'
    x,y = loadtxt(load_file, unpack=True, usecols=[0,1])
    plot1 = pl.plot(x,y,'r')
    pl.title('S(k) vs k')
    pl.xlabel('k')
    pl.ylabel('S(k)')
#    capt =  
#    pl.text(0.1, 0.1, capt)
    pl.savefig(save_file, bbox_inches=0)   

filenames = SaveFile(folder, curr_id)

code_packing        = './pack.out'
code_Sk             = './calc_Sk.out'

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
    
    f, axarr = pl.subplots(n_s, n_s+1)

    mpl.rcParams['axes.labelsize']= 'small'
    mpl.rcParams['axes.titlesize']= 'small'
    mpl.rcParams['xtick.labelsize'] = 'small'
    mpl.rcParams['ytick.labelsize'] = 'small'

    # calc Sk for each indiv species
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
            total_in.extend( [ x, y ] )
        proc.communicate('\n '.join(str(input) for input in inputs) )
        proc.wait()
        if plot_indiv:
            plotSk(out_Sk, filenames.Skfig, str(i))
        x,y = loadtxt(out_Sk, unpack=True, usecols=[0,1])
        axarr[i,0].plot(x,y,'r')
        axarr[i,0].set_title('Species ' + str(i))
        shutil.copyfile(out_Sk, filenames.Sk + str(i) + '.txt')
        print "inputs were: "+ ', '.join(str(input) for input in inputs)
        r_beg = r_beg + n_c + 6 #the 3 includes one dummy line and two empty lines
    print "total inputs were: "+ ', '.join(str(input) for input in total_in)

    # calc Sk for all species together
    proc = Popen(code_Sk, stdin = PIPE)
    proc.communicate('\n '.join(str(input) for input in total_in) )
    proc.wait()

    plt = plotSk(out_Sk, filenames.Skfig)
    x,y = loadtxt(out_Sk, unpack=True, usecols=[0,1])
    totfig = pl.subplot2grid((n_s,n_s+1),(0,1),colspan = n_s, rowspan = n_s)
    totfig.plot(x,y,'r')
    pl.title('all species')
    pl.xlabel('k')
    pl.ylabel('S(k)')
    f.savefig(filenames.Skfig+'all.png', bbox_inches=0)   
    shutil.copyfile(out_Sk, filenames.Sk + 'tot' + '.txt')
