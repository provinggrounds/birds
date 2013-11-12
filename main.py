#!/Users/ccl/Library/Enthought/Canopy_64bit/User/bin/python

from lib import tm, readcenters
from lib import birds
import shutil
import csv
import sys
import os

from subprocess import Popen, PIPE
from numpy import *

curr_id             = 'twospecies' #unique id to identify pattern and associated Sk

# parameters for packing code
num_species         = 2
num_stagegrowth     = 100
num_stagerelax      = 30
num_cells           = [100, 100]   #[250, 120, 95, 90, 60]
radii               = [0.01, 0.01]  #[0.012, 0.0095, 0.0085, 0.009, 0.008]
growthrate          = 1
transmod            = 0.0001
start_quench        = 0.58
start_config        = 0 # generate (0) or read (1)
num_MC              = 200
len_cell            = 0.08

# parameters for Sk code
maxx = 1
maxy = 1

class SaveFile():
    def __init__(self, curr_id):
        self.params = './dat/' + curr_id + '/' + curr_id + '_params.txt'
        self.centers= './dat/' + curr_id + '/' + curr_id + '_centers.txt'
        self.Sk     = './dat/' + curr_id + '/Sk/' + curr_id + '_Sk_'

filenames = SaveFile(curr_id)

code_pack        = './bin/pack.out'
code_Sk          = './bin/calc_Sk.out'

out_pack         = 'centers.txt'
out_Sk           = 'Sk_bin.txt'

#ensures that folder for file [f] exists
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

########  save parameters
def SaveParams():
    ensure_dir(filenames.params)
    args = ([num_species, num_stagegrowth, num_stagerelax] +
            num_cells + radii +
            [growthrate, transmod, start_quench, start_config, num_MC, len_cell])
    with open(filenames.params, 'w') as fout:
        for var in args:
            fout.write(str(var)+'\n')
        fout.close()

########  run packing code
def GetConfig():
    with open(filenames.params,'r') as fin:
        proc = Popen([code_pack], stdin = fin)
        proc.wait()
    shutil.move(out_pack, filenames.centers)

########  calculate Sk
def CalcSk():
    ensure_dir(filenames.Sk)

    data = list(csv.reader(open(filenames.centers,'rb'),delimiter='\t'))
    
    n_s =  int(float(data[0][0]))
    N = int(float(data[1][0]))
    r_beg = 3
    total_in = [maxx, maxy, N]

    # calc Sk for each species

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
            # print '(' + str(x) + ',' + str(y) + ')'
            inputs.extend( [ tmpx, tmpy ] )
            total_in.extend( [ tmpx, tmpy ] )
        
        tot_x.extend(x)
        tot_y.extend(y)
        
        proc.communicate('\n '.join(str(input) for input in inputs) )
        proc.wait()
        
        shutil.copyfile(out_Sk, filenames.Sk + str(i) + '.txt')
        
        # print "inputs were: "+ ', '.join(str(input) for input in inputs)
        
        r_beg = r_beg + n_c + 6 #the 3 includes one dummy line and two empty lines

    # print "total inputs were: "+ ', '.join(str(input) for input in total_in)

    # calc Sk for total
    proc = Popen(code_Sk, stdin = PIPE)
    proc.communicate('\n '.join(str(input) for input in total_in) )
    proc.wait()

    shutil.move(out_Sk, filenames.Sk + 'T' + '.txt')

def main():
    #SaveParams()
    #GetConfig()
    #CalcSk()
    birds.MakePlots(curr_id)
    #tm.Upload(2, curr_id, numbands=350) # 1, 2, 3 are for different resolutions.
    #tm.Download(1, curr_id)
    #tm.Analyze(1, curr_id)

main()