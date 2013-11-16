#!/Users/ccl/Library/Enthought/Canopy_64bit/User/bin/python

from lib import tm, readcenters, birds, Sk

import shutil
import ast
import sys
import os

from subprocess import Popen, PIPE
from numpy import *

num_args = len(sys.argv)

##### example of valid input:
# twospecies1 2 100 30 [200,100] [0,0] 1 0.0001 0.58 0 200 0.08

if num_args == 1:
    
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

elif num_args == 13:
    
    curr_id             = sys.argv[1]
    num_species         = int(sys.argv[2])
    num_stagegrowth     = int(sys.argv[3])
    num_stagerelax      = int(sys.argv[4])
    num_cells           = ast.literal_eval(sys.argv[5])
    radii               = ast.literal_eval(sys.argv[6])
    growthrate          = float(sys.argv[7])
    transmod            = float(sys.argv[8])
    start_quench        = float(sys.argv[9])
    start_config        = int(sys.argv[10])
    num_MC              = int(sys.argv[11])
    len_cell            = float(sys.argv[12])

else:
    
    sys.exit('incorrect args')

# parameters for Sk code
maxx = 1
maxy = 1

class SaveFile():
    def __init__(self, curr_id):
        self.params = './dat/' + curr_id + '/' + curr_id + '_params.txt'
        self.centers= './dat/' + curr_id + '/' + curr_id + '_centers.txt'

filenames = SaveFile(curr_id)

code_pack        = './bin/pack.out'
out_pack         = 'centers.txt'

########  save parameters
def SaveParams():
    
    readcenters.ensure_dir(filenames.params)

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

def main():
    
    os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
    
    newconfig = 0
    runSk = 0
    
    if(newconfig):
    
        SaveParams()
        GetConfig()

    if(runSk):

        Sk.CalcSk(curr_id)
        Sk.PlotSk(curr_id)


    birds.MakePlots(curr_id, 2)
    #tm.Upload(1, curr_id, numbands=350) # 1, 2, 3 are for different resolutions.
    #tm.Download(1, curr_id)
    #tm.Analyze(1, curr_id)

main()