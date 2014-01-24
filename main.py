#!/Users/ccl/Library/Enthought/Canopy_64bit/User/bin/python


from lib import mhudtools as mht
from lib import birds, Sk, NN, delanalyze, scripts, bandstructure

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
        
    sample_id             = 'two02' #unique id to identify pattern and associated Sk

    # parameters for packing code
    num_species         = 2
    num_stagegrowth     = 10000
    num_stagerelax      = 300
    num_cells           = [100, 100]   #[250, 120, 95, 90, 60]
    radii               = [0.02, 0.02]  #[0.012, 0.0095, 0.0085, 0.009, 0.008]
    growthrate          = 1
    transmod            = 0.0001
    start_quench        = 0.58
    start_config        = 0 # generate (0) or read (1)
    num_MC              = 200
    len_cell            = 0.08

elif num_args == 13:
    
    sample_id           = sys.argv[1]
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

else: sys.exit('incorrect args')

class SaveFile():
    def __init__(self, sample_id):
        self.params = './dat/' + sample_id + '/' + sample_id + '_params.txt'
        self.centers= './dat/' + sample_id + '/' + sample_id + '_centers.txt'
        self.centers2= './dat/' + sample_id + '/' + sample_id + '_centers_hard.txt'

filenames = SaveFile(sample_id)

code_pack        = './bin/pack.out'
out_pack         = 'centers.txt'
out_pack2        = 'centers_hard.txt'

########  save parameters
def SaveParams():
    
    mht.ensure_dir(filenames.params)

    args = ([num_species, num_stagegrowth, num_stagerelax] +
            num_cells + radii +
            [growthrate, transmod, start_quench, start_config, num_MC, len_cell])

    with open(filenames.params, 'w') as fout:
        for var in args:
            fout.write(str(var)+'\n')
        fout.close()

########  run packing code
def GetConfig():
    
    fout = './dat/' + sample_id + '/log/' + sample_id + '.log'

    mht.ensure_dir(fout)
    
    cmd = code_pack + ' < ' + filenames.params + ' | tee ' + fout
    
    os.system(cmd)
    
    shutil.move(out_pack, filenames.centers)
    shutil.move(out_pack2, filenames.centers2)

def RunRadDiel(sample_id, rad_res, diel_res, numbands, N_part):

    rad_min = 0.0
    rad_max = 0.5

    dielBack = 10.0
    diel_min = 1.0
    diel_max = 1.0

    print 'creating scripts for radius resolution {:d} and dielectric resolution {:d}'.format(rad_res, diel_res)
    
    [N, n_s, n_c, r_c, coords] = mht.read(sample_id)
    
    list_rad = mht.MakeRad(n_s, rad_res, rad_min, rad_max)
    
    list_dielCyl = mht.MakeDiel(n_s, diel_res, diel_min, diel_max, dielBack)
    
    fout_qsub = './dat/{0}/{0}_qsub.sh'.format(sample_id)
    mht.ensure_dir(fout_qsub)
    
    with open(fout_qsub,'w') as f:
        f.write('#!/bin/sh\n')
        for rad in list_rad:
            for dielCyl in list_dielCyl:
                cmd = scripts.Create(sample_id, n_s, rad, dielBack, dielCyl, N, n_c, coords, N_part, numbands)
                f.write(cmd)

# Runs ExtractBands on every file in filelist
def ExtractBandsList(filelist):
    for f in filelist:
        bandstructure.ExtractBandsFromNewFile(f)

def MakeFileList(sample_id, n_s, polarization):
    rad_res = 1
    diel_res = 1
    
    rad_min = 0.0
    rad_max = 0.5
    
    dielBack = 10.0
    diel_min = 1.0
    diel_max = 1.0

    list_rad = mht.MakeRad(n_s, rad_res, rad_min, rad_max)
    list_dielCyl = mht.MakeDiel(n_s, diel_res, diel_min, diel_max, dielBack)
    
    filelist = []

    for rad in list_rad:
        for dielCyl in list_dielCyl:
            str_rad = ''
            str_diel = 'dB{:0.2f}_'.format(dielBack)
            
            for i in range(0, n_s):
                str_rad += '{:0.4f}'.format(float(rad[i]))
                str_diel += '{:0.2f}'.format(float(dielCyl[i]))
                if(i < n_s - 1):
                    str_rad += '_'
                    str_diel += '_'
            
            name = './dat/{0}/{0}_{1}/out/{0}_{1}_{2}_{3}.out'.format(sample_id, str_rad, str_diel, polarization)
            filelist.append(name)
    return filelist

def main():
    
    #os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
    
    newconfig =     0
    runSk =         0
    
    #RunRadDiel(sample_id, rad_res = 1, diel_res = 1, numbands=450, N_part=200)
    
    filelist = MakeFileList('two02', 2, 'TM')
    #ExtractBandsList(filelist)
    bandstructure.GetBandsMin('two02', 0.01, 'TM')

    #NN.GetNNStatsMany('onespecies',1,23)
    #NN.GetNNStatsDelMany('twospecies',1,97)
    #birds.MakePlots(sample_id, 0)
    #delanalyze.MakePlots(sample_id)
    if(newconfig):
        SaveParams()
        GetConfig()
        birds.MakePlots(sample_id, 2)
    if(runSk):
        Sk.CalcSk(sample_id)
        Sk.PlotSk(sample_id)

main()