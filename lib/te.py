# tm.py
#

import readcenters
import numpy as np
import sys
import csv
import os

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def MakeRad(res):

    min = 0.0
    max = 0.15
    steps = 16
    
    return np.linspace(min, max, steps)


def CreateMakeFile(curr_id, r, N):
    print 'creating make...'
    
    fname = curr_id + '_TE_make.ctl'
    
    fout_folder = curr_id + '_{:0.3f}'.format(r)

    fout = './dat/' + curr_id + '/TEUpload/' + fout_folder + '/' + fname

    ensure_dir(fout)
    
    with open(fout,'w') as f:
        f.write('(set! geometry (list\n')
        for i in range(1,N+1):
            str_make ='(make cylinder (center CH_Centre_C_{:d}) (radius CylinderHorizontalRadius)  (height CH_Height_C_{:d}) (axis CH_Axis_C_{:d})  (material dielCylH))\n'.format(i,i,i)
            f.write(str_make)
        f.write('))')
        f.close()

    return fname

# creates define file
def CreateDefineFile(curr_id, N, r, walls):
    print 'creating define...',
    
    fname = curr_id + '_TE_define.ctl'

    fout_folder = curr_id + '_{:0.3f}'.format(r)

    fout = './dat/' + curr_id + '/TEUpload/' + fout_folder + '/' + fname
    
    ensure_dir(fout)

    i = 1

    with open(fout,'w') as f:
        
        for w in walls:
            print w
            xL = w[0] * (N**0.5)
            yL = w[1] * (N**0.5)

            xR = w[2] * (N**0.5)
            yR = w[3] * (N**0.5)

            str_CC_L = '(define CH_CC_C_{:d}_L  (vector3 {:0.6f} {:0.6f} 0.000000))\n'.format(i , xL, yL)
            str_CC_R = '(define CH_CC_C_{:d}_R  (vector3 {:0.6f} {:0.6f} 0.000000))\n'.format(i , xR, yR)
            str_CH_C = '(define CH_Centre_C_{:d}  (vector3*  (vector3+ CH_CC_C_{:d}_L   CH_CC_C_{:d}_R) 0.5) )\n'.format(i,i,i)
            str_CH_A = '(define CH_Axis_C_{:d}    (vector3*  (vector3- CH_CC_C_{:d}_L   CH_CC_C_{:d}_R) 0.5) )\n'.format(i,i,i)
            str_CH_H = '(define CH_Height_C_{:d}  (vector3-norm CH_Axis_C_{:d}))\n'.format(i,i)

            f.write( str_CC_L + str_CC_R + str_CH_C + str_CH_A + str_CH_H )
        
            i += 1
        
        f.close()

    return fname

# creates param file
def CreateFileParam(curr_id, numbands, N, r, fname_make, fname_define):

    fin_param = './lib/paramTE_DNT.ctl'
    
    fname_param = curr_id + '_TE_param_{:0.3f}'.format(r) + '.ctl'

    fout_folder = curr_id + '_{:0.3f}'.format(r)

    fout_param = './dat/' + curr_id + '/TEUpload/' + fout_folder + '/' + fname_param

    fname_append = '_' + curr_id + '_MHUDS'

    
    ensure_dir(fout_param)
    
    old0 = 'VAR_NUM_BANDS'
    new0 = str(numbands)

    old1 = 'VAR_NUM_PTS'
    new1 = str(N)

    old2 = 'VAR_MAKE_FILENAME'
    new2 = '\"' + fname_make + '\"'

    old3 = 'VAR_DEF_FILENAME'
    new3 = '\"' + fname_define + '\"'

    old4 = 'VAR_STRING_APPEND'
    new4 = '\"' + fname_append + '\"'
    
    old5 = 'VAR_RADIUS'
    new5 = str(r)

    #copies
    os.system("cp %s %s" % (fin_param, fout_param))

    #overrides line
    cmd0 = "sed -i~ -e 's/" + old0 + "/" + new0 + "/' " + fout_param
    cmd1 = "sed -i~ -e 's/" + old1 + "/" + new1 + "/' " + fout_param
    cmd2 = "sed -i~ -e 's/" + old2 + "/" + new2 + "/' " + fout_param
    cmd3 = "sed -i~ -e 's/" + old3 + "/" + new3 + "/' " + fout_param
    cmd4 = "sed -i~ -e 's/" + old4 + "/" + new4 + "/' " + fout_param
    cmd5 = "sed -i~ -e 's/" + old5 + "/" + new5 + "/' " + fout_param

    os.system(cmd0)
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)
    os.system(cmd5)

    os.system("rm %s" % (fout_param + '~'))

    return fname_param

# creates run file
def CreateFileRun(curr_id, numbands, r, fname_param):
    
    fin_run = './lib/runTE_DNT.ctl'
    
    fout_folder = curr_id + '_{:0.3f}'.format(r)
    
    fout_file = curr_id + '_TE_run_{:0.3f}'.format(r) + '.ctl'
        
    fout_run = './dat/' + curr_id + '/TEUpload/' + fout_folder + '/' + fout_file
    
    fname_out = './out/' + curr_id + '_TE_out_{:0.3f}'.format(r) + '.OUT'
    
    old0 = 'VAR_NUM_BANDS'
    new0 = str(numbands)

    old1 = 'VAR_PARAM_FILENAME'
    new1 = fname_param

    old2 = 'VAR_OUT_FILENAME'
    new2 = fname_out.replace('/','\/')

    old3 = 'VAR_CURR_ID'
    new3 = curr_id
    
    old4 = 'VAR_RAD'
    new4 = '{:0.3f}'.format(r)

    os.system("cp %s %s" % (fin_run, fout_run))

    cmd0 = "sed -i~ -e 's/" + old0 + "/" + new0 + "/' " + fout_run
    cmd1 = "sed -i~ -e 's/" + old1 + "/" + new1 + "/' " + fout_run
    cmd2 = "sed -i~ -e 's/" + old2 + "/" + new2 + "/' " + fout_run
    cmd3 = "sed -i~ -e 's/" + old3 + "/" + new3 + "/g' " + fout_run
    cmd4 = "sed -i~ -e 's/" + old4 + "/" + new4 + "/' " + fout_run

    os.system(cmd0)
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)

    os.system("rm %s" % (fout_run + '~'))


# makes script to upload files to della
def MakeQsub(curr_id, rad):
    
    fout = './dat/' + curr_id + '/TEUpload/' + curr_id + '_TE_qsub.sh'
    
    ensure_dir(fout)
    
    with open(fout,'w') as f:
        cmd = '#!/bin/sh\n'
        
        f.write(cmd)
        for r in rad:
            
            cur_file = curr_id + '_TE_run_{:0.3f}'.format(r) + '.ctl'
            cur_folder = curr_id + '_{:0.3f}'.format(r)
            cur_run = './' + cur_folder + '/' + cur_file
            
            cmd = 'qsub ' + cur_run + ';\n'
            
            f.write(cmd)
        
        f.close()


def readDel(curr_id):
    
    fin = './dat/' + curr_id + '/Del/' + curr_id + '_del_walls.txt'
    
    walls_in = list(csv.reader( open(fin, 'rb') , delimiter = '\t' ))
    
    walls = []
    
    N = len(walls_in)
    
    for w in walls_in:
        xL = float(w[0])
        yL = float(w[1])
        xR = float(w[2])
        yR = float(w[3])
        walls.append( [xL, yL, xR, yR] )

    return [N, walls]

def Upload(curr_id, numbands):
    print 'running TE.upload...'
    
    [N, walls] = readDel(curr_id)
    
    N1 = 584
    
    rad = MakeRad(1)
    
    for r in rad:
        fname_make = CreateMakeFile(curr_id, r, N)
        fname_define = CreateDefineFile(curr_id, N1, r, walls)
        fname_param = CreateFileParam(curr_id, numbands, N1, r, fname_make, fname_define)
        fname_run = CreateFileRun(curr_id, numbands, r, fname_param)
    
    MakeQsub(curr_id, rad)

    up_folder = './dat/' + curr_id + '/TEUpload/'
    cmd = 'scp -r ' + up_folder +' chaneyl@della.princeton.edu:/home/chaneyl/' + curr_id + '/TE/'

    os.system(cmd)

def Download(res, curr_id):
    
    print 'running TE.download...'

    down_folder = './dat/' + curr_id + '/Download/'
    
    ensure_dir(down_folder)
    
    cmd = 'scp ' + 'chaneyl@della.princeton.edu:/home/chaneyl/' + curr_id + '/out/* ' + down_folder

    os.system(cmd)


def Analyze(res, curr_id):
    print 'running tm.analyze...'
