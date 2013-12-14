# tm.py
#

import readcenters
import csv
import sys
import os

maxrad = 0.5

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

# assume box length of order 10 (number of points order 100)
def MakeRad(res, n_s):
    #print 'making radii for res {:d}...'.format(res),

    if res == 1:
        inc = 0.1
    elif res == 2:
        inc = 0.05
    elif res == 3:
        inc = 0.025
    else:
        sys.exit('not valid res')

    numper = int(maxrad / inc) + 1
    len = numper ** n_s

    rad = []
    for i in range(0,len):
        cur = []
        for j in range(0,n_s):
            cur.append( round((int(i / numper**(n_s - 1 - j)) % numper) * inc , 5))
        rad.append(cur)
    if res == 2:
        rad = [r for r in rad if r not in MakeRad(1,n_s)]
    elif res == 3:
        rad = [r for r in rad if r not in MakeRad(2,n_s)]

    return rad

def CreateMakeFile(curr_id, N, n_s, rad):
    print 'creating make...'
    
    fname = curr_id + '_TM_make.ctl'

    ffolder = curr_id

    for i in range(0, n_s):
        ffolder += '_r{:0.4f}'.format(float(rad[i]))
    
    fout = './dat/' + curr_id + '/TMUpload/' + ffolder + '/' + fname

    ensure_dir(fout)
    
    with open(fout,'w') as f:
        f.write('(set! geometry (list\n')
        for i in range(1,N+1):
            str_make = '(make cylinder (center CV_LC_C_{:d}) (radius (abs CV_R_{:d})) (height CV_H) (axis (CV_A theta phi)) (material dielCylV))\n'.format(i,i)
            f.write(str_make)
        f.write('))')
        f.close()

    return fname

# creates define file
def CreateDefineFile(curr_id, N, n_s, n_c, rad, coords, N_part):
    print 'creating define for rad ',
    print rad
    
    fname = curr_id + '_define'

    ffolder = curr_id

    for i in range(0, n_s):
        ffolder += '_r{:0.4f}'.format(float(rad[i]))
        fname += '_r{:0.4f}'.format(float(rad[i]))
    
    fname += '.ctl'

    fout = './dat/' + curr_id + '/TMUpload/' + ffolder + '/' + fname
    
    ensure_dir(fout)

    curr_index = 0

    with open(fout,'w') as f:
        for i in range(0, n_s):
            
            curr_rad = float(rad[i])
            
            for j in range(0, n_c[i]):
                
                curr_x = coords[i][j][0] * (N_part**0.5)
                curr_y = coords[i][j][1] * (N_part**0.5)
                
                curr_index += 1

                str_CC = '(define CV_CC_C_{:d} (vector3 {:0.6f} ' \
                    '{:0.6f} 0.000000))\n'.format(curr_index,curr_x, curr_y)
                str_LC  = '(define CV_LC_C_{:d} (cartesian->lattice '\
                    'CV_CC_C_{:d}))\n'.format(curr_index, curr_index)
                str_rad = '(define CV_R_{:d} {:0.6f})\n'.format(curr_index, curr_rad)

                f.write(str_CC + str_LC + str_rad)
    
        f.close()
    return fname

# creates param file
def CreateFileParam(curr_id, n_s, rad, numbands, N, fname_make, fname_define, N_part):

    fin_param = './lib/param_DNT.ctl'
    
    fname_param = curr_id + '_param'
    
    append_rad = ''
    for i in range(0, n_s):
        append_rad += '_r{:0.4f}'.format(float(rad[i]))
    
    fname_param = fname_param + append_rad + '.ctl'
    fname_append = 'MHUDS' + append_rad + '_'

    ffolder = curr_id + append_rad

    fout_param = './dat/' + curr_id + '/TMUpload/' + ffolder + '/' + fname_param
    
    ensure_dir(fout_param)
    old0 = 'VAR_NUM_BANDS'
    new0 = str(numbands)

    old1 = 'VAR_NUM_PTS'
    new1 = str(N_part)

    old2 = 'VAR_MAKE_FILENAME'
    new2 = '\"' + fname_make + '\"'

    old3 = 'VAR_DEF_FILENAME'
    new3 = '\"' + fname_define + '\"'

    old4 = 'VAR_STRING_APPEND'
    new4 = '\"' + fname_append + '\"'

    #copies
    os.system("cp %s %s" % (fin_param, fout_param))

    #overrides line
    cmd0 = "sed -i~ -e 's/" + old0 + "/" + new0 + "/' " + fout_param
    cmd1 = "sed -i~ -e 's/" + old1 + "/" + new1 + "/' " + fout_param
    cmd2 = "sed -i~ -e 's/" + old2 + "/" + new2 + "/' " + fout_param
    cmd3 = "sed -i~ -e 's/" + old3 + "/" + new3 + "/' " + fout_param
    cmd4 = "sed -i~ -e 's/" + old4 + "/" + new4 + "/' " + fout_param

    os.system(cmd0)
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)

    os.system("rm %s" % (fout_param + '~'))

    return fname_param

# creates run file
def CreateFileRun(curr_id, n_s, rad, numbands, fname_param):
    
    fin_run = './lib/run_DNT.ctl'

    append_rad = ''
    for i in range(0, n_s):
        append_rad += '_r{:0.4f}'.format(float(rad[i]))

    ffolder = curr_id + append_rad

    fout_run = './dat/' + curr_id + '/TMUpload/' + ffolder + '/' + curr_id + '_run' + append_rad + '.ctl'
    
    fname_out = './out/' + curr_id + '_out' + append_rad + '.OUT'
    
    old0 = 'VAR_NUM_BANDS'
    new0 = str(numbands)

    old1 = 'VAR_PARAM_FILENAME'
    new1 = fname_param

    old2 = 'VAR_OUT_FILENAME'
    new2 = fname_out.replace('/','\/')

    old3 = 'VAR_CURR_ID'
    new3 = curr_id

    old4 = 'VAR_RAD'
    new4 = append_rad

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
def MakeQsub(curr_id, res, n_s, rad):
    
    fout = './dat/' + curr_id + '/TMUpload/' + curr_id + '_qsub' + str(res) + '.sh'

    with open(fout,'w') as f:
        cmd = '#!/bin/sh\n'
        
        f.write(cmd)
        
        for r in rad:
            
            append_rad = ''
            for i in range(0, n_s):
                append_rad += '_r{:0.4f}'.format(float(r[i]))
            
            ffolder = curr_id + append_rad

            
            curr_run = curr_id + '_run' + append_rad + '.ctl'
            
            cmd = 'qsub ' + './' + ffolder + '/' + curr_run + ';\n'
            
            f.write(cmd)
    
        f.close()

def ReadDel(curr_id):
    fin = './dat/' + curr_id + '/Del/' + curr_id + '_del' + '_T.txt'

    cent_in = list(csv.reader( open(fin, 'rb') , delimiter = '\t' ))

    N = len(cent_in)

    centroids = []

    for c in cent_in:
        x = float(c[0])
        y = float(c[1])
        centroids.append([x,y])
    
    return [N, centroids]

def UploadDel(res, curr_id, numbands):
    print 'running tm.upload...'
    
    [N, centroids] = ReadDel(curr_id)
    
    rad = MakeRad(res, 1)
    
    fname_make = CreateMakeFile(curr_id, N)
    
    n_s = 1
    n_c = [N]
    
    for r in rad:
        fname_define = CreateDefineFile(curr_id, N, n_s, n_c, r, [centroids])
        fname_param = CreateFileParam(curr_id, n_s, r, numbands, N, fname_make, fname_define)
        fname_run = CreateFileRun(curr_id, n_s, r, numbands, fname_param)
    
    MakeQsub(curr_id, res, n_s, rad)
    
    up_folder = './dat/' + curr_id + '/Upload/'
    cmd = 'scp ' + up_folder + '*.ctl ' + up_folder + '*.sh chaneyl@della.princeton.edu:/home/chaneyl/' + curr_id + '/TM/'
    
    os.system(cmd)

def Upload(res, curr_id, numbands, N_part):
    print 'running tm.upload...'

    [n_s, n_c, r_c, coords] = readcenters.read(curr_id)
    
    N = 0
    for i in range(0,n_s):
        N += n_c[i]
    
    rad = MakeRad(res, n_s)
    
    for r in rad:
        fname_make = CreateMakeFile(curr_id, N, n_s, r)
        fname_define = CreateDefineFile(curr_id, N, n_s, n_c, r, coords, N_part)
        fname_param = CreateFileParam(curr_id, n_s, r, numbands, N, fname_make, fname_define, N_part)
        fname_run = CreateFileRun(curr_id, n_s, r, numbands, fname_param)

    MakeQsub(curr_id, res, n_s, rad)

#    up_folder = './dat/' + curr_id + '/TMUpload/'
#    cmd = 'scp -r ' + up_folder + ' chaneyl@della.princeton.edu:/home/chaneyl/' + curr_id + '/TM/'
#    os.system(cmd)

def Download(res, curr_id):
    
    print 'running tm.download...'

    [n_s, n_c, r_c, coords] = readcenters.read(curr_id)

    down_folder = './dat/' + curr_id + '/Download/'
    
    ensure_dir(down_folder)
    
    cmd = 'scp ' + 'chaneyl@della.princeton.edu:/home/chaneyl/' + curr_id + '/out/* ' + down_folder

    os.system(cmd)


def Analyze(res, curr_id):
    print 'running tm.analyze...'
