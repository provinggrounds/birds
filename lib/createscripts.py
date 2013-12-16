# createscripts.py
#

import readcenters
import sys
import os

maxrad = 0.5

class SaveFileName():
    def __init__(self, id, str_rad, str_diel):
        folder = './dat/{0}/{0}_{1}/'.format(id, str_rad)
        
        self.folder = folder
        self.id = id
        self.rad = str_rad
        
        self.TM_PARAM_TEMPLATE        = './lib/paramTM_DNT.ctl'
        self.TE_PARAM_TEMPLATE        = './lib/paramTE_DNT.ctl'
        
        self.TM_RUN_TEMPLATE        = './lib/runTM_DNT.ctl'
        self.TE_RUN_TEMPLATE        = './lib/runTE_DNT.ctl'
        
        str_append = '{0}_{1}_{2}'.format(id, str_rad, str_diel)
        
        self.NAME_MAKE_CYL_V    = '{0}_{1}'.format(str_append,'MakeCylV.ctl')
        self.NAME_DEFINE_CYL_V  = '{0}_{1}'.format(str_append,'DefineCylV.ctl')
        self.NAME_TM_PARAM      = '{0}_{1}'.format(str_append,'TMparam.ctl')
        self.NAME_TM_RUN        = '{0}_{1}'.format(str_append,'TMrun.ctl')
        self.NAME_TE_PARAM      = '{0}_{1}'.format(str_append,'TEparam.ctl')
        self.NAME_TE_RUN        = '{0}_{1}'.format(str_append,'TErun.ctl')

        
        self.PATH_MAKE_CYL_V    = '{0}{1}'.format(folder,self.NAME_MAKE_CYL_V)
        self.PATH_DEFINE_CYL_V  = '{0}{1}'.format(folder,self.NAME_DEFINE_CYL_V)
        self.PATH_TM_PARAM      = '{0}{1}'.format(folder,self.NAME_TM_PARAM)
        self.PATH_TM_RUN        = '{0}{1}'.format(folder,self.NAME_TM_RUN)
        self.PATH_TE_PARAM      = '{0}{1}'.format(folder,self.NAME_TE_PARAM)
        self.PATH_TE_RUN        = '{0}{1}'.format(folder,self.NAME_TE_RUN)

        self.TM_APPEND          = '{0}_{1}'.format(str_append,'TM_')
        self.TM_OUT             = './out/{0}_{1}'.format(str_append,'TM.out')
        self.TE_APPEND          = '{0}_{1}'.format(str_append,'TE_')
        self.TE_OUT             = './out/{0}_{1}'.format(str_append,'TE.out')


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d): os.makedirs(d)

# substitutes str_new for str_old in file
def substitute(str_old, str_new, file):
    
    cmd = "sed -i~ -e 's/{0}/{1}/' {2}".format(str_old, str_new, file)
    
    #print 'old: {0} \n new: {1} \n cmd: {2}'.format(str_old, str_new, cmd)
    
    os.system(cmd)
    os.system("rm %s" % (file + '~'))


# copies file to file_new
def copyfile(file, file_new):
    os.system("cp %s %s" % (file, file_new))

# makes list of radii
def MakeRad(res, n_s):
    
    if res == 1: inc = 0.1
    elif res == 2: inc = 0.05
    elif res == 3: inc = 0.025
    else: sys.exit('not valid res')
    
    numper = int(maxrad / inc) + 1
    len = numper ** n_s
    
    list_rad = []
    for i in range(0,len):
        cur = []
        for j in range(0,n_s):
            cur.append( round((int(i / numper**(n_s - 1 - j)) % numper) * inc , 5))
        list_rad.append(cur)
    if res == 2: list_rad = [r for r in rad if r not in MakeRad(1,n_s)]
    elif res == 3: list_rad = [r for r in rad if r not in MakeRad(2,n_s)]
    
    return list_rad

def CreateFileMakeCylV(files, N):

    fout = files.PATH_MAKE_CYL_V
    
    print 'creating {0}'.format(fout)
    
    with open(fout,'w') as f:
        f.write('(set! geometry (list\n')
        for i in range(1,N+1):
            str_make = '(make cylinder (center CV_LC_C_{:d}) (radius (abs CV_R_{:d})) (height CV_H) (axis (CV_A theta phi)) (material CV_DIELECTRIC_{:d}))\n'.format(i,i,i)
            f.write(str_make)
        f.write('))')
        f.close()

# creates define file
def CreateFileDefineCylV(files, n_s, n_c, coords, N_part):
    
    fout = files.PATH_DEFINE_CYL_V

    print 'creating {0}'.format(fout)

    curr_index = 0

    with open(fout,'w') as f:
        for i in range(0, n_s):
            for j in range(0, n_c[i]):
                
                curr_x = coords[i][j][0] * (N_part**0.5)
                curr_y = coords[i][j][1] * (N_part**0.5)
                
                curr_index += 1

                str_CC = '(define CV_CC_C_{:d} (vector3 {:0.6f} ' \
                    '{:0.6f} 0.000000))\n'.format(curr_index,curr_x, curr_y)
                str_LC  = '(define CV_LC_C_{:d} (cartesian->lattice '\
                    'CV_CC_C_{:d}))\n'.format(curr_index, curr_index)
                str_rad = '(define CV_R_{:d} radCyl_{:d})\n'.format(curr_index, i)
                str_diel= '(define CV_DIELECTRIC_{:d} (make dielectric (epsilon epsCyl_{:d})))\n'.format(curr_index, i)

                f.write(str_CC + str_LC + str_rad + str_diel)

        f.close()

def SubstituteParam(old, new, fout):
    len_old = len(old)
    len_new = len(new)
    if (len_old != len_new): sys.exit('different lengths')
    for i in range(len_old):
        substitute(old[i], new[i], fout)

# formats array of radii into formatted string that goes into param file
def FormatSchemeRad(rad):
    rad_scheme = ''
    for i in range(len(rad)):
        rad_scheme += '(define radCyl_{:d} {:0.6f})\n'.format(i, rad[i])
    return rad_scheme.replace('\n', '\\\n')

# formats dielectric constants into formatted string that goes into param file
def FormatSchemeDiel(dielBack, dielCyl):
    diel_scheme = '(define dielBack (make dielectric (epsilon {:0.2f})))\n'.format(dielBack)
    for i in range(len(dielCyl)):
        diel_scheme += '(define-param epsCyl_{:d} {:0.2f})\n'.format(i, dielCyl[i])
    return diel_scheme.replace('\n', '\\\n')

# creates param file
def CreateFileParam(files, numbands, rad, dielBack, dielCyl, N_part):
    #rad    {:0.6f}
    #diel
    #replaces variables
    old0 = 'VAR_NUM_BANDS'
    old1 = 'VAR_NUM_PTS'
    old2 = 'VAR_MAKE_FILENAME'
    old3 = 'VAR_DEF_FILENAME'
    old4 = 'VAR_RAD'
    old5 = 'VAR_DIEL'
    old6 = 'VAR_STRING_APPEND'
    
    old = [old0, old1, old2, old3, old4, old5, old6]
    
    new0 = str(numbands)
    new1 = str(N_part)
    new2 = '\"' + files.NAME_MAKE_CYL_V + '\"'
    new3 = '\"' + files.NAME_DEFINE_CYL_V + '\"'
    
    new4 = FormatSchemeRad(rad)
    new5 = FormatSchemeDiel(dielBack, dielCyl)
    
    new6_TM = '\"' + files.TM_APPEND + '\"'
    new6_TE = '\"' + files.TE_APPEND + '\"'
    
    new_TM = [new0, new1, new2, new3, new4, new5, new6_TM]
    new_TE = [new0, new1, new2, new3, new4, new5, new6_TE]

    fout_TM = files.PATH_TM_PARAM
    fout_TE = files.PATH_TE_PARAM

    copyfile(files.TM_PARAM_TEMPLATE, fout_TM)
    copyfile(files.TE_PARAM_TEMPLATE, fout_TE)
    
    SubstituteParam(old, new_TM, fout_TM)
    SubstituteParam(old, new_TE, fout_TE)

# creates run file
def CreateFileRun(files, numbands):
    
    #replaces variables
    old0 = 'VAR_NUM_BANDS'
    old1 = 'VAR_FOLDER'
    old2 = 'VAR_PARAM_FILENAME'
    old3 = 'VAR_OUT_FILENAME'

    old = [old0, old1, old2, old3]
    
    new0 = str(numbands)
    new1 = '{0}\/{0}_{1}'.format(files.id, files.rad)
    
    new2_TM = files.NAME_TM_PARAM
    new2_TE = files.NAME_TE_PARAM

    new3_TM = '{0}/'.format(files.TM_OUT).replace('/','\/')
    new3_TE = '{0}/'.format(files.TE_OUT).replace('/','\/')

    new_TM = [new0, new1, new2_TM, new3_TM]
    new_TE = [new0, new1, new2_TE, new3_TM]
    
    fout_TM = files.PATH_TM_RUN
    fout_TE = files.PATH_TE_RUN
    
    copyfile(files.TM_RUN_TEMPLATE, fout_TM)
    copyfile(files.TE_RUN_TEMPLATE, fout_TE)
    
    SubstituteParam(old, new_TM, fout_TM)
    SubstituteParam(old, new_TE, fout_TE)

def CreateScripts(id, n_s, rad, dielBack, dielCyl, N, n_c, coords, N_part, numbands):
    
    str_rad = ''
    
    str_diel = 'dB{:0.2f}_'.format(dielBack)
    
    for i in range(0, n_s):
        str_rad += '{:0.4f}'.format(float(rad[i]))
        str_diel += '{:0.2f}'.format(float(dielCyl[i]))
        if(i < n_s - 1):
            str_rad += '_'
            str_diel += '_'
    
    files = SaveFileName(id, str_rad, str_diel)
    
    ensure_dir(files.folder)
    
    # create make, define, param, run files for current set of radii
    CreateFileMakeCylV(files, N)
    CreateFileDefineCylV(files, n_s, n_c, coords, N_part)
    CreateFileParam(files, numbands, rad, dielBack, dielCyl, N_part)
    CreateFileRun(files, numbands)

    curr_run_TM = './{0}_{1}/'.format(id,str_rad,files.NAME_TM_RUN)
    curr_run_TE = './{0}_{1}/'.format(id,str_rad,files.NAME_TE_RUN)

    cmd_TM = 'qsub {0};\n'.format(curr_run_TM)
    cmd_TE = 'qsub {0};\n'.format(curr_run_TE)

    return cmd_TM + cmd_TE

def RunDifferentRadii(res, id, numbands, N_part):
    
    print 'creating scripts for radii resolution {:d}...'.format(res)
    
    [n_s, n_c, r_c, coords] = readcenters.read(id)
    
    dielCyl = []
    
    N = 0
    for i in range(0,n_s):
        N += n_c[i]
        dielCyl.append(1.)
    
    list_rad = MakeRad(res, n_s)
    
    dielBack = 11.
    
    fout_qsub = './dat/{0}/{0}_qsub{1}.sh'.format(id, str(res))
    
    ensure_dir(fout_qsub)
    
    with open(fout_qsub,'w') as f:
        
        f.write('#!/bin/sh\n')
        
        for rad in list_rad:

            # add run file to qsub
            cmd = CreateScripts(id, n_s, rad, dielBack, dielCyl, N, n_c, coords, N_part, numbands)
            f.write(cmd)