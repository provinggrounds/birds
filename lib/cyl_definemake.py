# change the format of the centers.txt to the format of the control file

# centers.txt contains
# num_species
# tot_cells

# then for each species,
# num cells
# radius
# interaction range
# [x,y]
# at the end of the [x,y] there is a linebreak, consisting of a bunch of #'s
#
#

import csv
import os
import sys

if (len(sys.argv)==1):

    curr_id     = 'twospeciesB'
    folder      = './'
    config_rad  = 0 # choose 0 for specifying radii, 1 for using radii elsewhere
    fin_name    = curr_id + '_centers.txt'
    out_folder  = './torun/'

    rad_sp      = 0
    radius      = 0.01
    ratio       = 1


elif (len(sys.argv)==5):

    curr_id     = sys.argv[1]
    folder      = './'
    config_rad  = 0 # choose 0 for specifying radii, 1 for using radii elsewhere
    fin_name    = curr_id + '_centers.txt'
    out_folder  = './torun/'

    rad_sp      = int(sys.argv[2])
    radius      = float(sys.argv[3])
    ratio       = float(sys.argv[4])

    print(rad_sp)
    print(radius)
    print(ratio)


##########

out_indiv   = 0 # if 0, doesn't save the define files for individual species

#########

num_bands   = 250

if(rad_sp == 0):
    rad = [radius, radius/ratio]
elif(rad_sp == 1):
    rad = [radius/ratio, radius]

data        = list(csv.reader(open(fin_name,'rb'),delimiter='\t'))

n_s         = int(float(data[0][0]))
N           = int(float(data[1][0]))
curr_index  = 0
r_beg       = 3

if(rad_sp == 0):
    fout_define   = curr_id + '_spT_rad{:0.4f}_rat{:0.2f}_N{:d}_CYL_DEF.ctl'.format(radius, ratio, N)
elif(rad_sp == 1):
    fout_define   = curr_id + '_spT_irad{:0.4f}_rat{:0.2f}_N{:d}_CYL_DEF.ctl'.format(radius, ratio, N)

fout_make  = curr_id + '_N' + str(N) + '_CYL_MAKE.ctl'

with open(out_folder + fout_define,'w') as f:

    print('num species ' + str(n_s) + '\n')

    for i in range(0,n_s):

        print ('species ' + str(i) + '...\n')
        
        n_c = int(float(data[r_beg][0]))
        
        if(config_rad==0):
            curr_rad = rad[i]
            print('using new radius ' + str(curr_rad) + '\n')
        else:
            curr_rad = float(data[r_beg+1][0])
            print('using old radius ' + str(curr_rad) + '\n')
        
        #species_define = curr_id + '_radius' + str(radius) + '_ratio' + str(ratio) + '_N' + str(n_c) + '_cyl_define.ctl'

        if(out_indiv):

            species_define = curr_id + '_rad' + str(curr_rad) + '_N' + str(n_c) + '_CYL_DEF.ctl'
            species_make = curr_id + '_N' + str(n_c) + '_CYL_MAKE.ctl'
            fsp_index = 1
            
            fsp_define = open(out_folder + species_define,'w')
            fsp_make = open(out_folder + species_make,'w')
            
            fsp_make.write('(set! geometry (list\n')
        
        for j in range(3,n_c+3):
            curr_x = float(data[r_beg+j][0])
            curr_y = float(data[r_beg+j][1])
            curr_index = curr_index + 1
            
            str_CC =    '(define CV_CC_C_{:d} (vector3 {:0.6f} ' \
                '{:0.6f} 0.000000))\n'.format(curr_index,curr_x, curr_y)
            str_LC =    '(define CV_LC_C_{:d} (cartesian->lattice '\
                'CV_CC_C_{:d}))\n'.format(curr_index, curr_index)
            str_rad =   '(define CV_R_{:d} {:0.6f})\n'.format(curr_index, curr_rad)
            
            f.write(str_CC + str_LC + str_rad)
            
            if(out_indiv):

                spstr_CC =    '(define CV_CC_C_{:d} (vector3 {:0.6f} ' \
                    '{:0.6f} 0.000000))\n'.format(fsp_index, curr_x, curr_y)
                spstr_LC =    '(define CV_LC_C_{:d} (cartesian->lattice '\
                    'CV_CC_C_{:d}))\n'.format(fsp_index, fsp_index)
                spstr_rad =   '(define CV_R_{:d} {:0.6f})\n'.format(fsp_index, curr_rad)
                
                fsp_define.write(spstr_CC + spstr_LC + spstr_rad)
                
                str_make = '(make cylinder (center CV_LC_C_{:d}) (radius  (abs CV_R_{:d})) (height CV_H) (axis (CV_A theta phi)) (material dielCylV))\n'.format(fsp_index, fsp_index)
                
                fsp_make.write(str_make)
                
                fsp_index = fsp_index + 1
    
        r_beg = r_beg + n_c + 6 # the 3 includes one dummy line and two empty lines
            
        if(out_indiv):
            fsp_define.close()
            fsp_make.write('))')
            fsp_make.close()

with open(out_folder + fout_make,'w') as f:
    f.write('(set! geometry (list\n')
    for i in range(1,N+1):
        str_make = '(make cylinder (center CV_LC_C_{:d}) (radius  (abs CV_R_{:d})) (height CV_H) (axis (CV_A theta phi)) (material dielCylV))\n'.format(i,i)
        f.write(str_make)
    f.write('))')

#########
# makes param.ctl file

old0 = 'VAR_NUM_BANDS'
new0 = str(num_bands)

old1 = 'VAR_NUM_PTS'
new1 = str(N)

old2 = 'VAR_DEF_FILENAME'
new2 = "\"" + fout_define + "\""

old3 = 'VAR_MAKE_FILENAME'
new3 = "\"" + fout_make + "\""

fin_param = 'param_DNT.ctl'

if(rad_sp == 0):
    fout_param = curr_id + '_spT_rad{:0.4f}_rat{:0.2f}_N{:d}_PARAM.ctl'.format(radius, ratio, N)
elif(rad_sp == 1):
    fout_param = curr_id + '_spT_irad{:0.4f}_rat{:0.2f}_N{:d}_PARAM.ctl'.format(radius, ratio, N)

#copies
os.system("cp %s %s" % (fin_param, out_folder + fout_param))

#overrides line
cmd0 = "sed -i~ -e 's/" + old0 + "/" + new0 + "/' " + out_folder + fout_param
cmd1 = "sed -i~ -e 's/" + old1 + "/" + new1 + "/' " + out_folder + fout_param
cmd2 = "sed -i~ -e 's/" + old2 + "/" + new2 + "/' " + out_folder + fout_param
cmd3 = "sed -i~ -e 's/" + old3 + "/" + new3 + "/' " + out_folder + fout_param

os.system(cmd0)
os.system(cmd1)
os.system(cmd2)
os.system(cmd3)

os.system("rm %s" % (out_folder + fout_param + '~'))

#########
# makes run.ctl file

old0 = 'VAR_NUM_BANDS'
new0 = str(num_bands)

old1 = 'VAR_PARAM_FILENAME'
new1 = fout_param

old2 = 'VAR_OUT_FILENAME'

if(rad_sp == 0):
    new2 = curr_id + '_spT_rad{:0.4f}_rat{:0.2f}_N{:d}.OUT'.format(radius, ratio, N)
elif(rad_sp == 1):
    new2 = curr_id + '_spT_irad{:0.4f}_rat{:0.2f}_N{:d}.OUT'.format(radius, ratio, N)

fin_run = 'run_DNT.ctl'

if(rad_sp == 0):
    fout_run = curr_id + '_spT_rad{:0.4f}_rat{:0.2f}_N{:d}_RUN.ctl'.format(radius, ratio, N)
elif(rad_sp == 1):
    fout_run = curr_id + '_spT_irad{:0.4f}_rat{:0.2f}_N{:d}_RUN.ctl'.format(radius, ratio, N)

os.system("cp %s %s" % (fin_run, out_folder + fout_run))

cmd0 = "sed -i~ -e 's/" + old0 + "/" + new0 + "/' " + out_folder + fout_run
cmd1 = "sed -i~ -e 's/" + old1 + "/" + new1 + "/' " + out_folder + fout_run
cmd2 = "sed -i~ -e 's/" + old2 + "/" + new2 + "/' " + out_folder + fout_run

os.system(cmd0)
os.system(cmd1)
os.system(cmd2)

os.system("rm %s" % (out_folder + fout_run + '~'))


