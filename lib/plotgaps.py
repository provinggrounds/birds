import pylab as plt
import csv
import numpy as np
import sys
import os

# twospeciesB_spT_irad0.0160_rat1.10_N200.OUT

fin = sys.argv[1]

fin_dat = fin.split('_')
#twospecies_spT_irad0.0130_rat1.50_N299.OUT
#['twospecies', 'spT', 'irad0.0130', 'rat1.50', 'N299.OUT']
print('parsing filename... ' + str(fin_dat))

fin_id = fin_dat[0]

#########
fin_rat = float(fin_dat[3][3:])

# specific to two-species, where i in irad signifies which species rad refers to

if (fin_dat[2][0]=='i'):
    rad_sp = 1
    fin_rad = float(fin_dat[2][4:])
    rad = [fin_rad/fin_rat, fin_rad]
elif (fin_dat[2][0]=='r'):
    rad_sp = 0
    fin_rad = float(fin_dat[2][3:])
    rad = [fin_rad, fin_rad/fin_rat]
else:
    print('invalid filename')
    sys.exit()
out_rad = '{:.4f}, {:.4f}'.format(rad[0],rad[1])

print('radius... ' + str(rad))

tmp_bands = 'tmp_bands.del'
cmd0 = 'cat ' + fin + ' | grep -i \"Range\" > tmp.del'
cmd1 = 'awk \'{print $4, $10}\' tmp.del > ' + tmp_bands
os.system(cmd0)
os.system(cmd1)


fig = plt.figure(figsize=(10,10))

ax = fig.gca()
ax.set_aspect(5)
ax.xaxis.set_ticks_position('none')
plt.tick_params(axis='y', which='major', labelsize=8)
plt.tick_params(axis='x', which='major', labelsize=0)
plt.tick_params(axis='both', which='minor', labelsize=0)
plt.setp(ax.get_xticklabels(), visible=False)

with open(tmp_bands,'r') as fin_bands:
    bands = [band.split() for band in fin_bands]
    numbands_f = len(bands)
    print 'num bands...' + str(numbands_f)
    
    prv_lo = 0.0
    prv_hi = 0.0

    for band in bands:
        cur_lo = float(band[0])
        cur_hi = float(band[1])

        if(cur_lo > prv_hi):
            plt.vlines(0, prv_hi, cur_lo, color='r', linewidths = 4)

        prv_lo = cur_lo
        prv_hi = cur_hi

plt.axis([-1,1,2.5,5])
plt.text(0,5.05,out_rad,size='8',horizontalalignment='center')
#plt.title(rad)
plt.savefig(fin + '.png', format='png', dpi=900, bbox_inches='tight')