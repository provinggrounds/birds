# bandstructure.py
#
# a file is stored in the sample_id file, which contains band structure for every sample of that file
#
# ExtractBands takes filename argument, and it will add it
#

import os
import csv
import mhudtools as mht

class SaveFileName():
    def __init__(self, sample_id, polarization):
        folder = './dat/{0}/'.format(sample_id)
        self.bands = '{0}{1}_{2}_bands.csv'.format(folder,sample_id, polarization)
        
class SaveFileNameMin():
    def __init__(self, sample_id, polarization, min):
        folder = './dat/{0}/'.format(sample_id)
        self.bands = '{0}{1}_{2}_bands.csv'.format(folder,sample_id, polarization)
        self.minbands = '{0}{1}_{2}_bands_'.format(folder,sample_id, polarization) + '{:0.2f}.csv'.format(min)

# extracts parameters from filename
# example filename:
#       ./two02_0.0000_0.1000/out/two02_0.0000_0.1000_dB10.00_1.00_1.00_TM.out
def NewFileParams(filepath):

    arr =  filepath.split('/') # [., [id_r1_r2], 'out', [id_r1...dB_dC..TM.out]]
    sample_id = arr[2]
    
    n_s = len(arr[3].split('_')) - 1

    filename = arr[5]
    filename_array = filename.split('_')

    dielB = float(filename_array[n_s+1][2:])
    rad = []
    dielC = []
    for i in range(n_s):
        rad.append(float(filename_array[i+1]))
        dielC.append(float(filename_array[i+n_s+2]))
    polarization = filename_array[ len(filename_array)-1].split('.')[0]

    out = [sample_id, n_s, rad, dielB, dielC, polarization]

    for i in out:
        print i

    return out


# if the file is being made, need to add first row which has headers
# example:
#         |  r0  |  r1  |  r2  |  dielB  |  d0  |  d1  |  d2  |
def MakeFirstRow(n_s):
    row = []
    for i in range(n_s):
        row.append('r{:d}'.format(i))
    row.append('dB')
    for i in range(n_s):
        row.append('d{:d}'.format(i))
    return row

# Gets existing bands data.
# If this is first time, it makes the header row
# Otherwise, gets all the data
# initializes: rows, cols, numbands_initial
def ReadBands(fout, n_s):
    if not os.path.exists(fout):
        firsttime = 1
        lines =[MakeFirstRow(n_s)]
        rows = 0
        cols = 0
        numbands_initial = 0
    else:
        firsttime = 0
        lines = list(csv.reader( open(fout, 'rb') , delimiter = ',' ))
        rows = len(lines)
        cols = len(lines[0])
        
        # Param = n_s + n_s + 1 = 2 n_s + 1 (because of radii, dielB, diel)
        # Cols - Params = 2 * Bands + (Bands - 1) = 3 * Bands - 1
        # Bands = (Cols - Params +  1) / 3
        #       = (cols - 2 n_s) / 3
        numbands_initial = (cols - 2 * n_s) / 3

    return [firsttime, lines, rows, cols, numbands_initial]

# Checks if file is already in bands file by comparing parameters with all
# existing configurations
def CheckAlreadyExist(lines, rad, dielB, diel, rows, cols, n_s):
    params = []
    for r in rad:
        params.append(r)
    params.append(dielB)
    for d in diel:
        params.append(d)

    for r in range(1,rows):
        check_counter = 0
        for c in range(2*n_s+1):
            file_param = float('{:0.4f}'.format(float(lines[r][c])))
            my_param = float('{:0.4f}'.format(float(params[c])))
            if(file_param != my_param): break
            else: check_counter += 1
            if(check_counter == 2*n_s+1):
                print 'already exists'
                return 1
    return 0

# outputs line containing all information from file
def ReadBandsFromNewFile(filename):
    tmp_bands = 'tmp_bands.del'
    tmp_bands2 = 'tmp.del'
    cmd0 = 'cat {0} | grep -i \"Range\" > {1}'.format(filename, tmp_bands2)
    cmd1 = 'awk \'{print $4, $10}\' ' + tmp_bands2 + ' > {0}'.format(tmp_bands)

    os.system(cmd0)
    os.system(cmd1)
    
    cur_line = []
    with open(tmp_bands,'r') as fin_bands:
        bands = [band.split() for band in fin_bands]
        numbands_final = len(bands)
        
        print 'number of bands {:d}'.format(numbands_final)
        
        prv_lo = 0.0
        prv_hi = 0.0
        
        for band in bands:
            cur_lo = float(band[0])
            cur_hi = float(band[1])
            if(cur_lo > prv_hi):
                gap = cur_lo - prv_hi
                midfreq = (cur_lo + prv_hi)/2
                normgap = gap/midfreq
            else:
                normgap = 0
            
            if(cur_lo>0):
                cur_line.append(normgap)
            cur_line.append(cur_lo)
            cur_line.append(cur_hi)
            
            prv_lo = cur_lo
            prv_hi = cur_hi

    os.system('rm {0}'.format(tmp_bands))
    os.system('rm {0}'.format(tmp_bands2))

    return [numbands_final, cur_line]

# for first time, will add bands to first row
# after, will add extra if needed
def AddToFirstRow(initiate, bands_initial, bands_final):
    line = []
    if(initiate):
        line.append('1L')
        line.append('1H')
        for i in range(bands_initial, bands_final):
            col1 = '{:d}to{:d}'.format(i-1,i)
            col2 = '{:d}L'.format(i)
            col3 = '{:d}H'.format(i)
            line.extend([col1,col2,col3])
    else:
        for i in range(numbands_i+1,numbands_f+1):
            col1 = '{:d}to{:d}'.format(i-1,i)
            col2 = '{:d}L'.format(i)
            col3 = '{:d}H'.format(i)
            line.extend([col1,col2,col3])
    return line

def ExtractBandsFromNewFile(filepath):
    
    print 'adding bands...'
    
    [sample_id, n_s, rad, dielB, dielC, polarization] = NewFileParams(filepath)
    
    files = SaveFileName(sample_id, polarization)
    
    fout = files.bands
    
    # read existing file
    [firsttime, lines, rows, cols, numbands_initial] = ReadBands(fout, n_s)
    
    if(CheckAlreadyExist(lines, rad, dielB, dielC, rows, cols, n_s)): return 1
    
    # read band data from filename
    
    # write all band lines to cur_line
    newband = rad
    newband.append(dielB)
    newband.extend(dielC)
    [numbands_final, cur_line] = ReadBandsFromNewFile(filepath)
    newband.extend(cur_line)
    
    lines.append(newband)
    
    # Add band numbers to topline if exceed existing

    if(firsttime):
        lines[0].extend(AddToFirstRow(1, 2, numbands_final + 1))
    elif(numbands_final > numbands_initial):
        lines[0].extend(AddToFirstRow(0, numbands_initial + 1, numbands_final + 1))

    # Finally, write all the data to the bands file
    with open(fout, 'w') as f:
        for line in lines:
            for el in line:
                f.write('%s,' % el)
            f.write('\n')
        f.close()

def GetBandsMin(sample_id, min, polarization):
    
    #print 'outputting all gaps for {0} greater than {:0.4f}'.format(sample_id, float(min))
    
    # we just need to get n_s
    [N, n_s, n_c, r_c, coords] = mht.read(sample_id)

    files = SaveFileNameMin(sample_id, polarization, min)
    
    fin = files.bands
    fout = files.minbands
    
    # read existing file
    [firsttime, lines, rows, cols, numbands_initial] = ReadBands(fin, n_s)

    y_ind = []
    
    # searches through all rows, from 1 to end.
    # c goes through indices of where the gaps appear
    # since there are 2 * n_s + 1 parameters, and the first band
    # is not indexed, the first column index is
    # 2 * n_s + 1 + 2
    
    min = float(min)
    
    for r in range(1, rows):
        for c in range(2 * n_s + 1 + 2, cols, 3):
            if c>=len(lines[r]) or lines[r][c] == "":
                break
            cur_gap = float(lines[r][c])
            if cur_gap > min:
                print 'comparing {:0.4f} to {:0.4f}'.format(cur_gap, min)
                y_ind.append(c)
    y_ind = list(set(y_ind))
    y_ind.sort()
    
    with open(fout,'w+') as f:
        for x in range(rows):
            for c in range(2*n_s + 1):
                f.write(lines[x][c] + ',')
                print lines[x][c],
                print '|',
            for y in y_ind:
                if y>=len(lines[x]):
                    break
                print lines[x][y-1],
                print lines[x][y+1],
                print lines[x][y],
                f.write(lines[x][y-1])
                f.write(',')
                f.write(lines[x][y+1])
                f.write(',')
                f.write(lines[x][y])
                f.write(',')
                print '|',
            f.write('\n')
            print '\n'

