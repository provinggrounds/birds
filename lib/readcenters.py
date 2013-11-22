# reads [curr_id]_centers.txt
# returns num species, num cells of each species, and (x,y)

import csv
import os

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def readdel(curr_id):

    print 'reading del for curr_id = ' + curr_id + '...'

    fin_path = './dat/' + curr_id + '/Del/' + curr_id + '_del' + '_T.txt'

    data = list(csv.reader(open(fin_path,'rb'),delimiter=' '))

    N = len(data)

    coords = []

    for i in range(N):
        x = float( data[i][0] )
        y = float( data[i][1] )
        coords.append([x,y])

    return [N, coords]

def read(curr_id):
    
    print 'reading curr_id = ' + curr_id + '...'
    
    fin_path = './dat/' + curr_id + '/' + curr_id + '_centers.txt'
    
    data = list(csv.reader(open(fin_path,'rb'),delimiter='\t'))
    
    n_s =  int(float(data[0][0]))
    N = int(float(data[1][0]))
    r_beg = 3
    
    n_c = []
    r_c = []
    coords = []

    for i in range(0, n_s):
        
        n_c.append( int(float(data[r_beg][0])) )
        r_c.append( float(data[r_beg+1][0]) )
        
        coords.append([])
        
        for j in range(3,n_c[i]+3):
            x = float(data[r_beg+j][0])
            y = float(data[r_beg+j][1])

            coords[i].append( [x,y] )
        
        r_beg = r_beg + n_c[i] + 6 #the 3 includes one dummy line and two empty lines

    print 'num species {:d}'.format(n_s)
    print 'num cells ',
    print n_c

    return [n_s, n_c, r_c, coords]