# reads [sample_id]_centers.txt
# returns num species, num cells of each species, and (x,y)

import csv
import os

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d): os.makedirs(d)

def readdel(sample_id):

    print 'reading del for sample_id = ' + sample_id + '...'

    fin_path = './dat/' + sample_id + '/Del/' + sample_id + '_del' + '_T.txt'

    data = list(csv.reader(open(fin_path,'rb'),delimiter=' '))

    N = len(data)

    coords = []

    for i in range(N):
        x = float( data[i][0] )
        y = float( data[i][1] )
        coords.append([x,y])

    return [N, coords]

def read(sample_id):
    
    print 'reading sample_id = ' + sample_id + '...'
    
    fin_path = './dat/' + sample_id + '/' + sample_id + '_centers.txt'
    
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

    print 'num species {:d}, '.format(n_s)
    print 'num cells ',
    print n_c

    return [N, n_s, n_c, r_c, coords]

# returns 1 if list contains something not equal to reference value
def ExistNonequal(list, reference):
    for elem in list:
        if elem != reference: return 1
    return 0

# makes list of radii
def MakeRad(n_s, res, min, max):
    
    if res == 1: inc = 0.1
    elif res == 2: inc = 0.05
    elif res == 3: inc = 0.025
    else: sys.exit('not valid res')
    
    numper = int( (max-min) / inc) + 1
    len = numper ** n_s
    
    list_rad = []
    for i in range(0,len):
        cur = []
        for j in range(0,n_s):
            cur.append( round((int(i / numper**(n_s - 1 - j)) % numper) * inc , 5) + min)
        # don't add if all zeros (maybe think of better way to implement this)
        if ExistNonequal(cur, 0): list_rad.append(cur)
    if res == 2: list_rad = [r for r in rad if r not in MakeRad(1,n_s)]
    elif res == 3: list_rad = [r for r in rad if r not in MakeRad(2,n_s)]
    
    return list_rad

def MakeDiel(n_s, res, min, max, dielBack):
    
    if res == 1: inc = 1.0
    elif res == 2: inc = 0.5
    elif res == 3: inc = 0.25
    else: sys.exit('not valid res')
    
    numper = int( (max - min) / inc) + 1
    len = numper ** n_s
    
    list_diel = []
    for i in range(0,len):
        cur = []
        for j in range(0,n_s):
            cur.append( round((int(i / numper**(n_s - 1 - j)) % numper) * inc , 2) + min)
        # don't add if all equal to background dielectric
        if ExistNonequal(cur, dielBack): list_diel.append(cur)
    if res == 2: list_diel = [e for e in eps if e not in MakeDiel(1,n_s)]
    elif res == 3: list_diel = [e for r in eps if e not in MakeDiel(2,n_s)]
    
    return list_diel