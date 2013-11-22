# NN.py
#
# gets nearest neighbor statistics

import numpy as np
from readcenters import ensure_dir, read, readdel

def GetDist(pos1, pos2):
    
    def min(a,b):
        if (a < b):
            return a
        else:
            return b
    
    boxlen = 1
    
    [x_1, y_1] = pos1
    [x_2, y_2] = pos2
    
    x_distsq = min( (x_1 - x_2)**2 , boxlen**2 - 2 * boxlen * (x_1 - x_2) + (x_1 - x_2)**2 )
    y_distsq = min( (y_1 - y_2)**2 , boxlen**2 - 2 * boxlen * (y_1 - y_2) + (y_1 - y_2)**2 )
    
    return (x_distsq + y_distsq)**0.5

def GetNNStatsDelMany(main_id, min, max):

    folder = './dat/' + main_id + '/'
    fname = main_id + 'stats_NN_del.csv'

    fout = folder + fname

    ensure_dir(fout)

    with open(fout,'w') as f:
        
        for i in range(min, max):
            
            curr_id = main_id + str(i)
            
            stats = GetNNStatsDelIndiv(curr_id)
            
            for s in stats:
                print s
                f.write(str(s))
                f.write(',')
            f.write('\n')
        
        f.close()

def GetNNStatsDelIndiv(curr_id):
    
    stats = [curr_id]

    [N, coords] = readdel(curr_id)


    NN_sp = []
    for i in range(N):
        curr_nbr = []
        pos1 = coords[i]
        for j in range(0 , N):
            if (i == j): continue
            pos2 = coords[j]
            dist = GetDist(pos1,pos2)
            curr_nbr.append(dist)
        tmp_a = np.array(curr_nbr)
        min = np.min(tmp_a)
        NN_sp.append(min)

    tmp_a = np.array(NN_sp)
    N = len(NN_sp)
    mean = np.mean(tmp_a)
    med = np.median(tmp_a)
    std = np.std(tmp_a)
    min = np.min(tmp_a)
    max = np.max(tmp_a)

    stats.extend([ N, mean, med, std, min, max])

    return stats

# goes through each
# returns radii, ratio (have option), n_c, [N, mean, median, stdev, min, max]
def GetNNStatsMany(main_id, min, max):

    folder = './dat/' + main_id + '/'
    fname = main_id + 'stats_NN.csv'

    fout = folder + fname
    
    ensure_dir(fout)

    with open(fout,'w') as f:
        
        for i in range(min, max):
            
            curr_id = main_id + str(i)
            
            [n_s, n_c, r_c, coords] = read(curr_id)

            stats = GetNNStatsIndiv(curr_id, n_s, n_c, r_c, coords)

            for s in stats:
                print s
                f.write(str(s))
                f.write(',')
            f.write('\n')

        f.close()

def GetNNStatsIndiv(curr_id, n_s, n_c, r_c, coords):
    
    print 'getting nearest neighbor information...'
    
    stats = [curr_id]
    stats.extend(r_c)
    
    for i in range(0,len(r_c)):
        stats.append(r_c[i]/ r_c[0])
    
    # do species first
    for i in range(0,n_s):
        NN_sp = []
        for j in range(0, n_c[i]):
            curr_nbr = []
            pos1 = coords[i][j]
            for l in range(0 , n_c[i]):
                if (l == j): continue
                pos2 = coords[i][l]
                dist = GetDist(pos1,pos2)
                curr_nbr.append(dist)
            tmp_a = np.array(curr_nbr)
            min = np.min(tmp_a)
            NN_sp.append(min)

        tmp_a = np.array(NN_sp)
        N = len(NN_sp)
        mean = np.mean(tmp_a)
        med = np.median(tmp_a)
        std = np.std(tmp_a)
        min = np.min(tmp_a)
        max = np.max(tmp_a)

        stats.extend([ N, mean, med, std, min, max])

    NN_all = []
    for i in range(0,n_s):
        for k in range(0, n_c[i]):
            curr_nbr = []
            pos1 = coords[i][k]
            for j in range(0,n_s):
                for l in range(0, n_c[j]):
                    if (i==j and k==l): continue
                    pos2 = coords[j][l]
                    dist = GetDist(pos1,pos2)
                    curr_nbr.append(dist)
            tmp_a = np.array(curr_nbr)
            min = np.min(tmp_a)
            NN_all.append(min)
    
    tmp_a = np.array(NN_all)
    N = len(NN_all)
    mean = np.mean(tmp_a)
    med = np.median(tmp_a)
    std = np.std(tmp_a)
    min = np.min(tmp_a)
    max = np.max(tmp_a)

    stats.extend([ N, mean, med, std, min, max])

    return stats