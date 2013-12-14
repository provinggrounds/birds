# delaunaymap.py
#

from scipy.spatial import Delaunay

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np


# parameters for plotting figures
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['font.size'] = 8.0
mpl.rcParams['figure.subplot.top'] = 0.8
mpl.rcParams['figure.subplot.wspace'] = 0.1
mpl.rcParams['figure.subplot.hspace'] = 0.1


# makes sure folder for file [f] exists
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def unique_rows(data):
    uniq = np.unique(data.view(data.dtype.descr * data.shape[1]))
    return uniq.view(data.dtype).reshape(-1, data.shape[1])

# draws circles on plots
def my_circle_scatter(axes, x_array, y_array, radius=0.5, **kwargs):
    for x, y in zip(x_array, y_array):
        circle = plt.Circle((x,y), radius=radius, **kwargs)
        axes.add_patch(circle)
    return True

def GetCentroids(tri):
    mapped_coords = []
    for t in tri:
        x = [vert[0] for vert in t]
        y = [vert[1] for vert in t]
        centroid = [sum(x) / len(t) , sum(y) / len(t)]
        mapped_coords.append( centroid )
    return mapped_coords

def GetWalls(tri, mapped_coords):
    
    # each triangle correspondings to a centroid, the indices match based
    # on construction of GetCentroids.
    # Go through each triangle, get neighbors. walls connect neighboring centroids
    # After done, remove duplicates
    
    indices = []
    walls = []
    
    # get neighboring indices
    
    num_tri = len(tri.vertices)
    print 'number of simplices = ' + str(num_tri)
    
    for i in range(num_tri):
        for j in tri.neighbors[i]:
            if (j != -1):
                if (i < j):
                    indices.append([i,j])
                else:
                    indices.append([j,i])

    # remove duplicates
    unique_indices = unique_rows(np.array(indices))

    N = len(unique_indices)

    print 'number of walls = ' + str(N)

    # get coords [for some reason, can't get the indices to work easily]

    for i in range(N):
        wall = [mapped_coords[unique_indices[i][0]], mapped_coords[unique_indices[i][1]]]
        walls.append(wall)

    return walls

def DelSaveCentroids(curr_id, mapped_coords):

    fout = './dat/' + curr_id + '/Del/' + curr_id + '_del' + '_T.txt'

    ensure_dir(fout)

    with open(fout,'w') as f:
        for p in mapped_coords:
            f.write(str(p[0]))
            f.write('\t')
            f.write(str(p[1]))
            f.write('\n')
        f.close()

def DelSaveWalls(curr_id, walls):

    fout = './dat/' + curr_id + '/Del/' + curr_id + '_del_walls.txt'

    N = len(walls)

    with open(fout,'w') as f:
        for w in walls:
            xL = str(w[0][0])
            yL = str(w[0][1])
            xR = str(w[1][0])
            yR = str(w[1][1])
            f.write(xL + '\t' + yL + '\t' + xR + '\t' + yR + '\n')
        f.close()


# computes bond-orientational-order parameter
# Q = \abs{ \frac{1}{N} \sum_{m=1}^N \frac{1}{N_b} \sum_{n=1}^{N_b} \exp( i 6 \theta_{mn}) }
def getWallsNonDel(n_s, coords):
    
    all_coords = []
    for i in range(n_s):
        all_coords.extend(coords[i])
    points = np.array(all_coords)
    tri = Delaunay(points)
    
    indices = []
    
    for t in tri.vertices:
        # append 0-1, 0-2, 1-2
        t0 = t[0]
        t1 = t[1]
        t2 = t[2]
        if(t0 < t1):
            indices.append([t0,t1])
        else:
            indices.append([t1,t0])
        if(t0 < t2):
            indices.append([t0,t2])
        else:
            indices.append([t2,t0])
        if(t1 < t2):
            indices.append([t1,t2])
        else:
            indices.append([t2,t1])
    
    unique_indices = unique_rows(np.array(indices))
    
    print points
    
    N = len(unique_indices)
    
    in_walls = []
    for i in range(N):
        wall = [points[unique_indices[i][0]], points[unique_indices[i][1]]]
        in_walls.append(wall)
    
    walls = []
    
    for w in in_walls:
        xL = w[0][0]
        yL = w[0][1]
        xR = w[1][0]
        yR = w[1][1]
        walls.append( [xL, yL, xR, yR] )
    
    return [N, walls]

def CalculateBOO(n_s, coords):
    
    [N, walls] = getWallsNonDel(n_s, coords)
    
    Q = 0.0 + 0.0j




def AddPlotStuff(ax):
    ax.grid()
    ax.set_aspect('equal')
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_xticklabels([],visible='false')
    ax.set_yticklabels([],visible='false')

####### runs delauney map, outputs centroids (and plots)
# first row is original
# second row shows simplices
# third row just shows centers
# radius is 0.01. no significance

def RunDelMap(curr_id, n_s, coords):
    
    print 'Running delaunay...'
    
    rad = 0.005
    
    # get color scheme
    jet = cm = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=0, vmax=n_s)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    
    fig = plt.figure(figsize=(12.0,12.0))
    
    axes_origcombined = fig.add_subplot(2,2,1)
    plt.title('all species, original')
    axes_tesselation = fig.add_subplot(2,2,2)
    plt.title('tesselation')
    axes_mapped = fig.add_subplot(2,2,3)
    plt.title('centroids')
    axes_walls = fig.add_subplot(2,2,4)
    plt.title('centroids, connected')
    AddPlotStuff(axes_origcombined)
    AddPlotStuff(axes_tesselation)
    AddPlotStuff(axes_mapped)
    AddPlotStuff(axes_walls)
    
    all_coords = []
    
    for i in range(n_s):
        
        all_coords.extend(coords[i])
    
        colorVal = scalarMap.to_rgba(i)
        
        for p in coords[i]:
            x = p[0]
            y = p[1]
            my_circle_scatter(axes_origcombined, [x], [y], radius=rad, alpha=0.5, color=colorVal)
        
    points = np.array(all_coords)
    
    tri = Delaunay(points)

    axes_tesselation.triplot(points[:,0], points[:,1], tri.vertices.copy())
    
    # repeating this line because I don't know how to change the ordering of things on same plot
    for i in range(n_s):
        
        colorVal = scalarMap.to_rgba(i)
        
        for p in coords[i]:
            x = p[0]
            y = p[1]
            my_circle_scatter(axes_tesselation, [x], [y], radius=rad, alpha=0.5, color=colorVal)

    Q = CalculateBOO(n_s, coords)

    mapped_coords = GetCentroids(points[tri.vertices])

    for p in mapped_coords:
        x = p[0]
        y = p[1]
        my_circle_scatter(axes_tesselation, [x], [y], radius=rad, alpha=0.5, color='r')
        my_circle_scatter(axes_mapped, [x], [y], radius=rad, alpha=0.5, color='r')
        my_circle_scatter(axes_walls, [x], [y], radius=rad, alpha=0.5, color='r')
    
    walls = GetWalls(tri, mapped_coords)

    print 'length of walls = ',
    print len(walls)
    for w in walls:
        x0 = w[0][0]
        y0 = w[0][1]
        x1 = w[1][0]
        y1 = w[1][1]
        plt.plot([x0,x1],[y0,y1], color = 'b')

    DelSaveCentroids(curr_id, mapped_coords)

    DelSaveWalls(curr_id, walls)

    fname = curr_id + '_DelCenters.eps'
    fout = './dat/' + curr_id + '/Del/' + fname
    
    ensure_dir(fout)
    
    plt.suptitle('Delaunay centroid map', fontsize=12)
    plt.savefig(fout, bbox_inches=0, dpi = 400)