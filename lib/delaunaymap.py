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
mpl.rcParams['figure.subplot.top'] = 0.9
mpl.rcParams['figure.subplot.wspace'] = 0.25
mpl.rcParams['figure.subplot.hspace'] = 0.25


# makes sure folder for file [f] exists
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

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

def DelSave(curr_id, n_s, i, mapped_coords):

    fout = './dat/' + curr_id + '/Del/' + curr_id + '_del' + '_'
    
    if i == n_s:
        fout = fout + 'T.txt'
    else:
        fout = fout + str(i) + '.txt'

    ensure_dir(fout)

    with open(fout,'w') as f:
        for p in mapped_coords:
            f.write(str(p[0]))
            f.write(' ')
            f.write(str(p[1]))
            f.write('\n')
        f.close()

def AddPlotStuff(ax):
    ax.grid()
    plt.axis('equal')
    plt.axis([0,1,0,1])
    plt.xlabel('x')
    plt.ylabel('y')

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
    
    fig = plt.figure(figsize=(12.0,4.0))
    
    axes_origcombined = plt.subplot(1, 3, 1)
    plt.title('all species, original')
    axes_tesselation = plt.subplot(1, 3, 2)
    plt.title('tesselation')
    axes_mapped = plt.subplot(1, 3, 3)
    plt.title('centroids')
    
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

    mapped_coords = GetCentroids(points[tri.vertices])

    for p in mapped_coords:
        x = p[0]
        y = p[1]
        my_circle_scatter(axes_tesselation, [x], [y], radius=rad, alpha=0.5, color='r')
        my_circle_scatter(axes_mapped, [x], [y], radius=rad, alpha=0.5, color='r')

    AddPlotStuff(axes_origcombined)
    AddPlotStuff(axes_tesselation)
    AddPlotStuff(axes_mapped)
    
    DelSave(curr_id, n_s, n_s, mapped_coords)

    fname = curr_id + '_DelCenters.eps'
    fout = './dat/' + curr_id + '/Del/' + fname
    
    ensure_dir(fout)
    
    plt.suptitle('Delaunay centroid map', fontsize=12)
    plt.savefig(fout, bbox_inches=0, dpi = 400)