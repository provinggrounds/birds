# poisson.py

import numpy as np
import random as rd
import shutil
import csv
import sys
import matplotlib
import pylab as pl
import matplotlib.pyplot as plt

from subprocess import Popen, PIPE
from numpy import *

x = []
y = []
N = 1000

filename = './calc_Sk.out'

inputs = [1, 1, N]

proc = Popen(filename, stdin = PIPE)

for i in range(0,N):
    tmpx = rd.random()
    tmpy = rd.random()
    x.append(tmpx)
    y.append(tmpy)
    inputs.extend([tmpx, tmpy])

proc.communicate('\n '.join(str(input) for input in inputs))
proc.wait()

outname = 'Sk_bin.txt'

k,Sk = loadtxt(outname, unpack=True, usecols=[0,1])

plt.figure(1)
plt.subplot(211)
plt.scatter(x,y)
plt.xlabel('x')
plt.ylabel('y')

plt.subplot(212)
plt.plot(k,Sk)
plt.xlabel('k')
plt.ylabel('S(k)')

plt.suptitle('Poisson')
plt.savefig('./data/Poisson.png', bbox_inches = 0)
plt.show()
