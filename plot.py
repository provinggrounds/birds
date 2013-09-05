# plot.py
import numpy as np
import matplotlib.pyplot as plt
import csv

folder = 'data'
curr_id = 'paper'

params = './' + folder + '/' + curr_id + '_params.txt'

data = list(csv.reader(open(params,'r'),delimiter = '\n'))

for i in data:
    for j in i:
        print j
print 
