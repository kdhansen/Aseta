#! /usr/bin/env python
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import brewer2mpl
from scipy.optimize import curve_fit
import time
import gatsp

#Parameters
velocity = 200
sec_pr_gen = 0.00031865
stat_name = '/home/kdh/Data/Code/new_aseta/statfiles/berlin52_1000/berlin52_seed_1.stat'

# Setup the algorithm
num_individuals = 100
mutate_rate = 0.3
crossover_rate = 0.9
generations = 10000
seed = 6

def func(k, Delta_m, C, s_star):
    return (2*Delta_m)/(k - 2*Delta_m*C) + s_star

def k_stop(Delta_m, C, v, epsilon):
    k = int(np.sqrt(2*Delta_m/epsilon/v) + 2*Delta_m*C)
    return k

stat_min = []
stat_i =[]
print 'reading file'
with open(stat_name, 'r') as csvfile:
    statreader = csv.reader(csvfile, delimiter=',')
    for i,row in enumerate(statreader):
        numbers = [float(x) for x in row]
        stat_min.append(min(numbers))
        stat_i.append(i)
stat_i = np.array(stat_i)
stat_min = np.array(stat_min)
print 'done!'

print 'Fitting'
popts = []
for e in xrange(10,len(stat_min),1):
    #Guesses
    s_star_guess = stat_min[e-1]
    C_guess = 1./(s_star_guess - stat_min[0])
    Delta_m_guess = -2./C_guess
    #Fitting
    sigma = np.linspace(1.0, 0.01, e, endpoint=False) # add ", sigma=sigma" to curve_fit if we want to use this.
    try:
        popt, pcov = curve_fit(func, stat_i[:e], stat_min[:e], p0=(Delta_m_guess, C_guess, s_star_guess), sigma=sigma)
        popts.append(popt)
    except RuntimeError:
        pass

#Plot
fig = plt.figure()
ax = fig.add_subplot(111)
set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
ax.set_color_cycle(set2)
ax.set_aspect('equal', 'box')
ax.set_ylim(top=len(stat_min))
k_terminate = []
for e,p in enumerate(popts):
    k_terminate.append(k_stop(p[0], p[1], velocity, sec_pr_gen))
ax.plot(np.linspace(10, len(stat_min), len(k_terminate)), k_terminate)
plt.show()
