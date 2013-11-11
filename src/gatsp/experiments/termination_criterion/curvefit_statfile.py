#! /usr/bin/env python
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import brewer2mpl
from scipy.optimize import curve_fit

stat_min = []
stat_i =[]
print 'reading file'
with open('/home/kdh/Code/new_aseta/statfiles/berlin52_1000/berlin52_seed_1.stat', 'r') as csvfile:
    statreader = csv.reader(csvfile, delimiter=',')
    for i,row in enumerate(statreader):
        numbers = [float(x) for x in row]
        stat_min.append(min(numbers))
        stat_i.append(i)
stat_i = np.array(stat_i)
stat_min = np.array(stat_min)
print 'done!'

def func(k, Delta_m, C, s_star):
    return (2*Delta_m)/(k + 2*Delta_m*C) + s_star

sigma = np.linspace(1.0, 0.01, len(stat_min)) # add ", sigma=sigma" to curve_fit if we want to use this.
print 'Fitting'
#Guesses
s_star_guess = stat_min[-1]
C_guess = -1./(s_star_guess - stat_min[0])
Delta_m_guess = 2./C_guess
print [Delta_m_guess, C_guess, s_star_guess]
popt, pcov = curve_fit(func, stat_i, stat_min, p0=(Delta_m_guess, C_guess, s_star_guess))
print popt
fitted_curve = func(stat_i, *popt)

fig_stats = plt.figure()
ax_stats = fig_stats.add_subplot(111)
set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
ax_stats.set_color_cycle(set2)

ax_stats.plot(stat_i, stat_min)
ax_stats.plot(stat_i, fitted_curve, ':')

plt.show()
