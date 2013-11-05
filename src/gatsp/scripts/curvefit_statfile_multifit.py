#! /usr/bin/env python
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import brewer2mpl
from scipy.optimize import curve_fit

#Parameters
n_fits = 10
velocity = 0.1
gen_pr_sec = 1./100.

def func(k, Delta_m, C, s_star):
    return (2*Delta_m)/(k + 2*Delta_m*C) + s_star

def k_stop(Delta_m, C, v, epsilon):
    print 'Delta_m: ' + str(Delta_m)
    print 'C: ' + str(C)
    print 'v: ' + str(v)
    print 'epsilon: ' + str(epsilon)
    k = int(np.sqrt(2*Delta_m*v/epsilon) + 2*Delta_m*C)
    print 'k: ' + str(k)
    return k

stat_min = []
stat_i =[]
print 'reading file'
with open('/home/kdh/Code/new_aseta/statfiles/berlin52_1000/berlin52_seed_10.stat', 'r') as csvfile:
    statreader = csv.reader(csvfile, delimiter=',')
    for i,row in enumerate(statreader):
        numbers = [float(x) for x in row]
        stat_min.append(min(numbers))
        stat_i.append(i)
stat_i = np.array(stat_i)
stat_min = np.array(stat_min)
print 'done!'

sigma = np.linspace(1.0, 0.01, len(stat_min), endpoint=False) # add ", sigma=sigma" to curve_fit if we want to use this.
print 'Fitting'

popts = []
ends = np.linspace(len(stat_min)/n_fits, len(stat_min), num=n_fits, endpoint=False).astype(int)
for e in ends:
    #Guesses
    s_star_guess = stat_min[e-1]
    C_guess = -1./(s_star_guess - stat_min[0])
    Delta_m_guess = 2./C_guess
    #Fitting
    popt, pcov = curve_fit(func, stat_i[:e], stat_min[:e], p0=(Delta_m_guess, C_guess, s_star_guess))
    popts.append(popt)

#Plot
fig_stats = plt.figure()
ax_stats = fig_stats.add_subplot(111)
set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
ax_stats.set_color_cycle(set2)
ax_stats.plot(stat_i, stat_min)
for e,p in enumerate(popts):
    fitted_curve = func(stat_i, *p)
    k_terminate = k_stop(p[0], p[1], velocity, gen_pr_sec)
    ax_stats.plot(stat_i, fitted_curve, ':')
    ax_stats.plot(k_terminate, fitted_curve[k_terminate], 'o')
    ax_stats.annotate(str(ends[e]), xy=(ends[e], fitted_curve[ends[e]]), xytext=(ends[e], fitted_curve[ends[e]]+1000), arrowprops=dict(arrowstyle='->'))

plt.show()
