#! /usr/bin/env python
from scipy.optimize import curve_fit
import brewer2mpl
import csv
import gatsp
import matplotlib
import numpy as np
import pickle
import sys
import time
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt

#Parameters
velocity = 200
sec_pr_gen = 0.00031865
stat_dir = '/home/kdh/Data/Code/new_aseta/statfiles/berlin52_1000/'

# Setup the algorithm
num_individuals = 100
mutate_rate = 0.3
crossover_rate = 0.9
generations = 10000
seed = 6

def stat_filename(problem_name, seed=None):
    name = problem_name
    if seed:
        name =  name + '_seed_' + str(seed)
    name = name + '.stat'
    return name

def func(k, Delta_m, C, s_star):
    return (2*Delta_m)/(k - 2*Delta_m*C) + s_star

k_results = []
def k_stop(Delta_m, C, v, epsilon):
    k = int(np.sqrt(2*Delta_m/epsilon/v) + 2*Delta_m*C)
    return k
for seed in xrange(1000):
    print '----'
    print seed
    stat_name = stat_dir + stat_filename('berlin52', seed=seed)
    stat_min = []
    stat_i =[]
    print 'Reading file'
    with open(stat_name, 'r') as csvfile:
        statreader = csv.reader(csvfile, delimiter=',')
        for i,row in enumerate(statreader):
            numbers = [float(x) for x in row]
            stat_min.append(min(numbers))
            stat_i.append(i)
    stat_i = np.array(stat_i)
    stat_min = np.array(stat_min)
    print 'Done!'

    print 'Fitting...'
    for e in xrange(10,len(stat_min),50):
        #Guesses
        s_star_guess = stat_min[e-1]
        div = (s_star_guess - stat_min[0])
        if div == 0.:
            div = - sys.float_info.epsilon
        C_guess = 1./div
        Delta_m_guess = -2./C_guess
        #Fitting
        sigma = np.linspace(1.0, 0.01, e, endpoint=False) # add ", sigma=sigma" to curve_fit if we want to use this.
        fit_error = False
        try:
            popt, pcov = curve_fit(func, stat_i[:e], stat_min[:e], p0=(Delta_m_guess, C_guess, s_star_guess), sigma=sigma)
            k_term = k_stop(popt[0], popt[1], velocity, sec_pr_gen)
            if k_term < e:
                break
        except RuntimeError:
            fit_error = True
    if not fit_error:
        sigma = np.linspace(1.0, 0.01, len(stat_min), endpoint=False)
        try:
            popt, pcov = curve_fit(func, stat_i, stat_min, p0=(Delta_m_guess, C_guess, s_star_guess), sigma=sigma)
            added_stats = stat_min/velocity + stat_i * sec_pr_gen
            k_min = added_stats.argmin()
            if (0 < k_term < len(stat_min)):
                T_ga = sec_pr_gen * k_term
                T_exe = stat_min[k_term] / velocity
                T_ga_min = sec_pr_gen * k_min
                T_exe_min = stat_min[k_min] / velocity
                T_wasted = T_ga + T_exe - T_ga_min - T_exe_min
                print k_term
                print k_min
                k_results.append({'k_term':k_term, 'k_min':k_min, 'T_exe':T_exe, 'T_ga_min':T_ga_min , 'T_exe_min':T_exe_min, 'T_wasted':T_wasted})
        except RuntimeError:
            pass


pickle.dump(k_results, open(stat_dir + 'k_results.pickle', 'wb'))

waste_factor = []
for r in k_results:
    waste_factor.append(r['T_wasted']/(r['T_ga_min'] + r['T_exe_min']))
pickle.dump(waste_factor, open(stat_dir + 'waste_factor.pickle', 'wb'))
plt.hist(waste_factor, bins=20)
plt.show()
