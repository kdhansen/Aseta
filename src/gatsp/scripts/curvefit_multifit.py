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
n_fits = 5
velocity = 1

# Setup the algorithm
num_individuals = 100
mutate_rate = 0.3
crossover_rate = 0.9
generations = 10000
seed = 6

stat_name = 'd1655.stat'

sec_pr_gen = 0.00031865

# Create a TSP from file
prob = gatsp.Euclidean3DProblem()
gatsp.read_tsp('/home/kdh/tsplib95/d1655.tsp', prob)
sol = prob.makeNearestNeighbor()

# # Run the GA 
# ga = gatsp.TraditionalGeneticAlgorithm(prob, sol, num_individuals, mutate_rate, crossover_rate, seed, stat_name)
# start_time = time.clock()
# ga.evolve(generations)
# run_time = time.clock() - start_time
# sec_pr_gen = run_time/generations
# print 'Run time: ' + str(run_time)
# print 'Seconds pr. generation: ' + str(sec_pr_gen)
# sol = ga.bestSolution()

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

# print 'Fitting'

# popts = []
# ends = np.linspace(len(stat_min)/n_fits, len(stat_min), num=n_fits, endpoint=False).astype(int)
# for e in ends:
#     #Guesses
#     s_star_guess = stat_min[e-1]
#     C_guess = 1./(s_star_guess - stat_min[0])
#     Delta_m_guess = -2./C_guess
#     print [s_star_guess, C_guess, Delta_m_guess]
#     #Fitting
#     sigma = np.linspace(1.0, 0.01, e, endpoint=False) # add ", sigma=sigma" to curve_fit if we want to use this.
#     try:
#         popt, pcov = curve_fit(func, stat_i[:e], stat_min[:e], p0=(Delta_m_guess, C_guess, s_star_guess), sigma=sigma)
#         popts.append(popt)
#     except RuntimeError:
#         pass

#Plot
fig_stats = plt.figure()
ax_stats = fig_stats.add_subplot(211)
set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
ax_stats.set_color_cycle(set2)
ax_stats.plot(stat_i, stat_min)
# for e,p in enumerate(popts):
#     fitted_curve = func(stat_i, *p)
#     k_terminate = k_stop(p[0], p[1], velocity, sec_pr_gen)
#     ax_stats.plot(stat_i, fitted_curve, ':')
#     print k_terminate
#     try:
#         ax_stats.plot(k_terminate, fitted_curve[k_terminate], 'o')
#     except IndexError:
#         pass
#     #ax_stats.annotate(str(ends[e]), xy=(ends[e], fitted_curve[ends[e]]), xytext=(ends[e], fitted_curve[ends[e]]+1000), arrowprops=dict(arrowstyle='->'))
route = prob.route(sol)
coords = np.zeros((len(route), 2))
for i,wp in enumerate(route):
    coords[i, :] = np.array([wp.x(), wp.y()])
ax_route = fig_stats.add_subplot(212)
ax_route.plot(coords[:,0], coords[:,1], ':o')
ax_route.set_aspect('equal')

plt.show()
