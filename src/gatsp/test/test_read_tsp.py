#! /usr/bin/env python
import gatsp
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import brewer2mpl
import timeit

# Create a TSP from file
prob = gatsp.Euclidean3DProblem()
gatsp.read_tsp('/home/kdh/tsplib95/d1655.tsp', prob)

# Run a GA
num_individuals = 100
mutate_rate = 0.3
crossover_rate = 0.9
seed = 4
ga = gatsp.TraditionalGeneticAlgorithm(prob, num_individuals, mutate_rate, crossover_rate, seed, "stat_file")
sol = ga.bestSolution()
route = prob.route(sol)
coords_before = np.zeros((len(route), 2))
for i,wp in enumerate(route):
    coords_before[i, :] = np.array([str(wp.x()), str(wp.y())])
ga.evolve(10000)
sol = ga.bestSolution()
route = prob.route(sol)
coords_after = np.zeros((len(route), 2))
for i,wp in enumerate(route):
    coords_after[i, :] = np.array([str(wp.x()), str(wp.y())])

# Plot the results
fig_route = plt.figure()
ax_route = fig_route.add_subplot(111)
ax_route.plot(coords_before[:,0], coords_before[:,1], ':o')
ax_route.plot(coords_after[:,0], coords_after[:,1], '-or')
ax_route.set_aspect('equal')

stat_max = []
stat_min = []
# stats = []
# stats_i =[]
ga.flushStatistics()
with open('stat_file', 'r') as csvfile:
    statreader = csv.reader(csvfile, delimiter=',')
    for i,row in enumerate(statreader):
        numbers = [float(x) for x in row]
        stat_max.append(max(numbers))
        stat_min.append(min(numbers))
        # for entry in row:
        #     stats.append(entry)
        #     stats_i.append(i)
fig_stats = plt.figure()
ax_stats = fig_stats.add_subplot(111)
set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
ax_stats.set_color_cycle(set2)
# ax_stats.plot(stats_i, stats, '.', markersize=5, color=set2[2])
ax_stats.plot(stat_min)
ax_stats.plot(stat_max)
plt.show()
