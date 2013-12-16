#! /usr/bin/env python
import gatsp
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import brewer2mpl
import time


# print '----'
# print 'Running a refueling problem GA on berlin52'
# # Create a TSP from file
# prob = gatsp.RefuelingProblem()
# gatsp.read_tsp('/home/kdh/tsplib95/berlin52.tsp', prob)
# # Run a GA
# num_individuals = 100
# mutate_rate = 0.3
# crossover_rate = 0.9
# seed = 4
# ga = gatsp.GeneticAlgorithm(
#     prob,
#     num_individuals,
#     mutate_rate,
#     crossover_rate,
#     seed,
#     "stat_file")
# berlin_generations = 5000
# start_berlin = time.clock()
# ga.evolve(berlin_generations)
# berlin_duration = time.clock() - start_berlin
# print 'Generations pr. second: ' + str(berlin_generations/berlin_duration)
# sol = ga.bestSolution()
# route = prob.route(sol)
# coords = np.zeros((len(route), 2))
# for i, wp in enumerate(route):
#     coords[i, :] = np.array([wp.x(), wp.y()])
# print 'Done!'
# print '----'

print '----'
print 'Running a small refueling problem'
# Create a TSP
print 'Creating problem'
# wps = np.array([[4.3, 9.5],
#                 [3.5, 2.4],
#                 [7.4, 8.2],
#                 [4.5, 2.9],
#                 [4.7, 8.4],
#                 [4.6, 6.5],
#                 [3.5, 3.3],
#                 [9.4, 4.2],
#                 [6.1, 5.3],
#                 [5.5, 2.4],
#                 [7.5, 3.4],
#                 [4.5, 1.9],
#                 [6.7, 3.4],
#                 [3.4, 5.5],
#                 [3.7, 6.3],
#                 [9.8, 6.2],
#                 [4.1, 6.3],
#                 [6.4, 5.7]])
wps = np.random.random((20, 2)) * 10
deps = np.array([[5.7, 6.3],
                 [6.2, 9.4]])
p = gatsp.Point(deps[0][0], deps[0][1], 0)
q = gatsp.Quaternion()
prob = gatsp.RefuelingProblem(gatsp.Waypoint(p, q))
for wp in wps:
    p = gatsp.Point(wp[0], wp[1], 0)
    q = gatsp.Quaternion()
    prob.addWaypoint(gatsp.Waypoint(p, q))
for dp in deps[1:]:
    p = gatsp.Point(dp[0], dp[1], 0)
    q = gatsp.Quaternion()
    prob.addDepot(gatsp.Waypoint(p, q))
prob.fuelCapacity(10.0)

# Run a GA
print 'Creating genetic algorithm'
num_individuals = 30
mutate_rate = 0.3
crossover_rate = 0.9
seed = 1
ga = gatsp.GeneticAlgorithm(
    prob,
    num_individuals,
    mutate_rate,
    crossover_rate,
    seed,
    "stat_file")
max_generations = 10000
start_berlin = time.clock()
print 'Starting evolution'
ga.evolve(max_generations)
duration = time.clock() - start_berlin
print 'Generations pr. second: ' + str(max_generations/(duration + 0.000000001))
sol = ga.bestSolution()
route = prob.route(sol)
coords = np.zeros((len(route), 2))
for i, wp in enumerate(route):
    coords[i, :] = np.array([wp.x(), wp.y()])
print 'Done!'
print '----'

fuel_left = prob.fuelLeft(sol)
for f in fuel_left:
    print f

# Plot the results
set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
fig_berlin = plt.figure()
ax_route = fig_berlin.add_subplot(211)
ax_route.plot(coords[:, 0], coords[:, 1], ':', color=set2[0])
ax_route.plot(wps[:, 0], wps[:, 1], 'o', color=set2[1])
ax_route.plot(deps[:, 0], deps[:, 1], 's', color=set2[2])
ax_route.set_aspect('equal')

stat_max = []
stat_min = []
stats = []
stats_i = []
ga.flushStatistics()
with open('stat_file', 'r') as csvfile:
    statreader = csv.reader(csvfile, delimiter=',')
    for i, row in enumerate(statreader):
        numbers = [float(x) for x in row]
        stat_max.append(max(numbers))
        stat_min.append(min(numbers))
        for entry in row:
            stats.append(entry)
            stats_i.append(i)
ax_stats = fig_berlin.add_subplot(212)
ax_stats.set_color_cycle(set2)
ax_stats.plot(stats_i, stats, '.', markersize=5, color=set2[2])
ax_stats.plot(stat_min)
ax_stats.plot(stat_max)
plt.show()
