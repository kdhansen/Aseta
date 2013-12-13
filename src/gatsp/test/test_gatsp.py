#! /usr/bin/env python
import gatsp
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import brewer2mpl
import time

print '----'
print 'Making a square 1x1:'
wp1 = gatsp.Waypoint(gatsp.Point(0.0, 0.0, 0.0), gatsp.Quaternion())
wp2 = gatsp.Waypoint(gatsp.Point(1.0, 0.0, 0.0), gatsp.Quaternion())
wp3 = gatsp.Waypoint(gatsp.Point(1.0, 1.0, 0.0), gatsp.Quaternion())
wp4 = gatsp.Waypoint(gatsp.Point(0.0, 1.0, 0.0), gatsp.Quaternion())
p4 = gatsp.Euclidean3DProblem()
p4.addWaypoint(wp1)
p4.addWaypoint(wp2)
p4.addWaypoint(wp3)
p4.addWaypoint(wp4)
s4 = p4.makeSolution()
print 'Cost should be 4. It is ' + str(p4.cost(s4))

print '----'
print 'Running a GA on berlin52'
# Create a TSP from file
prob = gatsp.Euclidean3DProblem()
gatsp.read_tsp('/home/kdh/tsplib95/berlin52.tsp', prob)
# Run a GA
num_individuals = 100
mutate_rate = 0.3
crossover_rate = 0.9
seed = 4
ga = gatsp.GeneticAlgorithm(
    prob,
    num_individuals,
    mutate_rate,
    crossover_rate,
    seed,
    "stat_file")
berlin_generations = 5000
start_berlin = time.clock()
ga.evolve(berlin_generations)
berlin_duration = time.clock() - start_berlin
print 'Generations pr. second: ' + str(berlin_generations/berlin_duration)
sol = ga.bestSolution()
route = prob.route(sol)
coords = np.zeros((len(route), 2))
for i, wp in enumerate(route):
    coords[i, :] = np.array([wp.x(), wp.y()])

# Plot the results
set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
fig_berlin = plt.figure()
ax_route = fig_berlin.add_subplot(211)
ax_route.plot(coords[:, 0], coords[:, 1], ':o', color=set2[0])
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

print '----'
print 'Running a refueling problem GA on berlin52'
# Create a TSP from file
prob = gatsp.RefuelingProblem()
gatsp.read_tsp('/home/kdh/tsplib95/berlin52.tsp', prob)
# Run a GA
num_individuals = 100
mutate_rate = 0.3
crossover_rate = 0.9
seed = 4
ga = gatsp.GeneticAlgorithm(
    prob,
    num_individuals,
    mutate_rate,
    crossover_rate,
    seed,
    "stat_file")
berlin_generations = 5000
start_berlin = time.clock()
ga.evolve(berlin_generations)
berlin_duration = time.clock() - start_berlin
print 'Generations pr. second: ' + str(berlin_generations/berlin_duration)
sol = ga.bestSolution()
route = prob.route(sol)
coords = np.zeros((len(route), 2))
for i, wp in enumerate(route):
    coords[i, :] = np.array([wp.x(), wp.y()])
print 'Done!'

# print '----'
# print 'Running a GA on d1655'
# # Create a TSP from file
# prob = gatsp.Euclidean3DProblem()
# gatsp.read_tsp('/home/kdh/tsplib95/d1655.tsp', prob)
# # Run a GA
# num_individuals = 100
# mutate_rate = 0.3
# crossover_rate = 0.9
# seed = 4
# ga = gatsp.GeneticAlgorithm(prob, num_individuals, mutate_rate, crossover_rate, seed)
# d1655_generations = 100
# start_d1655 = time.clock()
# ga.evolve(d1655_generations)
# d1655_duration = time.clock() - start_d1655
# print 'Generations pr. second: ' + str(d1655_generations/d1655_duration)
# sol = ga.bestSolution()
# route = prob.route(sol)
# coords = np.zeros((len(route), 2))
# for i,wp in enumerate(route):
#     coords[i, :] = np.array([wp.x(), wp.y()])

# # Plot the results
# fig_d1655 = plt.figure()
# ax_route = fig_d1655.add_subplot(111)
# ax_route.plot(coords[:,0], coords[:,1], ':o')
# ax_route.set_aspect('equal')

plt.show()
