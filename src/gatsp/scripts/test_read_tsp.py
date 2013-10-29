#! /usr/bin/env python
import gatsp
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np

# Create a TSP from file
prob = gatsp.Euclidean3DProblem()
gatsp.read_tsp('/home/kdh/tsplib95/berlin52.tsp', prob)
sol = prob.makeSolution()
route = prob.route(sol)
coords_before = np.zeros((len(route), 2))
for i,wp in enumerate(route):
    coords_before[i, :] = np.array([str(wp.x()), str(wp.y())])

# Run a GA
ga = gatsp.TraditionalGeneticAlgorithm(prob, 1, 1.0, 1.0, 1)
ga.evolve(1)
sol = ga.bestSolution()
route = prob.route(sol)
coords_after = np.zeros((len(route), 2))
for i,wp in enumerate(route):
    coords_after[i, :] = np.array([str(wp.x()), str(wp.y())])

# Plot the results
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(coords_before[:,0], coords_before[:,1], ':o')
ax.plot(coords_after[:,0], coords_after[:,1], '-or')
ax.set_aspect('equal')
plt.show()
