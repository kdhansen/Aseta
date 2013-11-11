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
gatsp.read_tsp('/home/kdh/tsplib95/berlin52.tsp', prob)
sol = prob.makeNearestNeighbor()
route = prob.route(sol)
coords = np.zeros((len(route), 2))
for i,wp in enumerate(route):
    coords[i, :] = np.array([wp.x(), wp.y()])

# Plot the results
fig_route = plt.figure()
ax_route = fig_route.add_subplot(111)
ax_route.plot(coords[:,0], coords[:,1], ':o')
ax_route.set_aspect('equal')

plt.show()
