#! /usr/bin/env python
import gatsp
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import brewer2mpl

# Create a TSP from file
prob = gatsp.Euclidean3DProblem()
gatsp.read_tsp('/home/kdh/tsplib95/berlin52.tsp', prob)

# Setup the algorithm
num_individuals = 100
mutate_rate = 0.3
crossover_rate = 0.9
generations = 10000;

iterations = 1000
prob_name = 'berlin52'

def stat_filename(problem_name, seed=None):
    name = problem_name
    if seed:
        name =  name + '_seed_' + str(seed)
    name = name + '.stat'
    return name

# Run the GA 1000 times
for seed in range(iterations):
    statfile_name = stat_filename(prob_name, seed=seed)
    ga = gatsp.TraditionalGeneticAlgorithm(prob, num_individuals, mutate_rate, crossover_rate, seed, statfile_name)
    ga.evolve(generations)
    ga.flushStatistics()
    print seed
