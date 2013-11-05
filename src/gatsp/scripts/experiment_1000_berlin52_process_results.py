#! /usr/bin/env python
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import brewer2mpl

# Setup-values
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

# Create a great big numpy array [generations x iterations]
# and compute mean along the second axis.
stat_min = np.zeros((generations, iterations))
for n_iter in range(iterations):
    statfile_name = stat_filename(prob_name, seed=n_iter)
    print statfile_name
    with open(statfile_name, 'r') as csvfile:
        statreader = csv.reader(csvfile, delimiter=',')
        for gen,row in enumerate(statreader):
            numbers = [float(x) for x in row]
            stat_min[gen, n_iter] = min(numbers)
stat_mean = np.mean(stat_min, axis=1)
stat_total_min = np.min(stat_min, axis=1)
stat_total_max = np.max(stat_min, axis=1)
# Write the results to a file
mean_file = stat_filename(prob_name + '_results')
with open(mean_file, 'w') as csvfile:
    statwriter = csv.DictWriter(csvfile, ['Mean', 'Minimum', 'Maximum'], delimiter=',')
    statwriter.writeheader()
    for i in xrange(len(stat_mean)):
        statwriter.writerow({'Mean':stat_mean[i], 'Minimum':stat_total_min[i], 'Maximum':stat_total_max[i]})

#Plot
fig_stats = plt.figure()
ax_stats = fig_stats.add_subplot(111)
set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
ax_stats.set_color_cycle(set2)
ax_stats.plot(stat_mean, color=set2[0]) # plot mean
ax_stats.plot(stat_total_max, color=set2[1]) # plot max/min
ax_stats.plot(stat_total_min, color=set2[1])
plt.show()
