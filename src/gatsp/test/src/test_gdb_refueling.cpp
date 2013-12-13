/*
Test to enable GDB profiling of the refueling problem.
Copyright (C) 2013 Karl D. Hansen (kdh@es.aau.dk)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include "gatsp/refueling_problem.h"

int main(int argc, char const *argv[])
{
    if (argc == 1)
    {
        std::cout << "Usage: test_gdb_refueling generations [dimension [population_size [mutation_rate [crossover_rate]]]]" << std::endl;
        return 0;
    }
    unsigned int generations;
    if (argc >= 2)
        std::istringstream(argv[1]) >> generations;
    unsigned int dimension = 100;
    if (argc >= 3)
        std::istringstream(argv[2]) >> dimension;
    unsigned int population_size = 100;
    if (argc >= 4)
        std::istringstream(argv[3]) >> population_size;
    double mutation_rate = 0.3;
    if (argc >= 5)
        std::istringstream(argv[4]) >> mutation_rate;
    double crossover_rate = 0.9;
    if (argc >= 6)
        std::istringstream(argv[5]) >> crossover_rate;
    std::cout << "Running GATSP test program." << std::endl;
    std::cout << "Generations: " << generations << std::endl;
    std::cout << "Dimension: " << dimension << std::endl;
    std::cout << "Population Size: " << population_size << std::endl;
    std::cout << "Mutation Rate: " << mutation_rate << std::endl;
    std::cout << "Crossover Rate: " << crossover_rate << std::endl;
    std::cout << "---" << std::endl;


    int seed = 0;
    // Make a set of waypoints and load them into a problem
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 10.0);
    RefuelingProblem prob(seed);
    for (size_t i = 0; i < dimension-2; ++i)
    {
        double x = distribution(generator);
        double y = distribution(generator);
        double z = distribution(generator);
        Waypoint wp(Point(x, y, z), Quaternion());
        prob.addWaypoint(wp);
    }
    for (size_t i = 0; i < 2; ++i)
    {
        double x = distribution(generator);
        double y = distribution(generator);
        double z = distribution(generator);
        Waypoint wp(Point(x, y, z), Quaternion());
        prob.addDepot(wp);
    }

    // Create a genetic algorithm
    GeneticAlgorithm ga(prob, population_size, mutation_rate, crossover_rate, seed, "");

    // Evolve and time it
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    ga.evolve(generations);
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "It took me " << time_span.count() << " seconds." << std::endl;
    std::cout << "That is " << generations / time_span.count() << " generations per second." << std::endl;

    // Get the solution
    std::unique_ptr<SolutionBase> solution(ga.bestSolution());

    // Get the route
    std::vector<Waypoint> route = prob.route(solution.get());

    // Get the cost
    double cost = prob.cost(solution.get());
    std::cout << "The cost became " << cost << std::endl;
    return 0;
}
