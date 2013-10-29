/*
Genetic Algorithm Base
Base class for the GATSP Genetic Algorithm
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

#include <algorithm>
#include <fstream>
#include <limits>
#include "gatsp/genetic_algorithm_base.h"

/// Constructor taking an already specified problem.
///
GeneticAlgorithmBase::GeneticAlgorithmBase(
    const ProblemBase& problem,
    int num_individuals = 100,
    double mutate_rate = 1.0,
    double crossover_rate = 0.0,
    unsigned int seed = 0,
    std::string statistics_file = ""
) :
    _problem(&problem),
    _num_individuals(num_individuals),
    _generations(0),
    _mutate_rate(mutate_rate),
    _crossover_rate(crossover_rate),
    _random_generator(std::mt19937(seed)),
    _best_cost(std::numeric_limits<double>::infinity())
{
    // Generate the population.
    _individuals.reserve(_num_individuals);
    for (size_t i = 0; i < _num_individuals; ++i)
    {
    	_individuals.push_back(problem.makeSolution());
    }
    evaluateIndividuals();
    // Create statistics file
    if (statistics_file != "")
    {
        _do_statistics = true;
        _statistics_file = statistics_file;
        createStatistics();
    }
}

/// Constructor taking a specified problem and a solution to seed the population.
///
GeneticAlgorithmBase::GeneticAlgorithmBase(
    const ProblemBase& problem,
    SolutionBase solution,
    int num_individuals = 100,
    double mutate_rate = 1.0,
    double crossover_rate = 0.0,
    unsigned int seed = 0,
    std::string statistics_file = ""
) :
    _problem(&problem),
    _num_individuals(num_individuals),
    _generations(0),
    _mutate_rate(mutate_rate),
    _crossover_rate(crossover_rate),
    _random_generator(std::mt19937(seed)),
    _best_cost(std::numeric_limits<double>::infinity())
{
    // Generate the population.
    _individuals.reserve(_num_individuals);
    for (size_t i = 0; i < _num_individuals; ++i)
    {
        _individuals.push_back(solution);
    }
    evaluateIndividuals();
    // Create statistics file
    if (statistics_file != "")
    {
        _do_statistics = true;
        _statistics_file = statistics_file;
        createStatistics();
    }
}

/// Copy constructor
///
/// Needs to be explicitly defined because of the use of mutex.
GeneticAlgorithmBase::GeneticAlgorithmBase(const GeneticAlgorithmBase & other):
    _problem(other._problem),
    _num_individuals(other._num_individuals),
    _generations(other._generations),
    _individuals(other._individuals),
    _costs(other._costs),
    _best_individual(other._best_individual),
    _best_cost(other._best_cost),
    _mutate_rate(other._mutate_rate),
    _crossover_rate(other._crossover_rate),
    _random_generator(other._random_generator),
    _isRunning(other._isRunning),
    _isDone(other._isDone)
{
    // The mutex is the only one that is not copied, it should be unique to each instantiation.
}

/// Assignment Operator
///
/// Needs to explicitly defined because of mutex. See copy constructor.
GeneticAlgorithmBase& GeneticAlgorithmBase::operator=(const GeneticAlgorithmBase& rhs)
{
    _problem = rhs._problem;
    _num_individuals = rhs._num_individuals;
    _generations = rhs._generations;
    _individuals = rhs._individuals;
    _costs = rhs._costs;
    _best_individual = rhs._best_individual;
    _best_cost = rhs._best_cost;
    _mutate_rate = rhs._mutate_rate;
    _crossover_rate = rhs._crossover_rate;
    _random_generator = rhs._random_generator;
    _isRunning = rhs._isRunning;
    _isDone = rhs._isDone;
    return *this;
}

/// Destructor
///
GeneticAlgorithmBase::~GeneticAlgorithmBase()
{}

/// Evolve the problem until the solution criterion is satisfied.
///
/// Simply step()s until the termination criterion returns true.
/// Should be ideal for threading (using boost::thread).
void GeneticAlgorithmBase::evolve(unsigned int max_generations)
{
    while (! terminate() && ! isRunning() && (_generations < max_generations))
    {
        try
        {
            step();
            boost::this_thread::interruption_point();
        }	
        catch (const boost::thread_interrupted&)
        {
            stop();
            break;
        }
    }

    // If we exited the while loop because the termination criterion
    // returned true, let's just note that we are done.
    if (terminate())
    {
        boost::lock_guard<boost::mutex> guard(_lock);
        _isDone = true;
    }
}

/// Stop the algorithm.
/// May be used both from the outside and inside of the algorithm.
void GeneticAlgorithmBase::stop()
{
    boost::lock_guard<boost::mutex> guard(_lock);
    _isRunning = false;
    return;
}

/// Test if the algorithm is running.
///
/// @returns true if the algorithm is running.
bool GeneticAlgorithmBase::isRunning()
{
    boost::lock_guard<boost::mutex> guard(_lock);
    return _isRunning;
}

/// Test whether the algorithm was stopped by the termination criterion.
///
/// @returns true if the algorithm was stopped by the termination criterion.
bool GeneticAlgorithmBase::isDone()
{
    boost::lock_guard<boost::mutex> guard(_lock);
    return _isDone;
}

/// Advance the algorithm one generation.
///
void GeneticAlgorithmBase::step()
{
    // Make babies, mutate the babies, and evaluate them.
    // Plus some statistics.
    crossover();
    mutate();
    evaluateIndividuals();
    elitism();

    _generations += 1;
    updateStatistics();
    return;
}

/// Get the best solution.
///
SolutionBase GeneticAlgorithmBase::bestSolution()
{
    return _best_individual;
}

/// Get the number of evolved generations.
unsigned int GeneticAlgorithmBase::generations() const
{
    return _generations;
}

/// Writes the statistics to the file
///
void GeneticAlgorithmBase::flushStatistics()
{
    std::ofstream myfile;
    myfile.open(_statistics_file, std::ios::app);
    for (auto g : _statistics)
    {
        std::string seperator("");
        for (auto c : g)
        {
            myfile << seperator << c;
            seperator = ", ";
        }
        myfile << "\n";
    }
    _statistics.clear();
    myfile.close();
}

/// Evaluate the individuals and update the all-time high.
///
void GeneticAlgorithmBase::evaluateIndividuals()
{
    std::vector<double> new_costs;
    for (auto ind : _individuals)
    {
        double cost = _problem->solutionCost(ind);
        new_costs.push_back(cost);
        if (cost < _best_cost)
        {
            _best_cost = cost;
            _best_individual = ind;
        }
    }
    _costs.swap(new_costs);
    return;
}

/// Replace worst solution with the best if there was no improvement.
///
void GeneticAlgorithmBase::elitism()
{
    auto best = std::min_element(_costs.begin(), _costs.end());
    if (*best > _best_cost)
    {        
        double worst_cost = 0;
        size_t worst_individual;
        for (size_t i = 0; i < _costs.size(); ++i)
        {
            if (_costs[i] > worst_cost)
            {
                worst_cost = _costs[i];
                worst_individual = i;
            }
        }
        _individuals[worst_individual] = _best_individual;
        evaluateIndividuals();
    }
}

/// Create the statistics objects and an empty file to dump it into.
///
void GeneticAlgorithmBase::createStatistics()
{
    _statistics.reserve(100);
    std::ofstream myfile;
    myfile.open(_statistics_file, std::ios::trunc);
    myfile.close();
}

/// Notes the costs of the individuals and dumps them into a file once in a while.
///
void GeneticAlgorithmBase::updateStatistics()
{
    if (!_do_statistics)
    {
        return;
    }

    _statistics.push_back(_costs);

    // Write chunk to file
    if ((_generations % 100) == 0)
    {
        flushStatistics();
    }
}