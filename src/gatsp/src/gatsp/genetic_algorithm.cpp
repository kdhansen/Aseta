/*
Genetic Algorithm
The GATSP Genetic Algorithm to solve instances of the ProblemBase.
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

#include <iostream>
#include <fstream>
#include "gatsp/genetic_algorithm.h"
#include "gatsp/refueling_problem.h"

/// Constructor
///
/// Takes a specified problem and possibly some algorithm specific
/// parameters.
///
/// @param problem The problem to solve.
/// @param num_individuals Number of individuals in the population.
/// @param mutate_rate The probability that an individual is mutated.
/// @param crossover_rate The probability that crossover will happen between 
///                       the mother and the father of a new individual.
/// @param seed The seed for the probabilistic engine.
/// @param statistics_file Name of the file that statistics will be dumped to.
///                        If no name is specified, no statistics will be recorded.
GeneticAlgorithm::GeneticAlgorithm(
    const ProblemBase& problem,
    int num_individuals,
    double mutate_rate,
    double crossover_rate,
    unsigned int seed,
    std::string statistics_file
) :
    _problem(&problem),
    _num_individuals(num_individuals),
    _generations(0),
    _mutate_rate(mutate_rate),
    _crossover_rate(crossover_rate),
    _random_generator(std::mt19937(seed)),
    _best_cost(std::numeric_limits<double>::infinity()),
    _isRunning(false),
    _isDone(false)
{
    // Generate the population.
    _population.individuals.reserve(_num_individuals);
    for (size_t i = 0; i < _num_individuals; ++i)
    {
        _population.individuals.push_back(problem.makeSolution());
    }
    
    // Make a dummy best individual becauses evaluate need it as a placeholder for the 
    // true best individual.
    _best_individual = _population.individuals[0]->clone();
    evaluateIndividuals();
    
    // Create statistics file
    if (statistics_file != "")
    {
        _do_statistics = true;
        _statistics_file = statistics_file;
        createStatistics();
    }
}

/// Evolve the problem until the termination criterion is satisfied.
///
/// Simply step()s until the termination criterion returns true.
/// Should be ideal for threading (using boost::thread).
void GeneticAlgorithm::evolve(unsigned int max_generations)
{
    _lock.lock();
    _isRunning = true;
    _lock.unlock();
    while (! terminate() && isRunning() && (_generations < max_generations))
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

/// Advance the algorithm one generation.
///
void GeneticAlgorithm::step()
{
    crossover();
    mutate();
    evaluateIndividuals();
    elitism();

    _generations += 1;
    updateStatistics();
    return;
}

/// Stop the algorithm.
/// May be used both from the outside and inside of the algorithm.
void GeneticAlgorithm::stop()
{
    boost::lock_guard<boost::mutex> guard(_lock);
    _isRunning = false;
    return;
}

/// Test if the algorithm is running.
///
/// @returns true if the algorithm is running.
bool GeneticAlgorithm::isRunning() const
{
    boost::lock_guard<boost::mutex> guard(_lock);
    return _isRunning;
}

/// Test whether the algorithm was stopped by the termination criterion.
///
/// @returns true if the algorithm was stopped by the termination criterion.
bool GeneticAlgorithm::isDone() const
{
    boost::lock_guard<boost::mutex> guard(_lock);
    return _isDone;
}

/// Get the best solution.
///
SolutionBase* GeneticAlgorithm::bestSolution() const
{
    boost::lock_guard<boost::mutex> guard(_lock);
    return _best_individual;
}

/// Get the number of evolved generations.
unsigned int GeneticAlgorithm::generations() const
{
    boost::lock_guard<boost::mutex> guard(_lock);
    return _generations;
}

/// Writes the statistics to the file
///
void GeneticAlgorithm::flushStatistics()
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

void GeneticAlgorithm::crossover()
{
    // Generator to decide whether to crossover
    std::uniform_real_distribution<double> probability(0,1);

    // Make a roulette wheel
    std::vector<double> roulette_wheel;
    roulette_wheel.reserve(_population.individuals.size());
    double max_score = *std::max_element(_costs.begin(), _costs.end());
    double total_score = 0.0;
    for (auto c : _costs)
    {
        total_score += -log(c/max_score);
        roulette_wheel.push_back(total_score);
    }
    std::uniform_real_distribution<double> roulette_ball(0, total_score);
    
    // Make children
    Population children;
    children.individuals.reserve(_population.individuals.size());
    for (size_t idx_child = 0; idx_child < _population.individuals.size(); ++idx_child)
    {
        // Choose dad
        size_t idx_dad = 0;
        double ball_value = roulette_ball(_random_generator);
        for (double r : roulette_wheel)
        {
            if (r >= ball_value)
                break;
            else
                ++idx_dad;
        }
        // We found the dad, now choose to crossover or else just copy him.
        if (probability(_random_generator) < _crossover_rate) 
        {
            // Choose mom
            size_t idx_mom = 0;
            ball_value = roulette_ball(_random_generator);
            for (double r : roulette_wheel)
            {
                if (r >= ball_value)
                    break;
                else
                    ++idx_mom;
            }
            SolutionBase* kid = _problem->crossover(_population.individuals[idx_dad], _population.individuals[idx_mom]);
            children.individuals.push_back(kid);
        }
        else
        {
            SolutionBase* clone = _population.individuals[idx_dad]->clone();
            children.individuals.push_back(clone);
        }
    }
    // Swap the children into place
    _population.individuals.swap(children.individuals);
}

void GeneticAlgorithm::mutate()
{
    // Make a probability generator to decide whether to mutate.
    std::uniform_real_distribution<double> probability(0, 1);
    // Mutate each individual.
    for (auto& sol : _population.individuals)
    {
        bool doMutate = probability(_random_generator) < _mutate_rate;
        if (doMutate)
        {
            _problem->mutate(sol);
        }
    }
    return;
}

/// Termination Criterion
///
bool GeneticAlgorithm::terminate()
{
        return false;
}

/// Evaluate the individuals and update the all-time high.
///
void GeneticAlgorithm::evaluateIndividuals()
{
    std::vector<double> new_costs;
    for (auto ind : _population.individuals)
    {
        double cost = _problem->cost(ind);
        new_costs.push_back(cost);
        if (cost < _best_cost)
        {
            _best_cost = cost;
            delete _best_individual;
            _best_individual = ind->clone();
        }
    }
    _costs.swap(new_costs);
    return;
}

/// Replace worst solution with the best if there was no improvement.
///
void GeneticAlgorithm::elitism()
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
        // Destroy the worst individual, he is no longer needed.
        delete _population.individuals[worst_individual];
        // Replace him with the best of the best.
        _population.individuals[worst_individual] = _best_individual->clone();
        _costs[worst_individual] = _best_cost;
    }
}

/// Create the statistics objects and an empty file to dump it into.
///
void GeneticAlgorithm::createStatistics()
{
    _statistics.reserve(100);
    std::ofstream myfile;
    myfile.open(_statistics_file, std::ios::trunc);
    myfile.close();
}

/// Notes the costs of the individuals and dumps them into a file once in a while.
///
void GeneticAlgorithm::updateStatistics()
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
