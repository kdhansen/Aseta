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

#ifndef GATSP_GENETIC_ALGORITHM_H
#define GATSP_GENETIC_ALGORITHM_H

#include <exception>
#include <random>
#include <boost/thread.hpp>

struct SolutionBase
{
    virtual ~SolutionBase() {};
    virtual SolutionBase* clone() = 0;
};

class InvalidSolution : public std::exception
{
public:
    InvalidSolution(std::string what = "") :
        _what(what)
    {}
    virtual const char* what() const throw()
    {
        if (_what == "")
            return "Solution invalid.";
        else
            return _what.c_str();
    }
private:
    std::string _what;
};

class ProblemBase
{
public:
    virtual ~ProblemBase() {};

    virtual SolutionBase* makeSolution() const = 0;
    virtual void mutate(SolutionBase* mutatee) const = 0;
    virtual SolutionBase* crossover(const SolutionBase* dad, const SolutionBase* mom) const = 0;
    virtual double cost(const SolutionBase* solution) const = 0;
};

struct Population
{
    ~Population()
    {
        for (SolutionBase* ind: individuals)
            delete ind;
    }
    std::vector<SolutionBase*> individuals;
};

class GeneticAlgorithm
{
public:
    GeneticAlgorithm(
        const ProblemBase& problem,
        int num_individuals = 100,
        double mutate_rate = 0.3,
        double crossover_rate = 0.9,
        unsigned int seed = 0,
        std::string statistics_file = ""
    );

    ~GeneticAlgorithm() {};

    void evolve(unsigned int num_generations);
    void step();
    void stop();
    bool isRunning() const;
    bool isDone() const;
    SolutionBase* bestSolution() const;
    unsigned int generations() const;
    void flushStatistics();

private:
    GeneticAlgorithm(const GeneticAlgorithm&);

    void mutate();
    void crossover();
    bool terminate();

    void evaluateIndividuals();
    void elitism();
    void createStatistics();
    void updateStatistics();

    const ProblemBase* _problem;
    
    int _num_individuals;
    unsigned int _generations;
    Population _population;
    std::vector<double> _costs;
    SolutionBase* _best_individual;
    double _best_cost;
    bool _do_statistics;
    std::string _statistics_file;
    std::vector<std::vector<double> > _statistics;

    double _mutate_rate;
    double _crossover_rate;

    std::mt19937 _random_generator;

    mutable boost::mutex _lock;
    bool _isRunning;
    bool _isDone;
};

#endif // GATSP_GENETIC_ALGORITHM_H
