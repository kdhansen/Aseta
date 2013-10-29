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

#ifndef GENETIC_ALGORITHM_BASE
#define GENETIC_ALGORITHM_BASE

#include <memory>
#include <random>
#include <string>
#include <boost/thread.hpp>
#include "gatsp/problem_base.h"
#include "gatsp/solution_base.h"

class GeneticAlgorithmBase
{
public:
	GeneticAlgorithmBase(
		const ProblemBase&,
	    int num_individuals, 
	    double mutate_rate, 
	    double crossover_rate, 
	    unsigned int seed,
	    std::string statistics_file
	);

	GeneticAlgorithmBase(
		const ProblemBase&,
		SolutionBase,
	    int num_individuals, 
	    double mutate_rate, 
	    double crossover_rate, 
	    unsigned int seed,
	    std::string statistics_file
	);

	GeneticAlgorithmBase(const GeneticAlgorithmBase&);

	GeneticAlgorithmBase& operator=(const GeneticAlgorithmBase&);

	virtual ~GeneticAlgorithmBase();

	void evolve(unsigned int);
	void stop();
	bool isRunning();
	bool isDone();
	virtual void step();
	SolutionBase bestSolution();
	unsigned int generations() const;
	void flushStatistics();

private:
	GeneticAlgorithmBase();

	virtual void mutate() = 0;
	virtual void crossover() = 0;
	virtual bool terminate() = 0;

	void evaluateIndividuals();
	void elitism();
	void createStatistics();
	void updateStatistics();

protected:
	const ProblemBase* _problem;
	
	int	_num_individuals;
	unsigned int _generations;
	std::vector<SolutionBase> _individuals;
	std::vector<double> _costs;
	SolutionBase _best_individual;
	double _best_cost;
	bool _do_statistics;
	std::string _statistics_file;
	std::vector<std::vector<double> > _statistics;

	double _mutate_rate;
	double _crossover_rate;

	std::mt19937 _random_generator;

	boost::mutex _lock;
	bool _isRunning;
	bool _isDone;
};

#endif // GENETIC_ALGORITHM_BASE