/*
Traditional Genetic Algorithm
GATSP Genetic Algorithm using inversion, displacement and exchange mutations.
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

#include "gatsp/genetic_algorithm_base.h"

class TraditionalGeneticAlgorithm : public GeneticAlgorithmBase
{
public:
    TraditionalGeneticAlgorithm(
        const ProblemBase&,
        int num_individuals = 100, 
        double mutate_rate = 1.0, 
        double crossover_rate = 1.0,
        unsigned int seed = 0,
        std::string statistics_file = ""
    );

    TraditionalGeneticAlgorithm(
        const ProblemBase&,
        SolutionBase,
        int num_individuals = 100, 
        double mutate_rate = 1.0, 
        double crossover_rate = 1.0, 
        unsigned int seed = 0,
        std::string statistics_file = ""
    );
    
    virtual ~TraditionalGeneticAlgorithm();

private:
    TraditionalGeneticAlgorithm();

    virtual void mutate();
    virtual void crossover();
    virtual bool terminate();

    void displacementMutation(SolutionBase&);
    void exchangeMutation(SolutionBase&);
    void inversionMutation(SolutionBase&);
};
