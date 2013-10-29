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

#include "gatsp/traditional_genetic_algorithm.h"

TraditionalGeneticAlgorithm::TraditionalGeneticAlgorithm(
    const ProblemBase& problem,
    int num_individuals, 
    double mutate_rate, 
    double crossover_rate, 
    unsigned int seed
) 
    : GeneticAlgorithmBase(problem, num_individuals, mutate_rate, crossover_rate, seed)
{}

TraditionalGeneticAlgorithm::TraditionalGeneticAlgorithm(
    const ProblemBase& problem,
    SolutionBase solution,
    int num_individuals, 
    double mutate_rate, 
    double crossover_rate, 
    unsigned int seed
) 
    : GeneticAlgorithmBase(problem, solution, num_individuals, mutate_rate, crossover_rate, seed)
{}

TraditionalGeneticAlgorithm::~TraditionalGeneticAlgorithm()
{}

void TraditionalGeneticAlgorithm::mutate()
{
    // Make a probability generator to decide whether to mutate.
    std::uniform_real_distribution<double> probability(0, 1);
    bool doMutate = (probability(_random_generator) < _mutate_rate);

    for (auto& sol : _individuals)
    {
        if ((2 < sol.size()) && doMutate)
        {
            // Choose mutation.
            std::uniform_int_distribution<unsigned int> which_mutation(0, 2);
            switch (which_mutation(_random_generator))
            {
                case 0:
                    displacementMutation(sol);
                    break;
                case 1:
                    exchangeMutation(sol);
                    break;
                case 2:
                    inversionMutation(sol);
                    break;
            }
        }
    }
    return;
}

void TraditionalGeneticAlgorithm::crossover()
{
    // Decide whether to crossover
    std::uniform_real_distribution<double> probability(0,1);
    bool doCrossover = (probability(_random_generator) < _crossover_rate);
    if(! doCrossover)
    {
        return;
    }

    // Make a roulette wheel
    std::vector<double> roulette_wheel;
    roulette_wheel.reserve(_individuals.size());
    double total_score = 0.0;
    for (auto c : _costs)
    {
        total_score += 1.0/c;
        roulette_wheel.push_back(total_score);
    }
    
    // Make children
    std::vector<SolutionBase> children;
    children.reserve(_individuals.size());
    for (size_t idx_child = 0; idx_child < _individuals.size(); ++idx_child)
    {
        std::uniform_real_distribution<double> roulette_ball(0, total_score);
        double ball_value = roulette_ball(_random_generator);
        size_t idx_parent = 0;
        for (double r : roulette_wheel)
        {
            if (r >= ball_value)
            {
                break;
            }
            else
            {
                ++idx_parent;
            }
        }
        SolutionBase s(_individuals[idx_parent]);
        children.push_back(s);
    }

    // Swap the children into place
    _individuals.swap(children);
}

/// Termination Criterion
///
bool TraditionalGeneticAlgorithm::terminate()
{
        return false;
}

/// Displacement mutation. Selects a range of indices and "slides" them
/// back or forth in the solution. One special case is where the 
/// insertion_point is where the range of indices already are located.
/// In that case, the vector is not mutated.
///
void TraditionalGeneticAlgorithm::displacementMutation(SolutionBase& s)
{
    size_t solution_length = s.size();
    std::uniform_int_distribution<size_t> dist_subtour_index(0, solution_length-1);
    size_t subtour_index = dist_subtour_index(_random_generator);
    std::uniform_int_distribution<size_t> dist_subtour_length(1, solution_length - subtour_index);
    size_t subtour_length = dist_subtour_length(_random_generator);
    std::uniform_int_distribution<size_t> dist_insertion_point(0, solution_length - subtour_length);
    size_t insertion_point = dist_insertion_point(_random_generator);

    auto subtour_begin = s.begin() + subtour_index;
    auto subtour_end = subtour_begin + subtour_length;      

    if (insertion_point > subtour_index)
    // Here the subtour is picked out in a temporary vector while the
    // affected region is moved.
    {
        SolutionBase subtour(subtour_begin, subtour_end);
        auto affected_portion_begin = subtour_end;
        auto affected_portion_end =  s.begin() 
                                   + insertion_point
                                   + subtour_length;
        auto to = subtour_begin;
        for (auto from = affected_portion_begin;
             from != affected_portion_end;
             ++from)
        {
            *to = *from;
            ++to;
        }
        for (auto from = subtour.begin();
             from != subtour.end();
             ++from)
        {
            *to = *from;
            ++to;
        }
    } 
    else if (insertion_point < subtour_index)
    // Here the affected region is picked out instead.
    {
        auto affected_portion_begin = s.begin() + insertion_point;
        auto affected_portion_end = subtour_begin;
        SolutionBase affected_portion(affected_portion_begin,
                                             affected_portion_end);

        auto to = affected_portion_begin;
        for (auto from = subtour_begin;
             from != subtour_end;
             ++from)
        {
            *to = *from;
            ++to;
        }
        for (auto from = affected_portion.begin();
             from != affected_portion.end();
             ++from)
        {
            *to = *from;
            ++to;
        }
    }
    return;
}

/// Exchange mutation. This mutation chooses two indices and swaps the
/// content.
///
void TraditionalGeneticAlgorithm::exchangeMutation(SolutionBase& s)
{
    size_t solution_length = s.size();
    std::uniform_int_distribution<size_t> dist(0, solution_length-1);
    size_t index1 = dist(_random_generator);
    size_t index2 = dist(_random_generator);
    while (index1 == index2)
    {
        index2 = dist(_random_generator);
    }
    auto section_begin = s.begin() + index1;
    auto section_end = s.begin() + index2;
    std::swap(*section_begin, *section_end);
    return;
}

/// Inversion mutation. This mutation chooses a range in the list and
/// reverses the sequence in that range.
///
void TraditionalGeneticAlgorithm::inversionMutation(SolutionBase& s)
{
    size_t solution_length = s.size();
    std::uniform_int_distribution<size_t> dist(0, solution_length-1);
    size_t index1 = dist(_random_generator);
    size_t index2 = dist(_random_generator);
    while (index1 == index2)
    {
        index2 = dist(_random_generator);
    }
    if (index1 > index2)
    {
        std::swap(index1, index2);
    }
    auto section_begin = s.begin() + index1;
    auto section_end = s.begin() + index2;
    size_t num_swaps = (index2 - index1) / 2 + 1;
    for (size_t i = 0; i < num_swaps; ++i)
    {
        std::swap(*section_begin, *section_end);
        ++section_begin;
        --section_end;
    }
    return;
}
