/*
Euclidean 3D Problem
Problem for the GATSP using euclidean distance measure as cost function.
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

#ifndef GATSP_PATH_REPRESENTATION_H
#define GATSP_PATH_REPRESENTATION_H


#include <random>
#include <vector>
#include <iostream>
#include "gatsp/genetic_algorithm.h"

template <class T>
struct PathSolution : public SolutionBase
{
    PathSolution() 
    {};

    PathSolution(
        const typename std::vector<T>::const_iterator begin, 
        const typename std::vector<T>::const_iterator end
    ) :
        genome(begin, end)
    {};

    virtual ~PathSolution() 
    {};

    PathSolution<T>* clone()
    {
        PathSolution<T>* clone = new PathSolution<T>;
        clone->genome = this->genome;
        return clone;
    };

    std::vector<T> genome;

    /// Exchange mutation. This mutation chooses two entries and swaps the
    /// content.
    ///
    /// @param random_generator A Mersenne Twister generator.
    void exchangeMutation(std::mt19937 & random_generator)
    {
        size_t solution_length = genome.size();
        std::uniform_int_distribution<size_t> dist(0, solution_length-1);
        size_t index1 = dist(random_generator);
        size_t index2 = dist(random_generator);
        while (index1 == index2)
        {
            index2 = dist(random_generator);
        }
        auto section_begin = genome.begin() + index1;
        auto section_end = genome.begin() + index2;
        // Invalidate costs
        section_begin->cost = -1;
        section_end->cost = -1;
        if (section_begin+1 != genome.end())
            (section_begin + 1)->cost = -1;
        else
            genome.begin()->cost = -1;
        if (section_end+1 != genome.end())
            (section_end + 1)->cost = -1;
        else
            genome.begin()->cost = -1;
        std::swap(*section_begin, *section_end);
        return;
    }

    /// Displacement mutation. Selects a range of indices and "slides" them
    /// back or forth in the solution. One special case is where the 
    /// insertion_point is where the range of indices already are located.
    /// In that case, the vector is not mutated.
    ///
    void displacementMutation(std::mt19937 & random_generator)
    {
        // You can think of the displacement mutation as 'sliding' a subtour
        // a number of indices along the original tour. But you can also think 
        // of it as exchanging two sections. This is what we do here. We exchange
        // the subtour and what we call the affected region.
        size_t solution_length = genome.size();
        std::uniform_int_distribution<size_t> dist_subtour_index(0, solution_length-1);
        size_t subtour_index = dist_subtour_index(random_generator);
        std::uniform_int_distribution<size_t> dist_subtour_length(1, solution_length - subtour_index);
        size_t subtour_length = dist_subtour_length(random_generator);
        std::uniform_int_distribution<size_t> dist_insertion_point(0, solution_length - subtour_length);
        size_t insertion_point = dist_insertion_point(random_generator);

        auto subtour_begin = genome.begin() + subtour_index;
        auto subtour_end = subtour_begin + subtour_length;

        // Invalidate subtour cut cost
        subtour_begin->cost = -1;
        if (subtour_end != genome.end())
            subtour_end->cost = -1;
        else
            genome.begin()->cost = -1;

        if (insertion_point > subtour_index)
        // Here the subtour is picked out in a temporary vector while the
        // affected region is moved.
        {
            PathSolution<T> subtour(subtour_begin, subtour_end);
            auto affected_portion_begin = subtour_end;
            auto affected_portion_end =  genome.begin() 
                                       + insertion_point
                                       + subtour_length;
            // Also invalidate affected portion's begining
            affected_portion_begin->cost = -1;
            auto to = subtour_begin;
            for (auto from = affected_portion_begin;
                 from != affected_portion_end;
                 ++from)
            {
                *to = *from;
                ++to;
            }
            for (auto from = subtour.genome.begin();
                 from != subtour.genome.end();
                 ++from)
            {
                *to = *from;
                ++to;
            }
        } 
        else if (insertion_point < subtour_index)
        // Here the affected region is picked out instead.
        {
            auto affected_portion_begin = genome.begin() + insertion_point;
            auto affected_portion_end = subtour_begin;
            // Also invalidate affected portion's begining
            affected_portion_begin->cost = -1;
            PathSolution<T> affected_portion(affected_portion_begin,
                                                 affected_portion_end);

            auto to = affected_portion_begin;
            for (auto from = subtour_begin;
                 from != subtour_end;
                 ++from)
            {
                *to = *from;
                ++to;
            }
            for (auto from = affected_portion.genome.begin();
                 from != affected_portion.genome.end();
                 ++from)
            {
                *to = *from;
                ++to;
            }
        }
        return;
    }

    /// Inversion mutation. This mutation chooses a range in the list and
    /// reverses the sequence in that range.
    ///
    void inversionMutation(std::mt19937 & random_generator)
    {
        // Coose two indices index2 is > than index1
        size_t solution_length = genome.size();
        std::uniform_int_distribution<size_t> dist(0, solution_length-1);
        size_t index1 = dist(random_generator);
        size_t index2 = dist(random_generator);
        while (index1 == index2)
        {
            index2 = dist(random_generator);
        }
        if (index1 > index2)
        {
            std::swap(index1, index2);
        }
        auto section_begin = genome.begin() + index1;
        auto section_end = genome.begin() + index2;

        // Invalidate the entry after the section
        if (section_end+1 != genome.end())
            (section_end + 1)->cost = -1;
        else
            genome.begin()->cost = -1;

        // Swap the entries in the section
        size_t num_swaps = (index2 - index1) / 2 + 1;
        for (size_t i = 0; i < num_swaps; ++i)
        {
            // Invalidate each entry in the section
            section_begin->cost = -1;
            section_end->cost = -1;

            std::swap(*section_begin, *section_end);
            ++section_begin;
            --section_end;
        }
        return;
    }

    /// Order crossover. Crosses two parents to an offspring, 
    /// trying to preserve the order of the cities.
    static PathSolution<T>* orderCrossover(const PathSolution<T>* dad, const PathSolution<T>* mom, std::mt19937 & random_generator)
    {
        // Create an empty kid
        PathSolution<T>* kid = new PathSolution<T>;
        size_t long_length = dad->genome.size();
        size_t short_length = mom->genome.size();
        if (short_length > long_length)
            std::swap(short_length, long_length);
        kid->genome.reserve(long_length);

        // Choose two cut points
        std::uniform_int_distribution<size_t> dist(0, short_length-1);
        size_t index1 = dist(random_generator);
        size_t index2 = dist(random_generator);
        while (index1 == index2)
            index2 = dist(random_generator);
        if (index1 > index2)
            std::swap(index1, index2);

        // Create a subtour from mom, and invalidate the first entry.
        PathSolution<T> subtour(mom->genome.begin()+index1, mom->genome.begin()+index2);
        subtour.genome.begin()->cost = -1;
        // Traverse the dad solution, copying entries != subtour.
        // and when the time is right, copy the subtour to kid in the
        // position it was cut.
        // To get some performance, we'll make a sorted subtour and use
        // binary search.
        PathSolution<T> sorted_subtour(mom->genome.begin()+index1, mom->genome.begin()+index2);
        std::sort(sorted_subtour.genome.begin(), sorted_subtour.genome.end());
        int i = 0;
        bool did_break = true;
        for (const auto& entry : dad->genome)
        {
            if (! std::binary_search(sorted_subtour.genome.begin(), sorted_subtour.genome.end(), entry))
            //if ( std::none_of(subtour.genome.begin(), subtour.genome.end(), [&](const Euclidean3DEntry& e){return e.index == entry.index;}) )
            {
                kid->genome.push_back(entry);
                // If the previous entry was broken, we'll invalidate the cost.
                if (did_break)
                {
                    kid->genome.rbegin()->cost = -1;
                    did_break = false;
                }
            }
            else
                did_break = true;
            if (i == index1)
            {
                kid->genome.insert(kid->genome.end(), subtour.genome.begin(), subtour.genome.end());
                did_break = true;
            }
            i += 1;
        }
        return kid;
    }

};

#endif