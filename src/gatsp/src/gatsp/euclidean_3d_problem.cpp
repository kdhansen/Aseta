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

#include <algorithm>
#include <cmath>
#include "gatsp/euclidean_3d_problem.h"


Euclidean3DProblem::Euclidean3DProblem(unsigned int seed) :
    _random_generator(std::mt19937(seed))
{}

Euclidean3DProblem::~Euclidean3DProblem()
{}

void Euclidean3DProblem::mutate(SolutionBase* mutatee) const
{
    // Cast as an Euclidean 3D Solution, we can't handle others
    Euclidean3DSolution* sol = static_cast<Euclidean3DSolution*>(mutatee);

    if (2 < sol->genome.size())
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
    return;
}

SolutionBase* Euclidean3DProblem::crossover(const SolutionBase* dad, const SolutionBase* mom) const
{
    // Cast as Euclidean 3D Solutions, we can't handle others
    const Euclidean3DSolution* euc_dad = static_cast<const Euclidean3DSolution*>(dad);
    const Euclidean3DSolution* euc_mom = static_cast<const Euclidean3DSolution*>(mom);
    return orderCrossover(euc_dad, euc_mom);
}

/// Get the cost of a solution
///
/// Advanced specialization that only computes invalidated costs.
///
/// @param solution The solution to check
/// @returns the Euclidean length of the solution.
double Euclidean3DProblem::cost(const SolutionBase* solution) const
{
    const Euclidean3DSolution* euc_sol = static_cast<const Euclidean3DSolution*>(solution);
    double dist = 0.0;
    if (2 <= euc_sol->genome.size()) // Cannot compute distance between less than two wps.
    {
        // Traverse the solution
        for (auto entry = euc_sol->genome.begin()+1; entry!= euc_sol->genome.end(); ++entry)
        {
            // If the step has been invalidated, well compute the distance.
            if (entry->cost < 0.0)
            {
                double d = euclideanDistance(_waypoints[(entry-1)->index], _waypoints[entry->index]);
                entry->cost = d;
                dist += d;
            }
            else
                dist += entry->cost;
        }

        // Close the loop
        auto from = euc_sol->genome.rbegin();
        auto to = euc_sol->genome.begin();
        if (to->cost < 0.0)
        {
            double d = euclideanDistance(_waypoints[from->index], _waypoints[to->index]);
            to->cost = d;
            dist += d;
        }
        else
            dist += to->cost;
    }
    return dist;
}


/// Make a solution for the genetic algorithm.
///
SolutionBase* Euclidean3DProblem::makeSolution() const
{
    return makeNearestNeighbor();
}


/// Displacement mutation. Selects a range of indices and "slides" them
/// back or forth in the solution. One special case is where the 
/// insertion_point is where the range of indices already are located.
/// In that case, the vector is not mutated.
///
void Euclidean3DProblem::displacementMutation(Euclidean3DSolution* s) const
{
    // You can think of the displacement mutation as 'sliding' a subtour
    // a number of indices along the original tour. But you can also think 
    // of it as exchanging two sections. This is what we do here. We exchange
    // the subtour and what we call the affected region.
    size_t solution_length = s->genome.size();
    std::uniform_int_distribution<size_t> dist_subtour_index(0, solution_length-1);
    size_t subtour_index = dist_subtour_index(_random_generator);
    std::uniform_int_distribution<size_t> dist_subtour_length(1, solution_length - subtour_index);
    size_t subtour_length = dist_subtour_length(_random_generator);
    std::uniform_int_distribution<size_t> dist_insertion_point(0, solution_length - subtour_length);
    size_t insertion_point = dist_insertion_point(_random_generator);

    auto subtour_begin = s->genome.begin() + subtour_index;
    auto subtour_end = subtour_begin + subtour_length;      

    // Invalidate subtour cut cost
    subtour_begin->cost = -1;
    if (subtour_end != s->genome.end())
        subtour_end->cost = -1;
    else
        s->genome.begin()->cost = -1;

    if (insertion_point > subtour_index)
    // Here the subtour is picked out in a temporary vector while the
    // affected region is moved.
    {
        Euclidean3DSolution subtour(subtour_begin, subtour_end);
        auto affected_portion_begin = subtour_end;
        auto affected_portion_end =  s->genome.begin() 
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
        auto affected_portion_begin = s->genome.begin() + insertion_point;
        auto affected_portion_end = subtour_begin;
        // Also invalidate affected portion's begining
        affected_portion_begin->cost = -1;
        Euclidean3DSolution affected_portion(affected_portion_begin,
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

/// Exchange mutation. This mutation chooses two indices and swaps the
/// content.
///
void Euclidean3DProblem::exchangeMutation(Euclidean3DSolution* s) const
{
    size_t solution_length = s->genome.size();
    std::uniform_int_distribution<size_t> dist(0, solution_length-1);
    size_t index1 = dist(_random_generator);
    size_t index2 = dist(_random_generator);
    while (index1 == index2)
    {
        index2 = dist(_random_generator);
    }
    auto section_begin = s->genome.begin() + index1;
    auto section_end = s->genome.begin() + index2;
    // Invalidate costs
    section_begin->cost = -1;
    section_end->cost = -1;
    if (section_begin+1 != s->genome.end())
        (section_begin + 1)->cost = -1;
    else
        s->genome.begin()->cost = -1;
    if (section_end+1 != s->genome.end())
        (section_end + 1)->cost = -1;
    else
        s->genome.begin()->cost = -1;
    std::swap(*section_begin, *section_end);
    return;
}

/// Inversion mutation. This mutation chooses a range in the list and
/// reverses the sequence in that range.
///
void Euclidean3DProblem::inversionMutation(Euclidean3DSolution* s) const
{
    // Coose two indices index2 is > than index1
    size_t solution_length = s->genome.size();
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
    auto section_begin = s->genome.begin() + index1;
    auto section_end = s->genome.begin() + index2;

    // Invalidate the entry after the section
    if (section_end+1 != s->genome.end())
        (section_end + 1)->cost = -1;
    else
        s->genome.begin()->cost = -1;

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
Euclidean3DSolution* Euclidean3DProblem::orderCrossover(const Euclidean3DSolution* dad, const Euclidean3DSolution* mom) const
{
    // Create an empty kid
    Euclidean3DSolution* kid = new Euclidean3DSolution;
    size_t solution_length = dad->genome.size();
    kid->genome.reserve(solution_length);

    // Choose two cut points
    std::uniform_int_distribution<size_t> dist(0, solution_length-1);
    size_t index1 = dist(_random_generator);
    size_t index2 = dist(_random_generator);
    while (index1 == index2)
        index2 = dist(_random_generator);
    if (index1 > index2)
        std::swap(index1, index2);

    // Create a subtour from mom, and invalidate the first entry.
    Euclidean3DSolution subtour(mom->genome.begin()+index1, mom->genome.begin()+index2);
    subtour.genome.begin()->cost = -1;
    // Traverse the dad solution, copying entries != subtour.
    // and when the time is right, copy the subtour to kid in the
    // position it was cut.
    // To get some performance, we'll make a sorted subtour and use
    // binary search.
    Euclidean3DSolution sorted_subtour(mom->genome.begin()+index1, mom->genome.begin()+index2);
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

/// Make a quick and feasible solution (but pretty bad probably).
///
Euclidean3DSolution* Euclidean3DProblem::makeTrivialSolution() const
{
    Euclidean3DSolution* solution = new Euclidean3DSolution;
    solution->genome.reserve(_waypoints.size());
    for(size_t i = 0; i < _waypoints.size(); ++i)
    {
        solution->genome.push_back(i);
    }
    return solution;
}

/// Make a nearest neighbor solution
///
Euclidean3DSolution* Euclidean3DProblem::makeNearestNeighbor() const
{
    // Prepare resulting solution
    Euclidean3DSolution* solution = new Euclidean3DSolution;
    solution->genome.reserve(_waypoints.size());
    // Get the unvisited points
    Euclidean3DSolution* unvisited = makeTrivialSolution();
    // Start with a random first point.
    std::uniform_int_distribution<size_t> dist(0, _waypoints.size()-1);
    size_t random_index = dist(_random_generator);
    auto from = unvisited->genome.begin() + random_index;
    solution->genome.push_back(*from);
    unvisited->genome.erase(from);
    while (unvisited->genome.size() > 0)
    {
        Euclidean3DEntry from = solution->genome.back();
        auto to = std::min_element(unvisited->genome.begin(), unvisited->genome.end(), 
            [&](Euclidean3DEntry i, Euclidean3DEntry j){return euclideanDistance(_waypoints[from.index], _waypoints[i.index]) < euclideanDistance(_waypoints[from.index], _waypoints[j.index]);}
        );
        solution->genome.push_back(*to);
        unvisited->genome.erase(to);
    }
    // Clean up
    delete unvisited;

    return solution;
}

/// Compile a route of waypoints from the given solution.
///
/// @param solution The solution that you want to use for the route.
/// @returns A feasible route
std::vector<Waypoint> Euclidean3DProblem::route(const SolutionBase* solution) throw(InvalidSolution)
{
    const Euclidean3DSolution* euc_solution = dynamic_cast<const Euclidean3DSolution*>(solution);
    if (0 == euc_solution)
        throw InvalidSolution("Solution not Euclidean3DSolution.");
    if (! isSolutionValid(*euc_solution))
        throw InvalidSolution();

    std::vector<Waypoint> output;
    output.reserve(_waypoints.size());
    for (const auto& sol_entry : euc_solution->genome)
    {
        output.push_back(_waypoints[sol_entry.index]);
    }
    return output;
}

/// Report the first waypoint in a sequence given a specific solution.
///
/// @param solution The solution that determines the first waypoint in the sequence.
/// @returns The first waypoint in the sequence.
Waypoint Euclidean3DProblem::firstWaypoint(const SolutionBase* solution) throw(InvalidSolution)
{
    const Euclidean3DSolution* euc_solution = dynamic_cast<const Euclidean3DSolution*>(solution);
    if (0 == euc_solution)
        throw InvalidSolution("Solution not Euclidean3DSolution.");

    if (euc_solution->genome.size() > 0 && isSolutionValid(*euc_solution))
        return _waypoints[euc_solution->genome[0].index];
    else
        throw InvalidSolution();
}

/// Remove the first waypoint in the sequence given a solution.
///
/// @param [in,out] solution The solution that determines the first waypoint, is also updated.
void Euclidean3DProblem::popWaypointFront(Euclidean3DSolution& solution) throw(InvalidSolution)
{   
    if (! isSolutionValid(solution))
        throw InvalidSolution();

    Euclidean3DEntry remove_entry = solution.genome[0];

    // Remove from the waypoint list
    _waypoints.erase(_waypoints.begin()+remove_entry.index);
    
    // Remove from the solution (erase-remove idiom)
    solution.genome.erase(std::remove(solution.genome.begin(), solution.genome.end(), remove_entry), solution.genome.end());

    // now the indexes of the entries > solution[0].index must be reduced by 1.
    for(auto& entry: solution.genome)
    {
        if (entry.index > remove_entry.index)
        {
            entry.index -= 1;
        }
    }

    // And invalidate the cost of the new front
    if (0 < solution.genome.size())
    solution.genome[0].cost = -1;
}

/// Add a waypoint to the problem.
///
/// @param waypoint Waypoint to add.
void Euclidean3DProblem::addWaypoint(const Waypoint& waypoint)
{
    _waypoints.push_back(waypoint);
}

// /// Add a waypoint to the problem and update related solutions
// ///
// /// @param waypoint Waypoint to add.
// /// @param [out] old_solution A solution that must be updated with the new waypoint.
// void Euclidean3DProblem::addWaypoint(const Waypoint& waypoint, SolutionBase& old_solution)
// {
//     addWaypoint(waypoint);
//     old_solution.push_back(_waypoints.size()-1);
// }

/// Check to see if solution is valid.
///
/// This is a simple check for validity, not wheter the solution is 
/// satisfactory.
///
/// @param solution The solution to check.
/// @returns True if the solution is valid.
bool Euclidean3DProblem::isSolutionValid(const Euclidean3DSolution& solution) const
{
    bool valid = true;
    for (const auto& sol_entry : solution.genome)
    {
        if (sol_entry.index > _waypoints.size() || sol_entry.index < 0)
        {
            valid = false;
            break;
        }
    }
    return valid;
}

// /// Attempt to repair solution.
// ///
// /// @param solution The solution to fix.
// /// @returns True if the solution could be repaired.
// bool Euclidean3DProblem::repairSolution(const SolutionBase& solution) const
// {
//     return true;
// }

/// Compute Euclidean distance between two waypoints.
///
/// @param from Where to go from.
/// @param to Where to go to.
/// @returns The distance between the waypoints.
double Euclidean3DProblem::euclideanDistance(const Waypoint& from, const Waypoint& to) const
{
    double dx = to.x() - from.x();
    double dy = to.y() - from.y();
    double dz = to.z() - from.z();
    return sqrt(dx*dx + dy*dy + dz*dz);
}
