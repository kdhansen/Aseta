/*
Euclidean 3D Problem with Refueling.
Problem for the GATSP using euclidean distance measure as cost function
with a vehicle using fuel, so it has to refuel once in a while.
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

#include "gatsp/refueling_problem.h"

RefuelingProblem::RefuelingProblem(unsigned int seed) :
    _random_generator(std::mt19937(seed))
{}

RefuelingProblem::~RefuelingProblem()
{}

/// Add a waypoint to the problem.
///
/// @param waypoint Waypoint to add.
void RefuelingProblem::addWaypoint(const Waypoint& waypoint)
{
    _waypoints.push_back(waypoint);
    _isDepot.push_back(false);
}

/// Add a depot to the problem.
///
/// @param depot Depot to add.
void RefuelingProblem::addDepot(const Waypoint& depot)
{
    _waypoints.push_back(depot);
    _isDepot.push_back(true);
    _depotIndices.push_back(_waypoints.size()-1);
}

/// Set the capacity of the fuel tank
///
/// @param metersPrTank How far can the helicopter fly on a tankfull.
void RefuelingProblem::fuelCapacity(double metersPrTank)
{
    _metersPrTank = metersPrTank;
}

/// Implementation of mutate needed by the genetic algorthm.
///
void RefuelingProblem::mutate(SolutionBase* mutatee) const
{
    // Cast as a RefuelingSolution, we can't handle others
    RefuelingSolution* sol = dynamic_cast<RefuelingSolution*>(mutatee);
    if (2 < sol->genome.size())
    {
        // Choose mutation.
        std::uniform_real_distribution<double> which_mutation(0.0, 1.0);
        double which = which_mutation(_random_generator);
        // if (which < 0.3)
        // {
        //     sol->displacementMutation(_random_generator);
        // }
        //else 
        // if (which < 0.9)
        // {
        //     sol->exchangeMutation(_random_generator);
        // }
        //else 
        if (which < 0.9)
        {
            sol->inversionMutation(_random_generator);
        }
        else if (which < 0.95)
        {
            injectDepot(sol);
        }
        else if (which <= 1.0)
        {
            extractDepot(sol);
        }
    }

    return;
}

/// Inject a fuel depot into the solution.
/// Problem specific mutation. It is complemented by extractDepot.
///
/// @param[in,out] mutatee Solution to mutate.
void RefuelingProblem::injectDepot(RefuelingSolution* mutatee) const
{
    size_t solution_length = mutatee->genome.size();
    std::uniform_int_distribution<size_t> dist(0, solution_length);
    auto inject_point = mutatee->genome.begin() + dist(_random_generator);

    size_t num_depots = _depotIndices.size();
    dist = std::uniform_int_distribution<size_t>(0, num_depots-1);
    size_t depot_index = _depotIndices[dist(_random_generator)];

    // Invalidate cost of nex entry
    if (inject_point != mutatee->genome.end())
        (inject_point)->cost = -1;
    else
        mutatee->genome.begin()->cost = -1;

    // Insert depot entry
    RefuelingEntry e(depot_index);
    mutatee->genome.insert(inject_point, e);

    return;
}

/// Extract a fuel depot into the solution.
/// Problem specific mutation. It is complemented by injectDepot.
/// It chooses an index and removes the first comming depot.
///
/// @param[in,out] mutatee Solution to mutate.
void RefuelingProblem::extractDepot(RefuelingSolution* mutatee) const
{
    size_t solution_length = mutatee->genome.size();
    std::uniform_int_distribution<size_t> dist(0, solution_length-1);
    auto extract_point = mutatee->genome.begin() + dist(_random_generator);
    // The stop_at trick will wrap around the vector starting from the 
    // extract_point (and also stopping at it.)
    auto stop_at = mutatee->genome.end();
    for (auto iter = extract_point; iter != stop_at; ++iter)
    {
        if (_isDepot[iter->index])
        {
            // Invalidate entry cost after depot
            if (iter + 1 != mutatee->genome.end())
                (iter + 1)->cost = -1;
            else
                mutatee->genome.begin()->cost = -1;

            // Remove depot
            mutatee->genome.erase(iter);
            break;
        }

        if (iter + 1 == mutatee->genome.end())
        {
            stop_at = extract_point;
            iter = mutatee->genome.begin() - 1;
        }
    }

    return;
}

/// Implementation of crossover needed by the genetic algorthm.
///
SolutionBase* RefuelingProblem::crossover(const SolutionBase* dad, const SolutionBase* mom) const
{
    // Cast as Euclidean 3D Solutions, we can't handle others
    const RefuelingSolution* ref_dad = static_cast<const RefuelingSolution*>(dad);
    const RefuelingSolution* ref_mom = static_cast<const RefuelingSolution*>(mom);

    // Create an empty kid
    RefuelingSolution* kid = new RefuelingSolution;
    size_t long_length = ref_dad->genome.size();
    size_t short_length = ref_mom->genome.size();
    if (short_length > long_length)
        std::swap(short_length, long_length);
    kid->genome.reserve(long_length);

    // Choose two cut points
    std::uniform_int_distribution<size_t> dist(0, short_length-1);
    size_t index1 = dist(_random_generator);
    size_t index2 = dist(_random_generator);
    while (index1 == index2)
        index2 = dist(_random_generator);
    if (index1 > index2)
        std::swap(index1, index2);

    // Create a subtour from mom, and invalidate the first entry.
    RefuelingSolution subtour(ref_mom->genome.begin()+index1, ref_mom->genome.begin()+index2);
    subtour.genome.begin()->cost = -1;
    // Traverse the dad solution, copying entries != subtour.
    // and when the time is right, copy the subtour to kid in the
    // position it was cut.
    // To get some performance, we'll make a sorted subtour and use
    // binary search.
    RefuelingSolution sorted_subtour(ref_mom->genome.begin()+index1, ref_mom->genome.begin()+index2);
    std::sort(sorted_subtour.genome.begin(), sorted_subtour.genome.end());
    int i = 0;
    bool did_break = true;
    for (const auto& entry : ref_dad->genome)
    {
        // We'll copy the entry into the kid if the entry was not in the subtour or it is a depot.
        if (! std::binary_search(sorted_subtour.genome.begin(), sorted_subtour.genome.end(), entry))
            // || _isDepot[entry.index])
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

/// Get the cost of a solution
///
/// Advanced specialization that only computes invalidated costs.
///
/// @param solution The solution to check
/// @returns the Euclidean length of the solution.
double RefuelingProblem::cost(const SolutionBase* solution) const
{
    const RefuelingSolution* sol = dynamic_cast<const RefuelingSolution*>(solution);

    double dist = 0.0;
    double maxFuel = 0.0; // find which entry that has the most fuelLeft.
    if (2 <= sol->genome.size()) // Cannot compute distance between less than two wps.
    {
        // Euclidean distance
        // Traverse the solution
        for (auto entry = sol->genome.begin()+1; entry != sol->genome.end(); ++entry)
        {
            // If the step has been invalidated, will compute the distance.
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
        auto from = sol->genome.rbegin();
        auto to = sol->genome.begin();
        if (to->cost < 0.0)
        {
            double d = euclideanDistance(_waypoints[from->index], _waypoints[to->index]);
            to->cost = d;
            dist += d;
        }
        else
            dist += to->cost;

        // Fuel usage
        // compute usage for each leg
        // cummulate fuel used between depots and set refuel amount
        auto last_depot = sol->genome.begin();
        double fuelUsed = 0.0;
        for (auto entry = sol->genome.begin() + 1; entry != sol->genome.end(); ++entry)
        {
            fuelUsed += entry->cost;
            if (_isDepot[entry->index])
            {
                last_depot->fuelLeft = fuelUsed;
                if (fuelUsed > maxFuel)
                    maxFuel = fuelUsed;
                fuelUsed = 0.0;
                last_depot = entry;
            }
            if (entry == sol->genome.end()-1) // we should count the first entry and refuel
            {
                fuelUsed += sol->genome.begin()->cost;
                last_depot->fuelLeft = fuelUsed;
                if (fuelUsed > maxFuel)
                    maxFuel = fuelUsed;
            }
        }
        // do another run to compute the fuel left
        for (auto entry = sol->genome.begin() + 1; entry != sol->genome.end(); ++entry)
        {
            if (! _isDepot[entry->index])
                entry->fuelLeft = (entry - 1)->fuelLeft - entry->cost;
        }
    }

    if (maxFuel > _metersPrTank)
        dist *= (maxFuel - _metersPrTank);

    return dist;
}

/// Compute Euclidean distance between two waypoints.
///
/// @param from Where to go from.
/// @param to Where to go to.
/// @returns The distance between the waypoints.
double RefuelingProblem::euclideanDistance(const Waypoint& from, const Waypoint& to) const
{
    double dx = to.x() - from.x();
    double dy = to.y() - from.y();
    double dz = to.z() - from.z();
    return sqrt(dx*dx + dy*dy + dz*dz);
}

/// Make a solution for the genetic algorithm.
/// It is possibly pretty bad. It is an ordered list of waypoints that
/// are not depots.
SolutionBase* RefuelingProblem::makeSolution() const
{
    RefuelingSolution* solution = new RefuelingSolution;
    solution->genome.reserve(_waypoints.size());
    for(size_t i = 0; i < _waypoints.size(); ++i)
    {
        if (! _isDepot[i])
            solution->genome.push_back(i);
    }
    return solution;
}

/// Get the waypoints of the problem.
///
/// This gets the waypoints (or targets) but not the depots.
/// See also RefuelingProblem::depots()
///
/// @returns a vector containing a copy of the waypoints.
std::vector<Waypoint> RefuelingProblem::waypoints()
{
    std::vector<Waypoint> result;
    for (size_t i = 0; i < _waypoints.size(); ++i)
    {
        if (! _isDepot[i])
            result.push_back(_waypoints[i]);
    }
    return result;
}

/// Get the depots of the problem.
///
/// This gets the refueling stations, not the targets.
/// See also RefuelingProblem::waypoints()
///
/// @returns a vector containing a copy of the depots.
std::vector<Waypoint> RefuelingProblem::depots()
{
    std::vector<Waypoint> result;
    for (size_t i = 0; i < _waypoints.size(); ++i)
    {
        if (_isDepot[i])
            result.push_back(_waypoints[i]);
    }
    return result;
}

/// Compile a route of waypoints from the given solution.
///
/// @param solution The solution that you want to use for the route.
/// @returns A feasible route
std::vector<Waypoint> RefuelingProblem::route(const SolutionBase* solution) throw(InvalidSolution)
{
    const RefuelingSolution* ref_solution = dynamic_cast<const RefuelingSolution*>(solution);
    if (0 == ref_solution)
        throw InvalidSolution("Solution not RefuelingSolution.");

    std::vector<Waypoint> output;
    output.reserve(_waypoints.size());
    for (const auto& sol_entry : ref_solution->genome)
    {
        output.push_back(_waypoints[sol_entry.index]);
    }
    return output;
}

/// Report the fuel left at each waypoints given the solution.
///
/// @param solution The solution that you want to use for the route.
/// @returns A fuel left
std::vector<double> RefuelingProblem::fuelLeft(const SolutionBase* solution) throw(InvalidSolution)
{
    const RefuelingSolution* ref_solution = dynamic_cast<const RefuelingSolution*>(solution);
    if (0 == ref_solution)
        throw InvalidSolution("Solution not RefuelingSolution.");

    std::vector<double> output;
    output.reserve(_waypoints.size());
    for (const auto& sol_entry : ref_solution->genome)
    {
        output.push_back(sol_entry.fuelLeft);
    }
    return output;
}
