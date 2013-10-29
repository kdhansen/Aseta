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

SolutionBase Euclidean3DProblem::makeSolution() const
{
    SolutionBase solution;
    solution.reserve(_waypoints.size());
    for(size_t i = 0; i < _waypoints.size(); ++i)
    {
        solution.push_back(i);
    }
    return solution;
}

/// Compile a route of waypoints from the given solution.
///
/// @param solution The solution that you want to use for the route.
/// @returns A feasible route
std::vector<Waypoint> Euclidean3DProblem::route(const SolutionBase& solution) throw(InvalidSolution)
{
    if (isSolutionValid(solution))
    {
        std::vector<Waypoint> output;
        output.reserve(_waypoints.size());
        for (auto sol_entry : solution)
        {
            output.push_back(_waypoints[sol_entry.index]);
        }
        return output;
    }
    else
    {
        throw InvalidSolution();
    }
}

/// Report the first waypoint in a sequence given a specific solution.
///
/// @param solution The solution that determines the first waypoint in the sequence.
/// @returns The first waypoint in the sequence.
Waypoint Euclidean3DProblem::firstWaypoint(const SolutionBase& solution) throw(InvalidSolution)
{
    if (isSolutionValid(solution))
    {
        if (solution.size() > 0)
            return _waypoints[solution[0].index];
    }
    else
    {
        throw InvalidSolution();
    }
}

/// Remove the first waypoint in the sequence given a solution.
///
/// @param [in,out] solution The solution that determines the first waypoint, is also updated.
void Euclidean3DProblem::popWaypointFront(SolutionBase& solution) throw(InvalidSolution)
{   
    if (! isSolutionValid(solution))
    {
        throw InvalidSolution();
    }

    SolutionEntryBase remove_entry = solution[0];

    // Remove from the waypoint list
    _waypoints.erase(_waypoints.begin()+remove_entry.index);
    
    // Remove from the solution (erase-remove idiom)
    solution.erase(std::remove(solution.begin(), solution.end(), remove_entry), solution.end());

    // now the indexes of the entries > solution[0].index must be reduced by 1.
    for(auto& entry: solution)
    {
        if (entry.index > remove_entry.index)
        {
            entry.index -= 1;
        }
    }
    return;
}

/// Add a waypoint to the problem.
///
/// @param waypoint Waypoint to add.
void Euclidean3DProblem::addWaypoint(const Waypoint& waypoint)
{
    _waypoints.push_back(waypoint);
    return;
}

/// Add a waypoint to the problem and update related solutions
///
/// @param waypoint Waypoint to add.
/// @param [out] old_solution A solution that must be updated with the new waypoint.
void Euclidean3DProblem::addWaypoint(const Waypoint& waypoint, SolutionBase& old_solution)
{
    addWaypoint(waypoint);
    old_solution.push_back(_waypoints.size()-1);
}

/// Get the cost of a solution
///
/// @param solution The solution to check
/// @returns the Euclidean length of the solution.
double Euclidean3DProblem::solutionCost(const SolutionBase& solution) const throw(InvalidSolution)
{
    if (!isSolutionValid(solution))
    {
        throw InvalidSolution();
    }

    double dist = 0.0;
    if (2 <= solution.size()) // Cannot compute distance between less than two wps.
    {
        // Traverse the solution
        for (auto s = solution.begin(); s!= solution.end()-1; ++s)
        {
            dist += euclideanDistance(_waypoints[s->index], _waypoints[(s+1)->index]);
        }

        // Close the loop
        auto from = _waypoints[(solution.end()-1)->index];
        auto to = _waypoints[0];
        dist += euclideanDistance(from, to);
    }

    return dist;
}

/// Check to see if solution is valid.
///
/// This is a simple check for validity, not wheter the solution is 
/// satisfactory.
///
/// @param solution The solution to check.
/// @returns True if the solution is valid.
bool Euclidean3DProblem::isSolutionValid(const SolutionBase& solution) const
{
    bool valid = true;
    for (auto sol_entry : solution)
    {
        if (sol_entry.index > _waypoints.size() || sol_entry.index < 0)
        {
            valid = false;
            break;
        }
    }
    return valid;
}

/// Attempt to repair solution.
///
/// @param solution The solution to fix.
/// @returns True if the solution could be repaired.
bool Euclidean3DProblem::repairSolution(const SolutionBase& solution) const
{
    return true;
}

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
