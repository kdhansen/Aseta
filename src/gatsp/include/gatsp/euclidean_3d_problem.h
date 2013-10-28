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

#ifndef EUCLIDEAN_3D_PROBLEM_H
#define EUCLIDEAN_3D_PROBLEM_H

#include "gatsp/problem_base.h"

class Euclidean3DProblem : public ProblemBase
{
public:
    virtual SolutionBase makeSolution();
    virtual std::vector<Waypoint> route(const SolutionBase&) throw(InvalidSolution);
    virtual Waypoint firstWaypoint(const SolutionBase&) throw(InvalidSolution);
    virtual void popWaypointFront(SolutionBase&) throw(InvalidSolution);
    virtual void addWaypoint(const Waypoint&);
    virtual void addWaypoint(const Waypoint&, SolutionBase&);

    virtual double solutionCost(const SolutionBase&) const throw(InvalidSolution);
    virtual bool isSolutionValid(const SolutionBase&) const;
    virtual bool repairSolution(const SolutionBase&) const;

private:
    std::vector<Waypoint> _waypoints;

    double euclideanDistance(const Waypoint&, const Waypoint&) const;
};

#endif // EUCLIDEAN_3D_PROBLEM_H