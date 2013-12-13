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

#ifndef REFULEING_PROBLEM_H
#define REFULEING_PROBLEM_H

#include "gatsp/euclidean_3d_problem.h"

struct RefuelingEntry
{
    RefuelingEntry(size_t i)
    {
        index = i;
        cost = -1;
        fuelLeft = 0;
    }

    bool operator==(const RefuelingEntry& other) const
    { return other.index == index; }

    bool operator<(const RefuelingEntry& other) const
    { return index < other.index; }

    size_t index;
    mutable double fuelLeft;
    mutable double cost;
};

typedef PathSolution<RefuelingEntry> RefuelingSolution;


class RefuelingProblem : public ProblemBase
{
public:
    RefuelingProblem(unsigned int seed = 0);
    virtual ~RefuelingProblem();

    virtual void addWaypoint(const Waypoint& waypoint);
    virtual void addDepot(const Waypoint& waypoint);
    void fuelCapacity(double metersPrTank);

    virtual void mutate(SolutionBase* mutatee) const;
    virtual SolutionBase* makeSolution() const;
    virtual SolutionBase* crossover(const SolutionBase* dad, const SolutionBase* mom) const;
    virtual double cost(const SolutionBase* solution) const;

    std::vector<Waypoint> route(const SolutionBase* solution) throw(InvalidSolution);
    std::vector<Waypoint> waypoints();
    std::vector<Waypoint> depots();
    std::vector<double> fuelLeft(const SolutionBase* solution) throw(InvalidSolution);

private:
    std::vector<Waypoint> _waypoints;
    std::vector<bool> _isDepot;
    std::vector<size_t> _depotIndices;
    double _metersPrTank;

    double euclideanDistance(const Waypoint& from, const Waypoint& to) const;
    void injectDepot(RefuelingSolution* mutatee) const;
    void extractDepot(RefuelingSolution* mutatee) const;

    mutable std::mt19937 _random_generator;
};

#endif