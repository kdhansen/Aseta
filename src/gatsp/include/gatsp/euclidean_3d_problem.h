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

#ifndef GATSP_EUCLIDEAN_3D_PROBLEM_H
#define GATSP_EUCLIDEAN_3D_PROBLEM_H

#include <random>
#include <vector>
#include "gatsp/genetic_algorithm.h"
#include "gatsp/waypoint.h"

struct Euclidean3DEntry
{
    Euclidean3DEntry(size_t i)
    {
        index = i;
        cost = -1;
    }
    bool operator==(const Euclidean3DEntry& other) const
    {
        return other.index == index;
    }
    bool operator<(const Euclidean3DEntry& other) const
    {
        return index < other.index;
    }
    size_t index;
    mutable double cost;
};

struct Euclidean3DSolution : public SolutionBase
{
    Euclidean3DSolution() {};
    Euclidean3DSolution(
        const std::vector<Euclidean3DEntry>::const_iterator begin, 
        const std::vector<Euclidean3DEntry>::const_iterator end
    );
    virtual ~Euclidean3DSolution() {};
    virtual Euclidean3DSolution* clone();

    std::vector<Euclidean3DEntry> genome;
};

class Euclidean3DProblem : public ProblemBase
{
public:
    Euclidean3DProblem(unsigned int seed = 0);
    virtual ~Euclidean3DProblem();

    virtual SolutionBase* makeSolution() const;
    virtual void mutate(SolutionBase* mutatee) const;
    virtual SolutionBase* crossover(const SolutionBase* dad, const SolutionBase* mom) const;
    virtual double cost(const SolutionBase* solution) const;

    virtual Euclidean3DSolution* makeNearestNeighbor() const;
    virtual Euclidean3DSolution* makeTrivialSolution() const;
    virtual std::vector<Waypoint> route(const SolutionBase*) throw(InvalidSolution);
    virtual Waypoint firstWaypoint(const SolutionBase*) throw(InvalidSolution);
    virtual void popWaypointFront(Euclidean3DSolution&) throw(InvalidSolution);
    virtual void addWaypoint(const Waypoint& waypoint);
    // virtual void addWaypoint(const Waypoint&, Euclidean3DSolution&);

    virtual bool isSolutionValid(const Euclidean3DSolution&) const;
    // virtual bool repairSolution(const Euclidean3DSolution&) const;

private:
    double euclideanDistance(const Waypoint&, const Waypoint&) const;

    void displacementMutation(Euclidean3DSolution* s) const;
    void exchangeMutation(Euclidean3DSolution* s) const;
    void inversionMutation(Euclidean3DSolution* s) const;

    Euclidean3DSolution* orderCrossover(const Euclidean3DSolution* dad, const Euclidean3DSolution* mom) const;

    std::vector<Waypoint> _waypoints;
    mutable std::mt19937 _random_generator;
};

#endif // GATSP_EUCLIDEAN_3D_PROBLEM_H