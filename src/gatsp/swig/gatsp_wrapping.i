%module gatsp_wrapping
%{
#include "gatsp/genetic_algorithm.h"

#include "gatsp/point.h"
#include "gatsp/quaternion.h"
#include "gatsp/waypoint.h"
#include "gatsp/euclidean_3d_problem.h"
#include "gatsp/refueling_problem.h"
%}

%include "std_vector.i"
%include "std_string.i"

%newobject ProblemBase::makeSolution();
%newobject ProblemBase::crossover();
%newobject GeneticAlgorithm::bestSolution();
%include "gatsp/genetic_algorithm.h"

%include "gatsp/point.h"
%include "gatsp/quaternion.h"
%include "gatsp/waypoint.h"
%template(WaypointVector) std::vector<Waypoint>;

%newobject Euclidean3DSolution::clone();
%newobject Euclidean3DProblem::makeSolution();
%newobject Euclidean3DProblem::crossover();
%include "gatsp/euclidean_3d_problem.h"

%newobject RefuelingProblem::makeSolution();
%newobject RefuelingProblem::crossover();
%include "gatsp/refueling_problem.h"
%template(DoubleVector) std::vector<double>;
