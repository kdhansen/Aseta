%module gatsp_wrapping
%{
#include "gatsp/point.h"
#include "gatsp/quaternion.h"
#include "gatsp/waypoint.h"
#include "gatsp/solution_base.h"
#include "gatsp/euclidean_3d_problem.h"
#include "gatsp/problem_base.h"
#include "gatsp/euclidean_3d_problem.h"
#include "gatsp/genetic_algorithm_base.h"
#include "gatsp/traditional_genetic_algorithm.h"
%}

%include "std_vector.i"
%include "gatsp/point.h"
%include "gatsp/quaternion.h"
%include "gatsp/waypoint.h"
%template(WaypointVector) std::vector<Waypoint>;
%include "gatsp/solution_base.h"
%template(SolutionVector) std::vector<size_t>;
%template(SolutionEntryBaseVector) std::vector<SolutionEntryBase>;
%include "gatsp/problem_base.h"
%include "gatsp/euclidean_3d_problem.h"
%include "gatsp/genetic_algorithm_base.h"
%include "gatsp/traditional_genetic_algorithm.h"
