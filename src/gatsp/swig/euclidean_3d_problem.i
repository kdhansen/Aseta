%module euclidean_3d_problem
%{
#include "gatsp/euclidean_3d_problem.h"
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
// %include "std_string.i"
