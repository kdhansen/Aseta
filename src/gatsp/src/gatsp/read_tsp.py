import re
import gatsp

def read_tsp(filename, tsp_problem):
    """Parse a .tsp file into a tsp_problem"""
    infile = open(filename)
    # Find dimension of problem and make sure the that a NODE_COORD_SECTION exists
    problem_dimension = None
    match = re.search(r'DIMENSION\s*:\s*(?P<dimension>\d*)', infile.read())
    if match:
        problem_dimension = int(match.group('dimension'))
    infile.seek(0)
    match = re.search(r'NODE_COORD_SECTION.*\n', infile.read())
    if match:
        has_coordinates = True
        section_start = match.end()
    else:
        has_coordinates = False

    # Read the coordinates
    if has_coordinates:
        infile.seek(section_start)
        for i in range(problem_dimension):
            line = infile.readline()
            match = re.match(r'\d*\s*(?P<x>[^\s]*)\s*(?P<y>[^\s]*)', line)
            if match:
                p = gatsp.Point(float(match.group('x')), float(match.group('y')), 0)
                q = gatsp.Quaternion()
                tsp_problem.addWaypoint(gatsp.Waypoint(p, q))
