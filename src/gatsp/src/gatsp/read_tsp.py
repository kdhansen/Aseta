import re

def read_tsp(filename):
    infile = open(filename)
    # find dimension of problem and make sure the that a NODE_COORD_SECTION exists
    problem_dimension = None
    match = re.search(r'DIMENSION\s*:\s*(\d*)', infile.read())
    if match:
        problem_dimension = int(match.group(1))
    infile.seek(0)
    match = re.search(r'NODE_COORD_SECTION.*\n', infile.read())
    if match:
        has_coordinates = True
        section_start = match.end()
    else:
        has_coordinates = False

    # read the coordinates
    if has_coordinates:
        coordinates = []
        infile.seek(section_start)
        for i in range(problem_dimension):
            line = infile.readline()
            match = re.match(r'\d*\s*(?P<x>[^\s]*)\s*(?P<y>[^\s]*)', line)
            if match:
                coordinates.append( (float(match.group('x')), float(match.group('y'))) )

        print coordinates
