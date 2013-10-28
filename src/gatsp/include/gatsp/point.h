// Header file for the Point structure of the GATSP project.
// 
// This file is part of the GATSP project.
// Copyright 2012-2013 Karl Damkjær Hansen - Aalborg University

#ifndef GATSP_POINT_H
#define GATSP_POINT_H

struct Point
{
	Point()
		: x(0), y(0), z(0)
	{
	}

	Point(double x, double y, double z)
		: x(x), y(y), z(z)
	{	
	}

	double x;
	double y;
	double z;
};

#endif // GATSP_POINT_H