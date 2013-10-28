// Header file for the Quaternion structure of the GATSP project.
// 
// This file is part of the GATSP project.
// Copyright 2012-2013 Karl Damkjær Hansen - Aalborg University

#ifndef GATSP_QUATERION_H
#define GATSP_QUATERION_H

struct Quaternion
{
	Quaternion()
		: x(0), y(0), z(0), w(0)
	{
	}

	Quaternion(double x, double y, double z, double w)
		: x(x), y(y), z(z), w(w)
	{
	}
	
	double x;
	double y;
	double z;
	double w;
};

#endif // GATSP_QUATERION_H