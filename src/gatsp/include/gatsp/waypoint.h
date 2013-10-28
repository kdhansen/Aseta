// Header file for the Waypoint class of the GATSP project.
// 
// This file is part of the GATSP project.
// Copyright 2012-2013 Karl Damkjær Hansen - Aalborg University

#ifndef GATSP_WAYPOINT_H
#define GATSP_WAYPOINT_H

#include "gatsp/point.h"
#include "gatsp/quaternion.h"

class Waypoint
{
public:
    Waypoint();
    Waypoint(Point, Quaternion);

    Point point() const;
    void point(Point);

    double x() const;
    void x(double);

    double y() const;
    void y(double);

    double z() const;
    void z(double);

    Quaternion quaternion() const;
    void quaternion(Quaternion);

private:
    Point _point;
    Quaternion _quaternion;
};

#endif // GATSP_WAYPOINT_H