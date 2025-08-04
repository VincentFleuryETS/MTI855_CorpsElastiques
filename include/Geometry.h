#pragma once

#include <Eigen/Dense>

// List of geometry type ids.
enum eGeometryType { kSphere, kPlane };

// Generic geometry interface.
//
class Geometry
{
public:
    virtual eGeometryType getType() const  = 0;
};

// Plane geometry.
//
// The infinite plane is defined by its perpendicular direction (normal),
// and we assume that the plane is centered at the body COM.
//
class Plane : public Geometry
{
public:
    Eigen::Vector3f p;          // Point on the plane.
    Eigen::Vector3f n;          // The plane normal.

    Plane(const Eigen::Vector3f& _p, const Eigen::Vector3f& _n) : n(_n), p(_p) {}
    virtual ~Plane() {}

    virtual eGeometryType getType() const override { return kPlane; }
};
