#pragma once

#include <Eigen/Dense>

struct Particle;

// Contact for non-interpenetration of soft bodies. 
// Penalty forces are used and applied to particles (nodes) in the tetrahedral mesh.
//
class Contact
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    static float k;                    // Contact stiffness (global)
    static float b;                    // Contact damping (global)

public:

    // Constructor with all parameters.
    Contact(Particle* _p, const Eigen::Vector3f& _n, float _phi);

    virtual ~Contact();

    // Computes a penalty contact force and adds it to p->fcontact.
    // The stiffness (k) and damping (b) are combined with the penetration depth
    // phi and particle velocity in a spring equation.
    //
    void computePenaltyForce();

    Particle* p;                // The tet node.
    Eigen::Vector3f n;          // The contact normal.
    float phi;                  // Signed penetration, should be negative for collisions.

protected:

    // Default constructor.
    explicit Contact();
};
