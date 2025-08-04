#pragma once

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <vector>
#include <array>

// Struct for storing mesh node (particle) data.
//
struct Particle
{
    Particle() : x(0,0,0), u(0,0,0), xdot(0,0,0), f(0,0,0), fmouse(0,0,0), fcontact(0,0,0), mass(1.0f), fixed(false), ind(0)
    {  }

    Eigen::Vector3f x;      // Current node position
    Eigen::Vector3f u;      // Undeformed node position
    Eigen::Vector3f xdot;   // Node velocity
    Eigen::Vector3f f;      // Elastic + external forces
    Eigen::Vector3f fmouse; // Mouse spring force
    Eigen::Vector3f fcontact; // Penalty-based contact forces
    Eigen::Vector3f dx;     // Auxiliary variable used for differential computation
    Eigen::Vector3f df;     // Force differentials
    float mass;             // Node mass.
    bool fixed;             // fixed flag
    unsigned int ind;     

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// Struct for storing mesh element data.
struct Tetrahedron
{
    std::array<unsigned int, 4> inds;   // The indices of the four nodes that comprise the tetrahedral element
    Eigen::Matrix3f Uinv;               // Linear basis matrix, where Uinv = inverse([ u1-u0, u2-u0, u3-u0 ])
    Eigen::Matrix3f F;                  // Deformation gradient, computed as F = [ x1-x0, x2-x0, x3-x0 ] * Uinv
    Eigen::Matrix3f Q;                  // Rotation matrix, computed by a polar decomposition of F (used only for co-rotated linear strain)
    float W;                            // The volume of the element

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// Surface triangles
struct Triangle
{
    std::array<unsigned int, 3> inds;   // The indices of the three nodes of the triangle
};


// Class to store a tetrahedral mesh.
//
class TetMesh
{
public:
    // Default constructor.
    TetMesh();

    // Read-only accessor for vertices
    const std::vector<Particle>& getParticles() const { return m_particles; }

    // Write accessor for vertices
    std::vector<Particle>& getParticles() { return m_particles; }

    // Read-only accessor for tetrahedra
    const std::vector<Tetrahedron>& getTetrahedra() const { return m_tets; }

    // Write accessor for tetrahedra
    std::vector<Tetrahedron>& getTetrahedra() { return m_tets; }

    // Read-only accessor for triangles
    const std::vector<Triangle>& getTriangles() const { return m_tris; }

    // Write accessor for triangles
    std::vector<Triangle>& getTriangles() { return m_tris; }

    // Initialize the tetrahedral mesh
    void init(const Eigen::Isometry3f& initTM = Eigen::Isometry3f::Identity());

    // Remove all elements and all nodes
    void clear();

    // Sets all elastic and external forces to zero
    void resetForces();

    // Computes elastic and external forces.
    void computeForces();

    // Computes the force differential for all elements based on values in Particle::dx
    void computeForceDifferentials();

    // Computes and caches the deformation gradient Tetrahedron::F for all elements
    void computeDeformationGradients();

    // mult() is used by the iterative conjugate gradient solver to compute
    // the product of A*v, where v is a current solution estimate and A is the
    // principal matrix:   A = M - dt * dt * K
    //
    void mult(float dt, const Eigen::VectorXf& v, Eigen::VectorXf& Av);

    // rhs() builds the right-hand side vector:
    //
    //   b = dt * f + dt * dt * K * xdot
    //
    void rhs(float dt, Eigen::VectorXf& b);

    float E;                            // Young's modulus
    float nu;                           // Poisson's ratio
    float density;                      // Masse volumique 

private:

    std::vector<Particle> m_particles;  // Mesh nodes
    std::vector<Tetrahedron> m_tets;    // Mesh elements
    std::vector<Triangle> m_tris;       // Surface triangles
};
