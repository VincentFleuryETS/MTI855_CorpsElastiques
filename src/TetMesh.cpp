/**
 * @file TetMesh.cpp
 *
 * @brief Tetrahedral mesh class.
 *
 * Nom: Vincent Fleury
 * Code permanent : FLEV81080005
 * Email : vincent.fleury.1@ens.etsmtl.ca
 *
 */

#include "TetMesh.h"

namespace
{
    // Compute the polar decomposition of F, returning the orthonormal rotation matrix in Q.
    static inline void polar(const Eigen::Matrix3f &F, Eigen::Matrix3f &Q)
    {
        Eigen::JacobiSVD<Eigen::Matrix3f> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const Eigen::Matrix3f& U = svd.matrixU();
        const Eigen::Matrix3f& V = svd.matrixV();
        Q = U * V.transpose();
    }

    // Compute the Lame parameters lambda and my from the Young's modulus (E) and Poisson's ratio (nu)
    static inline void computeLameParameters(float E, float nu, float& lambda, float& mu)
    {
        lambda = (E * nu) / ((1.0f + nu) * (1.0f - 2.0f*nu));
        mu = E / (2.0f * (1.0f + nu));
    }
}

TetMesh::TetMesh() : m_particles(), m_tets(), m_tris(), E(1e5f), nu(0.35f), density(1.0f)
{
    clear();
}

void TetMesh::clear()
{
    m_particles.clear();
    m_tets.clear();
    m_tris.clear();
}

void TetMesh::resetForces()
{
    // Reset all forces acting on the particles
    for(Particle& p : m_particles)
    {
        p.f.setZero();
        p.fmouse.setZero();
        p.fcontact.setZero();
    }
}

void TetMesh::init(const Eigen::Isometry3f& initTM)
{
    // Initialize the nodes (particles)
    //
    int count = 0;
    for(Particle& p : m_particles)
    {
        p.ind = count++;
        p.x = initTM * p.u;
        p.f.setZero();
        p.fmouse.setZero();
        p.xdot.setZero();
        p.dx.setZero();
        p.mass = 0.0f;
        p.df.setZero();
    }

    // TODO Initialize the tetrahedral elements
    for(Tetrahedron& t : m_tets)
    {
        // TODO Compute the linear shape matrix (undeformed).
        //      This is the matrix inverse(U) from the course notes.
        //      Store the result in t.Uinv
        //

        // We obtain the particles used by this Tetrahedron.
        Particle p0 = m_particles[t.inds[0]];
        Particle p1 = m_particles[t.inds[1]];
        Particle p2 = m_particles[t.inds[2]];
        Particle p3 = m_particles[t.inds[3]];

        Eigen::Matrix3f U;
        U.col(0) = (p1.u - p0.u);
        U.col(1) = (p2.u - p0.u);
        U.col(2) = (p3.u - p0.u);

        // We calculate and store the U inverse.
        t.Uinv = U.inverse();

        // Compute the volume of the tetrahedron as
        //         t.W = (1/6) * det(U)
        //
        //      https://en.wikipedia.org/wiki/Tetrahedron#Volume
        //
        t.W = (1/6.0f) * U.determinant();

        // Compute the individual mass of the tetrahedra at each node as
        //           individualMass = density * t.W / 4
        //
        const float individualMass = density * t.W / 4.0f;

        // Distribute the mass uniformly between the four nodes of the tetrahedron.
        // We add, and not set the mass, as the nodes can be part of multiple tetrahedron.
        //
        for(int i = 0; i < t.inds.size(); i ++){
            m_particles[t.inds[i]].mass += individualMass;
        }

        // Finally, we initialize the auxialiary variables for storing
        // the co-rotation (Q) and deformation gradient (F)
        //
        t.Q.setIdentity();
        t.F.setIdentity();
    }
}

void TetMesh::computeForces()
{
    // Apply gravitational forces.
    static const Eigen::Vector3f g(0.0f, -9.81f, 0.0f);
    for(Particle& p : m_particles)
    {
        p.f = p.mass * g;
    }

    float lambda, mu;
    computeLameParameters(E, nu, lambda, mu);

    // TODO For each tetrahedral element,
    // compute and apply elastic forces.
    //
    //  Here, you may assume that the deformation gradient, t.F,
    //  has already been computed.
    //
    for(const Tetrahedron& t : m_tets)
    {
        Eigen::Matrix3f strain, stress, H;

        // TODO Compute stress using the linear Green strain tensor
        //      and the 2nd Piola stress tensor.
        //
        //      Co-rotation should be used for this case, since otherwise the strain tensor
        //      is not invariant to the global rotation.
        //

        strain = ( 0.5f * ( t.F + t.F.transpose() ) ) - Eigen::Matrix3f::Identity();
        stress = lambda * strain.trace() * Eigen::Matrix3f::Identity() + 2 * mu * strain;

        // TODO: Compute the nodal force matrix H = -W * Q * stress * (U)^-T
        // and update Particle::f for each node in the element.
        //
        // Note: Q will be the identity matrix if co-rotation is not used,
        //       otherwise it will be the orthonormal rotation matrix computed
        //       using polar decomposition in the previous step.
        //
        H = -t.W * t.Q * stress * t.Uinv.transpose();

        // We add the elastic forces to each particle.
        // Equation (4.7) of https://viterbi-web.usc.edu/~jbarbic/femdefo/sifakis-courseNotes-TheoryAndDiscretization.pdf 

        m_particles[t.inds[1]].f += H.col(0);
        m_particles[t.inds[2]].f += H.col(1);
        m_particles[t.inds[3]].f += H.col(2);
        m_particles[t.inds[0]].f -= H.col(0) + H.col(1) + H.col(2);
    }
}


void TetMesh::computeForceDifferentials()
{
    // Compute Lame parameters.
    float lambda, mu;
    computeLameParameters(E, nu, lambda, mu);

    for(Particle& p : m_particles)
    {
        p.df = Eigen::Vector3f::Zero();
    }

    // TODO For each tetrahedral element,
    //      compute and sum the force differentials Particle::df
    //      based on the values of the nodal differentials Particle::dx
    for(Tetrahedron& t : m_tets)
    {
        // TODO Compute the differential of the deformation gradient as
        //           dF = dX * Uinv
        //      where dX is a 3x3 matrix that is assembled as
        //           dX = [ dx1 - dx0, dx2 - dx0,  dx3 - dx0 ]
        //
        Eigen::Matrix3f dF, dX, dH;

        Particle p0 = m_particles[t.inds[0]];
        Particle p1 = m_particles[t.inds[1]];
        Particle p2 = m_particles[t.inds[2]];
        Particle p3 = m_particles[t.inds[3]];

        dX.col(0) = (p1.dx - p0.dx);
        dX.col(1) = (p2.dx - p0.dx);
        dX.col(2) = (p3.dx - p0.dx);

        dF = dX * t.Uinv;

        // TODO Compute stress differentials.
        //
        // Compute the stress differential using the formula for
        // co-rotated linear Green strain tensor and
        // the 2nd Piola stress tensor.
        //
        Eigen::Matrix3f dstress, dstrain;
        dstrain = 0.5f * (dF.transpose() + dF);
        dstress = ( 2 * mu * dstrain + lambda * dstrain.trace() * Eigen::Matrix3f::Identity() );


        // TODO Compute the nodal force differential matrix dH = -W * Q * dstress * (U)^-T
        // and update Particle::df for each node in the element.
        //
        // Note: Q will be the identity matrix if co-rotation is not used.
        //

        dH = -t.W * t.Q * dstress * t.Uinv.transpose();

        m_particles[t.inds[1]].df += dH.col(0);
        m_particles[t.inds[2]].df += dH.col(1);
        m_particles[t.inds[3]].df += dH.col(2);
        m_particles[t.inds[0]].df -= dH.col(0) + dH.col(1) + dH.col(2);
    }
}

void TetMesh::computeDeformationGradients()
{
    for(Tetrahedron& t : m_tets)
    {
        // TODO Compute F = [ x1-x0, x2-x0, x3-x0 ] * U^-1
        //
        
        // We obtain the particles used by this Tetrahedron.
        Particle p0 = m_particles[t.inds[0]];
        Particle p1 = m_particles[t.inds[1]];
        Particle p2 = m_particles[t.inds[2]];
        Particle p3 = m_particles[t.inds[3]];

        Eigen::Matrix3f X;
        X.col(0) = (p1.x - p0.x);
        X.col(1) = (p2.x - p0.x);
        X.col(2) = (p3.x - p0.x);

        t.F = X * t.Uinv;

        // Update the co-rotation matrix.
        //
        polar(t.F, t.Q);
    }
}

void TetMesh::mult(float dt, const Eigen::VectorXf& deltav, Eigen::VectorXf& Av)
{
    // TODO Set dx for each particle to
    //      the corresponding values from the
    //      vector v.
    //
    for(Particle& p : m_particles)
    {
        p.dx = deltav.segment( 3 * p.ind, 3 );
    }

    // Next, we recompute the force differentials.
    // Effectively, this computes df = K*v
    //
    computeForceDifferentials();

    // TODO Compute the product of the principal matrix times the vector v.
    //   i.e.
    //       Av = (M - dt*dt*K) * v
    //  This is compute using a matrix free approach, per particle.
    //
    //  Hint : Compute this is in two steps.
    //    1. Compute Av = Mv.  Recall that M is a diagonal mass matrix
    //       and so computing Mv may be done using the mass stored at each particle.
    //
    //    2. Compute dt*dt*K*v and add it to the result from the first step.
    //       Recall that df = K*v has already been computed and stored for
    //       each particle by the call to computeForceDifferentials()
    //
    //
    Av.setZero( 3 * m_particles.size() );
    for(Particle& p : m_particles)
    {
        Av.segment( 3 * p.ind, 3 ) = p.mass * deltav.segment( 3 * p.ind, 3 ); // ** REMPLACER PAR p.dx?
        Av.segment( 3 * p.ind, 3 ) -= dt * dt * p.df;
    }
}

void TetMesh::rhs(float dt, Eigen::VectorXf& b)
{
    computeForces();

    // TODO Set Particle::dx to Particle::xdot for each particle,
    //      and set Particle::df to zero.
    for(Particle& p : m_particles)
    {
        p.dx = p.xdot;
    }

    // Next, we recompute the force differentials.
    // Effectively, this computes df = K*xdot
    //
    computeForceDifferentials();

    // TODO Compute the rhs vector and compute as:
    //    b = dt * (f + fmouse + fcontact)  + dt*dt*K*xdot
    //
    //   Recall that K*xdot has already been computed
    //   and stored in df for each particle.
    //
    b.setZero( 3 * m_particles.size() );
    for(Particle& p : m_particles)
    {
        b.segment( 3 * p.ind, 3 ) = dt * (p.f + p.fmouse + p.fcontact) + dt * dt * p.df;
    }
}
