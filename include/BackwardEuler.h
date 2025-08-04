#pragma once

#include "TetMesh.h"
#include <memory>

//
// Backward Euler integrator for tetrahedral meshes.
//
// This class uses a matrix-free Conjugate Gradient solve
// to compute the velocities of mesh nodes.

//  The functions TetMesh::computeForces() and TetMesh::computeForceDifferentials
//  are indirectly called by the functions mult() and rhs()
//
class BackwardEuler
{
public:
    BackwardEuler(TetMesh* _tetMesh) : m_tetMesh(_tetMesh), m_deltav(), maxIters(20)
    {
        init();
    }

    // Initialize the stepper
    void init()
    {
        auto& particles = m_tetMesh->getParticles();
        const unsigned int dim = 3 * particles.size();
        m_deltav.setZero(dim);
    }

    void step(float dt)
    {
        auto& particles = m_tetMesh->getParticles();
        const unsigned int dim = 3 * particles.size();

        m_tetMesh->computeDeformationGradients();

        // Build the rhs vector.
        //
        Eigen::VectorXf b;
        m_tetMesh->rhs(dt, b);

        // Solve for delta v, the velocity changes.
        //
        m_deltav.setZero(dim);
        cgsolve(dt, b, m_deltav);

        // Integrate the node positions
        //
        for(Particle& p : particles)
        {
            if( !p.fixed )
            {
                p.xdot += m_deltav.segment(3*p.ind, 3);
                p.x += dt * p.xdot;
            }
            else
            {
                p.xdot.setZero();
            }
        }

    }

    int maxIters;           // Maximum conjugate gradient iterations

private:

    // Conjugate gradient solve.
    //
    void cgsolve(float dt, const Eigen::VectorXf& b, Eigen::VectorXf& x)
    {
        Eigen::VectorXf Adv;
        m_tetMesh->mult(dt, x, Adv);

        Eigen::VectorXf r = b - Adv;
        Eigen::VectorXf p = r;
        float rTr = r.dot(r);
        for(int iter = 0; iter < maxIters; ++iter)
        {
            if(rTr < 1e-12f ) return;

            Eigen::VectorXf Ap;
            m_tetMesh->mult(dt, p, Ap);
            float alpha = rTr / p.dot(Ap);
            x += alpha * p;
            r = r - alpha * Ap;
            float rTrnew = r.dot(r);
            p = r + (rTrnew / rTr) * p;
            rTr = rTrnew;
        }

    }

    TetMesh* m_tetMesh;             // Pointer to the tetrahedral mesh.
    Eigen::VectorXf m_deltav;       // Solution vector of linear systems (M - dt*dt*K) deltav

};
