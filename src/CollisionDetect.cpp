/**
 * @file CollisionDetect.cpp
 *
 * Nom: Vincent Fleury
 * Code permanent : FLEV81080005
 * Email : vincent.fleury.1@ens.etsmtl.ca
 *
 */

#include "CollisionDetect.h"

#include "Contact.h"
#include "Geometry.h"
#include "TetMesh.h"

CollisionDetect::CollisionDetect() 
{

}

void CollisionDetect::detectCollisions(TetMesh* _tetMesh, const std::vector<Geometry*> _colliders, float dt)
{
    // Clear existing contacts.
    //
    clear();

    // Next, loop over all pairs of bodies and test for contacts.
    //
    auto& particles = _tetMesh->getParticles();
    for (Geometry* g : _colliders)
    {
        for (Particle& p : particles)
        {
            if ( p.fixed )
            {
                // Skip collision tests for 'fixed' particles.
                //
                continue;
            }
            else if (g->getType() == kPlane)
            {
                // Particle-plane collision. 
                //
                collisionDetectParticlePlane(&p, g);
            }

        }
    }

    // Compute all contact forces and add their contributions to Particle::fcontact
    //
    for (Contact* c : m_contacts)
    {
        c->computePenaltyForce();
    }
}

void CollisionDetect::clear()
{
    // Cleanup all contacts.
    //
    for(auto c : m_contacts)
    {
        delete c;
    }
    m_contacts.clear();
}

// Particle-plane collision test.
void CollisionDetect::collisionDetectParticlePlane(Particle* p, Geometry* g)
{
    Plane* plane = dynamic_cast<Plane*>(g);

    // TODO implement point-plane collision detection. 
    //
    // Compute the contact normal @a n and the signed penetration @a phi.
    //
    // If the collision exists, add it to the array @a m_contact.
    //
    //  v-------- Change
    //
    bool collision = false;
    Eigen::Vector3f n = plane->n;
    float phi = n.transpose() * (p->x - plane->p);

    collision = phi <= 0.0f;

    if ( collision )
    {
        m_contacts.push_back(new Contact(p, n, phi));
    }
}
