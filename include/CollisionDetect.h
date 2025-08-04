#pragma once

#include <Eigen/Dense>
#include <vector>

class Contact;
class TetMesh;
class Geometry;
struct Particle;

//
// The main collision detection class,
// it will test the geometries of each pair of rigid bodies
// and populate an array with contacts.
//
class CollisionDetect
{
public:

    // Constructor.
    //
    CollisionDetect();
    
    // Tests for collisions between all pairs of bodies in the rigid body system
    // and generates contacts for any intersecting geometries.
    // The array of generated contacts can be retrieved by calling getContacts().
    //
    void detectCollisions(TetMesh* _tetMesh, const std::vector<Geometry*> _colliders, float dt);

    void clear();

    // Accessors for the contacts
    const std::vector<Contact*>& getContacts() const { return m_contacts; }
    std::vector<Contact*>& getContacts() { return m_contacts; }

private:

    // Particle-plane collision test. Assumes that @a g is of type Plane.
    //
    void collisionDetectParticlePlane(Particle* p, Geometry* g);


    std::vector<Contact*> m_contacts;

};
