/**
 * @file Contact.cpp
 *
 * Nom: Vincent Fleury
 * Code permanent : FLEV81080005
 * Email : vincent.fleury.1@ens.etsmtl.ca
 *
 */

#include "Contact.h"
#include "TetMesh.h"

float Contact::k = 10000.0f;
float Contact::b = 40.0f;

Contact::Contact() : p(), n(), phi(0.0f)
{

}

Contact::Contact(Particle* _p, const Eigen::Vector3f& _n, float _phi) :
    p(_p), n(_n), phi(_phi)
{
}

Contact::~Contact()
{

}

void Contact::computePenaltyForce()
{
    p->fcontact += n * (-k*phi - n.dot(b * p->xdot));
}
