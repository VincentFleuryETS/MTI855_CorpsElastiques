#pragma once

#include <Eigen/Dense>
#include <vector>

namespace polyscope
{
    class VolumeMesh;
    class PointCloud;
    class CurveNetwork;
}

class BackwardEuler;
class CollisionDetect; 
class Geometry;
class TetMesh;
struct Particle;

// Viewer for a soft body simulation.
//
class SoftBodyViewer
{
public:
    SoftBodyViewer();
    virtual ~SoftBodyViewer();

    void start();

private:
    void loadTet();
    void loadCube();
    void loadBeam();
    void loadBunny();

    void initTetData();
    void updateTetData();
    void draw();
    void drawGUI();

    // The one and only tetrahedral model. When loading new meshes, call initTetData() and 
    // initialize with the nodes, elements, and triangles of the new mesh.
    std::unique_ptr<TetMesh> m_tetMesh;                 

    std::unique_ptr<BackwardEuler> m_stepper;           // Backward Euler integrator
    std::unique_ptr<CollisionDetect> m_collision;       // Backward Euler integrator

    polyscope::PointCloud* m_volPoints;
    polyscope::VolumeMesh* m_volMesh;
    polyscope::CurveNetwork* m_mouseCurve;

    // Simulation parameters
    float m_dt;                         // Time step parameter.
    bool m_paused;                      // Pause the simulation.
    bool m_stepOnce;                    // Advance the simulation by one frame and then stop.
    bool m_useExplicitIntegration;      // Enable/disable explicit integration.

    Particle* m_pickParticle;           // The picked particle for mouse spring interaction (null by default)
    float m_pickDepth;

    std::vector<Geometry*>  m_colliders;
};
