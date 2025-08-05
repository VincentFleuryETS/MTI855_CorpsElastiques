#include "SoftBodyViewer.h"

#include "polyscope/polyscope.h"
#include "polyscope/volume_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"
#include "polyscope/pick.h"
#include "polyscope/view.h"
#include "imgui.h"

#include <iostream>

using namespace std;

#include "BackwardEuler.h"
#include "CollisionDetect.h"
#include "Contact.h"
#include "Geometry.h"
#include "TetMeshLoader.h"
#include "TetMesh.h"

namespace polyscope
{
    namespace view
    {
        glm::vec3 screenCoordsToWorld(glm::vec2 screenCoords, float depth)
        {
            glm::mat4 view = getCameraViewMatrix();
            glm::mat4 proj = getCameraPerspectiveMatrix();
            glm::vec4 viewport = { 0., 0., view::windowWidth, view::windowHeight };

            glm::vec3 screenPos3{ screenCoords.x, view::windowHeight - screenCoords.y, depth };
            glm::vec3 worldPos = glm::unProject(screenPos3, view, proj, viewport);

            return worldPos;
        }

        glm::vec2 worldToScreenCoords(glm::vec3 worldPos, float* depth) 
        {
            glm::mat4 view = getCameraViewMatrix();
            glm::mat4 proj = getCameraPerspectiveMatrix();
            glm::vec4 viewport = { 0., 0., view::windowWidth, view::windowHeight };

            glm::vec3 screenPos3 = glm::project(worldPos, view, proj, viewport);
            if (depth) *depth = screenPos3.z;

            return glm::vec2(screenPos3.x, view::windowHeight - screenPos3.y);
        }
    }
}

namespace
{
    static const std::array<float, 3> pinColor = { 1.0f, 0.0f, 0.0f };
    static const std::array<float, 3> pointColor = { 1.0f, 1.0f, 0.0f };

    // Apply a mouse spring force to the selected particle.
    bool computeMouseSpringForce(Particle* _pickParticle, float pickDepth, glm::vec3& p0, glm::vec3& p1)
    {
        // Particle being dragged by the mouse spring.
        // Note: mouse spring force must be applied after
        // computeForces() is called in the integration loop.
        // 
        if (ImGui::IsMouseDown(0) && ImGui::GetIO().KeyCtrl && _pickParticle)
        {
            ImVec2 mouseP = ImGui::GetMousePos();
            glm::vec2 mouseVec = { mouseP.x, mouseP.y };
            p1 = polyscope::view::screenCoordsToWorld(mouseVec, pickDepth);
            p0 = { _pickParticle->x.x(), _pickParticle->x.y(), _pickParticle->x.z() };

            glm::vec3 f = { p1.x - p0.x, p1.y - p0.y, p1.z - p0.z };

            // Apply mouse spring force.
            {
                Eigen::Vector3f u(f.x, f.y, f.z);
                const float ulen = u.norm();

                if (ulen > 0.01f)
                {
                    u.normalize();
                    // Compute damping for stable spring forces
                    const float k = 5000.0f * _pickParticle->mass;
                    const float b = 50.0f * _pickParticle->mass;
                    _pickParticle->fmouse = (k * ulen - b * (_pickParticle->xdot.dot(u))) * u;
                }
            }
            return true;
        }

        return false;
    }

    static void updateContactPoints(polyscope::PointCloud* points, const std::vector<Contact*>& contacts)
    {
        const unsigned int numPoints = points->nPoints();
        const unsigned int numContacts = contacts.size();
        Eigen::MatrixXf contactN(numPoints, 3);

        contactN.setZero();

        for (unsigned int i = 0; i < numContacts; ++i)
        {
            contactN.row( contacts[i]->p->ind ) = contacts[i]->n.transpose();
        }

        points->addVectorQuantity("normal", contactN)->setVectorColor({ 1.0f, 0.0f, 0.0f })->setVectorLengthScale(0.1f)->setEnabled(true);
    }
}

SoftBodyViewer::SoftBodyViewer() :
    m_tetMesh(nullptr),
    m_dt(0.01f),
    m_paused(true),
    m_stepOnce(false),
    m_useExplicitIntegration(false),
    m_pickParticle(nullptr),
    m_pickDepth(0.0f),
    m_volMesh(nullptr), 
    m_volPoints(nullptr), 
    m_mouseCurve(nullptr)
{

}

SoftBodyViewer::~SoftBodyViewer()
{
}

void SoftBodyViewer::start()
{
    // Setup Polyscope
    polyscope::options::programName = "MTI855 Devoir 02 - Corps mous";
    polyscope::options::verbosity = 0;
    polyscope::options::usePrefsFile = false;
    polyscope::options::alwaysRedraw = true;
    polyscope::options::ssaaFactor = 2;
    polyscope::options::openImGuiWindowForUserCallback = true;
    polyscope::options::groundPlaneHeightFactor = 0.0f; // adjust the plane height
    polyscope::options::automaticallyComputeSceneExtents = false;
    polyscope::state::lengthScale = 4.;
    polyscope::state::boundingBox =
        std::tuple<glm::vec3, glm::vec3>{ {-1., -1., -1.}, {1., 1., 1.} };
    polyscope::options::buildGui = false;
    polyscope::options::maxFPS = 120;
    polyscope::view::windowWidth = 1920;
    polyscope::view::windowHeight = 1080;

    // initialize
    polyscope::init();

    // Specify the update callback
    polyscope::state::userCallback = std::bind(&SoftBodyViewer::draw, this);

    m_tetMesh = std::make_unique<TetMesh>();
    loadTet();

    m_stepper = std::make_unique<BackwardEuler>(m_tetMesh.get());
    m_collision = std::make_unique<CollisionDetect>();

    // Show the window
    polyscope::show();
}

void SoftBodyViewer::drawGUI()
{
    ImGui::Text("Simulation: ");
    ImGui::Checkbox("Pause", &m_paused);
    if (ImGui::Button("Step once"))
    {
        m_stepOnce = true;
    }

    ImGui::Text("Integration: ");
    ImGui::PushItemWidth(200);
    ImGui::SliderFloat("Time step", &m_dt, 0.0f, 0.1f, "%0.4f");
    ImGui::SliderInt("CG iterations", &m_stepper->maxIters, 1, 100, "%d");
    ImGui::Checkbox("Explicit integration", &m_useExplicitIntegration);
    ImGui::PopItemWidth();
    
    ImGui::Text("Material parameters: ");
    ImGui::PushItemWidth(200);
    ImGui::SliderFloat("Young's modulus", &m_tetMesh->E, 0.0f, 1e9f, "%6.0f", ImGuiSliderFlags_Logarithmic);
    ImGui::SliderFloat("Poisson ratio", &m_tetMesh->nu, 0.01f, 0.49f, "%0.2f", 0);
    ImGui::SliderFloat("Density", &m_tetMesh->density, 0.01f, 10.0f, "%2.1f", 0);
    ImGui::PopItemWidth();

    ImGui::Text("Contact parameters: ");
    ImGui::PushItemWidth(200);
    ImGui::SliderFloat("Stiffness", &Contact::k, 0.0f, 1e6f, "%6.0f", ImGuiSliderFlags_Logarithmic);
    ImGui::SliderFloat("Damping", &Contact::b, 0.0f, 1e5f, "%6.0f", ImGuiSliderFlags_Logarithmic);
    ImGui::PopItemWidth();

    ImGui::Text("Scenarios: ");
    ImGui::PushItemWidth(200);
    if (ImGui::Button("Single tet.")) 
    {
        loadTet();
    }
    else if (ImGui::Button("Cube")) 
    {
        loadCube();
    }
    else if (ImGui::Button("Beam")) 
    {
        loadBeam();
    }
    else if (ImGui::Button("Bunny")) 
    {
        loadBunny();
    }
    ImGui::PopItemWidth();

}

void SoftBodyViewer::initTetData()
{
    const auto& particles = m_tetMesh->getParticles();
    const auto& tets = m_tetMesh->getTetrahedra();
    const unsigned int numParticles = particles.size();
    const unsigned int numTets = tets.size();

    Eigen::MatrixXf meshV(numParticles, 3);
    Eigen::MatrixXi meshT(numTets, 4);

    for (int i = 0; i < numParticles; ++i)
    {
        meshV.row(i) << particles[i].x(0), particles[i].x(1), particles[i].x(2);
    }
    for (int i = 0; i < numTets; ++i)
    {
        meshT.row(i) << tets[i].inds[0], tets[i].inds[1], tets[i].inds[2], tets[i].inds[3];
    }

    // Register the mesh with Polyscope
    m_volMesh = polyscope::registerTetMesh("tets", meshV, meshT);
    m_volMesh->setMaterial("wax");
    m_volMesh->setTransparency(0.6f);
    m_volMesh->setEdgeWidth(1.0);
    // Register the particles point cloud with Polyscope
    m_volPoints = polyscope::registerPointCloud("particles", meshV);
    m_volPoints->setPointRadius(0.01);
    m_volPoints->setPointRenderMode(polyscope::PointRenderMode::Sphere);
    m_volPoints->addColorQuantity("colors", std::vector< std::array<float, 3> >(numParticles, pointColor))->setEnabled(true);

    // Curve for mouse spring
    std::vector<glm::vec3> points(2);
    m_mouseCurve = polyscope::registerCurveNetworkLine("mouseSpring", points);

}

void SoftBodyViewer::updateTetData()
{
    const auto& particles = m_tetMesh->getParticles();
    const auto& tets = m_tetMesh->getTetrahedra();
    const unsigned int numParticles = particles.size();
    const unsigned int numTets = tets.size();

    Eigen::MatrixXf meshV(numParticles, 3);

    for (int i = 0; i < numParticles; ++i)
    {
        meshV.row(i) << particles[i].x(0), particles[i].x(1), particles[i].x(2);
    }
    // Update tet mesh vertexes
    m_volMesh->updateVertexPositions(meshV);
    m_volPoints->updatePointPositions(meshV);

    // Polyscope requires updating *all* point quantities
     //
    std::vector< std::array<float, 3> > pointColors(numParticles, pointColor);
    for (int i = 0; i < numParticles; ++i)
    {
        if (particles[i].fixed)
            pointColors[i] = pinColor;
    }
    m_volPoints->addColorQuantity("colors", pointColors);
}

void SoftBodyViewer::draw()
{
    drawGUI();

    // Perform particle selection
    //
    if (ImGui::IsMouseDown(0) && ImGui::GetIO().KeyCtrl && m_pickParticle == nullptr)
    {
        const ImVec2 mouseP = ImGui::GetMousePos();
        const auto selection = polyscope::pick::evaluatePickQuery(mouseP.x, mouseP.y);

        if (m_volPoints == selection.first)
        {
            const unsigned int pickInd = selection.second;
            auto& particles = m_tetMesh->getParticles();
            m_pickParticle = &particles[pickInd];

            polyscope::view::worldToScreenCoords({ m_pickParticle->x.x(), m_pickParticle->x.y(), m_pickParticle->x.z() }, &m_pickDepth);
        }
    }
    // Perform particle pinning
    // 
    else if (ImGui::IsMouseClicked(1) && ImGui::GetIO().KeyCtrl && m_pickParticle == nullptr)
    {
        const ImVec2 mouseP = ImGui::GetMousePos();
        const auto selection = polyscope::pick::evaluatePickQuery(mouseP.x, mouseP.y);

        if (m_volPoints == selection.first)
        {
            const unsigned int pickInd = selection.second;
            auto& particles = m_tetMesh->getParticles();
            particles[pickInd].fixed = !(particles[pickInd].fixed);
        }
    }
    else if (ImGui::IsMouseReleased(0) && m_pickParticle)
    {
        m_pickParticle = nullptr;
    }

    // Reset forces for all particles (nodes)
    //
    m_tetMesh->resetForces();

    // Perform mouse spring interaction.
    //
    std::vector<glm::vec3> points(2);
    const bool picking = computeMouseSpringForce(m_pickParticle, m_pickDepth, points[0], points[1]);
    m_mouseCurve->updateNodePositions(points);
    m_mouseCurve->setEnabled(picking);

    // Simulation stepping
    //
    if (!m_paused || m_stepOnce)
    {
        m_collision->detectCollisions(m_tetMesh.get(), m_colliders, m_dt);

        if ( m_useExplicitIntegration )  // use explicit Euler integration
        {
            m_tetMesh->computeDeformationGradients();
            m_tetMesh->computeForces();

            for (Particle& p : m_tetMesh->getParticles())
            {
                if (!p.fixed)
                {
                    p.xdot += m_dt * (p.f + p.fmouse + p.fcontact) / p.mass;
                    p.x += m_dt * p.xdot;
                }
                else
                {
                    p.xdot.setZero();
                }
            }
        }
        else        // use Backward Euler integration
        {
            m_stepper->step(m_dt);
        }

        m_stepOnce = false;
    }

    // Draw contacts
    updateContactPoints(m_volPoints, m_collision->getContacts());

    // Update cloth mesh and point positions for rendering, 
    //  as well as cloth params.
    //
    updateTetData();
}

void SoftBodyViewer::loadTet()
{
    TetMeshLoader loader("resources/tet.mesh");
    loader.initMesh(*m_tetMesh);

    for (Particle& p : m_tetMesh->getParticles())
    {
        if (p.u.y() == 0.0f)
        {
            p.fixed = true;
        }
    }

    m_colliders.clear();
    m_colliders.push_back(new Plane({ 0,-1,0 }, { 0, 1, 0 }));

    initTetData();
}

void SoftBodyViewer::loadCube()
{

    Eigen::Isometry3f initTm;
    initTm.setIdentity();
    initTm.translation() = Eigen::Vector3f(0, 1.0f, 0);
    TetMeshLoader loader("resources/cube.mesh");
    loader.initMesh(*m_tetMesh, initTm);

    m_colliders.clear();
    m_colliders.push_back(new Plane({ 0,-1,0 }, { 0, 1, 0 }));

    initTetData();
}

void SoftBodyViewer::loadBeam()
{
    Eigen::Isometry3f initTm;
    initTm.setIdentity();
    initTm.translation() = Eigen::Vector3f(0, 1.0f, 0);

    TetMeshLoader loader("resources/beam.mesh");
    loader.initMesh(*m_tetMesh, initTm);

    for (Particle& p : m_tetMesh->getParticles())
    {
        if (p.u.x() <= -5.0f)
        {
            p.fixed = true;
        }
    }
    m_colliders.clear();
    m_colliders.push_back(new Plane({ 0,-1,0 }, { 0, 1, 0 }));

    initTetData();
}

void SoftBodyViewer::loadBunny()
{
    TetMeshLoader loader("resources/bunny.mesh");
    loader.initMesh(*m_tetMesh);

    m_colliders.clear();
    m_colliders.push_back(new Plane({ 0,-1,0 }, { 0, 1, 0 }));

    initTetData();
}
