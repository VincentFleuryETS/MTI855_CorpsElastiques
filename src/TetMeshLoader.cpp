#include "TetMeshLoader.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include "TetMesh.h"

//--------------------------------------------------------------------------------------------------
// Constructors / Destructors
TetMeshLoader::TetMeshLoader(const std::string& filename)
{
    readFile(filename);
}

TetMeshLoader::~TetMeshLoader()
{}

//--------------------------------------------------------------------------------------------------
// Load file
bool TetMeshLoader::readFile(const std::string& filename)
{
    // Open the input file
    std::ifstream file(filename.c_str(), std::ifstream::in);

    struct stat buffer;
    bool exists = (stat (filename.c_str(), &buffer) == 0);

    if (!file.is_open())
    {
        std::cout << "Error: Failed to open file " << filename << " for reading!" << std::endl;
        return false;
    }

    // Read file
    std::string line;
    while (std::getline(file, line))
    {

        // Remove trailing whitespace, just in case
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());

        if (line[0] == '#')
        {
            // Comments... just ignore the line
            continue;
        }
        else if (line == "Vertices")
        {
            // Vertex start
            std::getline(file, line);
            const unsigned int numVerts = std::stoi(line);

            m_verts.resize(numVerts);
            for(unsigned int i = 0; i < numVerts; ++i)
            {
                std::getline(file, line);
                std::stringstream ss(line);
                ss >> m_verts[i][0] >> m_verts[i][1] >> m_verts[i][2];
            }
        }
        else if (line == "Triangles")
        {
            // Triangles start
            std::getline(file, line);
            const unsigned int numTris = std::stoi(line);

            m_tris.resize(numTris);
            for(unsigned int i = 0; i < numTris; ++i)
            {
                std::getline(file, line);
                std::stringstream ss(line);
                unsigned int a, b, c;
                ss >> a >> b >> c;
                m_tris[i][0] = a-1;
                m_tris[i][1] = b-1;
                m_tris[i][2] = c-1;
            }
        }
        else if (line == "Tetrahedra")
        {
            // Tetrahedra start
            std::getline(file, line);
            const unsigned int numTets = std::stoi(line);

            m_tets.resize(numTets);
            for(unsigned int i = 0; i < numTets; ++i)
            {
                std::getline(file, line);
                std::stringstream ss(line);
                unsigned int a, b, c, d;
                ss >> a >> b >> c >> d;
                m_tets[i][0] = a-1;
                m_tets[i][1] = b-1;
                m_tets[i][2] = c-1;
                m_tets[i][3] = d-1;
            }
        }
    }

    // Close file
    file.close();

    return true;
}

void TetMeshLoader::initMesh(TetMesh& tetMesh, const Eigen::Isometry3f& initTM)
{
    tetMesh.clear();
    auto& verts = tetMesh.getParticles();
    auto& tets = tetMesh.getTetrahedra();
    auto& tris = tetMesh.getTriangles();

    const unsigned int numVerts = m_verts.size();
    verts.resize(numVerts);
    for(unsigned int i = 0; i < numVerts; ++i)
    {
        verts[i].u = Eigen::Vector3f(m_verts[i].data());
        verts[i].ind = i;
    }

    const unsigned int numTets = m_tets.size();
    tets.resize(numTets);
    for(unsigned int i = 0; i < numTets; ++i)
    {
        tets[i].inds = m_tets[i];
    }

    const unsigned int numTris = m_tris.size();
    tris.resize(numTris);
    for(unsigned int i = 0; i < numTris; ++i)
    {
        tris[i].inds = m_tris[i];
    }

    tetMesh.init(initTM);
}
