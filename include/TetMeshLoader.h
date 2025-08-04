#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <array>

class TetMesh;

// Class responsible for loading a tetrahedral mesh 
// from a TetGen .mesh file
//
class TetMeshLoader
{
public:
    TetMeshLoader(const std::string& filename);
    ~TetMeshLoader();

    typedef std::array<float, 3> Vertex;
    typedef std::array<unsigned int, 4> Tetrahedron;
    typedef std::array<unsigned int, 3> Triangle;

    const std::vector<Vertex>& getVerts() const { return m_verts; }
    const std::vector<Tetrahedron>& getTets() const { return m_tets; }
    const std::vector<Triangle>& getTris() const { return m_tris; }

    void initMesh(TetMesh& tetMesh, const Eigen::Isometry3f& initTM = Eigen::Isometry3f::Identity());

private:

    bool readFile(const std::string& filename);

    std::vector<Vertex>             m_verts;
    std::vector<Tetrahedron>        m_tets;
    std::vector<Triangle>           m_tris;

};
