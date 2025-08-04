cmake_minimum_required(VERSION 3.12)
project(MTI855-CorpsElastiques)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenGL REQUIRED)
		
if (APPLE)
  add_definitions( -DGL_SILENCE_DEPRECATION )
endif()

include(FetchContent)
FetchContent_Declare(
  polyscope
  GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
  GIT_TAG v2.4.0  
  GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(polyscope)
FetchContent_GetProperties(polyscope SOURCE_DIR polyscope_SRC_DIR BINARY_DIR polyscope_BIN_DIR)

FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE
  SOURCE_SUBDIR cmake
)
FetchContent_MakeAvailable(Eigen)
FetchContent_GetProperties(Eigen SOURCE_DIR Eigen_SRC_DIR BINARY_DIR Eigen_BIN_DIR)
		
set(softbody_HEADERS include/BackwardEuler.h
                include/CollisionDetect.h
				include/Contact.h
				include/Geometry.h
				include/SoftBodyViewer.h
				include/TetMesh.h
				include/TetMeshLoader.h )
set(softbody_SOURCE main.cpp
           src/CollisionDetect.cpp
		   src/Contact.cpp
           src/SoftBodyViewer.cpp
           src/TetMesh.cpp
           src/TetMeshLoader.cpp)
		 
add_executable (softBody main.cpp ${softbody_HEADERS} ${softbody_SOURCE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Eigen_SRC_DIR} ${COMMON_INCLUDES})
target_link_libraries(softBody OpenGL::GL polyscope)

source_group(src FILES ${softbody_SOURCE})
source_group(include FILES ${softbody_HEADERS})

if(MSVC) 
  set_property(TARGET softBody PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()
