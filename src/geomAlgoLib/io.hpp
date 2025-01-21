#pragma once

#include "types.hpp"

#include <iostream>
#include <fstream>

namespace geomAlgoLib
{

/// Read an OFF file and store the mesh in "mesh".
/// Returns true if it was loaded successfully.
bool readOFF(const std::string& filePath, Polyhedron& mesh);

/// Write a mesh at location "filePath"
void writeOFF(const Polyhedron& mesh, const std::string& filePath);

}