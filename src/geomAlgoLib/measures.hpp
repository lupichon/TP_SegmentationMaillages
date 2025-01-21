#pragma once
#include <iostream>
#include <io.hpp>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <types.hpp>

namespace geomAlgoLib
{
    Facet_double_map computeArea(const Polyhedron & mesh);
    void computeOFFFile(const Polyhedron & mesh, const Facet_double_map & faces, const std::string & filenameOFF);
}

