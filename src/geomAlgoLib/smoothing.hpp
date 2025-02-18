#pragma once
#include <iostream>
#include <io.hpp>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <types.hpp>

namespace geomAlgoLib{

    Polyhedron laplacianSmoothing(const Polyhedron & mesh, unsigned int iterations);
    Polyhedron gaussianSmoothing(const Polyhedron & mesh, unsigned int iterations, double lambda);
}