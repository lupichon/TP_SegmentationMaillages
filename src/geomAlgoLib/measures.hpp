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
    Facet_double_map computeSmallestAngle(const Polyhedron & mesh);
    void computeAngleNormalAxes(const Polyhedron & mesh, Facet_double_map & anglesNormalX, Facet_double_map & anglesNormalY, Facet_double_map & anglesNormalZ);
    Facet_double_map computeAltitudes(const Polyhedron & mesh);
    Facet_double_map computeRoughness(const Polyhedron & mesh);
    Facet_string_map computeSegmentation(const Facet_double_map & measures, float treshold, std::string class1, std::string class2);
    void computeColoredMesh(const Polyhedron & mesh, const Facet_string_map & segmentationMap, const std::string & filenameOFF, const Color & color1, const Color & color2, std::string class1, std::string class2);
    void computeColoredMeshMultiClass(const Polyhedron & mesh, const Facet_string_multimap & segmentationMultimap, const std::string & filenameOFF, const std::vector<Color> & colors, const std::vector<std::string> & classes);
}




