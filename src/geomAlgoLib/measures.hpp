/**
 * @file measure.hpp
 * 
 * This file contains the declaration of functions for mesh segmentation.
 * 
 * Functions include:
 * - computeArea: Calculates the area of each face.
 * - computeSmallestAngle: Computes the smallest angle between adjacent edges of each face.
 * - computeAngleNormalAxes: Computes the angles between the face normals and the X, Y, and Z axes.
 * - computeAltitudes: Computes the average altitude (Z-coordinate) of each face.
 * - computeRoughness: Computes the roughness of each facet based on the angle between adjacent facets.
 * - computeOFFFile: Exports the mesh to an OFF file with face properties such as area, with colors for visualization.
 * - computeSegmentation: Segments the mesh into distinct regions based on geometric properties or clustering.
 * - computeColoredMesh: Colors the mesh based on computed properties, like area or roughness, for better visualization.
 * - computeColoredMeshMultiClass: Applies multi-class color segmentation to the mesh, categorizing faces into multiple classes.
 */

#pragma once
#include <iostream>
#include <io.hpp>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <types.hpp>

namespace geomAlgoLib
{
    /**
     * @brief Computes the area of each facet in the given polyhedron mesh.
     * 
     * This function computes the area for each facet in the polyhedron mesh.
     * It returns a map where each facet is associated with its computed area.
     * 
     * @param mesh The polyhedron mesh to be processed.
     * @return A map of facets to their respective areas.
     */
    Facet_double_map computeArea(const Polyhedron & mesh);

    /**
     * @brief Exports the polyhedron mesh to an OFF file with facet properties.
     * 
     * This function exports the polyhedron mesh to an OFF file format, where each facet is associated with a
     * specific property, such as area or smallest angle, which is provided in the faces map.
     * 
     * @param mesh The polyhedron mesh to be exported.
     * @param faces A map of facets to their respective properties (e.g., area, angle).
     * @param filenameOFF The path to the output OFF file.
     */
    void computeOFFFile(const Polyhedron & mesh, const Facet_double_map & faces, const std::string & filenameOFF);

    /**
     * @brief Computes the smallest angle of each facet in the given polyhedron mesh.
     * 
     * This function calculates the smallest internal angle for each facet in the polyhedron mesh.
     * The result is returned as a map of facets to their respective smallest angles.
     * 
     * @param mesh The polyhedron mesh to be processed.
     * @return A map of facets to their respective smallest angles.
     */
    Facet_double_map computeSmallestAngle(const Polyhedron & mesh);

    /**
     * @brief Computes the angles of each facet's normal relative to the coordinate axes (X, Y, Z).
     * 
     * This function computes the angles between the normal vector of each facet and the X, Y, and Z axes.
     * The results are stored in three separate maps for the X, Y, and Z axes.
     * 
     * @param mesh The polyhedron mesh to be processed.
     * @param anglesNormalX A map to store the angles relative to the X axis.
     * @param anglesNormalY A map to store the angles relative to the Y axis.
     * @param anglesNormalZ A map to store the angles relative to the Z axis.
     */
    void computeAngleNormalAxes(const Polyhedron & mesh, Facet_double_map & anglesNormalX, Facet_double_map & anglesNormalY, Facet_double_map & anglesNormalZ);

    /**
     * @brief Computes the altitude (average Z-coordinate) of each facet in the given polyhedron mesh.
     * 
     * This function calculates the average Z-coordinate (altitude) for each facet in the polyhedron mesh.
     * The result is returned as a map of facets to their respective altitudes.
     * 
     * @param mesh The polyhedron mesh to be processed.
     * @return A map of facets to their respective altitudes.
     */
    Facet_double_map computeAltitudes(const Polyhedron & mesh);

    /**
     * @brief Computes the roughness (angle between adjacent facet normals) of each facet in the mesh.
     * 
     * This function computes the roughness for each facet in the polyhedron mesh by calculating the angle between
     * the normals of adjacent facets.
     * 
     * @param mesh The polyhedron mesh to be processed.
     * @return A map of facets to their respective roughness values (angles between normals).
     */
    Facet_double_map computeRoughness(const Polyhedron & mesh);

    /**
     * @brief Segments the polyhedron mesh based on a given threshold and assigns each facet to a class.
     * 
     * This function segments the polyhedron mesh into two classes based on a threshold value. Facets with values
     * greater than the threshold are assigned to `class1`, while those below the threshold are assigned to `class2`.
     * 
     * @param measures A map of facets to their respective measured values (e.g., area, angle, etc.).
     * @param threshold The threshold value used for segmentation.
     * @param class1 The label for the first class.
     * @param class2 The label for the second class.
     * @return A map of facets to their respective class labels.
     */
    Facet_string_map computeSegmentation(const Facet_double_map & measures, float threshold, std::string class1, std::string class2);

    /**
     * @brief Computes and exports a colored mesh to an OFF file based on a binary segmentation map.
     * 
     * This function generates an OFF file representation of the polyhedron mesh, where each facet is colored based
     * on its class from a binary segmentation map. Facets belonging to `class1` are assigned `color1`, while those
     * belonging to `class2` are assigned `color2`.
     * 
     * @param mesh The polyhedron mesh to be processed.
     * @param segmentationMap A map of facets to their respective class labels.
     * @param filenameOFF The path to the output OFF file where the colored mesh will be saved.
     * @param color1 The RGB color for `class1` facets.
     * @param color2 The RGB color for `class2` facets.
     * @param class1 The label for the first class.
     * @param class2 The label for the second class.
     */
    void computeColoredMesh(const Polyhedron & mesh, const Facet_string_map & segmentationMap, const std::string & filenameOFF, const Color & color1, const Color & color2, std::string class1, std::string class2);

    /**
     * @brief Computes and exports a multi-class colored mesh to an OFF file based on a multi-class segmentation map.
     * 
     * This function generates an OFF file representation of the polyhedron mesh, where each facet is colored based
     * on its associated class from a multi-class segmentation map. Facets that belong to multiple classes will have
     * an averaged color.
     * 
     * @param mesh The polyhedron mesh to be processed.
     * @param segmentationMultimap A multimap of facets to their respective class labels.
     * @param filenameOFF The path to the output OFF file where the colored mesh will be saved.
     * @param colors A vector of RGB color values for each class.
     * @param classes A vector of class labels.
     */
    void computeColoredMeshMultiClass(const Polyhedron & mesh, const Facet_string_multimap & segmentationMultimap, const std::string & filenameOFF, const std::vector<Color> & colors, const std::vector<std::string> & classes);
}



