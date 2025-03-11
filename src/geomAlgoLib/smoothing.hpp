/**
 * @file smoothing.hpp
 * 
 * This file contains the declarations of various mesh processing functions.
 * These functions include smoothing, bounding box calculations, influence computations, and free form deformation.
 * 
 * Functions include:
 * - laplacianSmoothing: Applies Laplacian smoothing to the given mesh over a specified number of iterations.
 * - gaussianSmoothing: Applies Gaussian smoothing to the given mesh with a specified lambda parameter.
 * - taubinSmoothing: Applies Taubin smoothing to the given mesh, with specified lambda and mu parameters.
 * - computeAABBCorner: Computes the 8 corners of the Axis-Aligned Bounding Box (AABB) for the given mesh.
 * - computeInfluences: Computes the influence of each vertex based on its proximity to the corners of the AABB.
 * - freeFormDeformation: Applies free-form deformation to the given mesh using a specified translation and corner index.
 */

#pragma once
#include <iostream>
#include <io.hpp>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <types.hpp>
#include <vector>
#include <limits>

namespace geomAlgoLib {
    
    /**
     * @brief Applies Laplacian smoothing to the given polyhedron mesh.
     * 
     * This function smooths the mesh using Laplacian smoothing
     * 
     * @param mesh The polyhedron mesh to be smoothed.
     * @param iterations The number of smoothing iterations to apply.
     * @return A new polyhedron mesh that has been smoothed.
     */
    Polyhedron laplacianSmoothing(const Polyhedron & mesh, unsigned int iterations);
    
    /**
     * @brief Applies Gaussian smoothing to the given polyhedron mesh.
     * 
     * This function smooths the mesh using Gaussian smoothing.
     * 
     * @param mesh The polyhedron mesh to be smoothed.
     * @param iterations The number of smoothing iterations to apply.
     * @param lambda The smoothing parameter that controls the amount of smoothing.
     * @return A new polyhedron mesh that has been smoothed.
     */
    Polyhedron gaussianSmoothing(const Polyhedron & mesh, unsigned int iterations, double lambda);
    
    /**
     * @brief Applies Taubin smoothing to the given polyhedron mesh.
     * 
     * This function smooths the mesh using Taubin's smoothing algorithm.
     * 
     * @param mesh The polyhedron mesh to be smoothed.
     * @param iterations The number of smoothing iterations to apply.
     * @param lambda The parameter that controls the positive smoothing iteration.
     * @param mu The parameter that controls the negative smoothing iteration.
     * @return A new polyhedron mesh that has been smoothed.
     */
    Polyhedron taubinSmoothing(const Polyhedron & mesh, unsigned int iterations, double lambda, double mu);
    
    /**
     * @brief Computes the 8 corners of the Axis-Aligned Bounding Box (AABB) of the mesh.
     * 
     * This function computes the 8 corners of the AABB for the given polyhedron mesh.
     * 
     * @param mesh The polyhedron mesh for which the AABB corners are computed.
     * @return A vector containing the 8 corners of the AABB as CGAL 3D points.
     */
    std::vector<CGAL::Point_3<Kernel>> computeAABBCorner(const Polyhedron & mesh);
    
    /**
     * @brief Computes the influence of each vertex in the mesh based on the AABB corners.
     * 
     * This function computes the influence of each vertex in the mesh. The influence is calculated 
     * for each corner, and the resulting influence values are stored in a map.
     * 
     * @param mesh The polyhedron mesh to be processed.
     * @return A vector of maps, each containing the influence of each vertex at a particular corner of the AABB.
     */
    std::vector<std::map<Polyhedron::Vertex_const_handle, double>> computeInfluences(const Polyhedron & mesh);
    
    /**
     * @brief Applies free form deformation (FFD) to the mesh.
     * 
     * This function applies free form deformation (FFD) to the polyhedron mesh by 
     * translating the mesh vertices based on a specified translation vector and a corner index.
     * 
     * @param mesh The polyhedron mesh to be deformed.
     * @param translation The translation vector to apply to the mesh.
     * @param cornerIndex The index of the corner that will be used as the anchor for deformation.
     * @return A new polyhedron mesh that has been deformed.
     */
    Polyhedron freeFormDeformation(Polyhedron & mesh, const CGAL::Point_3<Kernel>& translation, size_t cornerIndex);

}
