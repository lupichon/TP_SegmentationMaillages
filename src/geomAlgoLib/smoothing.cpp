/**
 * @file smoothing.cpp
 * 
 * This file contains the implementation of various mesh processing functions.
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

#include <smoothing.hpp>

namespace geomAlgoLib {
    
    Polyhedron laplacianSmoothing(const Polyhedron & mesh, unsigned int iterations) 
    {
        Polyhedron smoothed_mesh = mesh; 

        for (unsigned int iter = 0; iter < iterations; ++iter) 
        {
            std::vector<CGAL::Point_3<Kernel>> new_positions(smoothed_mesh.size_of_vertices(), CGAL::Point_3<Kernel>(0, 0, 0));
            
            unsigned int index = 0;
            for (auto v = smoothed_mesh.vertices_begin(); v != smoothed_mesh.vertices_end(); ++v, ++index) 
            {
                CGAL::Point_3<Kernel> centroid(0, 0, 0);
                Polyhedron::Halfedge_around_vertex_const_circulator h = v->vertex_begin(); 
                int num_neighbors = 0;
                int ctr = 0;
                do {
                    Vertex_iterator neighbor = h->opposite()->vertex();
                    centroid = centroid + (neighbor->point() - CGAL::Point_3<Kernel>(0, 0, 0));
                    ++num_neighbors;
                } while (++h != v->vertex_begin()); 
                
                if (num_neighbors > 0) 
                {
                    centroid = CGAL::Point_3<Kernel>(centroid.x() / num_neighbors, centroid.y() / num_neighbors, centroid.z() / num_neighbors);
                    new_positions[index] = centroid;
                } 
                else 
                {
                    new_positions[index] = v->point(); 
                }
            }

            index = 0;
            for (auto v = smoothed_mesh.vertices_begin(); v != smoothed_mesh.vertices_end(); ++v, ++index) 
            {
                v->point() = new_positions[index];
            }
        }
        return smoothed_mesh;
    }

    Polyhedron gaussianSmoothing(const Polyhedron & mesh, unsigned int iterations, double lambda)
    {
        Polyhedron smoothed_mesh = mesh; 

        for (unsigned int iter = 0; iter < iterations; ++iter) 
        {
            std::vector<CGAL::Point_3<Kernel>> new_positions(smoothed_mesh.size_of_vertices(), CGAL::Point_3<Kernel>(0, 0, 0));
            
            unsigned int index = 0;
            for (auto v = smoothed_mesh.vertices_begin(); v != smoothed_mesh.vertices_end(); ++v, ++index) 
            {
                CGAL::Point_3<Kernel> centroid(0, 0, 0);
                Polyhedron::Halfedge_around_vertex_const_circulator h = v->vertex_begin(); 
                int num_neighbors = 0;
                int ctr = 0;
                do {
                    Vertex_iterator neighbor = h->opposite()->vertex();
                    centroid = centroid + (neighbor->point() - CGAL::Point_3<Kernel>(0, 0, 0));
                    ++num_neighbors;
                } while (++h != v->vertex_begin()); 
                
                if (num_neighbors > 0) 
                {
                    centroid = CGAL::Point_3<Kernel>(centroid.x() / num_neighbors, centroid.y() / num_neighbors, centroid.z() / num_neighbors);
                    CGAL::Vector_3<Kernel> delta = centroid - v->point();  
                    new_positions[index] = v->point() + lambda * delta;

                } 
                else 
                {
                    new_positions[index] = v->point(); 
                }
            }

            index = 0;
            for (auto v = smoothed_mesh.vertices_begin(); v != smoothed_mesh.vertices_end(); ++v, ++index) 
            {
                v->point() = new_positions[index];
            }
        }
        return smoothed_mesh;
    }

    Polyhedron taubinSmoothing(const Polyhedron & mesh, unsigned int iterations, double lambda, double mu)
    {
        Polyhedron smoothed_mesh = mesh;

        for (unsigned int iter = 0; iter < iterations; ++iter)
        {
            std::vector<CGAL::Point_3<Kernel>> new_positions(smoothed_mesh.size_of_vertices(), CGAL::Point_3<Kernel>(0, 0, 0));
            
            unsigned int index = 0;
            for (auto v = smoothed_mesh.vertices_begin(); v != smoothed_mesh.vertices_end(); ++v, ++index)
            {
                CGAL::Point_3<Kernel> centroid(0, 0, 0);
                Polyhedron::Halfedge_around_vertex_const_circulator h = v->vertex_begin();
                int num_neighbors = 0;
                
                do {
                    Vertex_iterator neighbor = h->opposite()->vertex();
                    centroid = centroid + (neighbor->point() - CGAL::Point_3<Kernel>(0, 0, 0));
                    ++num_neighbors;
                } while (++h != v->vertex_begin());
                
                if (num_neighbors > 0)
                {
                    centroid = CGAL::Point_3<Kernel>(centroid.x() / num_neighbors, centroid.y() / num_neighbors, centroid.z() / num_neighbors);
                    CGAL::Vector_3<Kernel> delta = centroid - v->point();
                    new_positions[index] = v->point() + lambda * delta;
                }
                else
                {
                    new_positions[index] = v->point();
                }
            }

            index = 0;
            for (auto v = smoothed_mesh.vertices_begin(); v != smoothed_mesh.vertices_end(); ++v, ++index)
            {
                v->point() = new_positions[index];
            }

            std::vector<CGAL::Point_3<Kernel>> final_positions(smoothed_mesh.size_of_vertices(), CGAL::Point_3<Kernel>(0, 0, 0));

            index = 0;
            for (auto v = smoothed_mesh.vertices_begin(); v != smoothed_mesh.vertices_end(); ++v, ++index)
            {
                CGAL::Point_3<Kernel> centroid(0, 0, 0);
                Polyhedron::Halfedge_around_vertex_const_circulator h = v->vertex_begin();
                int num_neighbors = 0;
                
                do {
                    Vertex_iterator neighbor = h->opposite()->vertex();
                    centroid = centroid + (neighbor->point() - CGAL::Point_3<Kernel>(0, 0, 0));
                    ++num_neighbors;
                } while (++h != v->vertex_begin());

                if (num_neighbors > 0)
                {
                    centroid = CGAL::Point_3<Kernel>(centroid.x() / num_neighbors, centroid.y() / num_neighbors, centroid.z() / num_neighbors);
                    CGAL::Vector_3<Kernel> delta = centroid - v->point();
                    final_positions[index] = v->point() - mu * delta;  
                }
                else
                {
                    final_positions[index] = v->point();
                }
            }
            index = 0;
            for (auto v = smoothed_mesh.vertices_begin(); v != smoothed_mesh.vertices_end(); ++v, ++index)
            {
                v->point() = final_positions[index];
            }
        }
        return smoothed_mesh;
    }

    std::vector<CGAL::Point_3<Kernel>> computeAABBCorner(const Polyhedron & mesh)
    {
        std::vector<CGAL::Point_3<Kernel>> corners(8);

        CGAL::Point_3<Kernel> min_point = CGAL::Point_3<Kernel>(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
        CGAL::Point_3<Kernel> max_point = CGAL::Point_3<Kernel>(-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
        
        for (auto v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v)
        {
            const CGAL::Point_3<Kernel>& p = v->point();
            min_point = CGAL::Point_3<Kernel>(std::min(min_point.x(), p.x()), std::min(min_point.y(), p.y()), std::min(min_point.z(), p.z()));
            max_point = CGAL::Point_3<Kernel>(std::max(max_point.x(), p.x()), std::max(max_point.y(), p.y()), std::max(max_point.z(), p.z()));
        }

        corners[0] = CGAL::Point_3<Kernel>(min_point.x(), min_point.y(), min_point.z());
        corners[1] = CGAL::Point_3<Kernel>(min_point.x(), min_point.y(), max_point.z());
        corners[2] = CGAL::Point_3<Kernel>(min_point.x(), max_point.y(), min_point.z());
        corners[3] = CGAL::Point_3<Kernel>(min_point.x(), max_point.y(), max_point.z());
        corners[4] = CGAL::Point_3<Kernel>(max_point.x(), min_point.y(), min_point.z());
        corners[5] = CGAL::Point_3<Kernel>(max_point.x(), min_point.y(), max_point.z());
        corners[6] = CGAL::Point_3<Kernel>(max_point.x(), max_point.y(), min_point.z());
        corners[7] = CGAL::Point_3<Kernel>(max_point.x(), max_point.y(), max_point.z());

        return corners;
    }

    std::vector<std::map<Polyhedron::Vertex_const_handle, double>> computeInfluences(const Polyhedron& mesh)
    {
        std::vector<CGAL::Point_3<Kernel>> AABBCorners = computeAABBCorner(mesh);
        std::vector<std::map<Polyhedron::Vertex_const_handle, double>> influences(8);

        const CGAL::Point_3<Kernel>& min_point = AABBCorners[0];
        const CGAL::Point_3<Kernel>& max_point = AABBCorners[7];
        double V_total = (max_point.x() - min_point.x()) * (max_point.y() - min_point.y()) * (max_point.z() - min_point.z());

        for (auto v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v)
        {
            const CGAL::Point_3<Kernel>& P = v->point();

            double x_rel = (P.x() - min_point.x()) / (max_point.x() - min_point.x());
            double y_rel = (P.y() - min_point.y()) / (max_point.y() - min_point.y());
            double z_rel = (P.z() - min_point.z()) / (max_point.z() - min_point.z());

            for (size_t i = 0; i < AABBCorners.size(); ++i)
            {
                const CGAL::Point_3<Kernel>& corner = AABBCorners[i];

                double weight_x = 1.0 - std::abs(corner.x() - P.x()) / (max_point.x() - min_point.x());
                double weight_y = 1.0 - std::abs(corner.y() - P.y()) / (max_point.y() - min_point.y());
                double weight_z = 1.0 - std::abs(corner.z() - P.z()) / (max_point.z() - min_point.z());

                double influence = weight_x * weight_y * weight_z;

                influences[i][v] = influence;
            }
        }

        return influences;
    }

    Polyhedron freeFormDeformation(Polyhedron& mesh, const CGAL::Point_3<Kernel>& translation, size_t cornerIndex)
    {
        Polyhedron deformedMesh = mesh;

        std::vector<std::map<Polyhedron::Vertex_const_handle, double>> influences = computeInfluences(deformedMesh);

        std::vector<CGAL::Point_3<Kernel>> AABBCorners = computeAABBCorner(mesh);
        CGAL::Point_3<Kernel> targetCorner = AABBCorners[cornerIndex];
        CGAL::Vector_3<Kernel> translationVector = translation - targetCorner;
        
        AABBCorners[cornerIndex] = targetCorner + translationVector;

        for (auto v = deformedMesh.vertices_begin(); v != deformedMesh.vertices_end(); ++v)
        {
            CGAL::Point_3<Kernel> newPosition = v->point();
            CGAL::Vector_3<Kernel> deformation(0, 0, 0);

            for (size_t i = 0; i < 8; ++i)
            {
                double influence = influences[i][v];
                CGAL::Vector_3<Kernel> cornerToPoint = AABBCorners[i] - v->point();
                deformation = deformation + influence * cornerToPoint;
            }

            newPosition = v->point() + deformation;
            v->point() = newPosition;
        }

        return deformedMesh;
    }
}