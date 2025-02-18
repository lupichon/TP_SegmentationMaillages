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
}