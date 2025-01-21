#include "measures.hpp"


namespace geomAlgoLib
{

    Facet_double_map computeArea(const Polyhedron & mesh)
    {
        Facet_double_map areas;

        for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i) 
        {
            Halfedge_facet_circulator j = i->facet_begin();
            std::vector<CGAL::Point_3<Kernel>> face_points;
            do 
            {
                face_points.push_back(j->vertex()->point());
            }   while (++j != i->facet_begin());
            
            if(face_points.size() == 3)
            {
                areas.insert(std::make_pair(i, sqrt(CGAL::squared_area(face_points[0], face_points[1], face_points[2]))));
            }
            else
            {
                CGAL::Point_3<Kernel> p0 = face_points[0];
                double area = 0;
            
                for (long unsigned int k = 1; k < face_points.size() - 1; k++) {
                    CGAL::Point_3 p1 = face_points[k];
                    CGAL::Point_3 p2 = face_points[k + 1];
                    area += sqrt(CGAL::squared_area(p0, p1, p2));
                }
                areas.insert(std::make_pair(i, area));
            }
        }

        return areas;
    }

    void computeOFFFile(const Polyhedron & mesh, const Facet_double_map & faces, const std::string & filenameOFF)
    {
        
    }

}