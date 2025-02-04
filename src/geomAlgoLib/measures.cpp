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
        std::ofstream in_myfile;
        in_myfile.open(filenameOFF);

        CGAL::set_ascii_mode(in_myfile);

        in_myfile << "COFF" << std::endl << mesh.size_of_vertices() << ' ' << mesh.size_of_facets() << " 0" << std::endl; 

        std::copy(mesh.points_begin(), mesh.points_end(),
                std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));
        
        double min_area = std::numeric_limits<double>::max(); 
        double max_area = std::numeric_limits<double>::lowest(); 

        for (const auto &face : faces)
        {
            double area = face.second; 
            
            if (area < min_area)
            {
                min_area = area; 
            }
            
            if (area > max_area)
            {
                max_area = area;
            }
        }

        for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
        {
            Halfedge_facet_circulator j = i->facet_begin();

            CGAL_assertion(CGAL::circulator_size(j) >= 3);

            in_myfile << CGAL::circulator_size(j) << ' ';
            do
            {
                in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());

            } while (++j != i->facet_begin());

            double color = (faces.at(i) - min_area)/(max_area-min_area);

            in_myfile << std::fixed << std::setprecision(1) << ' ' << color << ' ' << color << ' ' << color;

            in_myfile << std::endl;
        }

        in_myfile.close();

        std::cout << "Colored mesh successfully exported at path: " << filenameOFF << " !" << std::endl;
    }


    Facet_double_map computeSmallestAngle(const Polyhedron &mesh)
    {
        Facet_double_map angles;

        for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i) 
        {
            Halfedge_facet_circulator h = i->facet_begin();
            double smallest_angle = std::numeric_limits<double>::infinity();  

            do {
                Kernel::Vector_3 edge1 = h->vertex()->point() - h->opposite()->vertex()->point();  
                Kernel::Vector_3 edge2 = h->next()->vertex()->point() - h->opposite()->vertex()->point(); 

                double angle = std::acos(edge1 * edge2 / (std::sqrt(edge1 * edge1) * std::sqrt(edge2 * edge2))) * 180/M_PI; 
                smallest_angle = std::min(smallest_angle, angle);

                ++h;  
            } while (h != i->facet_begin());  

            angles.insert(std::make_pair(i, smallest_angle));
        }

        return angles;
    }

    void computeAngleNormalAxes(const Polyhedron & mesh, Facet_double_map & anglesNormalX, Facet_double_map & anglesNormalY, Facet_double_map & anglesNormalZ)
    {
        Kernel::Vector_3 X_axis(1,0,0);
        Kernel::Vector_3 Y_axis(0,1,0);
        Kernel::Vector_3 Z_axis(0,0,1);

        for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i) 
        {
            CGAL::Point_3 p1 = i->facet_begin()->vertex()->point();
            CGAL::Point_3 p2 = i->facet_begin()->opposite()->vertex()->point();
            CGAL::Point_3 p3 = i->facet_begin()->next()->vertex()->point();

            CGAL::Vector_3 normal = CGAL::normal(p1,p2,p3);

            
            auto computeAngle = [](const Kernel::Vector_3 & v1, const Kernel::Vector_3 & v2) -> double {
                return std::acos(v1*v2/(std::sqrt(v1*v1)*std::sqrt(v2*v2))) * 180/M_PI ;
            };
            
            anglesNormalX.insert({i, computeAngle(normal, X_axis)});
            anglesNormalY.insert({i, computeAngle(normal, Y_axis)});
            anglesNormalZ.insert({i, computeAngle(normal, Z_axis)});
        }
    }

    Facet_string_map computeSegmentation(const Facet_double_map & measures, float treshold, std::string class1, std::string class2)
    {
        Facet_string_map segmentation; 

        for(const auto & measure : measures)
        {
            if(abs(measure.second) > treshold)
            {
                segmentation.insert({measure.first, class1});
            }
            else
            {
                segmentation.insert({measure.first, class2});
            }
        }
        return segmentation;
    }

    void computeColoredMesh(const Polyhedron & mesh, const Facet_string_map & segmentationMap, const std::string & filenameOFF, CGAL::Color color1, CGAL::Color color2, std::string class1, std::string class2)
    {
        std::ofstream in_myfile;
        in_myfile.open(filenameOFF);

        CGAL::set_ascii_mode(in_myfile);

        in_myfile << "COFF" << std::endl << mesh.size_of_vertices() << ' ' << mesh.size_of_facets() << " 0" << std::endl; 

        std::copy(mesh.points_begin(), mesh.points_end(),
                std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

        for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
        {
            Halfedge_facet_circulator j = i->facet_begin();

            CGAL_assertion(CGAL::circulator_size(j) >= 3);

            in_myfile << CGAL::circulator_size(j) << ' ';
            do
            {
                in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());

            } while (++j != i->facet_begin());

            if(segmentationMap.at(i) == class1)
            {
                in_myfile << std::fixed << std::setprecision(1) << ' ' << color1;
            }
            else
            {
                in_myfile << std::fixed << std::setprecision(1) << ' ' << color2;
            }

            in_myfile << std::endl;
        }

        in_myfile.close();

        std::cout << "Segmentation mesh successfully exported at path: " << filenameOFF << " !" << std::endl;
    }
}

