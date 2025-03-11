/**
 * @file measures.cpp
 * 
 * This file contains the implementation of functions for mesh segmentation.
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

#include "measures.hpp"

namespace geomAlgoLib
{
    /**
     * @brief Calculates the area of each face in a mesh.
     * 
     * This function iterates over each face in the mesh and computes its area. If the face is a triangle, the area is 
     * calculated using the `squared_area` function from CGAL. If the face has more than three points, it is treated as a 
     * non-triangular polygon, and its area is computed by dividing the face into triangles using the first point as a 
     * reference and iterating over the other points in the face.
     * 
     * @param mesh The 3D mesh for which the areas of the faces are to be calculated.
     * 
     * @return Facet_double_map A map associating each face in the mesh to its calculated area.
     */
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

    /**
     * @brief Exports a mesh to an OFF file with colors corresponding to face properties.
     * 
     * This function exports a 3D mesh to an OFF file format. It iterates over the faces of the mesh and associates a color 
     * with each face based on a property (e.g., area), normalizing the property value to a color scale. The face property 
     * values (such as area) are provided via a map, and the faces are colored according to their respective values, with 
     * the minimum and maximum face property values used for scaling the color.
     * 
     * The resulting OFF file includes the vertices of the mesh followed by the facets, each facet having its associated color.
     * 
     * @param mesh The 3D mesh to be exported.
     * @param faces A map associating each face of the mesh with a calculated property (e.g., area).
     * @param filenameOFF The name of the output OFF file where the mesh and its properties will be saved.
     * 
     */
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

    /**
     * @brief Computes the smallest angle between two adjacent edges for each face in the mesh.
     * 
     * This function calculates the smallest angle formed between any two adjacent edges for each face in the polyhedron. 
     * It iterates over all faces, calculates the angle between consecutive edges using the dot product, and keeps track 
     * of the smallest angle for each face.
     * 
     * The function uses the formula for the angle between two vectors, based on the dot product and vector magnitudes, 
     * to compute the angle between adjacent edges. The smallest angle is selected and stored for each face.
     * 
     * @param mesh The 3D polyhedron mesh for which the smallest angles between adjacent edges will be computed.
     * @return A map associating each face in the mesh with its smallest angle between two adjacent edges.
     * 
     * @note This function assumes that the faces are non-degenerate and are composed of at least two adjacent edges.
     */
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

    /**
     * @brief Computes the angle between the normal of each face and the X, Y, and Z axes.
     * 
     * This function calculates the angle between the normal vector of each face in the polyhedron mesh 
     * and the X, Y, and Z axes. It iterates through all the faces, computes the normal vector of each face, 
     * and then calculates the angle between the normal and each of the three coordinate axes (X, Y, and Z).
     * The result is returned in degrees.
     * 
     * @param mesh The 3D polyhedron mesh for which the angles between face normals and the coordinate axes are computed.
     * @param anglesNormalX A map where the key is a facet and the value is the angle between the normal and the X axis.
     * @param anglesNormalY A map where the key is a facet and the value is the angle between the normal and the Y axis.
     * @param anglesNormalZ A map where the key is a facet and the value is the angle between the normal and the Z axis.
     * 
     * @note This function computes the angles for each face using the normal vector defined by three vertices of the face.
     */
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


    /**
     * @brief Computes the average altitude (Z-coordinate) of each face in the polyhedron mesh.
     * 
     * This function calculates the average altitude of the vertices that define each face in the polyhedron mesh.
     * It iterates over all the faces and computes the average Z-coordinate (altitude) of the vertices for each face.
     * The altitude for each face is the mean of the Z-coordinates of the vertices that make up the face.
     * 
     * @param mesh The 3D polyhedron mesh for which the altitudes of the faces are computed.
     * @return A map where the key is a facet (face) and the value is the average altitude (Z-coordinate) of the face.
     * 
     */
    Facet_double_map computeAltitudes(const Polyhedron & mesh)
    {
        Facet_double_map altitudes;

        for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i) 
        {
            double altitude_sum = 0.0;
            int vertex_count = 0;

            Halfedge_facet_circulator h = i->facet_begin();
            Halfedge_facet_circulator end = h;
            
            do 
            {
                CGAL::Point_3 p = h->vertex()->point();
                altitude_sum += p.z();  
                ++vertex_count;
                ++h;
            } while (h != end);

            double altitude = altitude_sum / vertex_count;

            altitudes.insert({i, altitude});
        }

        return altitudes;
    }


    /**
     * @brief Computes the roughness of each facet in the polyhedron mesh based on the angle between adjacent facets.
     * 
     * This function calculates the roughness of each facet in the polyhedron mesh by measuring the angle between 
     * the normal vectors of adjacent facets. The angle between two adjacent facets is determined using the 
     * scalar product of their normal vectors. The roughness value is the angle between adjacent facets, 
     * which represents how much the adjacent faces differ in orientation.
     * 
     * The roughness for each facet is computed as the angle (in degrees) between its normal vector and the 
     * normal vector of the next adjacent facet in the mesh. If two facets are nearly aligned, the roughness 
     * will be small; if the facets are oriented differently, the roughness will be larger.
     * 
     * @param mesh The polyhedron mesh for which the roughness of each facet is computed.
     * @return A map where the key is a facet (face) and the value is the roughness of that facet in degrees.
     * 
     * @note The roughness between the first and last facets in the mesh is set to zero.
     */
    Facet_double_map computeRoughness(const Polyhedron & mesh)
    {
        Facet_double_map roughness;

        Facet_iterator i = mesh.facets_begin();

        CGAL::Point_3 p1 = i->facet_begin()->vertex()->point();
        CGAL::Point_3 p2 = i->facet_begin()->next()->vertex()->point();
        CGAL::Point_3 p3 = i->facet_begin()->next()->next()->vertex()->point();
        CGAL::Vector_3 normal1 = CGAL::normal(p1, p2, p3);

        for (; std::next(i) != mesh.facets_end(); ++i)
        {
            Facet_iterator next_face = std::next(i);

            p1 = next_face->facet_begin()->vertex()->point();
            p2 = next_face->facet_begin()->next()->vertex()->point();
            p3 = next_face->facet_begin()->next()->next()->vertex()->point();
            CGAL::Vector_3 normal2 = CGAL::normal(p1, p2, p3);

            double angle = std::acos(CGAL::scalar_product(normal1, normal2) /(std::sqrt(normal1.squared_length()) * std::sqrt(normal2.squared_length()))) * 180 / M_PI;

            roughness.insert({i, angle});

            normal1 = normal2;
        }

        auto last = std::prev(mesh.facets_end());
        p1 = last->facet_begin()->vertex()->point();
        p2 = last->facet_begin()->next()->vertex()->point();
        p3 = last->facet_begin()->next()->next()->vertex()->point();
        CGAL::Vector_3 normal2 = CGAL::normal(p1, p2, p3);

        roughness.insert({last, 0});

        return roughness;
    }


    /**
     * @brief Computes the segmentation of the facets in the polyhedron mesh based on a threshold value.
     * 
     * This function segments the facets of a polyhedron mesh into two classes based on a given threshold value.
     * Each facet is classified into one of two classes depending on whether its associated measure exceeds the threshold.
     * The function uses the provided `class1` and `class2` names to label the facets accordingly.
     * 
     * @param measures A map where each facet is associated with a measure value (e.g., area, angle, roughness, etc.).
     * @param treshold The threshold value used to classify the facets into one of two classes.
     * @param class1 The label to assign to facets whose measure exceeds the threshold.
     * @param class2 The label to assign to facets whose measure does not exceed the threshold.
     * @return A map where each facet is associated with a class label (`class1` or `class2`) based on the threshold comparison.
     * 
     */
    Facet_string_map computeSegmentation(const Facet_double_map & measures, float treshold, std::string class1, std::string class2)
    {
        Facet_string_map segmentation; 

        for(const auto & measure : measures)
        {
            if(measure.second > treshold)
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

    /**
     * @brief Computes and exports a colored mesh to an OFF file based on segmentation.
     * 
     * This function generates an OFF file representation of the polyhedron mesh, where each facet is colored
     * according to its assigned class from a segmentation map. Two colors are provided, and facets are colored 
     * based on whether they belong to `class1` or `class2`. The output mesh is saved in an OFF file format, 
     * with color information embedded for visualization.
     * 
     * @param mesh The polyhedron mesh to be processed.
     * @param segmentationMap A map that associates each facet with its corresponding class label (`class1` or `class2`).
     * @param filenameOFF The path to the output OFF file where the colored mesh will be saved.
     * @param color1 The RGB color to assign to facets belonging to `class1`.
     * @param color2 The RGB color to assign to facets belonging to `class2`.
     * @param class1 The label used for the first class of facets.
     * @param class2 The label used for the second class of facets.
     * 
     */
    void computeColoredMesh(const Polyhedron & mesh, const Facet_string_map & segmentationMap, const std::string & filenameOFF, const Color & color1, const Color & color2, std::string class1, std::string class2)
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
                in_myfile << std::fixed << std::setprecision(1) << ' ' << color1[0] << ' ' << color1[1] << ' ' << color1[2];
            }
            else
            {
                in_myfile << std::fixed << std::setprecision(1) << ' ' << color2[0] << ' ' << color2[1] << ' ' << color2[2];
            }

            in_myfile << std::endl;
        }

        in_myfile.close();

        std::cout << "Segmentation mesh successfully exported at path: " << filenameOFF << " !" << std::endl;
    }

    /**
     * @brief Computes and exports a multi-class colored mesh to an OFF file based on segmentation.
     * 
     * This function generates an OFF file representation of the polyhedron mesh, where each facet is colored
     * based on its associated class from a multi-class segmentation map. Multiple classes and their corresponding
     * colors are provided, and each facet can belong to one or more classes. The colors for each facet are averaged 
     * if it belongs to multiple classes. The output mesh is saved in an OFF file format with color information.
     * 
     * @param mesh The polyhedron mesh to be processed.
     * @param segmentationMultimap A multimap that associates each facet with one or more class labels.
     * @param filenameOFF The path to the output OFF file where the colored mesh will be saved.
     * @param colors A vector containing RGB color values for each class.
     * @param classes A vector containing the class labels.
     * 
     */
    void computeColoredMeshMultiClass(const Polyhedron & mesh, const Facet_string_multimap & segmentationMultimap, const std::string & filenameOFF, const std::vector<Color> & colors, const std::vector<std::string> & classes)
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


            auto values = segmentationMultimap.equal_range(i);
            float r=0, g=0, b=0;
            
            for (auto it = values.first; it != values.second; ++it) 
            {
                auto val = std::find(classes.begin(), classes.end(), it->second);

                if(val != classes.end())
                {
                    int index = std::distance(classes.begin(), val);
                    r += colors[index][0];
                    g += colors[index][1];
                    b += colors[index][2];
                }
            }
            r /= (classes.size()/4);
            g /= (classes.size()/4);
            b /= (classes.size()/4);
            in_myfile << std::fixed << std::setprecision(1) << ' ' << r << ' ' << g << ' ' << b;

            in_myfile << std::endl;
        }

        in_myfile.close();

        std::cout << "Segmentation mesh multiclass successfully exported at path: " << filenameOFF << " !" << std::endl;
    }
}