#include "io.hpp"

namespace geomAlgoLib
{

bool readOFF(const std::string& filePath, Polyhedron& mesh){

    std::ifstream input(filePath);

	if (!input || !(input >> mesh) || mesh.is_empty())
	{
		std::cerr << "Invalid .OFF file" << std::endl;
		return false;
	}

    return true;

}

void writeOFF(const Polyhedron& mesh, const std::string& filePath)
{
	std::ofstream in_myfile;
	in_myfile.open(filePath);

	CGAL::set_ascii_mode(in_myfile);

	in_myfile << "COFF" << std::endl // "COFF" makes the file support color information
			  << mesh.size_of_vertices() << ' ' 
			  << mesh.size_of_facets() << " 0" << std::endl; 
			  // nb of vertices, faces and edges (the latter is optional, thus 0)

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

		in_myfile << std::endl;
	}

	in_myfile.close();

	std::cout << "Successfully exported at path: " << filePath << " !" << std::endl;
}


}