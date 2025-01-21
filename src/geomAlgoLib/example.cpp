#include "example.hpp"

#include <iostream>


namespace geomAlgoLib
{
    
int computeGenus(const Polyhedron &mesh)
{
	unsigned int nbVerts = 0;
	for (Vertex_iterator i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		++nbVerts;
	}
	std::cout << "# Vertices : " << nbVerts << std::endl;

	unsigned int nbEdges = 0;
	for (Halfedge_iterator i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i)
	{
		++nbEdges;
	}
	nbEdges /= 2;
	std::cout << "# Edges: " << nbEdges << std::endl;

	unsigned int nbFaces = 0;
	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		++nbFaces;
	}
	std::cout << "# Faces: " << nbFaces << std::endl;

	unsigned int euler = nbVerts - nbEdges + nbFaces;
	unsigned int genus = (2 - euler) / 2;

	return genus;
}

}