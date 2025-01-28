#include <io.hpp>
#include <example.hpp>
#include <measures.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <types.hpp>

int main(int argc, char *argv[]){

    std::cout << "Hello !" << std::endl;

    if(argc < 2){
        throw std::invalid_argument("This program expects at least 1 argument (path to a mesh).");
    }

    const std::string meshPath = std::string{argv[1]};
    
    geomAlgoLib::Polyhedron myMesh;

    geomAlgoLib::readOFF(meshPath, myMesh);

    auto genus = geomAlgoLib::computeGenus(myMesh);
    std::cout << "The Genus of [" << meshPath << "] is = " << std::to_string(genus) << std::endl;

    geomAlgoLib::Facet_double_map areas = geomAlgoLib::computeArea(myMesh);
    geomAlgoLib::Facet_double_map angles = geomAlgoLib::computeSmallestAngle(myMesh);
    geomAlgoLib::Facet_doubleTab_map anglesNormal = geomAlgoLib::computeAngleNormalAxes(myMesh);

    auto itAreas = areas.begin();
    auto itAngles = angles.begin();
    auto itAnglesNormal = anglesNormal.begin();

    for(int i = 0; i < areas.size(); i++)
    {
        std::cout << "---Face " << i+1 << "---" << std::endl;
        std::cout << "Area = " << itAreas->second << std::endl;
        std::cout << "Smallest angle = " << itAngles->second << " (degrees)" << std::endl;
        std::cout << "Angles between normal and axes : X = " << itAnglesNormal->second[0] <<" Y = " << itAnglesNormal->second[1] << " Z = " << itAnglesNormal->second[2] << " degrees" << std::endl;
        std::cout << "------------" << std::endl;

        itAreas ++; 
        itAngles ++;
        itAnglesNormal ++;
    }

    geomAlgoLib::writeOFF(myMesh,"output.off");
    geomAlgoLib::computeSmallestAngle(myMesh);

    geomAlgoLib::computeOFFFile(myMesh, areas, "../data/coloredMesh.off");

    std::cout << "The end..." << std::endl;
    return 0;
}